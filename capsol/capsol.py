

import numpy as np
from numpy.lib.arraysetops import union1d
from scipy import linalg
from scipy import sparse
from scipy.sparse import linalg as la
from tqdm import tqdm
import capsol.newanalyzecapsol as nac
from datetime import datetime as dt

from dataclasses import dataclass

# R, Z = np.meshgrid(rs, zs) # z_size by r_size.


# 
def sphere(r, z, Rtip):
    """Check if grid points are in the tip sphere.
    
    Parameters:
    
    r: """
    return (r**2 + (z - Rtip)**2) <= Rtip**2

def cone(r, z, Rtip, theta, hcone):
    """Check if grid points are in the tip cone."""
    return np.where(np.logical_and(z >= Rtip, z<hcone), r <= (Rtip + (z - Rtip)*np.sin(theta)), False)

def body(r, z, hcone, dcant, rcant):
    return np.where((z>=hcone) * (z<(hcone+dcant)), r <= rcant, False)

def sample(z, n_sample):
    n, m = z.shape
    return np.s_[1:1+n_sample, :]

def epsilon_z(z, d, eps_r):
    """Returns epsilon_r(z) defined on the staggered grid.
    
    Ex: For d = 1.0 nm, a grid spaced every 0.5 nm, and eps_r = 3.0, we have

    z_ind    z                  epsz_ind     eps_z
    0       -2.0 -----------    
                                    0          3.0 
    1       -1.5 -- sample --
                                    1          3.0
    2       -1.0 ------------
                                    2          1.0
    3        -0.5 -- vaccuum --
                                    3          1.0
    4           0 -------------     
                        tip         
    

    Note that the array returned is always 1 shorter than the input array of z values.
    """
    tol = 1e-8 # Add a small margin for rounding error
    return np.where(z <= (-d+tol), eps_r, 1.0)[:-1]



def geometric_sum(r, a1, n):
    return sum(np.array([a1 * r ** i for i in range(1, n + 1)]))

def find_ratio(total, h0, terms):
    r = 1.00001
    s = 0
    while s < total:
        s = geometric_sum(r, h0, terms)
        r += 0.00001
    return r

def guni_grid(nuni, n, h0, grid_max):
    """Generates geometric grid which can be used for the radial and z directions. 
    
    nuni: number of uniform grid points 
    n: number of gird points
    h0: initial grid spacing
    grid_max: rho_max or z_max, maximum grid size"""
    
    
    r = np.zeros(n+1)
    dr = np.zeros(n) # Could also use n+1
    dr[:nuni] = h0
    r[1:nuni+1] = np.cumsum(dr[:nuni])
    r_left = grid_max-(nuni*h0)
    n_left=n-nuni
    ratio= find_ratio(r_left, h0, n_left)
    for i in range(n_left):
        dr[nuni+i]= h0*ratio**(i+1)
    r[1:] = np.cumsum(dr[:])
    return r, ratio


def gapsam_gridsize(h0, hsam, z):
    n_sample=int(np.round(hsam/h0))
    n_gap=int(np.round(z/h0))
    gapsam=n_sample + n_gap # No need for an extra point - use the bottom of the sample as the boundary
    return gapsam

def generate_gapsam_grid(h0, hsam, z):
    N=gapsam_gridsize(h0, hsam, z)
    z_minus=np.arange(N)*h0 + -(z+hsam)
    return z_minus


# def setup_probe(R, Z, Rtip, theta, hcone, dcant):
#     # Tip points if code, sphere, body
#     spm_tip = ( cs.cone(R, Z, Rtip, theta, hcone) + cs.sphere(R, Z, Rtip)
#                 + cs.body(R, Z, hcone, dcant, rcant) )
    
#     # Need to make sure there is one data point in the cantilever...
    
#     z_grid = Z[:, 0] # A single Z column 

#     b = np.ones_like(R) # set the boundaries
#     u = np.zeros_like(R)
#     u[spm_tip] = 1.0
# #     Ra = Rtip * (1 - np.sin(theta))
# #     Rc = Rtip * np.cos(theta)
#     probe_top = np.zeros_like(z_grid, dtype=int)


#     # -1 means that the probe is not directly above this data point at all
#     probe_bot = np.ones_like(r, dtype=int)*-1

#     for i, row in enumerate(spm_tip.T):
#         if any(row):
#             probe_bot[i] = np.argmax(row) # First z value where the probe is found
    

#     # Specify boundaries
#     b[0, :] = 0.0
#     b[-1, :] = 0.0
#     b[:, -1] = 0.0 

#     itop = np.argmax(z_grid > (hcone+dcant))
#     probe_top[:] = itop
#     return b, u, probe_top, probe_bot


def boundary_radial(Nr, Nz):
    """Returns points that are on the cylindrical coordinate boundary. As a 2D array u_ij,
    where i is the z-coordinate index and j is the radius index, this returns a boundary array where points 
    where i = 0 or Nz - 1 or j = Nr - 1 return True.
    """
    bc = np.zeros(Nr*Nz,dtype=bool)
    for i in range(Nz):
        for j in range(Nr):
            if i == 0 or i == (Nz-1) or j == (Nr-1):
                ind = i * Nr + j # This point
                bc[ind] = True
    return bc


def poisson_variable_spacing_radial(x, y):
    Nx = len(x)
    Ny = len(y)
    hx = np.diff(x)
    hy = np.diff(y)
    A = sparse.lil_matrix((Nx*Ny, Nx*Ny))
    for i in range(Ny):
        for j in range(Nx): # Radial
            ind = i * Nx + j # This point
            ixp = ind + 1    # +x 
            ixn = ind - 1    # -x
            iyp = (i+1)*Nx + j  # +y
            iyn = (i-1)*Nx + j  # -y
            
            
            Dx_plus = hx[j] if j < (Nx-1) else 0.0
            Dx_minus = hx[j-1] if j > 0 else hx[j]
            x0 = x[j]
            Dy_plus = hy[i] if i < (Ny-1) else 0.0
            Dy_minus = hy[i-1] if i > 0 else 0.0
            
            prefactor_x = 4/((Dx_plus+Dx_minus)*(Dx_plus**2 + Dx_minus**2))
            
            prefactor_y = 4/((Dy_plus+Dy_minus)*(Dy_plus**2 + Dy_minus**2))

            
            A[ind, ind] = (Dx_plus+Dx_minus) * prefactor_x + (Dy_plus+Dy_minus) * prefactor_y
            if j == 0:
                A[ind, ixp] = -2 * Dx_minus * prefactor_x # That's it, no radial derivative here...
            elif j < (Nx - 1):
                A[ind, ixp] = -1 * Dx_minus * prefactor_x + -1 / (x0 * (Dx_plus+Dx_minus))
            
                
            
            if j > 0:
                A[ind, ixn] = -1 * Dx_plus * prefactor_x + 1 / (x0 * (Dx_plus+Dx_minus))
            
            if j == (Nx - 1):
                A[ind, ind] += -1 / (x0 * (Dx_plus+Dx_minus)) # 1st order difference uses the grid point here...

            if i > 0:
                A[ind, iyn] = -1 * Dy_plus * prefactor_y
            if i < (Ny-1):
                A[ind, iyp] = -1 * Dy_minus * prefactor_y
    
    return sparse.csr_matrix(A) # Convert to better format for usage

class arrayBuilder:
    def __init__(self, estimated_size=None):
        self.rows = []
        self.cols = []
        self.data = []
    
    def __call__(self, row, col, val):
        self.rows.append(row)
        self.cols.append(col)
        self.data.append(val)

# This function should build the sparse matrix A faster, but currently seems
# to have an error in it? The matrix is different than the
# matrix produced by `poisson_variable_spacing_radial`
def poisson_variable_spacing_radial_faster(x, y):
    Nx = len(x)
    Ny = len(y)
    hx = np.diff(x)
    hy = np.diff(y)
    ab = arrayBuilder()
    for i in range(Ny):
        for j in range(Nx): # Radial
            ind = i * Nx + j # This point
            ixp = ind + 1    # +x 
            ixn = ind - 1    # -x
            iyp = (i+1)*Nx + j  # +y
            iyn = (i-1)*Nx + j  # -y
            
            
            Dx_plus = hx[j] if j < (Nx-1) else 0.0
            Dx_minus = hx[j-1] if j > 0 else hx[j]
            x0 = x[j]
            Dy_plus = hy[i] if i < (Ny-1) else 0.0
            Dy_minus = hy[i-1] if i > 0 else 0.0
            
            prefactor_x = 4/((Dx_plus+Dx_minus)*(Dx_plus**2 + Dx_minus**2))
            
            prefactor_y = 4/((Dy_plus+Dy_minus)*(Dy_plus**2 + Dy_minus**2))

            
            dia_ind = (Dx_plus+Dx_minus) * prefactor_x + (Dy_plus+Dy_minus) * prefactor_y
            if j == (Nx - 1):
                dia_ind += -1 / (x0 * (Dx_plus+Dx_minus)) # 1st order difference uses the g
            
            ab(ind, ind, dia_ind)

            if j == 0:
                ab(ind, ixp, -2 * Dx_minus * prefactor_x) # That's it, no radial derivative here...
            elif j < (Nx - 1):
                ab(ind, ixp, -1 * Dx_minus * prefactor_x + -1 / (x0 * (Dx_plus+Dx_minus)))
            
                
            
            if j > 0:
                ab(ind, ixn,  -1 * Dx_plus * prefactor_x + 1 / (x0 * (Dx_plus+Dx_minus)))

            if i > 0:
                ab(ind, iyn,  -1 * Dy_plus * prefactor_y)
            if i < (Ny-1):
                ab(ind, iyp,  -1 * Dy_minus * prefactor_y)
    
    return sparse.csr_matrix(sparse.coo_matrix((ab.data, (ab.rows, ab.cols)), shape=(Nx*Ny, Nx*Ny))) # Convert to better format for usage



def grid_area(r, z):
    """
    Parameters:
    
    r : 1d array of radii
    z : 1d array of z coordinates
    
    Returns the area of each grid cell (with z in rows, r in columns)."""
    # Area of each grid element is π * (R_outer^2 - R_inner^2) * (Δz) 
    dr2 = np.diff(r**2)
    return np.pi *  np.diff(z).reshape((-1, 1)) @ dr2.reshape((1, -1))

def E_field(u, r, z):
    """
    u: Voltage on the grid (by convention, rows are z, columns are r)
    r: 1D array of radii 
    z: 1D array of z coordinates
    
    Return: The electric field wrt to r and z, 
    Works for r, z (cylindrical coordinates), or x, y (Cartesian coordinates)."""
    Ny = len(z)
    Nx = len(r)
    Ey1 = np.diff(u.reshape((Ny, Nx)), axis=0) / np.diff(z).reshape((-1, 1))
    Ex1 = np.diff(u.reshape((Ny, Nx)), axis=1) / np.diff(r).reshape((1, -1))
    Ex = 0.5 * (Ex1[:-1, :] + Ex1[1:, :])
    Ey = 0.5 * (Ey1[:, :-1] + Ey1[:, 1:])
    return Ex + 1j*Ey
        


@dataclass
class Params:
    Rtip : float = 20.0
    theta_deg : float = 15.0
    Hcone : float = 15000.0
    Hcant : float = 500.0
    Rcant : float = 15000.0
    zMax : float = Rtip*1000.0
    rhoMax : float = Rtip*1000.0
    h0 : float = Rtip * 0.02
    d : float = Rtip
    Nuni : int = 50 # Uniformly spaced points in the r and z_plus directions
    Nr : int = 500
    Nz_plus : int = 500
    hsam : float = 0.0 # No sample thickness for now...

    
    @property
    def theta(self) -> float:
         return self.theta_deg * np.pi/180



class CapSol:
    def __init__(self, params: Params):
        self.params = params
        self.r, self.r_ratio = guni_grid(params.Nuni, params.Nr, params.h0, params.rhoMax)
        self.z_plus, self.z_ratio = guni_grid(params.Nuni, params.Nz_plus,
                                                params.h0, params.zMax)
        
        self.z_minus = generate_gapsam_grid(params.h0, params.hsam, params.d)

        # Make the final, overall, z grid:
        self.z = np.r_[self.z_minus, self.z_plus]

        self.R, self.Z = np.meshgrid(self.r, self.z)

        self.spm_tip = (sphere(self.R, self.Z, self.params.Rtip) +
                        cone(self.R, self.Z, params.Rtip, params.theta, params.Hcone) +
                        body(self.R, self.Z, params.Hcone, params.Hcant, params.Rcant)
                        )
        
        self.Nr = len(self.r)
        self.Nz = len(self.z)

        self.outer_boundary = boundary_radial(self.Nr, self.Nz)

        self.boundary = self.spm_tip.ravel() + self.outer_boundary

        self.u = np.zeros_like(self.R)
        self.u[self.spm_tip] = 1.0

    
    def setup_matrices(self):
        self.A = poisson_variable_spacing_radial(self.r, self.z)

        self.f = -self.A @ self.u.ravel()

        self.A_free = self.A[~self.boundary].T[~self.boundary].T

        self.f_free = self.f[~self.boundary]
    
    def solve(self):
        u_cut = la.spsolve(self.A_free, self.f_free)
        self.u = self.u.ravel()
        self.u[~self.boundary] = u_cut
        self.u = self.u.reshape((self.Nz, self.Nr))


    def process(self):
        self.dV = dV =  grid_area(self.r, self.z)

        self.energy = 0.5 * np.sum(dV * abs(E_field(self.u, self.r, self.z))**2) * 1e-9 * 8.854e-12

        self.energy_z = 0.5 * np.sum(dV * E_field(self.u, self.r, self.z).imag**2) * 1e-9 * 8.854e-12

        self.c=self.energy*2

        return self.c # In SI Units...

    def __repr__(self):
        return f"CapSol(params={repr(self.params)})"

class SphereTest:
    def __init__(self, params: Params):
        self.params = params
        self.r, self.r_ratio = guni_grid(params.Nuni, params.Nr, params.h0, params.rhoMax)
        self.z_plus, self.z_ratio = guni_grid(params.Nuni, params.Nz_plus,
                                                params.h0, params.zMax)
        
        self.z_minus = generate_gapsam_grid(params.h0, params.hsam, params.d)

        # Make the final, overall, z grid:
        self.z = np.r_[self.z_minus, self.z_plus]

        self.R, self.Z = np.meshgrid(self.r, self.z)

        self.spm_tip = sphere(self.R, self.Z, self.params.Rtip)
                        
        
        self.Nr = len(self.r)
        self.Nz = len(self.z)

        self.outer_boundary = boundary_radial(self.Nr, self.Nz)

        self.boundary = self.spm_tip.ravel() + self.outer_boundary

        self.u = np.zeros_like(self.R)
        self.u[self.spm_tip] = 1.0

    
    def setup_matrices(self):
        self.A = poisson_variable_spacing_radial(self.r, self.z)

        self.f = -self.A @ self.u.ravel()

        self.A_free = self.A[~self.boundary].T[~self.boundary].T

        self.f_free = self.f[~self.boundary]
    
    def solve(self):
        u_cut = la.spsolve(self.A_free, self.f_free)
        self.u = self.u.ravel()
        self.u[~self.boundary] = u_cut
        self.u = self.u.reshape((self.Nz, self.Nr))


    def process(self):
        self.dV = dV =  grid_area(self.r, self.z)

        self.energy = 0.5 * np.sum(dV * abs(E_field(self.u, self.r, self.z))**2) * 1e-9 * 8.854e-12

        self.energy_z = 0.5 * np.sum(dV * E_field(self.u, self.r, self.z).imag**2) * 1e-9 * 8.854e-12

        self.c=self.energy*2

        return self.c # In SI Units...

    def __repr__(self):
        return f"CapSol(params={repr(self.params)})"

def Totalsim(params, dmin, dmax, istep, fname, Test=0):
    capacitances=[]
    distances=np.arange(dmin, dmax, istep*params.h0)
    print(distances)
   
    for d in tqdm(distances):
        start_time= dt.now()
        params.d= d
        if Test==1:
            sim=SphereTest(params)
        else:
           sim=CapSol(params) 
        sim.setup_matrices()
        sim.solve()
        sim.process()
        capacitances.append(sim.c)
        end_time=dt.now()
        elapsed_time= end_time-start_time
        print(elapsed_time)
    np.savetxt(fname, np.c_[distances, capacitances], header='distance (nm) Capacitances(F)', footer=f'Totalsim(params={params}, dmin={dmin}, dmax={dmax}, istep={istep}, fname={fname}), Test={Test}')

    return distances, capacitances

def runnewcapsol(input_fname= "capsol.in", output_fname="C-Z.dat"):
    gp=nac.get_gridparameters(input_fname)
    params=Params(Rtip=gp["Rtip"], theta_deg=gp["half-angle"],Hcone=gp["HCone"], Hcant=gp["thickness_Cantilever"], Rcant=gp["RCantilever"], zMax=gp["z_max"],rhoMax=gp["rho_max"], h0=gp["h0"], d=gp["min"], Nuni=gp["Nuni"], Nr=gp["n"], Nz_plus=gp["m+"],hsam=gp["Thickness_sample"])
    totalsim=Totalsim(params, gp["min"], gp["max"], gp["istep"], output_fname, gp["Test"])
    return totalsim