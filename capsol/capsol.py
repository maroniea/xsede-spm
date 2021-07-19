

import numpy as np

# R, Z = np.meshgrid(rs, zs) # z_size by r_size.


# 
def sphere(r, z, Rtip):
    """Check if data points are in the tip sphere."""
    return (r**2 + (z - Rtip)**2) <= Rtip**2

def cone(r, z, Rtip, theta, hcone):
    return np.where(np.logical_and(z >= Rtip, z<hcone), r <= (Rtip + (z - Rtip)*np.sin(theta)), False)

def body(r, z, hcone, dcant, rcant):
    return np.where((z>=hcone) * (z<(hcone+dcant)), r <= rcant, False)

def sample(z, n_sample):
    n, m = z.shape
    return np.s_[1:1+n_sample, :]

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
    return r, dr


def gapsam_gridsize(h0, hsam, z):
    n_sample=int(np.round(hsam/h0))
    n_gap=int(np.round(z/h0))
    gapsam=n_sample + n_gap + 1 
    return gapsam

def generate_gapsam_grid(h0, hsam, z):
    N=gapsam_gridsize(h0, hsam, z)
    z_minus=np.arange(N)*h0 + -(z+hsam) + -h0
    return z_minus


def setup_probe(R, Z, Rtip, theta, hcone, dcant):
    # Tip points if code, sphere, body
    spm_tip = ( cs.cone(R, Z, Rtip, theta, hcone) + cs.sphere(R, Z, Rtip)
                + cs.body(R, Z, hcone, dcant, rcant) )
    
    # Need to make sure there is one data point in the cantilever...
    
    z_grid = Z[:, 0] # A single Z column 

    b = np.ones_like(R) # set the boundaries
    u = np.zeros_like(R)
    u[spm_tip] = 1.0
#     Ra = Rtip * (1 - np.sin(theta))
#     Rc = Rtip * np.cos(theta)
    probe_top = np.zeros_like(z_grid, dtype=int)


    # -1 means that the probe is not directly above this data point at all
    probe_bot = np.ones_like(r, dtype=int)*-1

    for i, row in enumerate(spm_tip.T):
        if any(row):
            probe_bot[i] = np.argmax(row) # First z value where the probe is found
    

    # Specify boundaries
    b[0, :] = 0.0
    b[-1, :] = 0.0
    b[:, -1] = 0.0 

    itop = np.argmax(z_grid > (hcone+dcant))
    probe_top[:] = itop
    return b, u, probe_top, probe_bot





        
