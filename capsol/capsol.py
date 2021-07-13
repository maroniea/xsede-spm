

import numpy as np

# R, Z = np.meshgrid(rs, zs) # z_size by r_size.


# 
# Check if a point is in the sphere
def sphere(r, z, Rtip):
    return (r**2 + (z - Rtip)**2) <= Rtip**2

def cone(r, z, Rtip, theta, hcone):
    return np.where(np.logical_and(z > Rtip, z<hcone), r < (Rtip + (z - Rtip)*np.sin(theta)), False)

def body(r, z, hcone, dcant, rcant):
    return np.where(np.logical_and(z>hcone, z<(hcone+dcant)), r < rcant, False)



sphere(2.0, 2.0, 20.0)

# def setup_probe(r, z_grid, Rtip, theta, hcone, dcant):

#     b = np.ones_like(R) # First the z grid, then the radial grid
#     Ra = Rtip * (1 - np.sin(theta))
#     Rc = Rtip * np.cos(theta)
#     probe_top = np.zeros_like(z_grid, dtype=int)

#     # Specify boundaries
#     b[0, :] = 0.0
#     b[-1, :] = 0.0
#     b[:, -1] = 0.0 

#     u = np.zeros_like(R)

#     itop = np.argmax(z_grid > (hcone+dcant))
#     probe_top[:] = itop

        
