
# This file tests that capsolpy works correctly

import numpy as np
import capsolpy




print(capsolpy.capsolcyl.__doc__)



params = dict(n=500,
m=500,
l_js=20,
h0=0.5,
rho_max=1e6,
z_max=1e6,
d_steps=20,
rtip=20.0,
theta_deg=15.0,
hcone=15000.0,
rcant=15000.0,
dcant=500.0,
eps_r=5.0,
hsam=10,
method_in='NOSOLVE',
test=0,
verbose=2
)


c_unitless, hess = capsolpy.capsolcyl(**params) # Pass a dictionary of parameters


print(hess.shape)