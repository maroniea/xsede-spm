
# This file tests that capsolpy works correctly

import numpy as np
from capsolpy import defcapsol
from copy import copy



make_hess = defcapsol.make_hess


print(defcapsol.generategrid.__doc__)
print(make_hess.__doc__)




params = dict(
n=500,
m=500,
l_js=20,
h0=0.5,
rho_max=1e6,
z_max=1e6,
z_steps=20, # Actual distance away
rtip=20.0,
theta_deg=15.0,
hcone=15000.0,
rcant=15000.0,
dcant=500.0,
eps_r=5.0,
hsam=10,
verbose=2
)

params['l'] = params['l_js'] + params['z_steps']
params['theta'] = params['theta_deg'] * 180 / np.pi

params_grid = copy(params)
del params_grid['l_js']
del params_grid['theta_deg']
del params_grid['z_steps']


# c_unitless, hess = make_hess(**params) # Pass a dictionary of parameters

hn,r,hm,zm = defcapsol.generategrid(**params_grid)

# print(hess.shape)