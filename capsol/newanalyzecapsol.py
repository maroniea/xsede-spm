

# To do
# Read in the .dat file (np.loadtxt)

# Function should take the name of the folder as an input, return the capacitance data?



# Function should take the name of the folder as an input, return the input paramters
# as a dictionary?
import numpy as np
import matplotlib as plt
from scipy import signal 


def get_gridparameters(f_name): 

    params = {}

    with open (f_name) as file_object:
        contents= file_object.readlines()
    line_1 = contents[0]
    line_2 = contents[1]
    line_3 = contents[2]
    line_4 = contents[3]
    line_5 = contents[4]
    line_6 = contents[5]
    line_7= contents[6]
    line_list1= line_1.split()
    line_list2= line_2.split()
    line_list3= line_3.split()
    line_list4= line_4.split()
    line_list5= line_5.split()
    line_list6= line_6.split()
    line_list7=line_7.split()
    params["n"]=(int(line_list1[0])) 
    params["m+"]=(int(line_list1[1]))
    params["m-"]=(int(line_list1[2]))
    params["h0"]=(float(line_list2[0]))
    params["rho_max"]=(float(line_list2[1]))
    params["z_max"]=(float(line_list2[2]))
    params["min"]= (float(line_list3[0]))
    params["max"]=(float(line_list3[1]))
    params["istep"]= (float(line_list3[2]))
    params["Rtip"]=(float(line_list4[0]))
    params["half-angle"]=(float(line_list4[1]))
    params["HCone"]=(float(line_list4[2]))
    params["RCantilever"]=(float(line_list4[3]))
    params["thickness_Cantilever"]=(float(line_list4[4]))
    params["eps_r"]= (float(line_list5[0]))
    params["Thickness_sample"]= (float(line_list5[1]))
    params["Solving Method"]=((line_list6[0]))
    params["Test"]=(int(line_list6[1]))
    params["Verbosity"]= (int(line_list6[2]))
    params["Nuni"]=(int(line_list7[0]))
    if line_list7[1]=="T":
        params["Equally spaced"]=True
    else:
        params["Equally spaced"]=False
    return(params) 
    #units in nm 



def process_data(params, data, smoothing= False, std=5*10**-9, fortran=True):
    if fortran:
        z_over_r= data[: , 0]
        r= params['Rtip']
        z= (z_over_r*r)


        c_over_pie0R= data[:, 2]
        e0= 8.854E-12 
        c= (c_over_pie0R*np.pi*e0*r)*10**-9
    else: 
        z=data[:, 0]
        c=data[:, 1]

    dz= (z[1]-z[0])*10**-9
    if smoothing:
         N_sigma = int(std/dz) #standard deviation in data points
         gaussian_smoothing= signal.gaussian(6*N_sigma, N_sigma)
         c_smoothed= np.convolve(c, gaussian_smoothing, mode='same')
         cz=np.gradient(c_smoothed, dz)
         czz= np.gradient(cz, dz)
         s= slice(3*N_sigma, -3*N_sigma)
         alpha=(2*(cz**2/c_smoothed)/czz)
         return dict(z=z[s] , c=c_smoothed[s] , cz=cz[s] , czz=czz[s] , alpha = alpha[s])
    else:
         cz=np.gradient(c, dz)
         czz= np.gradient(cz, dz)
         alpha=(2*(cz**2/c)/czz)
         return dict(z=z , c=c, cz=cz , czz=czz , alpha = alpha)
# import matplotlib.pyplot as plt
# plt.plot(z , cz)
# plt.ylabel("C'(m)")
# plt.xlabel("Z(m)")
# plt.show()
# fig, ax = plt.subplots()
# ax.plot(z, czz)
# ax.set_ylabel("C''(m)")
# ax.set_xlabel("Z(m)")
# fig.show()
# plt.show()

