from numpy.lib.index_tricks import AxisConcatenator
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
st.title("Connecting impedance model to Capsol Simulations")
k_cantilever = st.number_input(label="k_cantilever N/m" ,  value=3.5)
f0 = st.number_input(label= "Resonance Frequency (Fc) Hz" , value= 66500)
c=st.number_input(label= "C(F)" , value= 4E-14, format="%.3e")
Czz= st.number_input(label= 'C" (F/m^2)', value= 3.9E-2 , format="%.3e")
r_alpha= st.number_input(label = "Alpha", value= 0.2 )
c_sample=st.number_input(label= 'C sample', value=4E-14, format="%.3e")
r_sample= st.number_input(label= 'R sample', value=1E12, format='%.3e')
Czz_q= (1-r_alpha)*Czz
dCzz= r_alpha*Czz
Czz_q
dCzz
curvature_min= f0/(4*k_cantilever)*Czz_q
curvature_max=f0/(4*k_cantilever)*Czz
curvature_min
curvature_max

def Hfunc(f, c, c_sample, r_sample):
    s= 1j*f*(2*np.pi)
    z= 1/(1/r_sample + s*c_sample)
    h= (1/(s*c))/(z +(1/(s*c)))
    return h
def hbarfunc(c, c_sample, r_sample, fm, f0):
    hplus= Hfunc(fm+f0,c, c_sample, r_sample)
    hminus=Hfunc(fm-f0, c, c_sample, r_sample)  
    hbar=(1/2)*(hplus + hminus)
    return hbar    
f=np.logspace(0, 6)  
fm=np.logspace(0,6)  
h= Hfunc(f,c, c_sample, r_sample)
vm=1
fig,ax= plt.subplots()
ax.loglog(f,h)
st.pyplot(fig)
h_fm= Hfunc(fm, c, c_sample, r_sample)
h_fmsq=h_fm**2
hbar=hbarfunc(c, c_sample, r_sample, fm, f0)
df= -f0*(vm**2/(8*k_cantilever)*(Czz_q + (dCzz*hbar))*(h_fmsq))
df.real 
fig,ax=plt.subplots()
ax.semilogx((fm/f0), abs(df.real))
ax.semilogx((fm/f0), df.imag)
st.pyplot(fig)
