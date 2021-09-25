from dataclasses import dataclass, asdict
import streamlit as st
import numpy as np
import capsol.capsol as cap
import matplotlib.pyplot as plt

st.title("Connecting impedance model to Capsol Simulations")

# Add nice heading for cantilever

st.markdown("""## Cantilever Parameters""")
# st.image(r"output/IMG_9163.jpg", width=300)
k_cantilever = st.number_input(label="k_cantilever N/m" ,  value=3.5)
f0 = st.number_input(label= "Resonance Frequency (Fc) Hz" , value= 66500)

# Heading for tip model parameters...

def set_parameters(model, data_class_as_dict):
    for label, val in data_class_as_dict.items():
        setattr(model, label, st.number_input(label, value=val))


st.markdown("""## Tip and Sample Parameters""")
params_sample = cap.ParamsSample(Nr=900, Nz_plus=900)
set_parameters(params_sample, asdict(params_sample))

cp_sample = cap.CapSolSample(params_sample)

with st.form(key = "sample"):
    submit = st.form_submit_button()

if 'out' not in st.session_state:
    st.session_state.out = None

if submit:
    cp_sample.run()
    st.session_state.out =  cp_sample.c

st.markdown(f"{st.session_state.out} F")


c = 4e-14
# c=st.number_input(label= "C(F)" , value= 4E-14, format="%.3e")
# Czz= st.number_input(label= 'C" (F/m^2)', value= 3.9E-2 , format="%.3e")
# r_alpha= st.number_input(label = "Alpha", value= 0.2 )

# Sample model and parameters...
# Dictionary key: Model name
# Value: dict(parameters=dictionary_of_parameters, image=image, description=description)

@dataclass
class RCModel:
    """The parallel RC sample model."""
    Rsample : float = 1e12
    Csample : float = 4e-14

    def Z(self, f):
        s = 2j*np.pi*f
        return 1.0/(1.0/self.Rsample+s*self.Csample)




if "models" not in st.session_state:
    st.session_state.models = {'RC': dict(model=RCModel(),
                                        description="RC Model assumes the sample is a parallel resistor (Rsample) and capacitor (Csample).",
                                        image="image.png")}

selected_model_name = st.selectbox("Sample Model", list(st.session_state.models.keys()))

model_data = st.session_state.models[selected_model_name]
model = model_data['model']
# Show the picture associated with the model
# model_dict = models[model]

st.markdown(f"""### {selected_model_name}

{model_data['description']}
""")

for label, val in asdict(model).items():
    setattr(model, label, st.number_input(label, value=val, format="%.3e"))
            
st.markdown(
f"""### Model parameters
{asdict(model)}
""")


def H(f, sample_model, tip_model):
    Z = sample_model.Z(f)
    c = tip_model['C']
    s= 1j*f*(2*np.pi)
    return (1/(s*c))/(Z +(1/(s*c)))

f = st.number_input("f (Hz)", min_value=0.00, max_value=1e7,  value=1000.0)

H_at_f = H(f, model, dict(C=c))

st.markdown(f"H({f} Hz) = {H_at_f}")

# Inputs / outputs
# We select a model (selected_model)


# Czz_q= (1-r_alpha)*Czz
# dCzz= r_alpha*Czz

# Czz_q
# dCzz

# curvature_min= f0/(4*k_cantilever)*Czz_q
# curvature_max=f0/(4*k_cantilever)*Czz

# curvature_min
# curvature_max

# def Hfunc(f, c, c_sample, r_sample):
#     s= 1j*f*(2*np.pi)
#     z= 1/(1/r_sample + s*c_sample)
#     h= (1/(s*c))/(z +(1/(s*c)))
#     return h

# def hbarfunc(c, c_sample, r_sample, fm, f0):
#     hplus= Hfunc(fm+f0,c, c_sample, r_sample)
#     hminus=Hfunc(fm-f0, c, c_sample, r_sample)  
#     hbar=(1/2)*(hplus + hminus)
#     return hbar


# f=np.logspace(0, 6)  
# fm=np.logspace(0,6)

# h= Hfunc(f,c, c_sample, r_sample)
# vm=1
# fig,ax= plt.subplots()
# ax.loglog(f,h)
# st.pyplot(fig)

# h_fm= Hfunc(fm, c, c_sample, r_sample)
# h_fmsq=h_fm**2
# hbar=hbarfunc(c, c_sample, r_sample, fm, f0)
# df= -f0*(vm**2/(8*k_cantilever)*(Czz_q + (dCzz*hbar))*(h_fmsq))
# df.real 
# fig,ax=plt.subplots()
# ax.semilogx((fm/f0), abs(df.real))
# ax.semilogx((fm/f0), df.imag)
# st.pyplot(fig)
