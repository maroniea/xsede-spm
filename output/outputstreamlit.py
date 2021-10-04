from dataclasses import dataclass, asdict
import streamlit as st
import numpy as np
import capsol.capsol as cap
import matplotlib.pyplot as plt


# https://discuss.streamlit.io/t/cannot-print-the-terminal-output-in-streamlit/6602/9?u=ryanpdwyer
from contextlib import contextmanager
from io import StringIO
from streamlit.report_thread import REPORT_CONTEXT_ATTR_NAME
from threading import current_thread
import streamlit as st
import sys


@contextmanager
def st_redirect(src, dst):
    placeholder = st.empty()
    output_func = getattr(placeholder, dst)

    with StringIO() as buffer:
        old_write = src.write

        def new_write(b):
            if getattr(current_thread(), REPORT_CONTEXT_ATTR_NAME, None):
                buffer.write(b)
                output_func(buffer.getvalue())
            else:
                old_write(b)

        try:
            src.write = new_write
            yield
        finally:
            src.write = old_write


@contextmanager
def st_stdout(dst):
    with st_redirect(sys.stdout, dst):
        yield


@contextmanager
def st_stderr(dst):
    with st_redirect(sys.stderr, dst):
        yield



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
params_sample = cap.AllParams(Nr=500, Nz_plus=500, istep=2, dmin=5.0, dmax=8.0, h0=0.5)
with st.expander("CapSol Simulation Parameters"):
    params_to_set = asdict(params_sample)
    for param in ['pt', 'dmax']:
        params_to_set.pop(param)
    set_parameters(params_sample, params_to_set)

params_sample.dmax = params_sample.dmin + params_sample.h0 * params_sample.istep * 3

cp_sample = cap.CapSolAll(params_sample)


with st.form(key = "sample"):
    submit = st.form_submit_button()

st.markdown(f"""### Ratios

Should be 1.01 or less:

r_ratio: {cp_sample.r_ratio:.4f} (increase n to improve)

z_ratio: {cp_sample.z_ratio:.4f} (increase Nz_plus to improve)

""")

if 'out' not in st.session_state:
    st.session_state.out = None

if 'out_z' not in st.session_state:
    st.session_state.out_z = None

if 'out_zz' not in st.session_state:
    st.session_state.out_zz = None

if submit:
    st.markdown("Capsolpy output: ")
    with st_stdout("code"):
        cp_sample.run()
    st.session_state.out =  cp_sample.C
    st.session_state.out_z =  np.gradient(cp_sample.C) / (cp_sample.params.h0*cp_sample.params.istep*1e-9)
    st.session_state.out_zz =  np.gradient(st.session_state.out_z) / (cp_sample.params.h0*cp_sample.params.istep*1e-9)

if st.session_state.out is not None:
    C = st.session_state.out[1]
    Cz = st.session_state.out_z[1]
    Czz = st.session_state.out_zz[1]
    delta_Czz = 2*Cz**2 / C
    Czz_q = Czz - delta_Czz
    alpha_qosc = delta_Czz/Czz
    st.markdown("### CapSolPy results")
    st.markdown(f"$C$ = {C:.4e} F")
    st.markdown(f"$C'$ = {Cz:.4e} F/m")
    st.markdown(f"$C''$ = {Czz:.4e} F/m²")
    st.markdown(f"$C_q''$ = {Czz_q:.4e} F/m²")
    st.markdown(f"$\\Delta C''$ = {delta_Czz:.4e} F/m²")
    st.markdown(f"$\\alpha_\\mathrm{{q-osc}}$ = {alpha_qosc:.3f}")


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
