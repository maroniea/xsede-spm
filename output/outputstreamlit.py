
import streamlit as st
import output_streamlit_orig
import impedance
import impedance_kpfm


st.title("Kelvin Probe Force Microscopy Impedance Model")

apps = {"Electrochemical Impedance Fitting": impedance,
        "Capsolpy KPFM Simulation": output_streamlit_orig, 
        "KPFM Impedance Model": impedance_kpfm
}

app = st.selectbox("Choose an application:", list(apps.keys()))

apps[app].run()