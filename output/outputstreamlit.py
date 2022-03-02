
import streamlit as st
import output_streamlit_orig
import impedance


st.title("Kelvin Probe Force Microscopy Impedance Model")

apps = {"Electrochemical Impedance Fitting": impedance,
        "Original KPFM analysis": output_streamlit_orig, 
}

app = st.selectbox("Choose an application:", list(apps.keys()))

apps[app].run()