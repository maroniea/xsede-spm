
import streamlit as st
import output_streamlit_orig


st.title("Kelvin Probe Force Microscopy Impedance Model")

apps = {"Original KPFM analysis": output_streamlit_orig, 
}

app = st.selectbox("Choose an application:", list(apps.keys()))

apps[app].run()