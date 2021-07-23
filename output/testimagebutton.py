from typing import Dict
import streamlit as st
import numpy as np
html = """
  <style>
    /* Disable overlay (fullscreen mode) buttons */
    .overlayBtn {
      display: none;
    }

    /* Remove horizontal scroll */
    .button-container {
      width: auto !important;
    }

    .fullScreenFrame > div {
      width: auto !important;
    }

    /* 2nd thumbnail */
    .button-container:nth-child(4) {
      top: -266px;
      left: 350px;
    }

    /* 1st button */
    .button-container:nth-child(3) {
      left: 10px;
      top: -60px;
    }

    /* 2nd button */
    .button-container:nth-child(5) {
      left: 360px;
      top: -326px;
    }
  </style>
"""
if 'model' not in st.session_state: 
    st.session_state.model= None

st.title("TEST")
st.title("Sample Diagrams")
st.markdown(html, unsafe_allow_html=True)

st.markdown('<div class="button-container">', unsafe_allow_html=True)
st.image(r"IMG_9163.jpg", width=300)
RC=st.button("Show", key="1")

st.image("https://www.w3schools.com/howto/img_forest.jpg", width=300)
st.button("Show", key="2")
st.markdown('</div>', unsafe_allow_html= True)
if RC:
    st.session_state.model=RC
st.text("Circuit chosen:" + str(st.session_state.model))

from dataclasses import dataclass, field
@dataclass
class RCModel:
  name: str = "RC"
  image_fname: str =""
  params: dict=field(default_factory=lambda : dict(Rsample=1e12, Csample=1e-15))
  def Z(self, f):
    s= 1j*f*(2*np.pi)
    z= 1/(1/self.params["Rsample"] + s*self.params["Csample"])
    return z

@dataclass
class Series_ResistanceModel:
  name: str = "Series Resistance"
  image_fname: str =""
  params: dict=field(default_factory=lambda : dict(Rsample=1e12, Rseries=1e12 ,Csample=1e-15))
  def Z(self, f):
    s= 1j*f*(2*np.pi)
    z= 1/(1/self.params["Rsample"] + s*self.params["Csample"]) + self.params["Rseries"]
    return z

models=[RCModel(), Series_ResistanceModel()]
models[1].params

buttons= {}
for model in models:
  buttons[model.name]=(st.button(f'Select {model.name}'))

for button_name, button in buttons.items():
  if button:
    st.session_state.model= button_name
st.text("Circuit chosen:" + str(st.session_state.model))