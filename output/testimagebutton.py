import streamlit as st
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
