
import streamlit as st
import base64
import pandas as pd
import numpy as np
import io
import tempfile




# def create_file( suffix='.png'):
#     with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmpfile:
#         fig.savefig(tmpfile.name, format="png", dpi=300)
#     return tmpfile


# def download_figure(label, fig, default_filename, suffix='.png'):
#     filename = st.text_input("Filename:", default_filename)
#     download_figure = st.download_button(label, filename=filename+suffix,)




def process_file(f):
    data = None
    if f.name.endswith("csv"):
        data = pd.read_csv(f)
    elif f.name.endswith("xlsx") or f.name.endswith("xls"):
        data = pd.read_excel(f)
    elif f.name.endswith("Absorbance"):
        raw_data = np.loadtxt(f, skiprows=19, max_rows=2048)
        data = pd.DataFrame(raw_data, columns=["Wavelength (nm)", "Absorbance"])
    elif f.name.endswith("Transmittance"):
        raw_data = np.loadtxt(f, skiprows=19, max_rows=2048 )
        data = pd.DataFrame(raw_data, columns=["Wavelength (nm)", "Transmittance"])
    elif f.name.endswith("txt"):
        # UV-1800 Data Print
        data = pd.read_csv(f, skiprows=1)
    else:
        raise NotImplementedError(f"Data loading not supported for file {f.name}")
    return data


def find(val, array):
    return np.argmin(abs(array - val))


def write_excel(df, filename, label="Download Excel file"):
    towrite = io.BytesIO()
    downloaded_file = df.to_excel(towrite, encoding='utf-8', index=False, header=True)
    towrite.seek(0)  # reset pointer
    b64 = base64.b64encode(towrite.read()).decode()  # some strings
    linko = f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}" download="{filename}.xlsx">{label}</a>'
    st.markdown(linko, unsafe_allow_html=True)