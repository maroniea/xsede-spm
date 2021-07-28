# Make CapSol Simulations
# Use a Terminal / Powershell to run this file using streamlit with the command:
# streamlit run make_simulation.py

# Imports go at the top...
from re import T
import streamlit as st
import numpy as np
import os
import copy
import glob
import itertools

from streamlit.proto.Checkbox_pb2 import Checkbox
import SessionState 


# Helper functions used later in the file

def merge_dicts(dicts):
    """Combine multiple dictionaries.
    
    > merge_dicts([{"a": 2}, {"b": 3}])
    > {"a": 2, "b": 3}
    """
    super_dict = {}
    keys = []
    for d in dicts:
        for key, val in d.items():
            super_dict[key] = val
    return super_dict


def formatter(type):
    """Format floating point numbers in scientific notation
    (integers use normal integer formatting)."""
    if type == int:
        return "%d"
    elif type == float:
        return "%.3e"

def folder_selector(label="Select a file", folder_path='.'):
    filenames = glob.glob(folder_path+'/*/')
    return st.selectbox(label, filenames, format_func=os.path.relpath)


def write_inputfile(p):

    # if p{["Equally spaced"]}:
    #    es="Equally spaced"=True
    #    es==T
    # else:
    #     "Equally spaced" = False
    return f"""   {p['n']}  {p['m+']}   {p['m-']}                                # Number of grids: n, m+, m-
  {p['h0']}   {p['rho_max']}    {p['z_max']}                   # Resolution: h0, Box-size: Rho_max, Z_max    *** ALL LENGTHS IN NANOMETER *** 
     {p['min']}    {p['max']}   {p['istep']}                         # Tip-sample separation: min, max, istep (stepsize=istep*h0)
  {p['Rtip']}  {p['half-angle']}  {p['HCone']}  {p['RCantilever']}  {p['thickness_Cantilever']}  # Probe shape: Rtip,half-angle,HCone,RCantilever, thickness_Cantilever
   {p['eps_r']}      {p['Thickness_sample']}                             # Sample: eps_r, Thickness_sample
   LAPACK  {p['Test']}  {p['Verbosity']} 
   {p['Nuni']}  {str(p['Equally spaced'])[0]}"""

state = SessionState.get(folder=os.path.abspath('.'))

st.title("Create CapSol Simulation Input files")

st.markdown("""Combine scans to create a grid of simulation data points.
    For example, if you specify 3 scans with 4 simulations in each scan,
    4³ = 64 input files will be created.""")

# TO DO: Add the remaining variables and their default values...
default_params = {'n': 500,
 'm+': 500,
 'm-': 20,
 'h0': 0.5,
 'rho_max': 1000000.0,
 'z_max': 1000000.0,
 'min': 2.0,
 'max': 20.0,
 'istep': 2,
 'Rtip': 20.0,
 'half-angle': 15.0,
 'HCone': 15000.0,
 'RCantilever': 40000.0,
 'thickness_Cantilever': 500.0,
 'eps_r': 3.0,
 'Thickness_sample': 10.0,
 'Test': 0,
 'Verbosity': 0,
 'Nuni': 1,
 'Equally spaced': False}

# Infer from the default values
parameter_types = {key: type(val) for key, val in default_params.items()}

param_keys = list(default_params.keys())

st.markdown("## Update defaults")

st.markdown("""
This section contains default values that remain constant throughout the simulations - 
any parameters that you vary will be overwritten by your choices below.
""")

# parameter_types[key] gives int or float depending on the type of the default 

# parameter 

# We want to set step = 1 if the type is int, leave step = None otherwise 
updated_default_params = {}
for key, val in default_params.items():
    if parameter_types[key] == int:
        updated_default_params[key] = st.number_input(key,
                            value=val, step=1, 
                            format=formatter(
                                parameter_types[key]))
    elif parameter_types[key]==bool:
        updated_default_params[key] = st.checkbox(key, value=val)
    else:
        updated_default_params[key] = st.number_input(key,
                            value=val,
                            format=formatter(
                                parameter_types[key]))


st.markdown("## Setup Scans")

n_scans = st.number_input("Number of experimental scans:", 1, 5, 1, 1)

n_simulations = [] # The number of simulations for each scan
n_parameters = [] # The number of parameters varied in each scan
# all_params = [] # A list of list, containing the parameters varied in each scan
# all_values = [] # The values for each param

all_params_values = []

# all_params, all_values could be a dictionary

# Example:

# Scan 1 varying the tip size Rtip (7 possible values from 1 nm to 100 nm: 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0) and scan 2 varying the sample thickness sample_thickness (4 possible values: 1.0 nm, 3.0 nm, 10.0 nm, 30.0 nm)
# n_simulations = [7, 5]
# n_parameters = [1, 1]
# all_params = [['Rtip'], ['sample_thickness']]
# all_values = [[1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0], [1.0, 3.0, 10.0, 30.0]]

# Better is the list of dictionaries:
# all_params_values = [{'Rtip': [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]}, {'sample_thickness': [1.0, 3.0, 10.0, 30.0]}]

for i in range(n_scans):
    st.subheader(f"Scan {i+1}")
    n_sim = st.number_input(f"S{i+1} Number of simulations:", 1, 10, 5)
    n_simulations.append(n_sim)
    n_param = st.number_input(f"S{i+1} Parameters to vary at once:", 1, 5, 1)
    n_parameters.append(n_param)

    cols = st.beta_columns(n_param)

    params = []
    values = [] 

    for c, col in enumerate(cols):

        with col:
            params.append(st.selectbox(label=f"S{i+1} Parameter {c+1}", options=param_keys))

            values.append([st.number_input(f"S{i+1} {params[c]} #{j+1}:",
                            value=updated_default_params[params[c]],
                            format=formatter(
                                parameter_types[params[c]])
                                ) for j in range(n_sim)])

    all_params_values.append([{param: val for param, val in zip(params, val_row)} for val_row in zip(*values)])

st.header("Output information")

base_name = st.text_input("Output base filename")

st.markdown("### Output Folder")

st.markdown(f'**Selected folder:** {state.folder}')

st.write(f'Select another directory (up or down) below:')

left_col, right_col = st.beta_columns(2)
                    
with right_col:
    up_directory = st.button("↑📁")
    enter_directory = st.button("→")

if up_directory:
    state.folder = os.path.dirname(state.folder)
    st.experimental_rerun() # Restart the whole script to update everything

with left_col:
    save_location = folder_selector("Choose a folder",
                                    state.folder)
if enter_directory:
    if save_location:
        state.folder = os.path.join(state.folder, save_location.rstrip('/\\'))
        st.experimental_rerun() # Restart the whole script to update everything

# A little information about the number of simulation files that will be created.
n_total_simulations = np.prod(n_simulations)
n_simulation_strings = [str(n_sim) for n_sim in n_simulations]
simulation_str = "×".join(n_simulation_strings)

st.write(f"{simulation_str} = {n_total_simulations} total simulation files will be written.")

write_simulations = st.button("Write simulations")


# In theory, all the information you need to run a grid of simulations is in the parameters...
# Something like,
# {'n': [100, 200, 300, 400], }




if write_simulations:
    if not base_name:
        st.write("The Output base filename cannot be empty - fix above.")
    else:
        # Need a function to write the simulations...
        # Simulation string...
        # Need something like
        # new_params = copy.copy(default_params)
        # for i in range(n_scans):
        simulation_inputs = [merge_dicts(dicts) for dicts in itertools.product(*all_params_values)]
        st.write(simulation_inputs)
        # For loop here!
        for i, sim in enumerate(simulation_inputs):
            new_params = copy.copy(updated_default_params)
            new_params.update(sim) # Assign new parameter values...
            # Use the dictionary to write the simulation file...

            # TO DO: Write a function here that takes the parameter dictionary as an input
            #  and returns the properly formatted input file string as an output

            # Generate a filename
            output_filename = base_name+f"-{i:04d}"+'.in' # :04d formats the number as an integer with leading zeros so it is 4 digits long 
            
            # These two lines can be removed once we actually write the files  
            st.write(output_filename)        
            st.write(new_params) 
            st.text(write_inputfile(new_params))
            # Then once you have each string, open the file (by combining the folder and the filename, see os.path.join)
            # 
            with open (output_filename, 'w') as writer:
                writer.write(write_inputfile(new_params))