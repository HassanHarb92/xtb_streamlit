import streamlit as st
import pandas as pd
from fetch_smiles import fetch_smiles
from smiles_to_xyz import smiles_to_xyz
from generate_xtb_command_and_molecule import generate_xtb_command_and_molecule
from run_xtb_command import run_xtb_command
from filter_xtb_output import filter_xtb_output
from visualize_molecule import visualize_molecule

st.title("LLM-Powered xTB Energy Calculator")

if 'molecules_info' not in st.session_state:
    st.session_state.molecules_info = []
if 'selected_molecule' not in st.session_state:
    st.session_state.selected_molecule = None

user_input = st.text_area("Enter a command for multiple molecules (e.g., 'Calculate the energy of aspirin and paracetamol'):", height=150)


from filter_xtb_output_gpt import filter_xtb_output_gpt

# Modify this section in main.py:
if st.button("Submit"):
    if user_input:
        molecule_data = generate_xtb_command_and_molecule(user_input)

        if molecule_data:
            molecules_info = []
            for molecule_name, xtb_command in molecule_data:
                smiles = fetch_smiles(molecule_name)

                if smiles:
                    molecule_name = molecule_name.replace(" ", "_")
                    
                    if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
                        xtb_output_file = 'xtb_output.out'
                        xtb_command_with_output = f"{xtb_command} > {xtb_output_file}"
                        run_xtb_command(xtb_command_with_output)

                        # Use GPT to filter the xTB output
                        json_data = filter_xtb_output_gpt(xtb_output_file, molecule_name)
                        molecules_info.append(json_data)
                    
            st.session_state.molecules_info = molecules_info



if st.session_state.molecules_info:
    molecule_table = pd.DataFrame(st.session_state.molecules_info)
    st.table(molecule_table)

#### HH++

if st.session_state.molecules_info:
    molecule_choices = [info['molecule'] for info in st.session_state.molecules_info]
    selected_molecule = st.selectbox("Select a molecule to visualize:", molecule_choices, key='molecule_dropdown')


#### HH--

    if selected_molecule:
        st.write(f"Visualizing {selected_molecule}")
        visualize_molecule(f'{selected_molecule}.xyz')

