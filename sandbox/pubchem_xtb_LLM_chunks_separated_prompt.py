import openai
import pubchempy as pcp
import streamlit as st
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import os

# Load system prompt from external file (if applicable)
def load_system_prompt(filename):
    with open(filename, 'r') as file:
        return file.read()

# Fetch SMILES string from PubChem
def fetch_smiles(molecule_name):
    try:
        compound = pcp.get_compounds(molecule_name, 'name')
        if compound:
            return compound[0].isomeric_smiles
        else:
            st.error(f"Could not find molecule '{molecule_name}' in PubChem.")
            return None
    except Exception as e:
        st.error(f"Error fetching SMILES from PubChem: {e}")
        return None

# Convert SMILES to XYZ
def smiles_to_xyz(smiles, output_filename):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        
        # Create XYZ content
        conf = mol.GetConformer()
        xyz = f"{mol.GetNumAtoms()}\n\n"
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            xyz += f"{symbol} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
        
        # Save XYZ file
        with open(output_filename, "w") as f:
            f.write(xyz)
        return True
    except Exception as e:
        st.error(f"Error converting SMILES to XYZ: {e}")
        return False

# Function to generate xTB command and extract molecule name based on user input
def generate_xtb_command_and_molecule(prompt):
    # Load the system prompt from an external file
    system_prompt = load_system_prompt('xtb_system_prompt.txt')

    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Generate an xtb command based on the following prompt: {prompt}"}
        ],
        max_tokens=200
    )

    response_text = response['choices'][0]['message']['content'].strip()

    # Ensure the response contains the correct format and cleanly extract the molecule name and command
    try:
        lines = response_text.split('\n')
        molecule_name = lines[0].replace('Molecule: ', '').strip()
        xtb_command = lines[1].replace('xTB command: ', '').strip()
        return molecule_name, xtb_command
    except Exception as e:
        st.error(f"Error parsing response from LLM: {e}")
        return None, None

# Streamlit interface
st.title("LLM-Powered xTB Energy Calculator")

# Input from the user
user_input = st.text_area("Enter a command (e.g., 'Calculate the energy of aspirin' or 'Suggest a molecule similar to aspirin'):", height=150)

# Add a "Submit" button
if st.button("Submit"):
    if user_input:
        # Generate both the molecule name and the xTB command
        molecule_name, xtb_command = generate_xtb_command_and_molecule(user_input)

        if molecule_name and xtb_command:
            st.write(f"**Molecule:** {molecule_name}")
            st.write(f"**xTB Command:** `{xtb_command}`")

            # Fetch SMILES from PubChem
            smiles = fetch_smiles(molecule_name)

            if smiles:
                st.write(f"**SMILES for {molecule_name}:** {smiles}")
                
                # Convert SMILES to XYZ
                if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
                    # Run the xTB command
                    xtb_output_file = 'xtb_output.out'
                    xtb_command_with_output = f"{xtb_command} > {xtb_output_file}"
                    run_xtb_command(xtb_command_with_output)

                    # Extract and display the total energy
                    energy = extract_energy(xtb_output_file)
                    if energy:
                        st.write(f"**Total Energy:** {energy} Eh")

                    # Visualize the molecule
                    if os.path.exists(f'{molecule_name}.xyz'):
                        st.write("**Molecule Structure**:")
                        visualize_molecule(f'{molecule_name}.xyz')
                else:
                    st.error("Failed to generate XYZ file.")
            else:
                st.error(f"Could not fetch SMILES for {molecule_name}.")
    else:
        st.error("Please enter a valid command.")

