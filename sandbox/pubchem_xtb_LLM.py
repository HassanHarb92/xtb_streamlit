import openai
import os
import pubchempy as pcp
import subprocess
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit as st

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")

# Function to fetch SMILES string from PubChem
def fetch_smiles(molecule_name):
    try:
        compound = pcp.get_compounds(molecule_name, 'name')[0]
        return compound.isomeric_smiles
    except Exception as e:
        st.error(f"Error fetching SMILES from PubChem: {e}")
        return None

# Function to convert SMILES to XYZ
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
    system_prompt = """
    You are a computational chemistry assistant that generates xTB commands and identifies molecule names.
    Only use valid xTB flags such as:
    - '--chrg' for charge (followed by an integer)
    - '--opt' for geometry optimization
    - '--alpb' for implicit solvation (e.g., '--alpb water')
    - '--gfn1', '--gfn2', '--gfnff' for the method
    Do not use '--energy' or '--solvent' as they are not valid options.
    Extract the molecule name and generate the corresponding xTB command.
    Return the molecule name and xTB command in this format:
    Molecule: [molecule name]
    xTB command: [xtb command]
    """

    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Generate an xtb command based on the following prompt: {prompt}"}
        ],
        max_tokens=150
    )

    response_text = response['choices'][0]['message']['content'].strip()

    # Parsing the response to extract the molecule name and the xtb command
    lines = response_text.split('\n')
    molecule_name = lines[0].replace('Molecule: ', '').strip()
    xtb_command = lines[1].replace('xTB command: ', '').strip()

    return molecule_name, xtb_command

# Function to run the xTB command
def run_xtb_command(xtb_command):
    try:
        result = subprocess.run(xtb_command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            st.error(f"xTB command failed: {result.stderr}")
            return None
        return result.stdout
    except Exception as e:
        st.error(f"Error running xTB: {e}")
        return None

# Function to visualize XYZ file using py3Dmol
def visualize_molecule(xyz_file):
    try:
        with open(xyz_file, 'r') as f:
            xyz = f.read()
        
        viewer = py3Dmol.view(width=800, height=800)
        viewer.addModel(xyz, 'xyz')
        viewer.setStyle({'stick': {}, 'sphere': {'radius': 0.5}})
        viewer.zoomTo()
        viewer.show()
        st.components.v1.html(viewer._make_html(), width=800, height=800, scrolling=False)
    except Exception as e:
        st.error(f"Error visualizing molecule: {e}")

# Streamlit interface
st.title("LLM-Powered xTB Energy Calculator")

# Input from the user
user_input = st.text_input("Enter a command (e.g., 'Calculate the energy of aspirin'):")

if user_input:
    # Generate both the molecule name and the xTB command
    molecule_name, xtb_command = generate_xtb_command_and_molecule(user_input)

    st.write(f"**Molecule:** {molecule_name}")
    st.write(f"**xTB Command:** `{xtb_command}`")
    
    # Fetch SMILES from PubChem
    smiles = fetch_smiles(molecule_name)
    
    if smiles:
        st.write(f"**SMILES for {molecule_name}:** {smiles}")
        
        # Convert SMILES to XYZ
        if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
            # Run the xTB command
            xtb_output = run_xtb_command(xtb_command)
            
            if xtb_output:
                st.write("**xTB Output:**")
                st.text(xtb_output)
            
            # Visualize the molecule
            if os.path.exists(f'{molecule_name}.xyz'):
                st.write("**Molecule Structure**:")
                visualize_molecule(f'{molecule_name}.xyz')
        else:
            st.error("Failed to generate XYZ file.")
    else:
        st.error(f"Could not fetch SMILES for {molecule_name}.")

