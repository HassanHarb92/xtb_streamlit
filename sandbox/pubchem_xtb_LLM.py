import streamlit as st
import openai  # Assuming OpenAI API is used for LLM
import pubchempy as pcp
import subprocess
import os
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem


# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")




# Function to fetch SMILES from PubChem
def fetch_smiles(molecule_name):
    try:
        compound = pcp.get_compounds(molecule_name, 'name')[0]
        return compound.isomeric_smiles
    except Exception as e:
        st.error(f"Error fetching SMILES: {e}")
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

# Function to generate xTB command based on LLM output

def generate_xtb_command(prompt):
    response = openai.ChatCompletion.create(
        model="gpt-4",  # You can also use 'gpt-3.5-turbo'
        messages=[
            {"role": "system", "content": "You are an assistant that generates xTB commands."},
            {"role": "user", "content": f"Generate an xtb command based on the following prompt: {prompt}"}
        ],
        max_tokens=100
    )
    xtb_command = response['choices'][0]['message']['content'].strip()
    return xtb_command


# Function to run the generated xTB command
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

# Streamlit app interface
st.title("LLM-Powered xTB Energy Calculator")

# Input from the user
user_input = st.text_input("Enter a command (e.g., 'Calculate the energy of aspirin in water with a charge of -1'):")

if user_input:
    # Extract molecule name from user input
    if 'of' in user_input:
        molecule_name = user_input.split('of')[-1].strip().split()[0]
        
        # Fetch SMILES from PubChem
        smiles = fetch_smiles(molecule_name)
        
        if smiles:
            st.write(f"**SMILES for {molecule_name}:** {smiles}")
            
            # Convert SMILES to XYZ
            if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
                # Generate xTB command based on user input and molecule name
                xtb_command = generate_xtb_command(user_input, molecule_name, f'{molecule_name}.xyz')
                
                # Display the generated command
                st.write(f"**xTB Command:** `{xtb_command}`")
                
                # Run the generated xTB command
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
        st.error("Please specify a valid molecule in the format 'Calculate the energy of [molecule]'")










