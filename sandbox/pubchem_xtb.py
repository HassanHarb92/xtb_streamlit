import streamlit as st
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import pubchempy as pcp

# Function to fetch SMILES string from PubChem
def fetch_smiles_from_pubchem(molecule_name):
    try:
        compound = pcp.get_compounds(molecule_name, 'name')[0]  # Fetch the first compound match
        return compound.isomeric_smiles
    except Exception as e:
        st.error(f"Error fetching SMILES from PubChem: {e}")
        return None

# Function to calculate the number of electrons and determine multiplicity
def calculate_multiplicity(smiles, charge):
    mol = Chem.MolFromSmiles(smiles)
    num_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - charge
    multiplicity = 1 if num_electrons % 2 == 0 else 2  # Even electrons -> singlet (1), odd -> doublet (2)
    return num_electrons, multiplicity

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
        st.error(f"Error in converting SMILES to XYZ: {e}")
        return False

# Function to run xtb calculation
def run_xtb(input_filename, method, charge, multiplicity, solvent, solvent_name):
    try:
        uhf_flag = f"--uhf {multiplicity}" if multiplicity != 1 else ""
        solvent_flag = f"--{solvent} {solvent_name}" if solvent != 'none' else ""
        
        # Optimization command
        xtb_command = f"xtb {input_filename} --{method} --chrg {charge} {uhf_flag} --opt {solvent_flag} > xtb_output.out"
        result = subprocess.run(xtb_command, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            st.error(f"xtb optimization command failed: {result.stderr}")
            return False
        
        return True
    except Exception as e:
        st.error(f"Error running xtb calculation: {e}")
        return False

# Function to extract total energy from xtb output
def extract_results(output_file):
    try:
        energies = {}
        with open(output_file, 'r') as file:
            for line in file:
                if "TOTAL ENERGY" in line:
                    energies['Energy'] = line.split()[-3]
        return energies
    except FileNotFoundError:
        st.error(f"Output file {output_file} not found.")
    except Exception as e:
        st.error(f"Error extracting results: {e}")
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

# Function to create a downloadable results file
def create_results_file(energies, xyz_file):
    result_text = "Calculation Results:\n"
    for key, value in energies.items():
        result_text += f"{key}: {value} Eh\n"
    
    result_text += "\nOptimized XYZ Coordinates:\n"
    try:
        with open(xyz_file, 'r') as f:
            result_text += f.read()
    except FileNotFoundError:
        result_text += "XYZ file not found."

    return result_text

# Streamlit interface
st.title('xTB2Go - Energy Calculation App')

# Input options for the user
user_input = st.text_input("Enter calculation request (e.g., 'Calculate the energy of aspirin'):")

# Add a "Go" button right below the input box
if st.button("Go"):
    if user_input:
        # Extract the molecule name from user input (for simplicity, we'll just extract after "of")
        if 'of' in user_input:
            molecule_name = user_input.split('of')[-1].strip()
            
            # Fetch SMILES from PubChem
            smiles = fetch_smiles_from_pubchem(molecule_name)
            
            if smiles:
                st.write(f"**SMILES of {molecule_name}:** {smiles}")
                
                # Charge input with smaller box using columns
                charge_col, _ = st.columns([1, 4])  # Adjust column width to make input smaller
                with charge_col:
                    st.write("Charge")
                    charge = st.number_input("Charge:", value=0, step=1, format="%d", min_value=-10, max_value=10, label_visibility="collapsed")
                
                # Select method and solvent
                method = st.radio("Select method for xtb:", ['gfn1', 'gfn2', 'gfnff'], index=2)
                solvent = st.selectbox("Select solvation model:", ['none', 'alpb', 'gbsa'], index=0)
                
                solvent_name = ''
                if solvent != 'none':
                    solvent_options = {
                        'alpb': ['water', 'acetone', 'acetonitrile'],
                        'gbsa': ['water', 'acetone', 'acetonitrile']
                    }
                    solvent_name = st.selectbox("Select solvent:", solvent_options[solvent])
                
                # Immediately start the calculation without an additional button press
                num_electrons, multiplicity = calculate_multiplicity(smiles, charge)
                st.write(f"Number of Electrons: {num_electrons}, Multiplicity: {multiplicity}")

                # Convert SMILES to XYZ
                if smiles_to_xyz(smiles, 'input_molecule.xyz'):
                    # Run xtb calculation based on user inputs
                    if run_xtb('input_molecule.xyz', method, charge, multiplicity, solvent, solvent_name):
                        # Extract and display results
                        energies = extract_results('xtb_output.out')
                        
                        if energies:
                            st.write(f"**Total Energy**: {energies.get('Energy', 'N/A')} Eh")
                        
                        # Visualize optimized molecule
                        if os.path.exists('xtbopt.xyz'):
                            st.write("**Optimized Molecule Structure**:")
                            visualize_molecule('xtbopt.xyz')
                        
                        # Create downloadable results file
                        results_text = create_results_file(energies, 'xtbopt.xyz')
                        st.download_button("Download Results", data=results_text, file_name="results.txt", mime="text/plain")
                    else:
                        st.error("Failed to run calculation.")
                else:
                    st.error("Failed to convert SMILES to XYZ.")
        else:
            st.error("Please specify a valid molecule in the format 'Calculate the energy of [molecule]'")

