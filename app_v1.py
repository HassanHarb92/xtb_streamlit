import streamlit as st
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

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
        st.write(f"XYZ file created: {output_filename}")
        return True
    except Exception as e:
        st.error(f"Error in converting SMILES to XYZ: {e}")
        return False

# Function to run xtb optimization
def run_xtb(input_filename):
    try:
        xtb_command = f"xtb {input_filename} --opt > xtb_output.out"
        result = subprocess.run(xtb_command, shell=True, capture_output=True, text=True)
        st.write(f"Running xtb command: {xtb_command}")
        if result.returncode != 0:
            st.error(f"xtb command failed: {result.stderr}")
            return False
        else:
            st.write("xtb optimization completed successfully.")
            return True
    except Exception as e:
        st.error(f"Error running xtb: {e}")
        return False

# Updated function to extract total energy from xtb output
def extract_total_energy(output_file):
    try:
        with open(output_file, 'r') as file:
            for line in file:
                if "TOTAL ENERGY" in line:
                    # Split line by whitespace and get the energy value before 'Eh'
                    energy = line.split()[-3]  # Gets the energy value
                    st.write(f"Extracted total energy: {energy} Eh")
                    return energy
        st.error("TOTAL ENERGY not found in output file.")
    except FileNotFoundError:
        st.error(f"Output file {output_file} not found.")
    except Exception as e:
        st.error(f"Error extracting total energy: {e}")
    return None

# Function to visualize XYZ file using py3Dmol
def visualize_molecule(xyz_file):
    try:
        with open(xyz_file, 'r') as f:
            xyz = f.read()
        
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(xyz, 'xyz')
        viewer.setStyle({'stick': {}})
        viewer.zoomTo()
        viewer.show()
    except FileNotFoundError:
        st.error(f"Optimized XYZ file {xyz_file} not found.")
    except Exception as e:
        st.error(f"Error visualizing molecule: {e}")

# Streamlit interface
st.title('Molecule Optimization App')

# Input SMILES string
smiles = st.text_input("Enter SMILES string of a molecule:")

if smiles:
    # Convert SMILES to XYZ
    if smiles_to_xyz(smiles, 'input_molecule.xyz'):
        # Run xtb optimization
        if run_xtb('input_molecule.xyz'):
            # Extract total energy
            energy = extract_total_energy('xtb_output.out')
            
            if energy:
                # Display total energy
                st.write(f"**Total Energy**: {energy} Eh")
                
                # Visualize optimized molecule
                if os.path.exists('xtbopt.xyz'):
                    st.write("**Optimized Molecule Structure**:")
                    visualize_molecule('xtbopt.xyz')
                else:
                    st.error("Optimized XYZ file (xtbopt.xyz) not found after optimization.")
            else:
                st.error("Failed to extract total energy from xtb output.")
        else:
            st.error("Failed to run xtb optimization.")
    else:
        st.error("Failed to convert SMILES to XYZ.")

