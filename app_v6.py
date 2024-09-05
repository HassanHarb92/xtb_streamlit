import streamlit as st
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

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
def run_xtb(input_filename, method, charge, multiplicity, optimize, thermochemistry, solvent, solvent_name):
    try:
        uhf_flag = f"--uhf {multiplicity}" if multiplicity != 1 else ""
        solvent_flag = f"--{solvent} {solvent_name}" if solvent != 'none' else ""
        
        # Optimization command
        if optimize:
            xtb_command = f"xtb {input_filename} --{method} --chrg {charge} {uhf_flag} --opt {solvent_flag} > xtb_output.out"
            result = subprocess.run(xtb_command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                st.error(f"xtb optimization command failed: {result.stderr}")
                return False
            
            # Thermochemistry calculation
            if thermochemistry:
                optimized_filename = 'xtbopt.xyz'
                xtb_command_hess = f"xtb {optimized_filename} --{method} --chrg {charge} {uhf_flag} --hess {solvent_flag} > xtb_hess_output.out"
                result_hess = subprocess.run(xtb_command_hess, shell=True, capture_output=True, text=True)
                
                if result_hess.returncode != 0:
                    st.error(f"xtb thermochemistry command failed: {result_hess.stderr}")
                    return False
        else:
            # Single-point energy calculation
            xtb_command = f"xtb {input_filename} --{method} --chrg {charge} {uhf_flag} {solvent_flag} > xtb_output.out"
            result = subprocess.run(xtb_command, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                st.error(f"xtb command failed: {result.stderr}")
                return False
        
        return True
    except Exception as e:
        st.error(f"Error running xtb calculation: {e}")
        return False

# Function to extract total energy and thermochemistry results from xtb output
def extract_results(output_file):
    try:
        energies = {}
        with open(output_file, 'r') as file:
            for line in file:
                if "TOTAL ENERGY" in line:
                    energies['Energy'] = line.split()[-3]
                elif "TOTAL ENTHALPY" in line:
                    energies['Enthalpy'] = line.split()[-3]
                elif "TOTAL FREE ENERGY" in line:
                    energies['Free Energy'] = line.split()[-3]
        return energies
    except FileNotFoundError:
        st.error(f"Output file {output_file} not found.")
    except Exception as e:
        st.error(f"Error extracting results: {e}")
    return None

# Function to extract frequencies from g98.out file
def extract_frequencies(frequency_file):
    frequencies = []
    try:
        with open(frequency_file, 'r') as file:
            for line in file:
                if "Frequencies --" in line:
                    freqs = line.split()[2:]  # Get frequencies from the line
                    frequencies.extend(freqs)  # Add frequencies to the list
        return frequencies
    except FileNotFoundError:
        st.error(f"Frequency file {frequency_file} not found.")
    except Exception as e:
        st.error(f"Error extracting frequencies: {e}")
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
def create_results_file(energies, xyz_file, frequencies=None):
    result_text = "Calculation Results:\n"
    for key, value in energies.items():
        result_text += f"{key}: {value} Eh\n"
    
    result_text += "\nOptimized XYZ Coordinates:\n"
    try:
        with open(xyz_file, 'r') as f:
            result_text += f.read()
    except FileNotFoundError:
        result_text += "XYZ file not found."

    if frequencies:
        result_text += "\nFrequencies (cm^-1):\n"
        result_text += "\n".join(frequencies)

    return result_text

# Streamlit interface
st.title('xTB2Go')
st.markdown('## Molecule Optimization and Thermochemistry App')

# Input options for the user
smiles = st.text_input("Enter SMILES string of a molecule:")

# Charge input with smaller box using columns
charge_col, _ = st.columns([1, 4])  # Adjust column width to make input smaller
with charge_col:
    st.write("Charge")
    charge = st.number_input("Charge:", value=0, step=1, format="%d", min_value=-10, max_value=10, label_visibility="collapsed")

optimize = st.checkbox("Optimize geometry", value=True)
thermochemistry = False
if optimize:
    thermochemistry = st.checkbox("Perform Thermochemistry", value=False)

# Columns for method and solvent options
col1, col2, col3 = st.columns(3)

with col1:
    method = st.radio("Select method for xtb:", ['gfn1', 'gfn2', 'gfnff'], index=2)

with col2:
    # Dropdown for solvation model
    solvent = st.selectbox("Select solvation model:", ['none', 'alpb', 'gbsa'], index=0)

with col3:
    solvent_name = ''
    if solvent != 'none':
        solvent_options = {
            'alpb': ['water', 'acetone', 'acetonitrile', 'aniline', 'benzaldehyde', 'benzene', 'CH₂Cl₂', 'CHCl₃', 'CS₂', 'dioxane', 'DMF', 'DMSO', 'ether', 'ethylacetate', 'furane', 'hexadecane', 'hexane', 'methanol', 'nitromethane', 'octanol', 'octanol (wet)', 'phenol', 'toluene', 'THF'],
            'gbsa': ['water', 'acetone', 'acetonitrile', 'aniline', 'benzaldehyde', 'benzene', 'CH₂Cl₂', 'CHCl₃', 'CS₂', 'dioxane', 'DMF', 'DMSO', 'ether', 'ethylacetate', 'furane', 'hexadecane', 'hexane', 'methanol', 'nitromethane', 'octanol', 'octanol (wet)', 'phenol', 'toluene', 'THF']
        }
        solvent_name = st.selectbox("Select solvent:", solvent_options[solvent])

# "Run" button to execute the calculation
if st.button("Run"):
    if smiles:
        # Calculate the number of electrons and determine multiplicity
        num_electrons, multiplicity = calculate_multiplicity(smiles, charge)
        st.write(f"Number of Electrons: {num_electrons}, Multiplicity: {multiplicity}")

        # Convert SMILES to XYZ
        if smiles_to_xyz(smiles, 'input_molecule.xyz'):
            # Run xtb calculation based on user inputs
            if run_xtb('input_molecule.xyz', method, charge, multiplicity, optimize, thermochemistry, solvent, solvent_name):
                # Extract and display results
                if thermochemistry:
                    energies = extract_results('xtb_hess_output.out')
                    frequencies = extract_frequencies('g98.out')
                else:
                    energies = extract_results('xtb_output.out')
                    frequencies = None

                if energies:
                    st.write(f"**Total Energy**: {energies.get('Energy', 'N/A')} Eh")
                    if thermochemistry:
                        st.write(f"**Total Enthalpy**: {energies.get('Enthalpy', 'N/A')} Eh")
                        st.write(f"**Total Free Energy**: {energies.get('Free Energy', 'N/A')} Eh")

                # Display lowest five frequencies if available
                if frequencies:
                    frequencies.sort()  # Sort frequencies to find the lowest ones
                    lowest_frequencies = frequencies[:5]  # Get the lowest five frequencies
                    st.write("**Lowest Frequencies (cm^-1):**")
                    st.write(lowest_frequencies)
                # Visualize optimized molecule if optimization was selected
                if optimize and os.path.exists('xtbopt.xyz'):
                    st.write("**Optimized Molecule Structure**:")
                    visualize_molecule('xtbopt.xyz')
                elif not optimize:
                    st.write("No optimization was performed.")
                
                # Create downloadable results file
                results_text = create_results_file(energies, 'xtbopt.xyz', frequencies)
                st.download_button("Download Results", data=results_text, file_name="results.txt", mime="text/plain")
            else:
                st.error("Failed to run calculation.")
        else:
            st.error("Failed to convert SMILES to XYZ.")
    else:
        st.error("Please enter a valid SMILES string.")

