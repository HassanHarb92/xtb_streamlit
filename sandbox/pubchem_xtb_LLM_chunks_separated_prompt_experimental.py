## xtb command not working
## visualize molecule not working


import openai
import os
import pubchempy as pcp
import subprocess
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import streamlit as st
import pandas as pd
import json

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")

# Function to fetch SMILES string and properties from PubChem
def fetch_smiles_and_properties(molecule_name):
    try:
        compound = pcp.get_compounds(molecule_name, 'name')[0]
        smiles = compound.isomeric_smiles

        # Fetch physical properties (molecular weight, boiling point, etc.)
        molecular_weight = compound.molecular_weight if hasattr(compound, 'molecular_weight') else 'N/A'
        boiling_point = compound.boiling_point if hasattr(compound, 'boiling_point') else 'N/A'
        melting_point = compound.melting_point if hasattr(compound, 'melting_point') else 'N/A'
        exact_mass = compound.exact_mass if hasattr(compound, 'exact_mass') else 'N/A'
        
        properties = {
            "Molecular Weight": molecular_weight,
            "Boiling Point": boiling_point,
            "Melting Point": melting_point,
            "Exact Mass": exact_mass,
        }
        return smiles, properties
    except Exception as e:
        st.error(f"Error fetching data from PubChem for {molecule_name}: {e}")
        return None, None

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

# Function to read system prompt from a file
def read_system_prompt(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except Exception as e:
        st.error(f"Error reading system prompt from file: {e}")
        return None

# Function to generate xTB command and extract molecule names based on user input
def generate_xtb_command_and_molecule(prompt):
    system_prompt = read_system_prompt("xtb_system_prompt.txt")
    if system_prompt is None:
        return None

    response = openai.ChatCompletion.create(
        model="gpt-4o",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Generate xtb commands for multiple molecules based on the following prompt: {prompt}"}
        ],
        max_tokens=150
    )

    response_text = response['choices'][0]['message']['content'].strip()

    # Extract molecule names and xTB commands
    molecule_data = []
    for line in response_text.split('\n'):
        if "Molecule:" in line:
            molecule_name = line.replace('Molecule:', '').strip()
        if "xTB command:" in line:
            xtb_command = line.replace('xTB command:', '').strip()
            if molecule_name and xtb_command:
                molecule_data.append((molecule_name, xtb_command))
    
    return molecule_data

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

# Function to filter xTB output for relevant lines and return JSON format
def filter_xtb_output(output_file, molecule_name):
    relevant_lines = []
    keywords = ["TOTAL ENERGY", "ENTHALPY", "FREE ENERGY", "HOMO-LUMO GAP"]
    json_data = {}

    try:
        with open(output_file, 'r') as file:
            for line in file:
                if any(keyword in line for keyword in keywords):
                    relevant_lines.append(line)

        json_data['molecule'] = molecule_name
        json_data['energies'] = ''.join(relevant_lines)

        # Save to JSON
        with open(f'{molecule_name}.json', 'w') as json_file:
            json.dump(json_data, json_file)
        
        return json_data  # Return the JSON data for further use
    except Exception as e:
        st.error(f"Error filtering xTB output: {e}")
        return None

# Function to visualize XYZ file using py3Dmol

# Function to parse xTB output using GPT-4 and extract relevant information
def parse_xtb_output_with_llm(output_file, molecule_name):
    system_prompt = """
    You are a computational chemistry assistant. 
    Your task is to intelligently parse the xTB output file and extract all relevant chemical information.
    This may include:
    - Total energy
    - Enthalpy
    - Free energy
    - HOMO-LUMO gap
    - Dipole moment
    - Any other available properties
    Provide a clear, structured summary of the extracted information.
    Return the information in JSON format where possible.
    """

    try:
        with open(output_file, 'r') as file:
            xtb_output = file.read()

        # Chunk the xTB output if it's too long for the model
        chunk_size = 4000  # Adjust based on the model's context length
        parsed_output = ""

        for i in range(0, len(xtb_output), chunk_size):
            chunk = xtb_output[i:i + chunk_size]
            
            response = openai.ChatCompletion.create(
                model="gpt-4",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": f"Here is a portion of the xTB output for {molecule_name}:\n{chunk}"}
                ],
                max_tokens=300  # Adjust based on expected output
            )

            parsed_output += response['choices'][0]['message']['content'].strip()

        # Process the GPT-4 response to extract relevant data
        json_data = {"molecule": molecule_name}

        # Parse output as JSON if GPT-4 returns it in that format
        try:
            json_data.update(json.loads(parsed_output))
        except json.JSONDecodeError:
            # If not in JSON format, just store the raw parsed text
            json_data['parsed_data'] = parsed_output

        # Save to JSON file
        with open(f'{molecule_name}.json', 'w') as json_file:
            json.dump(json_data, json_file)

        return json_data  # Return the JSON data for further use

    except Exception as e:
        st.error(f"Error parsing xTB output with GPT-4: {e}")
        return None


# Streamlit interface
st.title("LLM-Powered xTB Energy Calculator")

# Initialize session state for storing molecule data and selections
if 'molecules_info' not in st.session_state:
    st.session_state.molecules_info = []
if 'selected_molecule' not in st.session_state:
    st.session_state.selected_molecule = None

# Input from the user
user_input = st.text_area("Enter a command for multiple molecules (e.g., 'Calculate the energy of aspirin and paracetamol'):", height=150)

# Main function where xTB output parsing happens
if st.button("Submit"):
    if user_input:
        # Generate molecule data and xTB commands
        molecule_data = generate_xtb_command_and_molecule(user_input)

        if molecule_data:
            molecules_info = []
            for molecule_name, xtb_command in molecule_data:
                # Ensure molecule_name is properly assigned
                if not molecule_name:
                    st.error("Molecule name could not be identified.")
                    continue

                # Fetch SMILES and properties from PubChem
                smiles, properties = fetch_smiles_and_properties(molecule_name)

                if smiles:
                    molecule_name = molecule_name.replace(" ", "_")                
                    
                    # Convert SMILES to XYZ
                    if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
                        # Run the xTB command
                        xtb_output_file = 'xtb_output.out'
                        xtb_command_with_output = f"{xtb_command} > {xtb_output_file}"
                        run_xtb_command(xtb_command_with_output)

                        # Parse xTB output using GPT-4 and extract relevant data
                        json_data = parse_xtb_output_with_llm(xtb_output_file, molecule_name)

                        # Append SMILES, energies, and PubChem properties to the molecule info
                        molecule_entry = {
                            'Molecule': molecule_name,
                            'SMILES': smiles,
                            'Molecular Weight': properties.get('Molecular Weight', 'N/A'),
                            'Boiling Point': properties.get('Boiling Point', 'N/A'),
                            'Melting Point': properties.get('Melting Point', 'N/A'),
                            'Exact Mass': properties.get('Exact Mass', 'N/A'),
                            'Parsed Data': json_data.get('parsed_data', 'N/A')  # Data from GPT-4
                        }
                        molecules_info.append(molecule_entry)

            # Store molecule information in session state
            st.session_state.molecules_info = molecules_info

# Display the table with molecule information
if st.session_state.molecules_info:
    molecule_table = pd.DataFrame(st.session_state.molecules_info)
    st.table(molecule_table)

# Dropdown to select and visualize the molecule
if st.session_state.molecules_info:
    molecule_choices = [info['Molecule'] for info in st.session_state.molecules_info]
    selected_molecule = st.selectbox("Select a molecule to visualize:", molecule_choices, key='molecule_dropdown')

    if selected_molecule:
        st.write(f"Visualizing {selected_molecule}")
        visualize_molecule(f'{selected_molecule}.xyz')

