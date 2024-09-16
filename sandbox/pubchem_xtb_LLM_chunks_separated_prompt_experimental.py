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

# Function to read system prompt from a file
def read_system_prompt(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except Exception as e:
        st.error(f"Error reading system prompt from file: {e}")
        return None

# Function to generate xTB command and extract molecule name based on user input
def generate_xtb_command_and_molecule(prompt):
    system_prompt = read_system_prompt("xtb_system_prompt_experimental.txt")
    if system_prompt is None:
        return None, None

    response = openai.ChatCompletion.create(
        model="gpt-4o",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"""
You are a computational chemistry expert with deep knowledge of molecular structures, chemical reactions, and xTB calculations. 

Please generate an xTB command or suggest molecules based on the following user request:

Task: {prompt}

Interpret the user’s request carefully and provide:
1. A **specific molecule name** if the user is requesting a molecule or a similar structure.
2. A valid **xTB command** if the user is asking for a calculation.
3. Any **modifications** (such as changes to molecular structures or charge) that may be necessary based on the user's request.

Always return the **molecule name** and the **xTB command** in this format:
Molecule: [molecule name]
xTB command: [xtb command]
"""
            }

            #{"role": "user", "content": f"Generate an xtb command based on the following prompt: {prompt}"}
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

# Function to filter xTB output for relevant lines
def filter_xtb_output(output_file):
    relevant_lines = []
    keywords = ["TOTAL ENERGY", "ENTHALPY", "FREE ENERGY", "HOMO-LUMO GAP"]

    try:
        with open(output_file, 'r') as file:
            for line in file:
                if any(keyword in line for keyword in keywords):
                    relevant_lines.append(line)

        return ''.join(relevant_lines)  # Return filtered content as a single string
    except Exception as e:
        st.error(f"Error filtering xTB output: {e}")
        return None

# Function to send xTB output to LLM in chunks
def parse_xtb_output_in_chunks(output_file):
    filtered_content = filter_xtb_output(output_file)
    
    if not filtered_content:
        return None

    chunk_size = 500  # Define the size of each chunk in characters
    parsed_output = ""

    system_prompt = """
    You are a computational chemistry assistant. 
    Your task is to extract all relevant energy values from an xTB output file.
    These may include:
    - Total energy
    - Enthalpy
    - Free energy
    - HOMO-LUMO gap
    Provide the extracted values clearly, with labels for each energy.
    """

    try:
        for i in range(0, len(filtered_content), chunk_size):
            chunk = filtered_content[i:i + chunk_size]
            response = openai.ChatCompletion.create(
                model="gpt-4o",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": f"Here is a portion of the xTB output:\n{chunk}"}
                ],
                max_tokens=200
            )

            parsed_output += response['choices'][0]['message']['content'].strip() + "\n"

        return parsed_output

    except Exception as e:
        st.error(f"Error parsing xTB output with LLM in chunks: {e}")
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
user_input = st.text_area("Enter a command (e.g., 'Calculate the energy of aspirin' or 'Suggest a molecule similar to aspirin'):", height=100)

# Add a "Submit" button
if st.button("Submit"):
    if user_input:
        # Generate both the molecule name and the xTB command
        molecule_name, xtb_command = generate_xtb_command_and_molecule(user_input)

        if molecule_name and xtb_command:
            st.write(f"**Molecule:** {molecule_name}")
            st.write(f"**xTB Command:** `{xtb_command}`")
#            molecule_name = molecule_name.replace(" ", "_")

            # Fetch SMILES from PubChem
            smiles = fetch_smiles(molecule_name)
            molecule_name = molecule_name.replace(" ", "_")
            if smiles:
                st.write(f"**SMILES for {molecule_name}:** {smiles}")
                
                # Convert SMILES to XYZ
                if smiles_to_xyz(smiles, f'{molecule_name}.xyz'):
                    # Run the xTB command
                    xtb_output_file = 'xtb_output.out'
                    xtb_command_with_output = f"{xtb_command} > {xtb_output_file}"
                    run_xtb_command(xtb_command_with_output)

                    # Send xTB output to LLM in chunks
                    parsed_output = parse_xtb_output_in_chunks(xtb_output_file)
                    if parsed_output:
                        st.write("**Extracted Energy Values from xTB Output:**")
                        st.write(parsed_output)
                    
                    # Visualize the molecule
                    copy_command = f'mv xtbopt.xyz {molecule_name}.xyz'
                    subprocess.run(copy_command, shell=True, capture_output=True, text=True)     
                    print (f"Molecule name {molecule_name}.xyz")
                    if os.path.exists(f'{molecule_name}.xyz'):
                        st.write("**Molecule Structure**:")
                        visualize_molecule(f'{molecule_name}.xyz')
                else:
                    st.error("Failed to generate XYZ file.")
            else:
                st.error(f"Could not fetch SMILES for {molecule_name}.")
    else:
        st.error("Please enter a valid command.")

