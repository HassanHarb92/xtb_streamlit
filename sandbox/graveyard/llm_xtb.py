import openai
import os

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")

# Function to generate xTB command and extract molecule name based on user input
def generate_xtb_command_and_molecule(prompt):
    # Provide a detailed system prompt for the LLM
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
        model="gpt-4",  # You can also use 'gpt-3.5-turbo'
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

# Example usage in Streamlit
#user_input = "Calculate the energy of aspirin in water with a charge of -1"
user_input = "Optimize the geometry of aspirin in water"
# Generate both the molecule name and the xTB command
molecule_name, xtb_command = generate_xtb_command_and_molecule(user_input)

print(f"Molecule: {molecule_name}")
print(f"xTB Command: {xtb_command}")



