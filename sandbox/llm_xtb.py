import openai
import os

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")

# Function to generate xTB command based on user input
def generate_xtb_command(prompt):
    # Provide a detailed system prompt for the LLM
    system_prompt = """
    You are a computational chemistry assistant that is expert in xTB and generates xTB commands.
    Only use valid xTB flags such as:
    - '--chrg' for charge (followed by an integer)
    - '--opt' for geometry optimization
    - '--alpb' for implicit solvation (e.g., '--alpb water')
    - '--gfn1', '--gfn2', '--gfnff' for the method
    - default method is --gfn2 so use that if no method is defined
    Do not use '--energy' or '--solvent' as they are not valid options.
    Generate xTB commands based on the user's input.
    """

    response = openai.ChatCompletion.create(
        model="gpt-4",  # You can also use 'gpt-3.5-turbo'
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Generate an xtb command based on the following prompt: {prompt}"}
        ],
        max_tokens=100
    )
    xtb_command = response['choices'][0]['message']['content'].strip()
    return xtb_command

# Example usage in Streamlit
user_input = "Calculate the energy of aspirin in water with a charge of -1"
xtb_command = generate_xtb_command(user_input)

print(xtb_command)  # This will print the generated xTB command

