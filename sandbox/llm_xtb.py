import openai
import os

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

if openai.api_key is None:
    raise ValueError("OpenAI API key not found. Please set the 'OPENAI_API_KEY' environment variable.")

# Function to generate xTB command based on user input
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

# Example usage in Streamlit
user_input = "Calculate the energy of aspirin in water with a charge of -1"
xtb_command = generate_xtb_command(user_input)

print(xtb_command)  # This will print the generated xTB command

