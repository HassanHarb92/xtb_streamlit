import openai
from read_system_prompt import read_system_prompt

def generate_xtb_command_and_molecule(prompt):
    system_prompt = read_system_prompt("xtb_system_prompt.txt")
    if system_prompt is None:
        return None, None

    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Generate xtb commands for multiple molecules based on the following prompt: {prompt}"}
        ],
        max_tokens=150
    )

    response_text = response['choices'][0]['message']['content'].strip()

    molecule_data = []
    for line in response_text.split('\n'):
        if "Molecule:" in line:
            molecule_name = line.replace('Molecule:', '').strip()
        if "xTB command:" in line:
            xtb_command = line.replace('xTB command:', '').strip()
            molecule_data.append((molecule_name, xtb_command))
    
    return molecule_data

