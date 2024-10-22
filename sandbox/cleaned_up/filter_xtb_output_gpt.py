import openai
import streamlit as st
import json

def filter_xtb_output_gpt(output_file, molecule_name):
    # Read the xTB output file
    try:
        with open(output_file, 'r') as file:
            xtb_output = file.read()
    except Exception as e:
        st.error(f"Error reading xTB output file: {e}")
        return None
    
    # Define the system prompt for GPT to parse xTB output
    system_prompt = (
        "You are an expert in computational chemistry. Extract key information "
        "from xTB output such as Total Energy, Enthalpy, Free Energy, and HOMO-LUMO Gap."
    )
    
    # Send xTB output to GPT for parsing
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": f"Here is the xTB output for {molecule_name}:\n{xtb_output}"}
            ],
            max_tokens=200
        )

        if 'choices' not in response or len(response['choices']) == 0:
            st.error("GPT response is empty or invalid.")
            return None

        parsed_output = response['choices'][0]['message']['content'].strip()
        if not parsed_output:
            st.error("Parsed output from GPT is empty.")
            return None
        
        # Save the parsed output to a JSON file
        json_data = {
            "molecule": molecule_name,
            "parsed_xtb_output": parsed_output
        }
        
        try:
            with open(f'{molecule_name}_parsed.json', 'w') as json_file:
                json.dump(json_data, json_file)
        except Exception as e:
            st.error(f"Error saving parsed output to JSON file: {e}")
            return None
        
        return json_data

    except Exception as e:
        st.error(f"Error processing xTB output with GPT: {e}")
        return None

