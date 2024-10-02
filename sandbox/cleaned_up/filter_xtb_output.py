import json
import streamlit as st

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

        with open(f'{molecule_name}.json', 'w') as json_file:
            json.dump(json_data, json_file)
        
        return json_data
    except Exception as e:
        st.error(f"Error filtering xTB output: {e}")
        return None



