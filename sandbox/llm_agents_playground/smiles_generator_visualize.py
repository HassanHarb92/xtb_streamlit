import streamlit as st
import openai
import os
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

# Initialize the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Streamlit App Title
st.title("AI-Powered SMILES Generator")

# Sidebar Configuration
st.sidebar.markdown("### Settings")
selected_temperature = st.sidebar.slider("Generation Creativity (Temperature)", min_value=0.0, max_value=1.0, value=0.7)

# Input Sections for SMILES and Generation Count
st.text("Provide a list of SMILES strings and the desired number of new SMILES:")
smiles_input = st.text_area("Input SMILES (separated by commas):")
generation_count = st.number_input("Number of New SMILES to Generate:", min_value=1, max_value=100, value=10)

# Function to generate similar SMILES
def generate_similar_smiles(smiles_list, num_generated):
    try:
        # Format input for the model
        formatted_smiles = ", ".join(smiles_list)
        
        # Prompt to guide the LLM in generating similar SMILES
        smiles_prompt = f"""
        You are an expert in cheminformatics and molecular design. The user has provided the following SMILES strings:
        {formatted_smiles}

        Your task is to generate {num_generated} new SMILES strings that are structurally similar to these input molecules.
        Ensure that:
        - The generated SMILES are chemically valid.
        - They maintain core features of the provided molecules while introducing slight variations.
        - They are diverse yet recognizable in relation to the input SMILES.
        
        Provide the generated SMILES as a comma-separated list.
        """
        
        # Query the OpenAI model
        response = openai.ChatCompletion.create(
            model="gpt-4o",
            messages=[{"role": "system", "content": smiles_prompt}],
            max_tokens=1000,
            temperature=selected_temperature,
        )

        # Extract and return the response
        return response["choices"][0]["message"]["content"].split(",")
    except Exception as e:
        return [f"An error occurred: {e}"]

# Function to visualize SMILES
def visualize_smiles(smiles_list):
    try:
        smiles_list = smiles_list[:10]  # Limit visualization to 10 molecules
        mols = [Chem.MolFromSmiles(s.strip()) for s in smiles_list if Chem.MolFromSmiles(s.strip())]
        img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200))
        return img
    except Exception as e:
        return f"An error occurred while generating the visualization: {e}"

# Handle SMILES Generation
if st.button("Generate Similar SMILES"):
    if smiles_input:
        smiles_list = [s.strip() for s in smiles_input.split(",")]
        generated_smiles = generate_similar_smiles(smiles_list, generation_count)
        
        st.markdown("### Generated SMILES")
        for smi in generated_smiles:
            st.markdown(f"- `{smi.strip()}`")
        
        # Generate visualization
        st.markdown("### SMILES Visualization")
        img = visualize_smiles(generated_smiles)
        if isinstance(img, Image.Image):
            st.image(img, caption="Generated SMILES Structures", use_column_width=True)
        else:
            st.error(img)
    else:
        st.error("Please provide at least one SMILES string.")

# Export Generated SMILES
if st.button("Export SMILES"):
    if smiles_input:
        try:
            with open("generated_smiles.txt", "w") as f:
                f.write("\n".join(generated_smiles))
            st.success("SMILES exported successfully!")
        except Exception as e:
            st.error(f"Failed to export SMILES: {e}")
    else:
        st.error("No SMILES available to export. Please generate some first.")

