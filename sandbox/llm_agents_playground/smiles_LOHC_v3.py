import streamlit as st
import openai
import os
import PyPDF2
import pandas as pd
from rdkit import Chem

# Initialize the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Streamlit App Title
st.title("AI-Powered LOHC SMILES Generator")

# Sidebar Configuration
st.sidebar.markdown("### Settings")
selected_temperature = st.sidebar.slider("Generation Creativity (Temperature)", min_value=0.0, max_value=1.0, value=0.7)

# Input Sections for SMILES and Generation Count
st.text("Provide a list of known LOHC SMILES strings and the desired number of new SMILES:")
smiles_input = st.text_area("Input LOHC SMILES (separated by commas):")

# Path to the LOHC knowledge document
lohc_pdf_path = "LOHC_DD.pdf"

def extract_pdf_text(pdf_path):
    """Extracts and summarizes key content from a given PDF file."""
    try:
        with open(pdf_path, "rb") as file:
            reader = PyPDF2.PdfReader(file)
            text = " ".join([page.extract_text() for page in reader.pages if page.extract_text()])
        return text[:5000]  # Limit input to avoid excessive token usage
    except Exception as e:
        return f"Error reading PDF: {e}"

lohc_knowledge = extract_pdf_text(lohc_pdf_path)

def is_valid_smiles(smiles):
    """Checks if a SMILES string is chemically valid using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def generate_unique_smiles(smiles_list, knowledge_text):
    """Generates unique LOHC SMILES strings from each given SMILES."""
    try:
        all_generated_smiles = set()
        
        for smiles in smiles_list:
            attempts = 0
            valid_smiles = []
            
            while not valid_smiles and attempts < 3:  # Retry up to 3 times if no valid SMILES are generated
                lohc_prompt = f"""
                You are an expert in molecular design, specializing in Liquid Organic Hydrogen Carriers (LOHCs).
                
                The user provided this known LOHC SMILES:
                {smiles}
                
                Based on the knowledge from the LOHC document:
                {knowledge_text}
                
                Your task is to generate 10 novel LOHC SMILES strings that:
                - Are chemically valid.
                - Retain key structural features relevant to hydrogen storage.
                - Are distinct from the provided SMILES.
                - Exhibit potential for practical application.
                
                Provide the new SMILES in a comma-separated list. Ensure that the output contains only valid SMILES strings.
                """
                
                response = openai.ChatCompletion.create(
                    model="gpt-4o",
                    messages=[{"role": "system", "content": lohc_prompt}],
                    max_tokens=1000,
                    temperature=selected_temperature,
                    stop=["\n"]  # Ensure clear separation of outputs
                )
                
                raw_output = response["choices"][0]["message"]["content"]
                smiles_candidates = [s.strip() for s in raw_output.split(",")]
                
                valid_smiles = [s for s in smiles_candidates if is_valid_smiles(s)]
                attempts += 1
            
            if valid_smiles:
                all_generated_smiles.update(valid_smiles)
            else:
                st.error(f"No valid SMILES generated for input: {smiles} after {attempts} attempts.")
        
        return list(all_generated_smiles)  # Ensure uniqueness
    except Exception as e:
        return [f"An error occurred: {e}"]

# Handle SMILES Generation
if st.button("Generate New LOHC SMILES"):
    if smiles_input:
        smiles_list = [s.strip() for s in smiles_input.split(",")]
        generated_smiles = generate_unique_smiles(smiles_list, lohc_knowledge)
        
        if generated_smiles:
            st.markdown("### Generated LOHC SMILES")
            for smi in generated_smiles:
                st.markdown(f"- `{smi}`")
        else:
            st.error("No valid SMILES generated. Try adjusting parameters or input examples.")
    else:
        st.error("Please provide at least one known LOHC SMILES string.")

# Export Generated SMILES to CSV
if st.button("Export SMILES to CSV"):
    if smiles_input:
        try:
            df = pd.DataFrame({"Generated LOHC SMILES": generated_smiles})
            df.to_csv("generated_lohc_smiles.csv", index=False)
            st.success("LOHC SMILES exported successfully as CSV!")
        except Exception as e:
            st.error(f"Failed to export SMILES: {e}")
    else:
        st.error("No SMILES available to export. Please generate some first.")

