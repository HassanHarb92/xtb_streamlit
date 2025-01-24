import streamlit as st
import openai
import os

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
generation_count = st.number_input("Number of New LOHC SMILES to Generate:", min_value=1, max_value=100, value=10)

# Path to the LOHC knowledge document
lohc_pdf_path = "LOHC_DD.pdf"

# Function to generate new LOHC SMILES
def generate_lohc_smiles(smiles_list, num_generated):
    try:
        # Format input for the model
        formatted_smiles = ", ".join(smiles_list)
        
        # Prompt to guide the LLM in generating novel LOHC SMILES
        lohc_prompt = f"""
        You are an expert in molecular design with a specialization in Liquid Organic Hydrogen Carriers (LOHCs). 
        The user has provided a set of known LOHC SMILES strings:
        {formatted_smiles}
        
        Additionally, you have access to the document `{lohc_pdf_path}`, which contains detailed knowledge about LOHC molecular design principles. From that document, you need to learn the design rules from them as well.        
        Your task is to generate {num_generated} new LOHC SMILES strings that:
        - Are chemically valid.
        - Maintain key structural features relevant to efficient hydrogen storage and release.
        - Are not present in the provided set of known LOHC SMILES.
        - Exhibit potential for practical application based on insights from the LOHC_DD.pdf document.
        
        Provide the generated SMILES as a comma-separated list.
        """
        
        # Query the OpenAI model
        response = openai.ChatCompletion.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": f"Refer to `{lohc_pdf_path}` for additional LOHC knowledge."},
                {"role": "system", "content": lohc_prompt}
            ],
            max_tokens=1000,
            temperature=selected_temperature,
        )

        # Extract and return the response
        return response["choices"][0]["message"]["content"].split(",")
    except Exception as e:
        return [f"An error occurred: {e}"]

# Handle SMILES Generation
if st.button("Generate New LOHC SMILES"):
    if smiles_input:
        smiles_list = [s.strip() for s in smiles_input.split(",")]
        generated_smiles = generate_lohc_smiles(smiles_list, generation_count)
        
        st.markdown("### Generated LOHC SMILES")
        for smi in generated_smiles:
            st.markdown(f"- `{smi.strip()}`")
    else:
        st.error("Please provide at least one known LOHC SMILES string.")

# Export Generated SMILES
if st.button("Export SMILES"):
    if smiles_input:
        try:
            with open("generated_lohc_smiles.txt", "w") as f:
                f.write("\n".join(generated_smiles))
            st.success("LOHC SMILES exported successfully!")
        except Exception as e:
            st.error(f"Failed to export SMILES: {e}")
    else:
        st.error("No SMILES available to export. Please generate some first.")

