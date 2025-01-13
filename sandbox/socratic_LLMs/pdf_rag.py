import streamlit as st
import openai
import os

# Initialize the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Streamlit App Title
st.title("AI-Powered Scientific Paper Introduction Generator")

# Sidebar Configuration
st.sidebar.markdown("### Settings")
selected_temperature = st.sidebar.slider("Response Creativity (Temperature)", min_value=0.0, max_value=1.0, value=0.7)

# Input Section
st.text("Provide guidelines for your scientific paper's introduction:")
user_input = st.text_area("Your Guidelines:")

# Path to the introduction guide PDF
intro_guide_path = "/Users/hassan/Downloads/intro_paper.pdf"

# Function to read and process the PDF guide
def process_intro_guide(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()
        return content
    except Exception as e:
        return f"An error occurred while reading the guide: {e}"

# RAG Function to Generate Introductions Based on the Guide
def generate_intro_with_rag(guidelines):
    try:
        # Load the content of the guide
        guide_content = process_intro_guide(intro_guide_path)

        # Enhanced RAG prompt for stronger introductions
        rag_prompt = f"""
        You are an expert in scientific writing with a focus on creating impactful introductions for research papers. Below is a comprehensive guide on writing strong introductions:

        {guide_content}

        Using this guide, construct a compelling introduction based on the user's guidelines. Follow these instructions:

        1. Begin with a powerful opening sentence that captures the significance of the topic.
        2. Contextualize the research area by referencing key challenges, existing gaps, or opportunities in the field.
        3. Clearly articulate the purpose of the study and its broader relevance to science, technology, or society.
        4. Expand on the user's provided guidelines, ensuring logical flow and elaboration of key points.
        5. Maintain clarity, specificity, and adherence to scientific writing conventions while avoiding generic statements.

        User Guidelines:
        {guidelines}

        Ensure the introduction is well-structured, logically developed, and exceeds the quality expected in high-impact journals.
        """

        # Query the OpenAI model
        response = openai.ChatCompletion.create(
            model="gpt-4o",
            messages=[{"role": "system", "content": rag_prompt}],
            max_tokens=5000,
            temperature=selected_temperature,
        )

        # Extract and return the response
        return response["choices"][0]["message"]["content"]
    except Exception as e:
        return f"An error occurred: {e}"

# Handle Submission
if st.button("Generate Introduction"):
    if user_input:
        introduction = generate_intro_with_rag(user_input)
        st.markdown("### Generated Introduction")
        st.markdown(introduction)

# Export Generated Introduction
if st.button("Export Introduction"):
    try:
        with open("generated_introduction.md", "w") as f:
            f.write(introduction)
        st.success("Introduction exported successfully!")
    except Exception as e:
        st.error(f"Failed to export the introduction: {e}")

