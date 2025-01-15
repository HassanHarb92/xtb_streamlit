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

# Input Sections for Big Picture and Paper-Specific Details
st.text("Provide the big picture and one-sentence summary for your scientific paper:")
big_picture_input = st.text_area("Big Picture of the Science (e.g., research on liquid organic hydrogen carriers):")
paper_summary_input = st.text_area("One Sentence About the Paper (e.g., In this study, we employ cheminformatics methods and quantum chemical calculations to explore...):")

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
def generate_intro_with_rag(big_picture, paper_summary):
    try:
        # Load the content of the guide
        guide_content = process_intro_guide(intro_guide_path)

        # Enhanced RAG prompt for stronger introductions
        rag_prompt = f"""
        You are an expert in scientific writing with a focus on creating impactful introductions for research papers. Below is a comprehensive guide on writing strong introductions:

        {guide_content}

        Using this guide, construct a compelling and comprehensive introduction based on the user's inputs. Follow these detailed instructions inspired by the book "Trees, Maps, and Theorems":

        1. **Context**: Begin broadly by defining the research territory. Introduce the general importance of the topic and provide sufficient background information to orient readers and establish relevance. Use the provided big picture for this step.
        2. **Need**: Narrow the focus by presenting the gap in current knowledge or limitations in the existing literature. Emphasize the contrast between the current state ("what we have") and the desired state ("what we want").
        3. **Task**: Clearly state what actions were undertaken to address the identified need. Use precise and action-oriented language to describe your contributions, incorporating the one-sentence summary provided.
        4. **Object**: Preview the structure and key findings of the paper. Prepare readers for what is to come without going into excessive detail.

        Ensure the introduction is:
        - Well-structured, logically flowing from general to specific.
        - Adhering to the principles and structure outlined in the guide.
        - Comprehensive, with a word count of no less than 800 words.
        - Written in a professional tone suitable for high-impact journals.

        Big Picture of the Science:
        {big_picture}

        One-Sentence Summary of the Paper:
        {paper_summary}

        The final introduction must effectively engage the reader, motivate the research, and provide a clear roadmap for the rest of the paper.
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
    if big_picture_input and paper_summary_input:
        introduction = generate_intro_with_rag(big_picture_input, paper_summary_input)
        st.markdown("### Generated Introduction")
        st.markdown(introduction)
    else:
        st.error("Please fill in both the Big Picture and One-Sentence Summary fields.")

# Export Generated Introduction
if st.button("Export Introduction"):
    try:
        with open("generated_introduction.md", "w") as f:
            f.write(introduction)
        st.success("Introduction exported successfully!")
    except Exception as e:
        st.error(f"Failed to export the introduction: {e}")


