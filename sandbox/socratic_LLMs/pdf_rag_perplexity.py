import streamlit as st
import openai
import os
import requests

# Initialize the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')
perplexity_api_key = os.getenv('PERPLEXITY_API_KEY')

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

# Function to query Perplexity API for citations
def get_citations(text):
    try:
        headers = {
            "Authorization": f"Bearer {perplexity_api_key}",
            "Content-Type": "application/json",
        }
        payload = {"query": text, "num_results": 5}
        response = requests.post("https://api.perplexity.ai/v1/search", headers=headers, json=payload)

        if response.status_code == 200:
            results = response.json().get("results", [])
            citations = [result.get("url", "No URL Found") for result in results]
            return citations
        else:
            return [f"Error: {response.status_code} - {response.text}"]
    except Exception as e:
        return [f"An error occurred while fetching citations: {e}"]

# RAG Function to Generate Introductions Based on the Guide
def generate_intro_with_rag(guidelines):
    try:
        # Load the content of the guide
        guide_content = process_intro_guide(intro_guide_path)

        # Enhanced RAG prompt for stronger introductions
        rag_prompt = f"""
        You are an expert in scientific writing with a focus on creating impactful introductions for research papers. Below is a comprehensive guide on writing strong introductions:

        {guide_content}

        Using this guide, construct a compelling and comprehensive introduction based on the user's guidelines. Follow these instructions:

        1. Begin with a powerful and general opening sentence that captures the significance of the topic and motivates the reader to continue.
        2. Gradually narrow down from the general context to the specific focus of the research, providing a logical and compelling flow.
        3. Contextualize the research area by referencing key challenges, existing gaps, or opportunities in the field.
        4. Clearly articulate the purpose of the study and its broader relevance to science, technology, or society.
        5. Expand extensively on the user's provided guidelines, ensuring a word count of no less than 800 words.
        6. Adhere strictly to the structure, flow, and principles outlined in the provided guide.
        7. Maintain clarity, specificity, and adherence to scientific writing conventions while avoiding vague or generic statements.

        User Guidelines:
        {guidelines}

        Ensure the introduction is comprehensive, logically structured, and exceeds the quality expected in high-impact journals. The final output must strictly adhere to the principles and structure described in the guide.
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

        # Generate citations for each sentence in the introduction
        st.markdown("### Citations")
        sentences = introduction.split(".")
        for i, sentence in enumerate(sentences):
            if sentence.strip():
                citations = get_citations(sentence)
                st.markdown(f"**Sentence {i+1}:** {sentence.strip()}")
                st.markdown("Citations:")
                for citation in citations:
                    st.markdown(f"- {citation}")

# Export Generated Introduction
if st.button("Export Introduction"):
    try:
        with open("generated_introduction.md", "w") as f:
            f.write(introduction)
        st.success("Introduction exported successfully!")
    except Exception as e:
        st.error(f"Failed to export the introduction: {e}")

