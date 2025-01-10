import streamlit as st
import openai
import os

# Set the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Load Socratic principles from Socratic_Prompts.md
def load_principles(file_path="Socratic_Prompts.md"):
    with open(file_path, "r") as file:
        return file.read()

socratic_prompts_md = load_principles()

# App title
st.title("AI-Powered Socratic Prompt Engineering with RAG")

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = [
        {"role": "system", "content": "You are a Socratic assistant specializing in chemistry and materials discovery. Your role is to transform user prompts into Socratic prompts based on an appropriate principle and generate follow-up questions."}
    ]

# Function to query the LLM
def get_response(user_prompt):
    st.session_state.messages.append({"role": "user", "content": user_prompt})
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=st.session_state.messages,
        max_tokens=1000,
        temperature=0.7,
    )
    reply = response["choices"][0]["message"]["content"]
    st.session_state.messages.append({"role": "assistant", "content": reply})
    return reply

# User input
st.text("Enter your chemistry or materials-related prompt below:")
user_input = st.text_input("Your Prompt:")

if st.button("Submit"):
    if user_input:
        # Construct the RAG-enhanced prompt
        rag_prompt = f"""
        The following is a reference document containing information about Socratic principles and their applications:

        {socratic_prompts_md}

        Based on the user's input, determine the most appropriate Socratic principle from the document above and justify your choice. Then, transform the user's prompt into a Socratic prompt using the selected principle and generate three follow-up questions.

        **Output Format**:
        Selected Principle: [Principle Name]
        Justification: [Why this principle was chosen]
        Socratic Prompt: [Transformed Socratic prompt]
        Follow-Up Questions:
        1. [First follow-up question]
        2. [Second follow-up question]
        3. [Third follow-up question]

        User's Prompt: {user_input}
        """

        # Get the LLM response
        response = get_response(rag_prompt)

        # Parse and display the response
        if "Selected Principle:" in response:
            try:
                principle = response.split("Selected Principle:")[1].split("Justification:")[0].strip()
                justification = response.split("Justification:")[1].split("Socratic Prompt:")[0].strip()
                socratic_prompt = response.split("Socratic Prompt:")[1].split("Follow-Up Questions:")[0].strip()
                follow_ups = response.split("Follow-Up Questions:")[1].strip()
            except IndexError:
                principle = "Error parsing selected principle."
                justification = "Error parsing justification."
                socratic_prompt = "Error parsing Socratic prompt."
                follow_ups = "Error parsing follow-up questions."
        else:
            principle = "No principle selected."
            justification = "No justification provided."
            socratic_prompt = "Error generating Socratic prompt."
            follow_ups = "No follow-up questions generated."

        # Display results
        st.markdown(f"### Selected Principle: {principle}")
        st.markdown(f"**Justification:** {justification}")
        st.markdown("### Socratic Prompt")
        st.markdown(socratic_prompt)
        st.markdown("### Follow-Up Questions")
        st.markdown(follow_ups)

# Optional: Allow user to upload a custom Socratic_Prompts.md file
st.sidebar.markdown("### Upload Custom Socratic Prompts File")
uploaded_file = st.sidebar.file_uploader("Upload your Socratic_Prompts.md file", type=["md"])
if uploaded_file:
    socratic_prompts_md = uploaded_file.read().decode("utf-8")
    st.sidebar.success("Custom file loaded successfully!")

