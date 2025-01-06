import streamlit as st
import openai
import os

# Set the OpenAI API key from environment variables
openai.api_key = os.getenv("OPENAI_API_KEY")

# App title
st.title("💬 Socratic Chatbot for Chemistry & Materials")

# Sidebar with options
with st.sidebar:
    # "Download Script" button for conversation history
    script = "\n".join([f'{msg["role"]}: {msg["content"]}' for msg in st.session_state.get("messages", [])])
    st.download_button(
        label="Download Script",
        data=script,
        file_name="socratic_chat_script.txt",
        mime="text/plain",
    )

    # "Reset Chat" button
    if st.button("Reset Chat"):
        st.session_state["messages"] = [{"role": "assistant", "content": "How can I help you?"}]
        st.experimental_rerun()

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state["messages"] = [{"role": "assistant", "content": "How can I help you?"}]

# Function for Socratic prompt generation
def generate_socratic_prompt(user_prompt):
    socratic_system_prompt = """
    You are a Socratic assistant specializing in chemistry and materials discovery. Your role is to transform user prompts into Socratic prompts that encourage critical thinking, exploration, and refinement of ideas. You must also generate three follow-up questions to deepen inquiry.

    When transforming user input:
    1. Reformulate the prompt into a Socratic prompt by clarifying terms, encouraging exploration, and guiding iterative reasoning.
    2. Ensure the Socratic prompt is focused, precise, and adheres to the principles of the Socratic method (clarification, exploration, and critical thinking).
    3. Generate three follow-up questions related to the same topic that align with Socratic principles and drive further exploration.

    **Output Format**:
    Socratic Prompt:
    [Your transformed Socratic prompt here]

    Follow-Up Questions:
    1. [First follow-up question]
    2. [Second follow-up question]
    3. [Third follow-up question]
    """
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": socratic_system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        max_tokens=500,
        temperature=0.7,
    )
    return response["choices"][0]["message"]["content"]

# Function for querying the main LLM
def query_llm(prompt):
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[{"role": "user", "content": prompt}],
        max_tokens=1000,
        temperature=0.7,
    )
    return response["choices"][0]["message"]["content"]

# Chat interface
for msg in st.session_state.messages:
    st.chat_message(msg["role"]).write(msg["content"])

if prompt := st.chat_input("Enter your prompt"):
    st.session_state["messages"].append({"role": "user", "content": prompt})

    # Convert the user's input to a Socratic prompt and generate follow-ups
    socratic_response = generate_socratic_prompt(prompt)

    # Extract Socratic prompt and follow-up questions
    if "Socratic Prompt:" in socratic_response and "Follow-Up Questions:" in socratic_response:
        try:
            socratic_prompt = socratic_response.split("Socratic Prompt:")[1].split("Follow-Up Questions:")[0].strip()
            follow_ups = socratic_response.split("Follow-Up Questions:")[1].strip().split("\n")
        except IndexError:
            socratic_prompt = "Error parsing Socratic prompt. Please try again."
            follow_ups = ["Error parsing follow-up questions. Please try again."]
    else:
        socratic_prompt = "Error: Socratic prompt not generated. Please refine your input."
        follow_ups = ["No follow-up questions generated."]

    # Add the Socratic prompt to the chat
    st.session_state.messages.append({"role": "assistant", "content": f"Socratic Prompt: {socratic_prompt}"})

    # Query the main LLM with the Socratic prompt
    llm_response = query_llm(socratic_prompt)
    st.session_state.messages.append({"role": "assistant", "content": llm_response})

    # Add follow-up questions to the chat
    for i, question in enumerate(follow_ups, 1):
        question = question.strip()
        if question:
            follow_up_button = st.button(f"Follow-Up {i}: {question}")
            if follow_up_button:
                st.session_state["messages"].append({"role": "user", "content": question})
                st.experimental_rerun()  # Treat the follow-up as a new prompt

    # Display everything in the chat
    st.chat_message("assistant").write(f"Socratic Prompt: {socratic_prompt}")
    st.chat_message("assistant").write(llm_response)

