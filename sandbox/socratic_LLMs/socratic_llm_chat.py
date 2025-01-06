import streamlit as st
import openai
import os

# App Sidebar
with st.sidebar:
    openai_api_key = os.getenv('OPENAI_API_KEY') # st.text_input("OpenAI API Key", type="password")
    st.markdown("[Get an OpenAI API key](https://platform.openai.com/account/api-keys)")
    st.markdown("[View the source code](https://github.com/streamlit/llm-examples/blob/main/Chatbot.py)")

# App Title
st.title("🧠 Socratic Chemistry Chatbot")

# Initialize session state
if "messages" not in st.session_state:
    st.session_state["messages"] = []
if "socratic_prompt" not in st.session_state:
    st.session_state["socratic_prompt"] = ""
if "follow_up_questions" not in st.session_state:
    st.session_state["follow_up_questions"] = []

# Functions
def get_socratic_prompt(user_prompt):
    """First agent to convert user prompt into a Socratic prompt with follow-up questions."""
    system_prompt = """
    You are a Socratic assistant specializing in chemistry and materials discovery. 
    Your role is to:
    1. Convert user input into a Socratic prompt that encourages critical thinking and exploration.
    2. Generate three follow-up questions based on the topic for deeper exploration.

    Output format:
    Socratic Prompt:
    [Your Socratic prompt here]

    Follow-Up Questions:
    1. [First follow-up question]
    2. [Second follow-up question]
    3. [Third follow-up question]
    """
    messages = [{"role": "system", "content": system_prompt}, {"role": "user", "content": user_prompt}]
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=messages,
        max_tokens=500,
        temperature=0.7,
    )
    return response["choices"][0]["message"]["content"]

def get_model_response(prompt):
    """Second agent to get the LLM's response based on the Socratic prompt."""
    messages = [{"role": "user", "content": prompt}]
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=messages,
        max_tokens=1000,
        temperature=0.7,
    )
    return response["choices"][0]["message"]["content"]

# User input
if prompt := st.chat_input():
    if not openai_api_key:
        st.info("Please add your OpenAI API key to continue.")
        st.stop()

    # Get Socratic Prompt and Follow-Up Questions
    socratic_output = get_socratic_prompt(prompt)
    try:
        socratic_prompt = socratic_output.split("Socratic Prompt:")[1].split("Follow-Up Questions:")[0].strip()
        follow_up_questions = socratic_output.split("Follow-Up Questions:")[1].strip()
    except IndexError:
        socratic_prompt = "Error parsing Socratic prompt. Please try again."
        follow_up_questions = "Error parsing follow-up questions. Please try again."

    # Save Socratic prompt and follow-ups
    st.session_state.socratic_prompt = socratic_prompt
    st.session_state.follow_up_questions = follow_up_questions.split("\n")

    # Get response to the Socratic prompt
    response = get_model_response(socratic_prompt)

    # Save messages
    st.session_state.messages.append({"role": "user", "content": prompt})
    st.session_state.messages.append({"role": "assistant", "content": response})

# Display chat messages
for msg in st.session_state.messages:
    st.chat_message(msg["role"]).write(msg["content"])

# Display follow-up questions
if st.session_state.follow_up_questions:
    st.markdown("### Follow-Up Questions:")
    for i, question in enumerate(st.session_state.follow_up_questions, start=1):
        if st.button(f"Ask: {question}", key=f"follow_up_{i}"):
            follow_up_response = get_model_response(question)
            st.session_state.messages.append({"role": "user", "content": question})
            st.session_state.messages.append({"role": "assistant", "content": follow_up_response})
            st.chat_message("user").write(question)
            st.chat_message("assistant").write(follow_up_response)

# Download button
if st.button("Download Script"):
    chat_history = "\n".join([f"{msg['role'].capitalize()}: {msg['content']}" for msg in st.session_state.messages])
    st.download_button(
        label="Download Conversation",
        data=chat_history,
        file_name="socratic_chemistry_chat.txt",
        mime="text/plain",
    )

