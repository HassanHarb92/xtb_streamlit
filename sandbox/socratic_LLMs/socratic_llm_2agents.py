import streamlit as st
import openai
import os

# Set the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# App title
st.title("Socratic Prompt Engineering for Chemistry")

# Refined System Prompt for the first LLM
system_prompt = """
You are a Socratic assistant specializing in chemistry and materials discovery. Your role is to transform user prompts into Socratic prompts that encourage critical thinking, exploration, and refinement of ideas. You must also generate three follow-up questions to deepen inquiry.

When transforming user input:
1. Reformulate the prompt into a Socratic prompt by clarifying terms, encouraging exploration, and guiding iterative reasoning.
2. Ensure the Socratic prompt is focused, precise, and adheres to the principles of the Socratic method (clarification, exploration, and critical thinking).
3. Generate three follow-up questions related to the same topic that align with Socratic principles and drive further exploration.

**Output Format** (always follow this structure):
Socratic Prompt:
[Your transformed Socratic prompt here]

Follow-Up Questions:
1. [First follow-up question]
2. [Second follow-up question]
3. [Third follow-up question]

Always adhere to this structure, even if the input is vague or unclear. If the input is unclear, reformulate it based on possible interpretations relevant to chemistry and materials discovery.
"""

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = [{"role": "system", "content": system_prompt}]

# Function to get a response from the first LLM
def get_socratic_response(user_prompt):
    st.session_state.messages.append({"role": "user", "content": user_prompt})
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=st.session_state.messages,
        max_tokens=500,
        temperature=0.7,
    )
    reply = response["choices"][0]["message"]["content"]
    st.session_state.messages.append({"role": "assistant", "content": reply})
    return reply

# Function to get the second LLM's response based on the Socratic prompt
def get_second_llm_response(socratic_prompt):
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[{"role": "user", "content": socratic_prompt}],
        max_tokens=1000,
        temperature=0.7,
    )
    return response["choices"][0]["message"]["content"]

# User input
st.text("Enter your chemistry or materials-related prompt below:")
user_input = st.text_input("Your Prompt:")

if st.button("Submit"):
    if user_input:
        socratic_response = get_socratic_response(user_input)

        # Parse the response for Socratic Prompt
        if "Socratic Prompt:" in socratic_response:
            try:
                socratic_prompt = socratic_response.split("Socratic Prompt:")[1].split("Follow-Up Questions:")[0].strip()
                follow_ups = socratic_response.split("Follow-Up Questions:")[1].strip()
            except IndexError:
                socratic_prompt = "Error parsing Socratic prompt. Please try again."
                follow_ups = "Error parsing follow-up questions. Please try again."
        else:
            socratic_prompt = "The response does not include a Socratic prompt. Please refine your input."
            follow_ups = "No follow-up questions were generated. Please try again."

        # Display the results from the first LLM
        st.markdown("### Socratic Prompt")
        st.markdown(socratic_prompt)
        st.markdown("### Suggested Follow-Up Questions")
        st.markdown(follow_ups)

        # Use the second LLM to generate a detailed response
        if "Error" not in socratic_prompt:
            second_llm_output = get_second_llm_response(socratic_prompt)
            st.markdown("### Second LLM Output")
            st.markdown(second_llm_output)

