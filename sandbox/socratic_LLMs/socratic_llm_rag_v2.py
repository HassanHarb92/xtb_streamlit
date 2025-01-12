import streamlit as st
import openai
import os
from graphviz import Digraph
from langchain.vectorstores import FAISS
from langchain.embeddings import OpenAIEmbeddings

# Set the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Load Socratic principles from Socratic_Prompts.md
def load_principles(file_path="Socratic_Prompts.md"):
    with open(file_path, "r") as file:
        return file.read()

socratic_prompts_md = load_principles()

# App title
st.title("AI-Powered Socratic Prompt Engineering with RAG")

# Initialize session state
if "messages" not in st.session_state:
    st.session_state.messages = [
        {"role": "system", "content": "You are a Socratic assistant specializing in chemistry and materials discovery. Your role is to transform user prompts into Socratic prompts based on an appropriate principle and generate follow-up questions."}
    ]

if "dialogue_history" not in st.session_state:
    st.session_state.dialogue_history = []

# Function to query the LLM
def get_response(user_prompt, temperature=0.7):
    st.session_state.messages.append({"role": "user", "content": user_prompt})
    response = openai.ChatCompletion.create(
        model="gpt-4",
        messages=st.session_state.messages,
        max_tokens=1000,
        temperature=temperature,
    )
    reply = response["choices"][0]["message"]["content"]
    st.session_state.messages.append({"role": "assistant", "content": reply})
    return reply

# Function to visualize dialogue
def visualize_dialogue(prompt, follow_ups):
    dot = Digraph()
    dot.node("User Input", prompt)
    for i, question in enumerate(follow_ups.split("\n")):
        dot.node(f"Q{i+1}", question)
        dot.edge("User Input", f"Q{i+1}")
    return dot

# Load vectorstore for RAG
embeddings = OpenAIEmbeddings()
vectorstore = FAISS.from_texts([socratic_prompts_md], embeddings)

# Sidebar
st.sidebar.markdown("### Settings")
selected_temperature = st.sidebar.slider("Response Creativity (Temperature)", min_value=0.0, max_value=1.0, value=0.7)

domain = st.sidebar.selectbox("Choose Inquiry Domain:", ["General", "Chemistry", "AI Ethics"])
st.sidebar.markdown("### Socratic Principles")
for principle, description in {
    "Definition": "Clarifies terms and concepts to ensure shared understanding.",
    "Generalization": "Derives broad principles from observed patterns or theories.",
    "Induction": "Forms hypotheses based on evidence, recognizing inherent uncertainties.",
    "Elenchus": "Tests the consistency of beliefs through cross-examination.",
    "Hypothesis Elimination": "Challenges assumptions by identifying counterexamples.",
    "Maieutics": "Encourages self-reflection to uncover existing knowledge.",
    "Dialectic": "Explores opposing perspectives to gain insights.",
    "Recollection": "Draws out latent knowledge through questioning.",
    "Irony": "Highlights gaps in understanding to provoke thought.",
    "Analogy": "Uses comparisons to clarify complex concepts."
}.items():
    st.sidebar.markdown(f"**{principle}:** {description}")

# User input
st.text("Enter your prompt below:")
user_input = st.text_input("Your Prompt:")

if st.button("Submit"):
    if user_input:
        # RAG-enhanced prompt
        relevant_context = vectorstore.similarity_search(user_input, k=1)
        rag_prompt = f"""
        The following is a reference document containing information about Socratic principles and their applications:

        {relevant_context[0].page_content}

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

        # Get response
        response = get_response(rag_prompt, temperature=selected_temperature)

        # Parse response
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

        # Display results
        st.markdown(f"### Selected Principle: {principle}")
        st.markdown(f"**Justification:** {justification}")
        st.markdown("### Socratic Prompt")
        st.markdown(socratic_prompt)
        st.markdown("### Follow-Up Questions")
        st.markdown(follow_ups)

        # Auto-prompt and answer
        if socratic_prompt:
            auto_prompt = f"""
            Based on the following Socratic prompt, provide a thoughtful and critical answer:
            {socratic_prompt}
            """
            auto_response = get_response(auto_prompt, temperature=selected_temperature)
            st.markdown("### Answer to the Socratic Prompt")
            st.markdown(auto_response)

        # Visualize dialogue
        st.markdown("### Dialogue Visualization")
        st.graphviz_chart(visualize_dialogue(user_input, follow_ups))

# Save conversation
if st.button("Export Conversation"):
    with open("conversation.md", "w") as f:
        for message in st.session_state.messages:
            f.write(f"{message['role'].capitalize()}: {message['content']}\n")
    st.success("Conversation exported successfully!")

