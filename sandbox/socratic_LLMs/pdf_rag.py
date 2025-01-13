import streamlit as st
import openai
import os

# Initialize the OpenAI API key
openai.api_key = os.getenv('OPENAI_API_KEY')

# Streamlit App Title
st.title("AI-Powered Book Q&A with RAG")

# Sidebar Configuration
st.sidebar.markdown("### Settings")
selected_temperature = st.sidebar.slider("Response Creativity (Temperature)", min_value=0.0, max_value=1.0, value=0.7)

# Input Section
st.text("Ask a question about the Socratic book:")
user_input = st.text_input("Your Question:")

# Define RAG Function to Answer Questions Based on the Book
def answer_question_with_rag(question):
    try:
        # RAG-enhanced prompt
        rag_prompt = f"""
        You are an assistant who answers questions specifically based on the content of a book about Socratic questioning. Use only the content of the book to answer the user's query.

        Question: {question}

        Provide a clear, concise answer based on the book's content.
        """

        # Query the OpenAI model
        response = openai.ChatCompletion.create(
            model="gpt-4",
            messages=[{"role": "system", "content": rag_prompt}],
            max_tokens=500,
            temperature=selected_temperature,
        )

        # Extract and return the response
        return response["choices"][0]["message"]["content"]
    except Exception as e:
        return f"An error occurred: {e}"

# Handle Submission
if st.button("Submit"):
    if user_input:
        answer = answer_question_with_rag(user_input)
        st.markdown("### Answer")
        st.markdown(answer)

# Export Conversation
if st.button("Export Conversation"):
    try:
        with open("conversation.md", "w") as f:
            for message in st.session_state.get("messages", []):
                f.write(f"{message['role'].capitalize()}: {message['content']}\n")
        st.success("Conversation exported successfully!")
    except Exception as e:
        st.error(f"Failed to export conversation: {e}")

