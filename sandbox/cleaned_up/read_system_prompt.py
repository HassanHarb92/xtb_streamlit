
import streamlit as st

def read_system_prompt(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except Exception as e:
        st.error(f"Error reading system prompt from file: {e}")
        return None

