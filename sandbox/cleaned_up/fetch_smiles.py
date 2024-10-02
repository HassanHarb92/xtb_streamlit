import pubchempy as pcp
import time
import streamlit as st

def fetch_smiles(molecule_name, retries=3, delay=5):
    attempt = 0
    while attempt < retries:
        try:
            compound = pcp.get_compounds(molecule_name, 'name')[0]
            return compound.isomeric_smiles
        except IndexError:
            st.error(f"No compound found for {molecule_name}")
            return None
        except Exception as e:
            st.error(f"Error fetching SMILES from PubChem: {e}")
            attempt += 1
            time.sleep(delay)  # Wait before retrying
    st.error("PubChem service is unavailable. Please try again later.")
    return None

