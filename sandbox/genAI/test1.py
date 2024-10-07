import pandas as pd
import json
import os
import openai

# Initialize OpenAI API using environment variable
openai.api_key = os.getenv('OPENAI_API_KEY')

# Step 1: Load the QM9 dataset and extract SMILES

def load_smiles_from_csv(csv_file, column_name='unsat_SMILE'):
    """
    Load SMILES strings from the specified column of a CSV file.
    
    :param csv_file: Path to the CSV file containing SMILES strings.
    :param column_name: Name of the column containing SMILES.
    
    :return: List of SMILES strings.
    """
    df = pd.read_csv(csv_file)
    smiles_list = df[column_name].dropna().tolist()  # Drop any rows with NaN
    return smiles_list

# Step 2: Format SMILES for fine-tuning in JSONL format

def format_smiles_for_finetuning(smiles_list, output_file='smiles_fine_tune_data.jsonl'):
    """
    Format SMILES strings as prompt-completion pairs for fine-tuning and save to a JSONL file.
    
    :param smiles_list: List of SMILES strings.
    :param output_file: Output file in JSONL format.
    """
    with open(output_file, 'w') as outfile:
        for smiles in smiles_list:
            for i in range(1, len(smiles)):
                prompt = smiles[:i]
                completion = smiles[i]
                entry = {
                    "prompt": prompt,
                    "completion": completion
                }
                json.dump(entry, outfile)
                outfile.write("\n")
    print(f"Data has been formatted and saved to {output_file}")

# Step 3: Fine-Tuning the model (use CLI for this step)

# You will use OpenAI CLI commands to fine-tune the model, which we assume you're familiar with:
# Run this in your terminal after preparing the data:
# openai tools fine_tunes.prepare_data -f smiles_fine_tune_data.jsonl
# openai api fine_tunes.create -t <your-file-id> -m "gpt-3.5-turbo"

# Step 4: Generate new SMILES strings using the fine-tuned model

def generate_smiles(prompt, model_name, num_samples=5, max_length=50, temperature=0.7):
    """
    Generate new SMILES strings based on the input prompt using a fine-tuned model.
    
    :param prompt: Input SMILES string as a seed for generation.
    :param model_name: Fine-tuned model ID.
    :param num_samples: Number of SMILES strings to generate.
    :param max_length: Maximum length of generated sequences.
    :param temperature: Sampling temperature (0-1) for controlling randomness.
    
    :return: List of generated SMILES strings.
    """
    response = openai.Completion.create(
        model=model_name,
        prompt=prompt,
        max_tokens=max_length,
        n=num_samples,
        temperature=temperature
    )
    
    # Extract generated SMILES
    generated_smiles = [choice.text.strip() for choice in response.choices]
    return generated_smiles

# Step 5: Validate generated SMILES using RDKit

from rdkit import Chem

def validate_smiles(smiles_list):
    """
    Validate the generated SMILES strings using RDKit.
    
    :param smiles_list: List of SMILES strings to validate.
    
    :return: List of valid SMILES strings.
    """
    valid_smiles = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                valid_smiles.append(smi)
        except:
            continue
    return valid_smiles

# Example usage of the full pipeline

if __name__ == "__main__":
    # Step 1: Load SMILES from the CSV file
    csv_file = 'QM9_G4MP2_all.csv'  # Path to your CSV file
    smiles_list = load_smiles_from_csv(csv_file, column_name='unsat_SMILE')
    
    # Step 2: Format SMILES for fine-tuning
    format_smiles_for_finetuning(smiles_list, output_file='smiles_fine_tune_data.jsonl')
    
    # Step 3: Fine-tuning is done externally via OpenAI CLI as mentioned above
    
    # Step 4: Generate SMILES using the fine-tuned model (assuming the fine-tuning is complete)
    fine_tuned_model_name = "your-fine-tuned-model-id"  # Replace with your model ID from OpenAI
    
    # Generate new SMILES similar to an existing one
    seed_smiles = "C=C"  # Example seed SMILES from the dataset
    generated_smiles = generate_smiles(seed_smiles, model_name=fine_tuned_model_name, num_samples=5)
    
    print("Generated SMILES strings:")
    for idx, smi in enumerate(generated_smiles, 1):
        print(f"SMILES {idx}: {smi}")
    
    # Step 5: Validate the generated SMILES
    valid_smiles = validate_smiles(generated_smiles)
    print("\nValid SMILES strings:")
    for smi in valid_smiles:
        print(smi)

