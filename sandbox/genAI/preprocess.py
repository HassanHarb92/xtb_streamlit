import pandas as pd
import json

# Load the QM9 dataset and extract SMILES strings
csv_file = 'QM9_G4MP2_all.csv'  # Path to your dataset
df = pd.read_csv(csv_file)

# Extract the 'unsat_SMILE' column
smiles_list = df['unsat_SMILE'].dropna().tolist()

# Convert SMILES to JSONL format for fine-tuning
output_file = 'smiles_fine_tune_data.jsonl'

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

