# xTB2go

**Stable Version**: v6

## Capabilities

- **SMILES to 3D Structure Conversion**  
  Converts user-provided SMILES strings to 3D XYZ molecular structures using RDKit, enabling visualization and computational chemistry calculations.

- **Geometry Optimization**  
  Performs geometry optimization of molecular structures using the xTB tight binding program, providing optimized geometries and energy data.

- **Thermochemical Calculations**  
  Computes thermodynamic properties, including total energy, enthalpy, and free energy, using Hessian calculations. The results are extracted and presented clearly.

- **Frequency Analysis**  
  Extracts vibrational frequencies from xTB's output files to identify the lowest frequencies, which are crucial for understanding molecular vibrations and confirming optimized geometries.

- **Multiplicity and Spin State Control**  
  Automatically determines the appropriate multiplicity based on the number of electrons, allowing for accurate calculations of open-shell and closed-shell molecules.

- **Solvent Effects Modeling**  
  Supports implicit solvation models (e.g., ALPB and GBSA) to study molecules in different solvent environments, enhancing the relevance of calculations to experimental conditions.

- **Customizable Calculation Parameters**  
  Provides users with options to customize charge, method, and solvent model settings, allowing for tailored quantum chemical calculations.

- **Dynamic Visualization of Molecules**  
  Integrates with py3Dmol for interactive 3D visualization of molecular structures directly in the web app interface.

- **Downloadable Calculation Results**  
  Offers easy download of calculation results, including energies, optimized geometries, and frequencies, in a user-friendly `.txt` format.

- **User-Friendly Streamlit Interface**  
  Utilizes a simple, interactive web interface built with Streamlit, making advanced quantum chemical calculations accessible to users without a steep learning curve.

- **Error Handling and User Guidance**  
  Provides real-time feedback and error messages to guide users in correcting input errors or adjusting calculation parameters for successful runs.

- **Advanced Computational Chemistry Tools Integration**  
  Seamlessly integrates with xTB for advanced semi-empirical methods, supporting a range of calculations from basic energy evaluation to complex thermochemical properties and conformational analyses.

## Dependencies

- `streamlit`
- `subprocess`
- `os`
- `rdkit`
- `py3Dmol`
- `openai`

## In the works:

xTB calculator powered by large language model. Current version uses OpenAI API.

- **Version 1**: One LLM assistant --> gets xTB command line from prompt & SMILES from PubChem
- **Version 2**: Two LLM assistants --> Second one parses the xTB output file to get energies and properties

![LLM xTB Demo](video/LLM_xtb_demo.gif)



## To-Do List

1. **Conformer Search with CREST**  
   Implement conformer search functionality using CREST for more accurate conformations.

2. **Charge / Multiplicity Combinations**  
   Add support for specifying multiplicity in addition to charge for various molecular species.

3. **Hessian Calculation**  
   Perform Hessian calculations to obtain and print out thermochemistry data. **Status**: *Done*

4. **Ionization Potentials and Electron Affinities**  
   Introduce an option to calculate ionization potentials and electron affinities for molecules.

5. **Additional Outputs**  
   Add options to print:
   - **HOMO-LUMO Gap**
   - **Thermochemistry Data** (if available)

6. Run multiple jobs (from .csv file) and allow downloading csv output file

7. Connect to PubChem API – Similarity search

