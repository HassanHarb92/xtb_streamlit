# xTB2go

**Stable Version**: v6

## Dependencies

- `streamlit`
- `subprocess`
- `os`
- `rdkit`
- `py3Dmol`

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
