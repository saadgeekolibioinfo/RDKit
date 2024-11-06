#!/usr/bin/env python3

import sys
import os
from rdkit import Chem

def converter(directory):
    # Open the output file for writing SMILES strings
    with open('output_smiles.txt', "w") as out_file:
        # Loop through all files in the directory
        for file_name in os.listdir(directory):
            if file_name.endswith(".sdf"):  # Process only .sdf files
                sdf_file = os.path.join(directory, file_name)
                # Read the SDF file
                sdf_supplier = Chem.SDMolSupplier(sdf_file)
                
                for mol in sdf_supplier:
                    if mol is not None:  # Avoid compounds that cannot be loaded
                        smi = Chem.MolToSmiles(mol)
                        # Write SMILES to the output file
                        out_file.write(f"{file_name}\t{smi}\n")
                    else:
                        print(f"Warning: Could not process molecule in {file_name}")
    print("Conversion complete! SMILES saved to 'output_smiles.txt'")

if __name__ == "__main__":
    # Get the directory containing the .sdf files from the command-line argument
    directory = sys.argv[1]
    converter(directory)
