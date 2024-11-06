#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

def calculate_properties(smiles):
    """Calculate various properties of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None

    # Calculating properties using RDKit
    mol_weight = Descriptors.MolWt(mol)
    rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    logP = Crippen.MolLogP(mol)
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    return {
        "Molecular Weight": mol_weight,
        "Rotatable Bonds": rotatable_bonds,
        "LogP": logP,
        "H-Bond Donors": h_donors,
        "H-Bond Acceptors": h_acceptors,
        "TPSA": tpsa
    }

def process_smiles_file(input_file='output_smiles.txt', output_file='output_properties.txt'):
    """Read SMILES from an input file (with compound IDs) and calculate properties, writing them to an output file."""
    with open(input_file, 'r') as in_file:
        with open(output_file, 'w') as out_file:
            # Write header to output file
            out_file.write("Compound ID\tSMILES\tMolecular Weight\tRotatable Bonds\tLogP\tH-Bond Donors\tH-Bond Acceptors\tTPSA\n")
            
            # Process each line in the input file (Compound ID + SMILES)
            for line in in_file:
                parts = line.strip().split("\t")  # Assuming compound ID and SMILES are tab-separated
                if len(parts) < 2:
                    continue  # Skip lines that don't contain a compound ID and SMILES
                
                compound_id, smiles = parts[0], parts[1]  # Split compound ID and SMILES

                # Calculate properties for this SMILES
                properties = calculate_properties(smiles)
                
                if properties:
                    # Write the results to the output file
                    out_file.write(f"{compound_id}\t{smiles}\t{properties['Molecular Weight']}\t{properties['Rotatable Bonds']}\t"
                                   f"{properties['LogP']}\t{properties['H-Bond Donors']}\t"
                                   f"{properties['H-Bond Acceptors']}\t{properties['TPSA']}\n")
                else:
                    print(f"Warning: Invalid SMILES string '{smiles}'")

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) > 2:
        # If command-line arguments are passed, use them
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        process_smiles_file(input_file, output_file)
    else:
        # Otherwise, use default input and output filenames
        process_smiles_file()
