import os
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

def add_hydrogens_and_generate_conformers(smiles, num_conformers=1):
    """
    Adds missing hydrogens, generates conformers, and returns the molecule with conformers.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)  # Generate conformers
    return mol

def optimize_geometry(mol):
    """
    Performs geometry optimization on the conformers and returns the optimized molecule.
    """
    for conf_id in range(mol.GetNumConformers()):
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)  # Perform geometry optimization using UFF
    return mol

def save_conformers(mol, name, output_dir, suffix=""):
    """
    Saves the conformers of a molecule in SDF format.
    """
    file_path = os.path.join(output_dir, f"{name}{suffix}.sdf")
    writer = Chem.SDWriter(file_path)
    for conf_id in range(mol.GetNumConformers()):
        mol.SetProp("_Name", f"{name}_conf{conf_id + 1}")  # Name each conformer
        writer.write(mol, confId=conf_id)
    writer.close()
    return file_path

def main(input_csv, output_dir, num_conformers=1):
    # Read CSV
    df = pd.read_csv(input_csv)
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process each SMILES and its corresponding name
    for index, row in df.iterrows():
        smiles = row[0]  # SMILES in the first column
        name = row[1]    # Name in the second column
        
        try:
            print(f"Processing: {name} - {smiles}")
            
            # Generate conformers
            mol_with_conformers = add_hydrogens_and_generate_conformers(smiles, num_conformers=num_conformers)
            
            # Save initial conformers
            save_conformers(mol_with_conformers, name, output_dir)

            # Perform geometry optimization
            optimized_mol = optimize_geometry(mol_with_conformers)

            # Save optimized conformers
            save_conformers(optimized_mol, name, output_dir, suffix="_optimized")
            
            print(f"Saved conformers for {name} to {output_dir}")
        
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add hydrogens, generate conformers, and perform geometry optimization on SMILES from a CSV file.")
    parser.add_argument("input_csv", help="Input CSV file with SMILES and names")
    parser.add_argument("output_dir", help="Output directory to save conformers")
    parser.add_argument("--num_conformers", type=int, default=1, help="Number of conformers to generate (default: 1)")

    args = parser.parse_args()

    main(args.input_csv, args.output_dir, args.num_conformers)

