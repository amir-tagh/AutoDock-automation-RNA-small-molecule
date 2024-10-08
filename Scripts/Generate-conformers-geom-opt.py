import argparse
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel

def load_smiles_from_csv(input_csv):
    """Load SMILES strings from the CSV file."""
    df = pd.read_csv(input_csv)
    if 'smiles' not in df.columns:
        raise ValueError("The CSV file must contain a 'smiles' column.")
    return df['smiles'].tolist()

def add_hydrogens_and_generate_conformers(smiles, num_conformers=1):
    """
    Add hydrogens, generate 3D conformers for a molecule, and return the molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)  # Generate conformers
    return mol

def save_molecule_as_sdf(mol, file_path):
    """Save the RDKit molecule as an SDF file."""
    writer = Chem.SDWriter(file_path)
    writer.write(mol)
    writer.close()

def optimize_geometry(input_sdf, output_sdf):
    """
    Perform geometry optimization using Open Babel.
    Convert the molecule to 3D, optimize the geometry, and save as an SDF file.
    """
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("sdf")
    obConversion.SetOutFormat("sdf")

    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, input_sdf):
        raise ValueError(f"Could not read SDF file: {input_sdf}")

    # Perform geometry optimization
    forcefield = openbabel.OBForceField.FindForceField("mmff94")
    forcefield.Setup(mol)
    forcefield.ConjugateGradients(500, 1.0e-6)  # Optimize using conjugate gradients
    forcefield.GetCoordinates(mol)

    if not obConversion.WriteFile(mol, output_sdf):
        raise ValueError(f"Could not write optimized SDF file: {output_sdf}")

def process_smiles(smiles_list, output_dir, num_conformers):
    """
    Process a list of SMILES strings: add hydrogens, generate conformers, optimize geometry.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for idx, smiles in enumerate(smiles_list):
        try:
            print(f"Processing molecule {idx + 1}/{len(smiles_list)}: {smiles}")
            # Add hydrogens and generate conformers
            mol = add_hydrogens_and_generate_conformers(smiles, num_conformers)

            # Save to SDF
            sdf_file = os.path.join(output_dir, f"molecule_{idx + 1}.sdf")
            save_molecule_as_sdf(mol, sdf_file)

            # Perform geometry optimization
            optimized_sdf = os.path.join(output_dir, f"molecule_{idx + 1}_optimized.sdf")
            optimize_geometry(sdf_file, optimized_sdf)

            print(f"Molecule {idx + 1} optimized and saved to {optimized_sdf}")
        except Exception as e:
            print(f"Error processing SMILES: {smiles}. Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Process a list of SMILES: Add hydrogens, generate conformers, and perform geometry optimization.")
    parser.add_argument("input_csv", help="CSV file containing SMILES strings")
    parser.add_argument("output_dir", help="Directory to save the 3D structures and optimized geometries")
    parser.add_argument("--num_conformers", type=int, default=1, help="Number of conformers to generate for each molecule")

    args = parser.parse_args()

    # Load SMILES from CSV
    smiles_list = load_smiles_from_csv(args.input_csv)

    # Process each SMILES string
    process_smiles(smiles_list, args.output_dir, args.num_conformers)

if __name__ == "__main__":
    main()

