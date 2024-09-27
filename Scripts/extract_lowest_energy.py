##################################################################################
#extract the lowest energy structure from *.dlg file
#convert the pdbqt-->pdb
##run the script as:
#python extract_and_convert_lowest_model.py docking_results.dlg output_directory
#Amirhossein Taghavi
#UF Scripps
#07/23/2024
##################################################################################



#!/usr/bin/env python

import os
import re
import argparse
from openbabel import openbabel

def extract_lowest_energy_structure(dlg_file, output_pdbqt):
    """
    Extract the lowest energy structure from the DLG file and save it as a PDBQT file,
    removing the 'DOCKED:' prefix from each line.
    """
    lowest_energy = float('inf')
    lowest_energy_structure = None

    with open(dlg_file, 'r') as f:
        model = []
        in_conformer = False
        for line in f:
            if line.startswith("DOCKED: MODEL"):
                model = [line.strip()[7:].strip()]  # Remove 'DOCKED:' prefix
                in_conformer = True
            elif line.startswith("DOCKED: ENDMDL"):
                if in_conformer:
                    model.append(line.strip()[7:].strip())  # Remove 'DOCKED:' prefix
                    in_conformer = False
                    energy_line = next((l for l in model if "Estimated Free Energy of Binding" in l), None)
                    if energy_line:
                        energy = float(re.search(r"[-+]?\d*\.\d+|\d+", energy_line).group())
                        if energy < lowest_energy:
                            lowest_energy = energy
                            lowest_energy_structure = model
            elif in_conformer:
                model.append(line.strip()[7:].strip())  # Remove 'DOCKED:' prefix

    if lowest_energy_structure is None:
        raise ValueError("No conformer found in the DLG file.")

    with open(output_pdbqt, 'w') as out_f:
        for line in lowest_energy_structure:
            out_f.write(line + '\n')

    return output_pdbqt

def convert_pdbqt_to_pdb(input_pdbqt, output_pdb):
    """
    Convert a PDBQT file to PDB format using Open Babel.
    """
    # Check if the PDBQT file exists and has content
    if not os.path.isfile(input_pdbqt):
        raise FileNotFoundError(f"The PDBQT file '{input_pdbqt}' does not exist.")
    
    if os.path.getsize(input_pdbqt) == 0:
        raise ValueError(f"The PDBQT file '{input_pdbqt}' is empty.")

    # Print the content of the PDBQT file for debugging
    with open(input_pdbqt, 'r') as file:
        content = file.read()
        print(f"Content of '{input_pdbqt}':\n{content}")

    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdbqt")
    obConversion.SetOutFormat("pdb")

    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, input_pdbqt):
        raise ValueError(f"Could not read the PDBQT file '{input_pdbqt}' using Open Babel.")
    
    if not obConversion.WriteFile(mol, output_pdb):
        raise ValueError(f"Could not write the PDB file '{output_pdb}' using Open Babel.")
    
    # Print the content of the PDB file for debugging
    with open(output_pdb, 'r') as file:
        pdb_content = file.read()
        print(f"Content of '{output_pdb}':\n{pdb_content}")

def main():
    parser = argparse.ArgumentParser(description="Extract the lowest energy structure from a DLG file and convert it to PDB format.")
    parser.add_argument("dlg_file", help="DLG file from AutoDock")
    parser.add_argument("output_pdb", help="Output PDB file")

    args = parser.parse_args()

    # Extract the lowest energy structure to a temporary PDBQT file
    temp_pdbqt = "temp_lowest_energy.pdbqt"
    extract_lowest_energy_structure(args.dlg_file, temp_pdbqt)

    # Convert the PDBQT file to PDB format
    convert_pdbqt_to_pdb(temp_pdbqt, args.output_pdb)

    # Clean up the temporary file
    os.remove(temp_pdbqt)

if __name__ == "__main__":
    main()

