############################################################################################################
#Extract the Lowest Energy Model Only.
#The script only keeps track of the lowest energy model and saves that one model to a PDBQT file.
#Save the Lowest Energy Model:
#The extract_lowest_energy_model function only saves the model with the lowest energy.
#File Naming:
#The PDBQT and PDB files are saved with a fixed name lowest_energy_model to reflect the lowest energy model.
#Running the Script
#python extract_and_convert_lowest_model.py docking_results.dlg output_directory
#Amirhossein Taghavi
#UF Scripps
#07/22/2024
##########################################################################################################

#!/usr/bin/env python

import os
import re
import argparse
from openbabel import openbabel

def extract_lowest_energy_model(dlg_file, output_dir):
    """
    Extract the lowest energy model from the DLG file and save it as a PDBQT file.
    Return the model number and the file name.
    """
    lowest_energy = float('inf')
    lowest_energy_model = None
    lowest_energy_lines = None

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    current_model = None
    model_lines = []

    with open(dlg_file, 'r') as f:
        for line in f:
            if line.startswith("DOCKED: MODEL"):
                if current_model is not None:
                    # Check if this model is the lowest energy
                    energy_match = re.search(r"Estimated Free Energy of Binding\s+=\s*(-?\d+\.\d+)", ''.join(model_lines))
                    if energy_match:
                        energy = float(energy_match.group(1))
                        if energy < lowest_energy:
                            lowest_energy = energy
                            lowest_energy_model = current_model
                            lowest_energy_lines = model_lines

                # Start a new model
                current_model = line.split()[2]
                model_lines = [line.strip()[7:].strip()]  # Remove 'DOCKED:' prefix
            elif line.startswith("DOCKED: ENDMDL"):
                if current_model is not None:
                    model_lines.append(line.strip()[7:].strip())  # Remove 'DOCKED:' prefix
                    # Save the last model if it is the lowest energy
                    energy_match = re.search(r"Estimated Free Energy of Binding\s+=\s*(-?\d+\.\d+)", ''.join(model_lines))
                    if energy_match:
                        energy = float(energy_match.group(1))
                        if energy < lowest_energy:
                            lowest_energy = energy
                            lowest_energy_model = current_model
                            lowest_energy_lines = model_lines

                    current_model = None
                    model_lines = []
            elif current_model is not None:
                model_lines.append(line.strip()[7:].strip())  # Remove 'DOCKED:' prefix

    if lowest_energy_model is None:
        raise ValueError("No models found in the DLG file.")

    # Save the lowest energy model
    pdbqt_file = os.path.join(output_dir, f"lowest_energy_model.pdbqt")
    with open(pdbqt_file, 'w') as out_f:
        for line in lowest_energy_lines:
            out_f.write(line + '\n')

    return lowest_energy_model, pdbqt_file

def convert_pdbqt_to_pdb(input_pdbqt, output_pdb):
    """
    Convert a PDBQT file to PDB format using Open Babel.
    """
    # Check if the PDBQT file exists and has content
    if not os.path.isfile(input_pdbqt):
        raise FileNotFoundError(f"The PDBQT file '{input_pdbqt}' does not exist.")
    
    if os.path.getsize(input_pdbqt) == 0:
        raise ValueError(f"The PDBQT file '{input_pdbqt}' is empty.")

    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdbqt")
    obConversion.SetOutFormat("pdb")

    mol = openbabel.OBMol()
    if not obConversion.ReadFile(mol, input_pdbqt):
        raise ValueError(f"Could not read the PDBQT file '{input_pdbqt}' using Open Babel.")
    
    if not obConversion.WriteFile(mol, output_pdb):
        raise ValueError(f"Could not write the PDB file '{output_pdb}' using Open Babel.")

def main():
    parser = argparse.ArgumentParser(description="Extract the lowest energy model from a DLG file, save it as PDBQT, and convert to PDB format.")
    parser.add_argument("dlg_file", help="DLG file from AutoDock")
    parser.add_argument("output_dir", help="Directory to save output files")

    args = parser.parse_args()

    # Extract the lowest energy model and save as PDBQT file
    _, pdbqt_file = extract_lowest_energy_model(args.dlg_file, args.output_dir)

    # Convert the PDBQT file to PDB format
    pdb_file = os.path.join(args.output_dir, "lowest_energy_model.pdb")
    convert_pdbqt_to_pdb(pdbqt_file, pdb_file)
    print(f"Lowest energy model saved as {pdbqt_file}")
    print(f"Converted to {pdb_file}")

if __name__ == "__main__":
    main()

