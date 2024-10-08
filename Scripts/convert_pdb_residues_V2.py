import os

def convert_nucleotide_residue(line):
    """
    Converts non-standard nucleotide residues to standard ones, while maintaining the PDB formatting.
    """
    residue = line[17:20].strip()

    conversion_map = {
        'RG5': 'G',   # Guanine
        'RA':  'A',   # Adenine
        'RC':  'C',   # Cytosine
        'RG':  'G',   # Guanine
        'RU':  'U',   # Uracil
        'RC3': 'C'    # Cytosine
    }

    if residue in conversion_map:
        converted_residue = conversion_map[residue].rjust(3)
        line = line[:17] + converted_residue + line[20:]

    return line

def convert_pdb_residues(input_pdb, output_pdb):
    """
    Reads a PDB file, converts the residues, and writes the updated content to a new PDB file.
    """
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                line = convert_nucleotide_residue(line)
            outfile.write(line)

def generate_output_filename(input_pdb):
    """
    Generates a new filename by appending '_corrected' before the file extension.
    """
    base, ext = os.path.splitext(input_pdb)
    return f"{base}_corrected{ext}"

def process_directory(directory):
    """
    Processes all PDB files in the directory, sorts them alphabetically,
    converts the residues, and saves the corrected versions.
    """
    # Get all PDB files in the directory and sort them alphabetically
    pdb_files = sorted([f for f in os.listdir(directory) if f.lower().endswith('.pdb')])

    if not pdb_files:
        print(f"No PDB files found in directory: {directory}")
        return

    print(f"Processing {len(pdb_files)} PDB files in directory: {directory}")

    for pdb_file in pdb_files:
        input_pdb = os.path.join(directory, pdb_file)
        output_pdb = generate_output_filename(input_pdb)

        print(f"Processing file: {input_pdb}")
        convert_pdb_residues(input_pdb, output_pdb)
        print(f"Saved corrected PDB file as: {output_pdb}")

def main(input_pdb=None, directory=None):
    """
    Main function to process either a single PDB file or all PDB files in a directory.
    """
    if directory:
        if not os.path.isdir(directory):
            print(f"Error: Directory {directory} does not exist.")
            return
        process_directory(directory)
    elif input_pdb:
        if not os.path.exists(input_pdb):
            print(f"Error: File {input_pdb} does not exist.")
            return

        output_pdb = generate_output_filename(input_pdb)
        convert_pdb_residues(input_pdb, output_pdb)
        print(f"Converted PDB file saved as {output_pdb}")
    else:
        print("Error: You must provide either a PDB file or a directory containing PDB files.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert specific nucleotide residues in PDB files to standard ones.')
    parser.add_argument('-f', '--file', type=str, help='Input PDB file')
    parser.add_argument('-d', '--directory', type=str, help='Directory containing PDB files')
    
    args = parser.parse_args()

    if args.file:
        main(input_pdb=args.file)
    elif args.directory:
        main(directory=args.directory)
    else:
        print("Error: You must specify either a PDB file or a directory with the -f or -d option.")

