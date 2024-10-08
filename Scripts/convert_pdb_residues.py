import os

def convert_nucleotide_residue(line):
    """
    Converts non-standard nucleotide residues to standard ones, while maintaining the PDB formatting.
    """
    # Extract the residue name from the PDB line (columns 18-20, 1-based index)
    residue = line[17:20].strip()

    # Define mapping from non-standard to standard nucleotides
    conversion_map = {
        'RG5': 'G',   # Guanine
        'RA':  'A',   # Adenine
        'RC':  'C',   # Cytosine
        'RG':  'G',   # Guanine
        'RU':  'U',   # Uracil
        'RC3': 'C'    # Cytosine
    }

    # Convert the residue if it's in the map, otherwise leave it unchanged
    if residue in conversion_map:
        # Ensure we preserve the same line structure by using fixed-width formatting
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
                # Convert residues for ATOM and HETATM records
                line = convert_nucleotide_residue(line)
            # Write the line to the output file
            outfile.write(line)

def generate_output_filename(input_pdb):
    """
    Generates a new filename by appending '_corrected' before the file extension.
    """
    base, ext = os.path.splitext(input_pdb)
    return f"{base}_corrected{ext}"

def main(input_pdb, output_pdb=None):
    if not os.path.exists(input_pdb):
        print(f"Error: File {input_pdb} does not exist.")
        return

    # Automatically generate the output filename if not provided
    if output_pdb is None:
        output_pdb = generate_output_filename(input_pdb)

    convert_pdb_residues(input_pdb, output_pdb)
    print(f"Converted PDB file saved as {output_pdb}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert specific nucleotide residues in PDB files to standard ones.')
    parser.add_argument('input_pdb', type=str, help='Input PDB file')
    parser.add_argument('-o', '--output_pdb', type=str, default=None, help='Output PDB file (optional)')
    
    args = parser.parse_args()
    main(args.input_pdb, args.output_pdb)

