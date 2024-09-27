import argparse
from Bio import PDB

def extract_rna_sequence(pdb_file):
    """
    Extracts the RNA sequence(s) from the given PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        list: A list of RNA sequences, one per chain.
    """
    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("RNA_structure", pdb_file)

    rna_sequences = []

    # Loop through all chains in the structure
    for model in structure:
        for chain in model:
            chain_sequence = []
            for residue in chain:
                # Only extract RNA residues
                if residue.resname in ['A', 'G', 'C', 'U']:
                    nucleotide = three_letter_to_one_letter(residue.resname)
                    if nucleotide:
                        chain_sequence.append(nucleotide)
            if chain_sequence:
                rna_sequences.append("".join(chain_sequence))

    return rna_sequences

def three_letter_to_one_letter(res_name):
    """
    Converts a 3-letter RNA residue name to a 1-letter code.

    Args:
        res_name (str): The 3-letter residue name.

    Returns:
        str: Corresponding 1-letter code if valid, otherwise None.
    """
    rna_codes = {
        'A': 'A',  # Adenine
        'G': 'G',  # Guanine
        'C': 'C',  # Cytosine
        'U': 'U',  # Uracil
    }
    return rna_codes.get(res_name, None)

def save_sequence_to_file(rna_sequences, output_file):
    """
    Saves RNA sequences to a specified output file.

    Args:
        rna_sequences (list): List of RNA sequences.
        output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as f:
        for i, seq in enumerate(rna_sequences):
            f.write(f">Chain_{i+1}_RNA_Sequence\n")
            f.write(seq + "\n")
    print(f"RNA sequences saved to {output_file}")

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Extract RNA sequences from a PDB file and save them to a file.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("-o", "--output", default="rna_sequences.txt", help="Output file to save RNA sequences (default: rna_sequences.txt).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose mode to print sequences to the console.")

    # Parse arguments
    args = parser.parse_args()

    # Extract RNA sequences
    rna_sequences = extract_rna_sequence(args.pdb_file)

    # Save sequences to file
    save_sequence_to_file(rna_sequences, args.output)

    # Verbose mode: print sequences to console
    if args.verbose:
        for i, seq in enumerate(rna_sequences):
            print(f"Chain {i+1} RNA Sequence: {seq}")

if __name__ == "__main__":
    main()

