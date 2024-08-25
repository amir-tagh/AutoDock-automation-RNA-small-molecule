import sys

def read_mol2_file(input_file):
    with open(input_file, 'r') as f:
        content = f.read()
    return content

def extract_molecules(content, start, end):
    molecules = content.split('@<TRIPOS>MOLECULE')
    
    if start < 1:
        start = 1
    if end > len(molecules) - 1:
        end = len(molecules) - 1

    selected_molecules = molecules[start:end + 1]
    
    # Add back the '@<TRIPOS>MOLECULE' identifier to each molecule
    selected_molecules = ['@<TRIPOS>MOLECULE' + mol for mol in selected_molecules]
    return ''.join(selected_molecules)

def write_mol2_file(output_file, selected_molecules):
    with open(output_file, 'w') as f:
        f.write(selected_molecules)

def main(input_file, output_file, start, end):
    content = read_mol2_file(input_file)
    selected_molecules = extract_molecules(content, start, end)
    write_mol2_file(output_file, selected_molecules)
    print(f"Extracted molecules from {start} to {end} and saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} <input_file> <output_file> <start> <end>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])

    main(input_file, output_file, start, end)

