import os
import subprocess
import argparse
from Bio import PDB

def extract_first_model(pdb_file, output_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('rna', pdb_file)
    first_model = structure[0]
    
    io = PDB.PDBIO()
    io.set_structure(first_model)
    io.save(output_file)
    
    return output_file

def remove_water_and_ions(pdb_file, filtered_pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('rna', pdb_file)
    
    io = PDB.PDBIO()
    io.set_structure(structure)
    
    class NonWaterAndNonIonSelect(PDB.Select):
        def accept_residue(self, residue):
            return not (residue.get_resname() in ['HOH', 'WAT'] or 
                        residue.get_resname() in ['NA', 'CL', 'MG', 'CA', 'ZN','WO2'])
    
    io.save(filtered_pdb_file, NonWaterAndNonIonSelect())
    
    return filtered_pdb_file

def extract_ligand_info(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('rna', pdb_file)
    ligand_chain = None
    ligand_residue = None

    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()
                if isinstance(res_id, tuple) and res_id[0].startswith('H_'):
                    ligand_chain = chain.id
                    ligand_residue = res_id[1]
                    residue_name = residue.get_resname()
                    print(f"Found HETATM: Chain {ligand_chain}, Residue {ligand_residue}, Residue Name {residue_name}")

    if ligand_chain and ligand_residue:
        return ligand_chain, ligand_residue
    else:
        raise ValueError('No ligand found in the PDB file. Please ensure that the PDB file contains HETATM records for the ligand.')

def extract_ligand_only(pdb_file, output_pdb, ligand_chain, ligand_residue):
    pymol_script = f"""
    load {pdb_file}
    remove solvent
    remove (resn hoh or resn wat)
    select ligand, chain {ligand_chain} and resi {ligand_residue}
    cmd.do('print("Ligand selection count: ", cmd.count_atoms("ligand"))')
    save {output_pdb}, ligand
    cmd.do('print("Ligand saved to: {output_pdb}")')
    quit
    """

    with open('extract_ligand.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract ligand: {pdb_file}")
    subprocess.run(['pymol', '-cq', 'extract_ligand.pml'], check=True)
    os.remove('extract_ligand.pml')

    if os.path.getsize(output_pdb) == 0:
        print(f"Ligand extraction failed or empty for: {pdb_file}")
        return None

    return output_pdb

def extract_rna_only(pdb_file, output_pdb):
    pymol_script = f"""
    load {pdb_file}
    remove solvent
    remove (resn hoh or resn wat)
    select rna, polymer
    cmd.do('print("RNA selection count: ", cmd.count_atoms("rna"))')
    save {output_pdb}, rna
    cmd.do('print("RNA saved to: {output_pdb}")')
    quit
    """

    with open('extract_rna.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract RNA: {pdb_file}")
    subprocess.run(['pymol', '-cq', 'extract_rna.pml'], check=True)
    os.remove('extract_rna.pml')

    if os.path.getsize(output_pdb) == 0:
        print(f"RNA extraction failed or empty for: {pdb_file}")
        return None

    return output_pdb

def extract_binding_pocket(pdb_file, output_pdb, ligand_chain, ligand_residue):
    radius = 6.0

    pymol_script = f"""
    load {pdb_file}
    remove solvent
    remove (resn hoh or resn wat)
    select ligand, chain {ligand_chain} and resi {ligand_residue}
    cmd.do('print("Ligand selection count: ", cmd.count_atoms("ligand"))')
    select pocket, (all within {radius} of ligand and not resn hoh and not resn wat)
    cmd.do('print("Pocket selection count: ", cmd.count_atoms("pocket"))')
    save {output_pdb}, pocket
    cmd.do('print("Pocket saved to: {output_pdb}")')
    quit
    """

    with open('extract_pocket.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract binding pocket: {pdb_file}")
    subprocess.run(['pymol', '-cq', 'extract_pocket.pml'], check=True)
    os.remove('extract_pocket.pml')

def get_binding_box_coordinates(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('rna', pdb_file)

    min_coords = [float('inf')] * 3
    max_coords = [-float('inf')] * 3

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    for i in range(3):
                        if coord[i] < min_coords[i]:
                            min_coords[i] = coord[i]
                        if coord[i] > max_coords[i]:
                            max_coords[i] = coord[i]

    center_coords = [(min_coords[i] + max_coords[i]) / 2 for i in range(3)]
    dimensions = [max_coords[i] - min_coords[i] for i in range(3)]

    return center_coords, dimensions

def write_gpf_file(output_file, center_coords, dimensions):
    with open(output_file, 'w') as f:
        f.write(f"REMARK  Grid parameter file generated by script\n")
        f.write(f"REMARK  Binding box center: {center_coords}\n")
        f.write(f"REMARK  Binding box dimensions: {dimensions}\n")
        f.write("GRIDFILE  grid.dx\n")
        f.write("CENTER  {:.3f} {:.3f} {:.3f}\n".format(*center_coords))
        f.write("SPACING  0.375\n")
        f.write("XMIN     {:.3f}\n".format(center_coords[0] - dimensions[0] / 2))
        f.write("XMAX     {:.3f}\n".format(center_coords[0] + dimensions[0] / 2))
        f.write("YMIN     {:.3f}\n".format(center_coords[1] - dimensions[1] / 2))
        f.write("YMAX     {:.3f}\n".format(center_coords[1] + dimensions[1] / 2))
        f.write("ZMIN     {:.3f}\n".format(center_coords[2] - dimensions[2] / 2))
        f.write("ZMAX     {:.3f}\n".format(center_coords[2] + dimensions[2] / 2))

def prepare_ligand(pdb_file, ligand_pdb, ligand_pdbqt):
    # Convert ligand PDB to PDBQT using obabel
    command = [
        'obabel', ligand_pdb,
        '-O', ligand_pdbqt,
        '--partialcharge', 'gasteiger'
    ]
    subprocess.run(command, check=True)
    print(f"Ligand PDBQT file generated: {ligand_pdbqt}")

def prepare_receptor(pdb_file, receptor_pdb, receptor_pdbqt):
    # Convert receptor PDB to PDBQT using obabel
    command = [
        'obabel', receptor_pdb,
        '-O', receptor_pdbqt,
        '--partialcharge', 'gasteiger'
    ]
    subprocess.run(command, check=True)
    print(f"Receptor PDBQT file generated: {receptor_pdbqt}")

def process_pdb_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.lower().endswith('.pdb'):
            pdb_file = os.path.join(input_dir, filename)
            pdb_id = filename.split('.')[0].upper()

            pdb_output_dir = os.path.join(output_dir, pdb_id)
            if not os.path.exists(pdb_output_dir):
                os.makedirs(pdb_output_dir)

            first_model_pdb = os.path.join(pdb_output_dir, f"{pdb_id}_first_model.pdb")
            filtered_pdb_file = os.path.join(pdb_output_dir, f"{pdb_id}_filtered.pdb")
            ligand_pdb = os.path.join(pdb_output_dir, f"{pdb_id}_ligand.pdb")
            rna_pdb = os.path.join(pdb_output_dir, f"{pdb_id}_rna.pdb")
            binding_pocket_pdb = os.path.join(pdb_output_dir, f"{pdb_id}_pocket.pdb")
            gpf_file = os.path.join(pdb_output_dir, f"{pdb_id}.gpf")
            ligand_pdbqt = os.path.join(pdb_output_dir, f"{pdb_id}_ligand.pdbqt")
            receptor_pdbqt = os.path.join(pdb_output_dir, f"{pdb_id}_receptor.pdbqt")

            print(f"Processing file: {pdb_file}")

            extract_first_model(pdb_file, first_model_pdb)
            remove_water_and_ions(first_model_pdb, filtered_pdb_file)

            try:
                ligand_chain, ligand_residue = extract_ligand_info(filtered_pdb_file)
                ligand_pdb_file = extract_ligand_only(filtered_pdb_file, ligand_pdb, ligand_chain, ligand_residue)

                if ligand_pdb_file:
                    extract_binding_pocket(filtered_pdb_file, binding_pocket_pdb, ligand_chain, ligand_residue)
                    center_coords, dimensions = get_binding_box_coordinates(binding_pocket_pdb)
                    write_gpf_file(gpf_file, center_coords, dimensions)
                    
                    prepare_ligand(filtered_pdb_file, ligand_pdb, ligand_pdbqt)
                    prepare_receptor(filtered_pdb_file, filtered_pdb_file, receptor_pdbqt)
                else:
                    print(f"Ligand extraction failed or empty for: {pdb_file}")

                rna_pdb_file = extract_rna_only(filtered_pdb_file, rna_pdb)
                if not rna_pdb_file:
                    print(f"RNA extraction failed or empty for: {pdb_file}")

            except ValueError as e:
                print(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files to extract binding pockets, RNA, and generate GPF files for docking.")
    parser.add_argument('input_dir', type=str, help='Directory containing input PDB files.')
    parser.add_argument('output_dir', type=str, help='Directory to save the processed PDB files and GPF files.')

    args = parser.parse_args()

    process_pdb_files(args.input_dir, args.output_dir)

