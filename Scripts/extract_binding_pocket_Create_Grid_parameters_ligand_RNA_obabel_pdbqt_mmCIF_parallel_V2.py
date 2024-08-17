import os
import subprocess
import argparse
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from multiprocessing import Pool

def extract_first_model(cif_file, output_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('rna', cif_file)
    first_model = structure[0]
    
    io = MMCIFIO()
    io.set_structure(first_model)
    io.save(output_file)
    
    return output_file

def remove_water_and_ions(cif_file, filtered_cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('rna', cif_file)
    
    io = MMCIFIO()
    io.set_structure(structure)
    
    class NonWaterAndNonIonSelect(Select):
        def accept_residue(self, residue):
            return not (residue.get_resname() in ['HOH', 'WAT'] or 
                        residue.get_resname() in ['NA', 'CL', 'MG', 'CA', 'ZN','WO2','CD','K','ACA','PHA'])
    
    io.save(filtered_cif_file, NonWaterAndNonIonSelect())
    
    return filtered_cif_file

def extract_ligand_info(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('rna', cif_file)
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
        raise ValueError('No ligand found in the mmCIF file. Please ensure that the mmCIF file contains HETATM records for the ligand.')

def extract_ligand_only(cif_file, output_cif, ligand_chain, ligand_residue):
    pymol_script = f"""
    load {cif_file}
    remove solvent
    remove (resn hoh or resn wat)
    select ligand, chain {ligand_chain} and resi {ligand_residue}
    cmd.do('print("Ligand selection count: ", cmd.count_atoms("ligand"))')
    save {output_cif}, ligand
    cmd.do('print("Ligand saved to: {output_cif}")')
    quit
    """

    with open('extract_ligand.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract ligand: {cif_file}")
    subprocess.run(['pymol', '-cq', 'extract_ligand.pml'], check=True)
    os.remove('extract_ligand.pml')

    if os.path.getsize(output_cif) == 0:
        print(f"Ligand extraction failed or empty for: {cif_file}")
        return None

    return output_cif

def extract_rna_only(cif_file, output_cif):
    pymol_script = f"""
    load {cif_file}
    remove solvent
    remove (resn hoh or resn wat)
    # Select only RNA by including common RNA residue names
    select rna, (polymer and (resn A or resn C or resn G or resn U))
    cmd.do('print("RNA selection count: ", cmd.count_atoms("rna"))')
    save {output_cif}, rna
    cmd.do('print("RNA saved to: {output_cif}")')
    quit
    """

    with open('extract_rna.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract RNA: {cif_file}")
    subprocess.run(['pymol', '-cq', 'extract_rna.pml'], check=True)
    os.remove('extract_rna.pml')

    if os.path.getsize(output_cif) == 0:
        print(f"RNA extraction failed or empty for: {cif_file}")
        return None

    return output_cif

def extract_binding_pocket(cif_file, output_cif, ligand_chain, ligand_residue):
    radius = 10.0

    pymol_script = f"""
    load {cif_file}
    remove solvent
    remove (resn hoh or resn wat)
    select ligand, chain {ligand_chain} and resi {ligand_residue}
    cmd.do('print("Ligand selection count: ", cmd.count_atoms("ligand"))')
    select pocket, (all within {radius} of ligand and not resn hoh and not resn wat)
    cmd.do('print("Pocket selection count: ", cmd.count_atoms("pocket"))')
    save {output_cif}, pocket
    cmd.do('print("Pocket saved to: {output_cif}")')
    quit
    """

    with open('extract_pocket.pml', 'w') as f:
        f.write(pymol_script)

    print(f"Running PyMOL script to extract binding pocket: {cif_file}")
    subprocess.run(['pymol', '-cq', 'extract_pocket.pml'], check=True)
    os.remove('extract_pocket.pml')

def get_binding_box_coordinates(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('rna', cif_file)

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
        f.write("gridcenter  {:.3f} {:.3f} {:.3f}\n".format(*center_coords))
        f.write("npts 60 60 60\n")
        #f.write(f"REMARK  Grid parameter file generated by script\n")
        #f.write(f"REMARK  Binding box center: {center_coords}\n")
        #f.write(f"REMARK  Binding box dimensions: {dimensions}\n")
        #f.write("GRIDFILE  grid.dx\n")
        #f.write("CENTER  {:.3f} {:.3f} {:.3f}\n".format(*center_coords))
        #f.write("SPACING  0.375\n")
        #f.write("XMIN     {:.3f}\n".format(center_coords[0] - dimensions[0] / 2))
        #f.write("XMAX     {:.3f}\n".format(center_coords[0] + dimensions[0] / 2))
        #f.write("YMIN     {:.3f}\n".format(center_coords[1] - dimensions[1] / 2))
        #f.write("YMAX     {:.3f}\n".format(center_coords[1] + dimensions[1] / 2))
        #f.write("ZMIN     {:.3f}\n".format(center_coords[2] - dimensions[2] / 2))
        #f.write("ZMAX     {:.3f}\n".format(center_coords[2] + dimensions[2] / 2))

def prepare_ligand(cif_file, ligand_cif, ligand_pdbqt):
    # Convert ligand CIF to PDBQT using obabel
    command = [
        'obabel', ligand_cif,
        '-O', ligand_pdbqt,
        '--partialcharge', 'gasteiger'
    ]
    subprocess.run(command, check=True)
    print(f"Ligand PDBQT file generated: {ligand_pdbqt}")

def prepare_receptor(cif_file, receptor_cif, receptor_pdbqt):
    # Convert receptor CIF to PDBQT using obabel
    command = [
        'obabel', receptor_cif,
        '-O', receptor_pdbqt,
        '--partialcharge', 'gasteiger'
    ]
    subprocess.run(command, check=True)
    print(f"Receptor PDBQT file generated: {receptor_pdbqt}")

def process_single_cif(cif_file, output_dir):
    pdb_id = os.path.basename(cif_file).split('.')[0].upper()
    cif_output_dir = os.path.join(output_dir, pdb_id)
    if not os.path.exists(cif_output_dir):
        os.makedirs(cif_output_dir)

    first_model_cif = os.path.join(cif_output_dir, f'{pdb_id}_first_model.cif')
    filtered_cif = os.path.join(cif_output_dir, f'{pdb_id}_filtered.cif')
    ligand_cif = os.path.join(cif_output_dir, f'{pdb_id}_ligand.cif')
    rna_cif = os.path.join(cif_output_dir, f'{pdb_id}_rna.cif')
    pocket_cif = os.path.join(cif_output_dir, f'{pdb_id}_pocket.cif')
    gpf_file = os.path.join(cif_output_dir, f'{pdb_id}.gpf')
    ligand_pdbqt = os.path.join(cif_output_dir, f'{pdb_id}_ligand.pdbqt')
    receptor_pdbqt = os.path.join(cif_output_dir, f'{pdb_id}_receptor.pdbqt')

    print(f"Processing file: {cif_file}")

    extract_first_model(cif_file, first_model_cif)
    remove_water_and_ions(first_model_cif, filtered_cif)
    ligand_chain, ligand_residue = extract_ligand_info(filtered_cif)
    extract_ligand_only(filtered_cif, ligand_cif, ligand_chain, ligand_residue)
    extract_rna_only(filtered_cif, rna_cif)
    extract_binding_pocket(filtered_cif, pocket_cif, ligand_chain, ligand_residue)
    center_coords, dimensions = get_binding_box_coordinates(pocket_cif)
    write_gpf_file(gpf_file, center_coords, dimensions)
    prepare_ligand(filtered_cif, ligand_cif, ligand_pdbqt)
    prepare_receptor(filtered_cif, rna_cif, receptor_pdbqt)

def process_cif_files(input_dir, output_dir, num_cpus):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cif_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.lower().endswith('.cif')]
    
    with Pool(processes=num_cpus) as pool:
        pool.starmap(process_single_cif, [(cif_file, output_dir) for cif_file in cif_files])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process mmCIF files to extract RNA, ligands, and binding pocket information.')
    parser.add_argument('input_dir', type=str, help='Directory containing mmCIF files')
    parser.add_argument('output_dir', type=str, help='Directory to save processed files')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use for parallel processing')

    args = parser.parse_args()
    
    process_cif_files(args.input_dir, args.output_dir, args.cpus)

