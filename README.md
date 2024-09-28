# Download mmCIF files for RNA-small molecule complexes from the Protein Data Bank (PDB), follow these steps:

1. Read PDB IDs from the CSV File\
   1.1 download the CSV file from HARIBOSS database
2. Make sure you have the pdb_ids.csv file in the same directory as your script.
3. Run the script. It will read the PDB IDs from the CSV file and download the corresponding mmCIF files into the specified output directory.
This script automates the download process by using the PDB IDs listed in your CSV file.
# Prepare the pdbqt files for RNA and the corresponding small molecule (in progress)
1. extract the first model from ensembles (NMR with more than one structure)
2. Remove water, ions, WO2,ACA and PHA
3. Extract the ligand (pymol)
4. Extract the RNA only
5. Extract the binding pocket (residues within 10 Angstrom from the ligand)
6. Get the binding box coordinates
7. write the gpf file for AutoDock
8. Prepare the pdbqt file for ligand (obabel)
9. prepare the pdbqt file for the receptor (obabel)
# Multiprocessing to parallelize tasks across multiple CPU cores
Explanation:
cpu_count():

Automatically detects the number of available CPUs on your machine.
Pool(processes=num_cpus):

Creates a pool of worker processes equal to the number of available CPUs.
pool.starmap():

Distributes the download_mmcif function across the pool of worker processes, allowing multiple PDB IDs to be downloaded concurrently.

# Generate conformers from SMILES and geometry optimization

python Generate-conformers-geom-opt.py -h\
usage: Generate-conformers-geom-opt-V2.py [-h] [--num_conformers NUM_CONFORMERS] input_csv output_dir\
Add hydrogens, generate conformers, and perform geometry optimization on SMILES from a CSV file.\
positional arguments:\
  input_csv             Input CSV file with SMILES and names\
  output_dir            Output directory to save conformers\
optional arguments:\
  -h, --help            show this help message and exit\
  --num_conformers NUM_CONFORMERS (Number of conformers to generate (default: 1))

