# Download mmCIF files for RNA-small molecule complexes from the Protein Data Bank (PDB), follow these steps:

1. Read PDB IDs from the CSV File\
   1.1 download the CSV file from HARIBOSS database
2. Make sure you have the pdb_ids.csv file in the same directory as your script.
3. Run the script. It will read the PDB IDs from the CSV file and download the corresponding mmCIF files into the specified output directory.
This script automates the download process by using the PDB IDs listed in your CSV file.
# Prepare the pdbqt files for RNA and the corresponding small molecule (in progress)
1. extract the first model from ensembles (NMR with more than one structure)
2. Remove water, ions, WO2,ACA and PHA
3. 
