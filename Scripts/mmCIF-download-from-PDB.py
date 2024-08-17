import os
import requests
import pandas as pd

def download_mmcif(pdb_id, output_dir):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(os.path.join(output_dir, f"{pdb_id}.cif"), 'wb') as file:
            file.write(response.content)
        print(f"Downloaded: {pdb_id}.cif")
    else:
        print(f"Failed to download: {pdb_id}")

def download_mmcif_files(pdb_ids, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for pdb_id in pdb_ids:
        download_mmcif(pdb_id, output_dir)

def get_pdb_ids_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    return df['id'].tolist()

# CSV file containing PDB IDs
csv_file = 'pdb_ids.csv'

# Output directory for downloaded mmCIF files
output_dir = 'mmCIF_files'

# Get the PDB IDs from the CSV file
pdb_ids = get_pdb_ids_from_csv(csv_file)

# Download the files
download_mmcif_files(pdb_ids, output_dir)

