import requests
import os
import re

cif_dir = 'cifs'
os.makedirs(cif_dir, exist_ok=True)

pdb_ids = ['7TRE', '7TRD', '7TRF', '7TRC', '2I46']

for pid in pdb_ids:
    cif_path = os.path.join(cif_dir, f'{pid}.cif')
    if not os.path.exists(cif_path):
        url = f'https://files.rcsb.org/download/{pid}.cif'
        print(f'Downloading {url}...')
        r = requests.get(url)
        with open(cif_path, 'wb') as f:
            f.write(r.content)
    print(f'{pid}.cif ready ({os.path.getsize(cif_path)} bytes)')

print('All CIF files downloaded.')
