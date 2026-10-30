import urllib.request
import re

def inspect_cif_residues(pdb_id, target_auth_chain):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    print(f"Fetching {pdb_id}.cif header and atoms...")
    
    atom_site_cols = []
    in_atom_site_loop = False
    
    seen_residues = {}
    with urllib.request.urlopen(req) as resp:
        for raw_line in resp:
            line = raw_line.decode('utf-8', errors='ignore').strip()
            if line.startswith('_atom_site.'):
                in_atom_site_loop = True
                col_name = line.split('.')[1]
                atom_site_cols.append(col_name)
            elif in_atom_site_loop:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    parts = line.split()
                    row = dict(zip(atom_site_cols, parts[:len(atom_site_cols)]))
                    auth_asym = row.get('auth_asym_id')
                    label_asym = row.get('label_asym_id')
                    if auth_asym == target_auth_chain:
                        auth_seq = row.get('auth_seq_id')
                        label_seq = row.get('label_seq_id')
                        comp_id = row.get('label_comp_id')
                        if auth_seq not in seen_residues:
                            seen_residues[auth_seq] = (label_seq, comp_id, auth_asym, label_asym)
                elif line.startswith('#') or line.startswith('loop_'):
                    if atom_site_cols:
                        break

    print(f"=== {pdb_id} (auth_chain {target_auth_chain}) ===")
    keys = list(seen_residues.keys())
    print(f"Total modeled residues: {len(keys)}")
    print(f"First 5 modeled: {[(k, seen_residues[k]) for k in keys[:5]]}")
    print(f"Last 5 modeled: {[(k, seen_residues[k]) for k in keys[-5:]]}")
    # Search for Glu 168-171, Leu 183, Leu 212
    for query in ['168', '169', '170', '171', '180', '183', '212']:
        if query in seen_residues:
            print(f"  auth_seq {query}: {seen_residues[query]}")
        else:
            print(f"  auth_seq {query}: NOT FOUND in auth_seq_id")

inspect_cif_residues('7TRE', 'C') # TPP1 in 7TRE
inspect_cif_residues('2I46', 'A') # TPP1 in 2I46
