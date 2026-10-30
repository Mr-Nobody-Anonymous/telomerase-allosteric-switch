import urllib.request

def inspect_cif_tert(pdb_id, target_auth_chain):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    print(f"Fetching {pdb_id}.cif header and atoms for TERT...")
    
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
                if line.startswith('ATOM'):
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

    print(f"=== {pdb_id} TERT (auth_chain {target_auth_chain}) ===")
    keys = list(seen_residues.keys())
    print(f"Total modeled residues: {len(keys)}")
    print(f"First 5 modeled: {[(k, seen_residues[k]) for k in keys[:5]]}")
    print(f"Last 5 modeled: {[(k, seen_residues[k]) for k in keys[-5:]]}")
    # Inspect catalytic triad / active site (e.g. Asp 712, Asp 868, Asp 869 in hTERT)
    for q in ['1', '10', '42', '100', '199', '330', '540', '712', '868', '869', '1120']:
        if q in seen_residues:
            print(f"  auth_seq {q}: {seen_residues[q]}")
        else:
            print(f"  auth_seq {q}: NOT MODELED")

inspect_cif_tert('7TRE', 'A')
inspect_cif_tert('7TRD', 'A')
inspect_cif_tert('7TRF', 'A')
