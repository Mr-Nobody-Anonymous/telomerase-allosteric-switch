import math

# Parse atom coordinates for 7TRE
# We want TERT (label_asym_id 'B', auth_asym_id 'A') and TPP1 (label_asym_id 'D', auth_asym_id 'C')

def parse_atoms(cif_path):
    atoms = []
    with open(cif_path, 'r', encoding='utf-8') as f:
        in_atoms = False
        headers = []
        for line in f:
            line_str = line.strip()
            if line_str == 'loop_':
                in_atoms = False
                headers = []
                continue
            if line_str.startswith('_atom_site.'):
                in_atoms = True
                headers.append(line_str)
                continue
            if in_atoms:
                if line_str.startswith('#') or line_str.startswith('_') or line_str.startswith('loop_'):
                    in_atoms = False
                    continue
                parts = line_str.split()
                if len(parts) >= len(headers):
                    row = dict(zip(headers, parts[:len(headers)]))
                    if row.get('_atom_site.group_PDB') == 'ATOM':
                        atoms.append(row)
    return atoms

atoms = parse_atoms('cifs/7TRE.cif')
print(f'Total ATOM records in 7TRE: {len(atoms)}')

tert_atoms = [a for a in atoms if a.get('_atom_site.label_asym_id') == 'B']
tpp1_atoms = [a for a in atoms if a.get('_atom_site.label_asym_id') == 'D']
print(f'TERT atoms (label B / auth A): {len(tert_atoms)}')
print(f'TPP1 atoms (label D / auth C): {len(tpp1_atoms)}')

# Let's find contacting pairs within 4.5 A
contacts = []
for t1 in tpp1_atoms:
    x1, y1, z1 = float(t1['_atom_site.Cartn_x']), float(t1['_atom_site.Cartn_y']), float(t1['_atom_site.Cartn_z'])
    t1_seq = t1['_atom_site.auth_seq_id']
    t1_res = t1['_atom_site.auth_comp_id']
    t1_atom = t1['_atom_site.auth_atom_id']
    
    for t2 in tert_atoms:
        x2, y2, z2 = float(t2['_atom_site.Cartn_x']), float(t2['_atom_site.Cartn_y']), float(t2['_atom_site.Cartn_z'])
        # quick box filter
        if abs(x1 - x2) > 4.5 or abs(y1 - y2) > 4.5 or abs(z1 - z2) > 4.5:
            continue
        d = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        if d <= 4.5:
            t2_seq = t2['_atom_site.auth_seq_id']
            t2_res = t2['_atom_site.auth_comp_id']
            t2_atom = t2['_atom_site.auth_atom_id']
            contacts.append((t1_res, t1_seq, t1_atom, t2_res, t2_seq, t2_atom, d))

print(f'Total atomic contacts <= 4.5 A: {len(contacts)}')

# Summarize residue pairs
res_pairs = {}
for c in contacts:
    pair = (f"{c[0]}{c[1]}", f"{c[3]}{c[4]}")
    if pair not in res_pairs or c[6] < res_pairs[pair]:
        res_pairs[pair] = c[6]

print(f'Unique residue-residue contact pairs <= 4.5 A: {len(res_pairs)}')
print('Residue-residue contacts (TPP1 -> TERT):')
for pair, mindist in sorted(res_pairs.items(), key=lambda x: (int(x[0][0][3:]), x[1])):
    print(f"  TPP1 {pair[0]:>7} -- TERT {pair[1]:>7} (min dist: {mindist:.2f} A)")
