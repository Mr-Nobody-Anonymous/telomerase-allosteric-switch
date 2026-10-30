import urllib.request
import numpy as np

url = "https://files.rcsb.org/download/7TRE.cif"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
print("Parsing 7TRE atoms...")

atom_cols = []
in_loop = False
tert_atoms = [] # list of (auth_seq, resname, atom_name, x, y, z)
tpp1_atoms = []

with urllib.request.urlopen(req) as resp:
    for raw in resp:
        line = raw.decode('utf-8', errors='ignore').strip()
        if line.startswith('_atom_site.'):
            in_loop = True
            atom_cols.append(line.split('.')[1])
        elif in_loop:
            if line.startswith('ATOM'):
                parts = line.split()
                row = dict(zip(atom_cols, parts[:len(atom_cols)]))
                auth_asym = row.get('auth_asym_id')
                atom_type = row.get('type_symbol')
                if atom_type == 'H': continue # ignore hydrogens
                auth_seq = row.get('auth_seq_id')
                resname = row.get('label_comp_id')
                atom_name = row.get('auth_atom_id')
                x = float(row.get('Cartn_x'))
                y = float(row.get('Cartn_y'))
                z = float(row.get('Cartn_z'))
                
                if auth_asym == 'A': # TERT
                    tert_atoms.append((auth_seq, resname, atom_name, x, y, z))
                elif auth_asym == 'C': # TPP1
                    tpp1_atoms.append((auth_seq, resname, atom_name, x, y, z))
            elif line.startswith('#') or line.startswith('loop_'):
                if atom_cols: break

print(f"Loaded {len(tert_atoms)} TERT heavy atoms, {len(tpp1_atoms)} TPP1 heavy atoms.")

tert_coords = np.array([[a[3], a[4], a[5]] for a in tert_atoms], dtype=np.float32)
tpp1_coords = np.array([[a[3], a[4], a[5]] for a in tpp1_atoms], dtype=np.float32)

# Fast grid hashing neighbor search in pure NumPy
def find_contacts(coords_a, coords_b, cutoff):
    # Spatial cell hashing
    cell_size = cutoff
    grid = {}
    for idx_b, pt in enumerate(coords_b):
        cell = (int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size))
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx_b)
    
    contacts = []
    cutoff_sq = cutoff * cutoff
    for idx_a, pt in enumerate(coords_a):
        cx, cy, cz = int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor_cell = (cx + dx, cy + dy, cz + dz)
                    if neighbor_cell in grid:
                        for idx_b in grid[neighbor_cell]:
                            d_sq = np.sum((pt - coords_b[idx_b])**2)
                            if d_sq <= cutoff_sq:
                                contacts.append((idx_a, idx_b, np.sqrt(d_sq)))
    return contacts

for cutoff in [3.5, 4.0, 4.5]:
    atom_pairs = find_contacts(tert_coords, tpp1_coords, cutoff)
    res_pairs = set()
    for idx_a, idx_b, dist in atom_pairs:
        t = tert_atoms[idx_a]
        p = tpp1_atoms[idx_b]
        res_pairs.add((f"TERT_{t[1]}{t[0]}", f"TPP1_{p[1]}{p[0]}"))
    print(f"\n--- Cutoff {cutoff} Å ---")
    print(f"Total atom-atom contact pairs: {len(atom_pairs)}")
    print(f"Total residue-residue contact pairs: {len(res_pairs)}")
    tel_contacts = [rp for rp in res_pairs if any(k in rp[1] for k in ['168', '169', '170', '171', '180', '183', '212'])]
    print(f"TEL-patch residue contacts ({len(tel_contacts)}): {sorted(tel_contacts)}")
