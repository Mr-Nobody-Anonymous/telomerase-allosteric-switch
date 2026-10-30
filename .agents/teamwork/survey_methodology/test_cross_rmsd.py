import urllib.request
import numpy as np

def parse_ca_atoms(pdb_id, auth_chain):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    print(f"Fetching {pdb_id}.cif...")
    
    atom_cols = []
    in_loop = False
    ca_dict = {} # auth_seq_id (int) -> (resname, [x, y, z])
    
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
                    if row.get('auth_asym_id') == auth_chain and row.get('auth_atom_id') == 'CA':
                        auth_seq = int(row.get('auth_seq_id'))
                        resname = row.get('label_comp_id')
                        x = float(row.get('Cartn_x'))
                        y = float(row.get('Cartn_y'))
                        z = float(row.get('Cartn_z'))
                        ca_dict[auth_seq] = (resname, np.array([x, y, z], dtype=float))
                elif line.startswith('#') or line.startswith('loop_'):
                    if atom_cols: break
    return ca_dict

def kabsch(p, q):
    centroid_p = np.mean(p, axis=0)
    centroid_q = np.mean(q, axis=0)
    p_c = p - centroid_p
    q_c = q - centroid_q
    h = np.dot(p_c.T, q_c)
    u, s, vt = np.linalg.svd(h)
    d = np.linalg.det(np.dot(vt.T, u.T))
    e = np.eye(3)
    if d < 0: e[2, 2] = -1.0
    r = np.dot(vt.T, np.dot(e, u.T))
    t = centroid_q - np.dot(r, centroid_p)
    p_rot = np.dot(p, r.T) + t
    rmsd = np.sqrt(np.mean(np.sum((p_rot - q)**2, axis=1)))
    return rmsd, r, t, p_rot

ca_7tre = parse_ca_atoms('7TRE', 'A')
ca_7trd = parse_ca_atoms('7TRD', 'A')

common_res = sorted(set(ca_7tre.keys()) & set(ca_7trd.keys()))
print(f"Common CA residues between 7TRE and 7TRD: {len(common_res)}")

coords_tre = np.array([ca_7tre[r][1] for r in common_res])
coords_trd = np.array([ca_7trd[r][1] for r in common_res])

# 1. Global Superposition
rmsd_global, r_glob, t_glob, tre_aligned_glob = kabsch(coords_tre, coords_trd)
print(f"Global CA RMSD (all {len(common_res)} residues): {rmsd_global:.3f} A")

# 2. Catalytic Core Superposition (res >= 326)
core_idx = [i for i, r in enumerate(common_res) if r >= 326]
coords_tre_core = coords_tre[core_idx]
coords_trd_core = coords_trd[core_idx]
rmsd_core, r_core, t_core, _ = kabsch(coords_tre_core, coords_trd_core)
print(f"Catalytic Core CA RMSD (res >= 326, n={len(core_idx)}): {rmsd_core:.3f} A")

# Transform all 7TRE coordinates using catalytic core rotation/translation
tre_aligned_by_core = np.dot(coords_tre, r_core.T) + t_core

# Evaluate TEN domain RMSD when aligned by catalytic core (rigid body shift!)
ten_idx = [i for i, r in enumerate(common_res) if r < 200]
coords_trd_ten = coords_trd[ten_idx]
tre_ten_on_core = tre_aligned_by_core[ten_idx]
disp_ten_on_core = np.sqrt(np.mean(np.sum((tre_ten_on_core - coords_trd_ten)**2, axis=1)))
print(f"TEN domain displacement when aligned on Catalytic Core (n={len(ten_idx)}): {disp_ten_on_core:.3f} A")

# TEN domain internal RMSD (superimposed on TEN alone)
rmsd_ten_internal, _, _, _ = kabsch(coords_tre[ten_idx], coords_trd[ten_idx])
print(f"TEN domain internal CA RMSD (superimposed on TEN alone): {rmsd_ten_internal:.3f} A")
