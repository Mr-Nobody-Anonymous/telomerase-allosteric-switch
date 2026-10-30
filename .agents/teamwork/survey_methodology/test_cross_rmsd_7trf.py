import urllib.request
import numpy as np
from test_cross_rmsd import parse_ca_atoms, kabsch

ca_7tre = parse_ca_atoms('7TRE', 'A')
ca_7trf = parse_ca_atoms('7TRF', 'A')

common_res = sorted(set(ca_7tre.keys()) & set(ca_7trf.keys()))
print(f"Common CA residues between 7TRE and 7TRF: {len(common_res)}")

coords_tre = np.array([ca_7tre[r][1] for r in common_res])
coords_trf = np.array([ca_7trf[r][1] for r in common_res])

# 1. Global Superposition
rmsd_global, r_glob, t_glob, tre_aligned_glob = kabsch(coords_tre, coords_trf)
print(f"Global CA RMSD (all {len(common_res)} residues): {rmsd_global:.3f} A")

# 2. Catalytic Core Superposition (res >= 326)
core_idx = [i for i, r in enumerate(common_res) if r >= 326]
coords_tre_core = coords_tre[core_idx]
coords_trf_core = coords_trf[core_idx]
rmsd_core, r_core, t_core, _ = kabsch(coords_tre_core, coords_trf_core)
print(f"Catalytic Core CA RMSD (res >= 326, n={len(core_idx)}): {rmsd_core:.3f} A")

# Transform all 7TRE coordinates using catalytic core rotation/translation
tre_aligned_by_core = np.dot(coords_tre, r_core.T) + t_core

# Evaluate TEN domain RMSD when aligned by catalytic core (rigid body shift!)
ten_idx = [i for i, r in enumerate(common_res) if r < 200]
coords_trf_ten = coords_trf[ten_idx]
tre_ten_on_core = tre_aligned_by_core[ten_idx]
disp_ten_on_core = np.sqrt(np.mean(np.sum((tre_ten_on_core - coords_trf_ten)**2, axis=1)))
print(f"TEN domain displacement when aligned on Catalytic Core (n={len(ten_idx)}): {disp_ten_on_core:.3f} A")

# TEN domain internal RMSD (superimposed on TEN alone)
rmsd_ten_internal, _, _, _ = kabsch(coords_tre[ten_idx], coords_trf[ten_idx])
print(f"TEN domain internal CA RMSD (superimposed on TEN alone): {rmsd_ten_internal:.3f} A")
