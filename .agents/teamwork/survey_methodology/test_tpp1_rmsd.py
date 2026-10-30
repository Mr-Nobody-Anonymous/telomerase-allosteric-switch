import urllib.request
import numpy as np
from test_cross_rmsd import parse_ca_atoms, kabsch

ca_7tre_tpp1 = parse_ca_atoms('7TRE', 'C')
ca_2i46_tpp1 = parse_ca_atoms('2I46', 'A')

print(f"7TRE TPP1 CA count: {len(ca_7tre_tpp1)}, 2I46 TPP1 CA count: {len(ca_2i46_tpp1)}")

common_res = sorted(set(ca_7tre_tpp1.keys()) & set(ca_2i46_tpp1.keys()))
print(f"Common CA residues between 7TRE and 2I46 TPP1: {len(common_res)}")

coords_tre = np.array([ca_7tre_tpp1[r][1] for r in common_res])
coords_2i46 = np.array([ca_2i46_tpp1[r][1] for r in common_res])

rmsd_global, r, t, tre_aligned = kabsch(coords_tre, coords_2i46)
print(f"Global CA RMSD between bound 7TRE and unbound 2I46 TPP1: {rmsd_global:.3f} A")

# Check displacement at TEL patch residues specifically
tel_patch = [168, 169, 170, 171, 180, 183, 212]
print("\nDisplacement at TEL-patch residues after global superposition:")
for res in tel_patch:
    if res in common_res:
        idx = common_res.index(res)
        dist = np.linalg.norm(tre_aligned[idx] - coords_2i46[idx])
        print(f"  Residue {ca_7tre_tpp1[res][0]}{res}: {dist:.3f} A")
