import gemmi
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree

# Load 7TRE.cif
cif_path = Path("data/raw/7TRE.cif")
doc = gemmi.cif.read_file(str(cif_path))
st = gemmi.make_structure_from_block(doc.sole_block())

tert_chain = st[0]["A"]
tpp1_chain = st[0]["C"]

# Extract all heavy atoms
tert_atoms = []
tert_coords = []
for res in tert_chain:
    for at in res:
        if at.element.name != "H":
            tert_atoms.append({
                "resnum": res.seqid.num,
                "resname": res.name,
                "atom_name": at.name,
            })
            tert_coords.append([at.pos.x, at.pos.y, at.pos.z])

tpp1_atoms = []
tpp1_coords = []
for res in tpp1_chain:
    for at in res:
        if at.element.name != "H":
            tpp1_atoms.append({
                "resnum": res.seqid.num,
                "resname": res.name,
                "atom_name": at.name,
            })
            tpp1_coords.append([at.pos.x, at.pos.y, at.pos.z])

tert_coords = np.array(tert_coords)
tpp1_coords = np.array(tpp1_coords)

tree_tpp1 = cKDTree(tpp1_coords)

# Find all atom-atom pairs <= 4.5 A
pairs_35 = 0
pairs_40 = 0
pairs_45 = 0

pair_min_dist = {}

for i, pt in enumerate(tert_coords):
    matches = tree_tpp1.query_ball_point(pt, r=4.5)
    for j in matches:
        dist = np.linalg.norm(pt - tpp1_coords[j])
        pairs_45 += 1
        if dist <= 4.0:
            pairs_40 += 1
            if dist <= 3.5:
                pairs_35 += 1
        
        res_pair = (tert_atoms[i]["resnum"], tert_atoms[i]["resname"], tpp1_atoms[j]["resnum"], tpp1_atoms[j]["resname"])
        if res_pair not in pair_min_dist or dist < pair_min_dist[res_pair]:
            pair_min_dist[res_pair] = dist

print(f"Total atom contacts: 3.5A={pairs_35}, 4.0A={pairs_40}, 4.5A={pairs_45}")
print(f"Total residue pairs <= 4.5A: {len(pair_min_dist)}")

df_sens = pd.read_csv("results/interface/sensitivity_analysis.csv")
print(f"sensitivity_analysis.csv rows: {len(df_sens)}")

# Check each residue pair in df_sens against independently computed distances
mismatches = []
for _, row in df_sens.iterrows():
    pair = (int(row["resnum_tert"]), row["resname_tert"], int(row["resnum_tpp1"]), row["resname_tpp1"])
    if pair not in pair_min_dist:
        mismatches.append(f"Pair {pair} in CSV but not found in independent calculation <= 4.5A")
        continue
    indep_min_d = pair_min_dist[pair]
    csv_min_d = float(row["min_distance_angstrom"])
    if abs(indep_min_d - csv_min_d) > 0.01:
        mismatches.append(f"Pair {pair}: CSV min dist={csv_min_d}, Indep={indep_min_d:.3f}")
    
    # Check persistence class
    expected_class = None
    if indep_min_d <= 3.5:
        expected_class = "Core Persistent"
    elif indep_min_d <= 4.0:
        expected_class = "Intermediate Boundary"
    else:
        expected_class = "Peripheral / Sensitive"
    
    if row["persistence_class"] != expected_class:
        mismatches.append(f"Pair {pair}: CSV class='{row['persistence_class']}', expected='{expected_class}'")

print(f"Discrepancies found: {len(mismatches)}")
for m in mismatches:
    print("  ", m)
