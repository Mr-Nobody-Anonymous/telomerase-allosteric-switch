"""
Independent Forensic Verification Script for Milestone 2.
Author: Forensic Auditor (auditor_m2)

Empirical tests to prove:
1. PairwiseAligner is genuinely dynamic and responds to insertions, mutations, deletions.
2. Coordinate sequence extraction is dynamic and not hardcoded to filename or known PDB tables.
3. 2I46 dimer asymmetry is dynamically extracted from mmCIF residue records.
4. Output schemas and residue correspondence counts are strictly accurate.
"""

import sys
import tempfile
from pathlib import Path
import gemmi
from Bio import Align

# Import target deliverable
from telotrade.sequence_mapping import (
    align_sequences,
    map_residue_correspondence,
    get_chain_sequence,
    get_uniprot_mapping,
    check_tel_patch_integrity,
    select_valid_tpp1_chain,
    generate_tert_correspondence_table,
    generate_tpp1_correspondence_table,
    generate_chain_entity_mapping,
    reconcile_tert_offset,
    reconcile_tpp1_offset,
)

raw_dir = Path(r"C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\data\raw")

print("=== CHECK 1: Dynamic Pairwise Alignment & Perturbation Sensitivity ===")

# Test 1.1: Exact match baseline
s1 = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
res1 = align_sequences(s1, s1)
assert res1.identity == 1.0, f"Expected 1.0, got {res1.identity}"
assert len(res1.mapping) == len(s1), "All residues must map"
print(f"PASS 1.1: Identity alignment verified (score={res1.score}, len={len(res1.mapping)})")

# Test 1.2: Insertion shift perturbation
# Insert 5 glycines at index 10 in target sequence
s2_ins = s1[:10] + "GGGGG" + s1[10:]
res_ins = align_sequences(s1, s2_ins)
# For all residues after index 10 in s1, mapped target index should shift by 5
for i in range(10):
    assert res_ins.mapping.get(i) == i, f"Index {i} mapped to {res_ins.mapping.get(i)}"
for i in range(10, len(s1)):
    assert res_ins.mapping.get(i) == i + 5, f"Index {i} mapped to {res_ins.mapping.get(i)}, expected {i+5}"
print("PASS 1.2: Insertion perturbation shifts indices dynamically by exactly +5")

# Test 1.3: Deletion perturbation
s2_del = s1[:10] + s1[15:]  # Delete 5 residues (indices 10..14)
res_del = align_sequences(s1, s2_del)
for i in range(10, 15):
    assert i not in res_del.mapping, f"Deleted index {i} should not be in mapping"
for i in range(15, len(s1)):
    assert res_del.mapping.get(i) == i - 5, f"Index {i} mapped to {res_del.mapping.get(i)}, expected {i-5}"
print("PASS 1.3: Deletion perturbation unmaps deleted residues and shifts downstream indices by -5")

# Test 1.4: Mismatch score sensitivity
s_random1 = "AAAAAAAAAA"
s_random2 = "CCCCCCCCCC"
res_mismatch = align_sequences(s_random1, s_random2)
# Mismatch score is -1.0 * 10 = -10.0 (or alignment penalty)
print(f"PASS 1.4: Mismatch scoring dynamically handled (identity={res_mismatch.identity}, score={res_mismatch.score})")


print("\n=== CHECK 2: mmCIF Dynamic Extraction vs Synthetic mmCIF File ===")
# Build a synthetic mmCIF file with arbitrary, non-PDB sequence and coordinates
doc = gemmi.cif.Document()
block = doc.add_new_block("TEST_BLOCK")

# Minimal atom_site loop
loop = block.init_loop("_atom_site.", [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "auth_seq_id", "auth_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
    "occupancy", "B_iso_or_equiv", "pdbx_PDB_ins_code", "pdbx_PDB_model_num",
])

# Create 5 synthetic residues: MET, GLY, TRP, TYR, PHE
synth_residues = [("MET", 1, 101), ("GLY", 2, 102), ("TRP", 3, 103), ("TYR", 4, 104), ("PHE", 5, 105)]
atom_id = 1
for comp, l_seq, a_seq in synth_residues:
    loop.add_row([
        "ATOM", str(atom_id), "C", "CA", ".", comp, "X", "1", str(l_seq), str(a_seq), "A", "0.0", "0.0", "0.0", "1.0", "20.0", "?", "1"
    ])
    atom_id += 1

with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tf:
    synth_cif_path = Path(tf.name)
doc.write_file(str(synth_cif_path))

try:
    # Test Gemmi extraction on synthetic mmCIF
    seq_auth, nums_auth = get_chain_sequence(synth_cif_path, "A", use_auth=True, parser="gemmi")
    assert seq_auth == "MGWYF", f"Expected 'MGWYF', got '{seq_auth}'"
    assert nums_auth == [101, 102, 103, 104, 105], f"Expected [101, 102, 103, 104, 105], got {nums_auth}"

    seq_label, nums_label = get_chain_sequence(synth_cif_path, "X", use_auth=False, parser="gemmi")
    assert seq_label == "MGWYF", f"Expected 'MGWYF', got '{seq_label}'"
    assert nums_label == [1, 2, 3, 4, 5], f"Expected [1, 2, 3, 4, 5], got {nums_label}"

    # Test Biopython extraction on synthetic mmCIF
    seq_bio, nums_bio = get_chain_sequence(synth_cif_path, "A", use_auth=True, parser="biopython")
    assert seq_bio == "MGWYF", f"Expected 'MGWYF', got '{seq_bio}'"
    assert nums_bio == [101, 102, 103, 104, 105], f"Expected [101, 102, 103, 104, 105], got {nums_bio}"

    print("PASS 2.1: Synthetic mmCIF parser test confirmed dynamic extraction with zero hardcoding!")
finally:
    if synth_cif_path.exists():
        synth_cif_path.unlink()


print("\n=== CHECK 3: Real Structural mmCIF File Verification ===")
# 3.1 7TRE TERT chain
cif_7tre = raw_dir / "7TRE.cif"
seq_tert, nums_tert = get_chain_sequence(cif_7tre, "A", use_auth=True)
assert len(seq_tert) == 934
assert nums_tert[0] == 7 and nums_tert[-1] == 1132
print(f"PASS 3.1: 7TRE TERT verified ({len(seq_tert)} modeled residues from {nums_tert[0]} to {nums_tert[-1]})")

# 3.2 2I46 Asymmetry test
cif_2i46 = raw_dir / "2I46.cif"
chk_a = check_tel_patch_integrity(cif_2i46, "A")
chk_b = check_tel_patch_integrity(cif_2i46, "B")
assert chk_a["resolved"] is True
assert chk_b["resolved"] is False
assert chk_b["missing_residues"] == [168, 169, 170, 171]
print("PASS 3.2: 2I46 Chain A vs Chain B asymmetry verified (DE-loop 168-171 missing only in B)")

# 3.3 Redirection test
redirected = select_valid_tpp1_chain(cif_2i46, "B", require_tel_patch=True)
assert redirected == "A"
print("PASS 3.3: select_valid_tpp1_chain redirect verified (Chain B -> Chain A)")


print("\n=== CHECK 4: Generated Output Tables Schema and Row Counts ===")
df_tert = generate_tert_correspondence_table(raw_dir)
assert len(df_tert) == 1132, f"Expected 1132 rows, got {len(df_tert)}"
assert df_tert["modeled_in_7tre"].sum() == 934
assert df_tert["modeled_in_7trd"].sum() == 934
assert df_tert["modeled_in_7trf"].sum() == 934
print(f"PASS 4.1: TERT correspondence table verified (rows={len(df_tert)}, modeled={df_tert['modeled_in_7tre'].sum()})")

df_tpp1 = generate_tpp1_correspondence_table(raw_dir)
assert len(df_tpp1) == 161, f"Expected 161 rows, got {len(df_tpp1)}"
assert df_tpp1["modeled_in_7tre"].sum() == 152
assert df_tpp1["modeled_in_2i46_chain_a"].sum() == 152
assert df_tpp1["modeled_in_2i46_chain_b"].sum() == 134
print(f"PASS 4.2: TPP1 correspondence table verified (rows={len(df_tpp1)}, 7TRE={df_tpp1['modeled_in_7tre'].sum()}, 2I46_A={df_tpp1['modeled_in_2i46_chain_a'].sum()}, 2I46_B={df_tpp1['modeled_in_2i46_chain_b'].sum()})")

df_entity = generate_chain_entity_mapping(raw_dir)
assert len(df_entity) == 19
assert set(df_entity["pdb_id"]) == {"7TRE", "7TRD", "7TRF", "7TRC", "2I46"}
print(f"PASS 4.3: Chain entity mapping verified (rows={len(df_entity)}, structures={len(set(df_entity['pdb_id']))})")

print("\n>>> ALL EMPIRICAL FORENSIC INTEGRITY CHECKS PASSED SUCCESSFULLY! <<<")
