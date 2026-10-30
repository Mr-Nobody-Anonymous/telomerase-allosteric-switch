import gemmi
from pathlib import Path

root = Path(r"C:\Users\hp\.gemini\antigravity\scratch\telomerase-research")
cif_7tre = root / "data" / "raw" / "7TRE.cif"
cif_2i46 = root / "data" / "raw" / "2I46.cif"

st_7tre = gemmi.read_structure(str(cif_7tre))
model_7tre = st_7tre[0]

print("=== 7TRE Chain C (TPP1) Residues 160-175 ===")
for ch in model_7tre:
    if ch.name == "C": # author chain C
        for r in ch:
            if 160 <= r.seqid.num <= 175:
                print(f"auth_seq_id: {r.seqid.num}, label_seq: {r.label_seq}, name: {r.name}, subchain: {r.subchain}")

st_2i46 = gemmi.read_structure(str(cif_2i46))
model_2i46 = st_2i46[0]

print("\n=== 2I46 Chain A Residues 160-175 ===")
for ch in model_2i46:
    if ch.name == "A":
        for r in ch:
            if 160 <= r.seqid.num <= 175:
                print(f"auth_seq_id: {r.seqid.num}, label_seq: {r.label_seq}, name: {r.name}, subchain: {r.subchain}")

print("\n=== 7TRE Entity 4 Sequence ===")
doc = gemmi.cif.read_file(str(cif_7tre))
block = doc.sole_block()
poly = block.find("_entity_poly_seq.", ["entity_id", "num", "mon_id"])
for row in poly:
    if row[0] == "4" and 70 <= int(row[1]) <= 95:
        print(f"entity 4 num {row[1]}: mon_id {row[2]}")
