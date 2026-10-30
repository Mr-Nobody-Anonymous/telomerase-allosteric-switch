import gemmi

doc = gemmi.cif.read_file('C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw/7TRE.cif')
block = doc.sole_block()

print("=== CHECKING AUTH_ASYM_ID vs LABEL_ASYM_ID IN 7TRE ===")
seen = set()
for row in block.find('_atom_site.', ['label_asym_id', 'auth_asym_id', 'label_entity_id']):
    key = (row[0], row[1], row[2])
    if key not in seen:
        seen.add(key)
        print(f"label_asym_id: {row[0]} | auth_asym_id: {row[1]} | entity_id: {row[2]}")

st = gemmi.make_structure_from_block(block)
for chain in st[0]:
    print(f"Chain {chain.name}: {len(chain)} residues. First: {chain[0].name} {chain[0].seqid.num}, Last: {chain[-1].name} {chain[-1].seqid.num}")
