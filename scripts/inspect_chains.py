import gemmi

doc = gemmi.cif.read_file('C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw/7TRE.cif')
block = doc.sole_block()

print("=== ENTITIES IN 7TRE ===")
entities = {}
for row in block.find('_entity.', ['id', 'type', 'pdbx_description']):
    eid, etype, edesc = row[0], row[1], row[2]
    entities[eid] = (etype, edesc)
    print(f"Entity {eid}: {etype} | {edesc}")

print("\n=== CHAINS / STRUCT_ASYM IN 7TRE ===")
for row in block.find('_struct_asym.', ['id', 'entity_id']):
    cid, eid = row[0], row[1]
    etype, edesc = entities.get(eid, ("unknown", "unknown"))
    print(f"Chain {cid} -> Entity {eid} ({edesc})")
