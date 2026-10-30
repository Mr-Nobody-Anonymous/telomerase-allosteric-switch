import json

with open(r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\survey_raw.json', 'r') as f:
    data = json.load(f)

for pdb_id, val in data.items():
    print(f"=== {pdb_id} ===")
    entry = val.get("entry", {})
    title = entry.get("struct", {}).get("title", "")
    print(f"Title: {title}")
    entities = val.get("entities", {})
    if isinstance(entities, dict):
        entities = entities.values()
    for ent in entities:
        ent_id = ent.get("rcsb_id")
        desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        seq = ent.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        asym_ids = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("asym_ids", [])
        auth_asym_ids = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])
        print(f"  Entity {ent_id}: {desc}")
        print(f"    label_asym: {asym_ids}, auth_asym: {auth_asym_ids}, seq_len: {len(seq)}")
    print()
