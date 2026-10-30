import json

with open(r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\survey_raw.json', 'r') as f:
    data = json.load(f)

for pdb_id in ['7TRE', '2I46']:
    val = data[pdb_id]
    print(f"=== {pdb_id} ===")
    entities = val.get("entities", {})
    if isinstance(entities, dict):
        for eid, ent in entities.items():
            desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
            if "Adrenocortical" in desc or "TPP1" in desc:
                alignments = ent.get("rcsb_polymer_entity_align", [])
                for aln in alignments:
                    print("  Reference DB:", aln.get("reference_database_name"), aln.get("reference_database_accession"))
                    for reg in aln.get("aligned_regions", []):
                        print("    entity_beg:", reg.get("entity_beg_seq_id"), "ref_beg:", reg.get("ref_beg_seq_id"), "length:", reg.get("length"))
