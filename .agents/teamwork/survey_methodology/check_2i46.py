import json
import urllib.request

# Check 2I46 mmCIF or PDB sequence numbering in survey_raw.json
with open(r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\survey_raw.json', 'r') as f:
    data = json.load(f)

print("2I46 entities:")
for eid, ent in data['2I46']['entities'].items():
    print("Entity:", eid)
    c_ids = ent.get('rcsb_polymer_entity_container_identifiers', {})
    print("  ref_seq_ids:", c_ids.get('reference_sequence_identifiers'))

# Check literature: Wang et al 2007 (Structure of TPP1 OB fold, PDB 2I46)
# What was the residue numbering in 2I46?
# Let's fetch the first few lines of 2I46 cif from RCSB or inspect RCSB entry
url_2i46 = "https://data.rcsb.org/rest/v1/core/polymer_entity_instance/2I46/A"
req = urllib.request.Request(url_2i46, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as resp:
        inst_data = json.load(resp)
        print("2I46 instance A keys:", inst_data.keys())
        rcsb_inst = inst_data.get('rcsb_polymer_entity_instance_container_identifiers', {})
        print("  auth_asym_id:", rcsb_inst.get('auth_asym_id'))
        print("  asym_id:", rcsb_inst.get('asym_id'))
except Exception as e:
    print("Error fetching instance:", e)
