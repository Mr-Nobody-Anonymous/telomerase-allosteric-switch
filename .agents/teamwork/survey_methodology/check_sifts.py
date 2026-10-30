import urllib.request
import json

# Check 2I46 polymer_entity 1 SIFTS / residue mapping
url = "https://data.rcsb.org/rest/v1/core/polymer_entity/2I46/1"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    pe_2i46 = json.load(resp)

print("2I46 entity alignment:")
for aln in pe_2i46.get('rcsb_polymer_entity_align', []):
    print("  DB:", aln.get('reference_database_name'), aln.get('reference_database_accession'))
    for reg in aln.get('aligned_regions', []):
        print(f"    entity_beg={reg.get('entity_beg_seq_id')}, ref_beg={reg.get('ref_beg_seq_id')}, len={reg.get('length')}")

# Also check 7TRE polymer_entity 4 (TPP1)
url_7tre = "https://data.rcsb.org/rest/v1/core/polymer_entity/7TRE/4"
req = urllib.request.Request(url_7tre, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as resp:
    pe_7tre = json.load(resp)

print("\n7TRE entity 4 alignment:")
for aln in pe_7tre.get('rcsb_polymer_entity_align', []):
    print("  DB:", aln.get('reference_database_name'), aln.get('reference_database_accession'))
    for reg in aln.get('aligned_regions', []):
        print(f"    entity_beg={reg.get('entity_beg_seq_id')}, ref_beg={reg.get('ref_beg_seq_id')}, len={reg.get('length')}")
