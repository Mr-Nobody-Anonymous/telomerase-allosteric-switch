import json

with open('survey_raw.json') as f:
    data = json.load(f)

for pid in ['7TRE', '7TRD', '7TRF', '7TRC', '2I46']:
    print('==============================================')
    print(f'PDB ID: {pid}')
    e = data[pid]['entry']
    pub = e.get('citation', [{}])[0]
    print(f"Title: {e.get('struct', {}).get('title')}")
    print(f"Citation: {pub.get('title')} | Journal: {pub.get('journal_abbrev')} ({pub.get('year')}) | DOI: {pub.get('pdbx_database_id_doi')}")
    print(f"Method: {e.get('exptl', [{}])[0].get('method')}")
    print(f"Resolution: {e.get('rcsb_entry_info', {}).get('resolution_combined')}")
    print(f"EMDB: {e.get('rcsb_entry_container_identifiers', {}).get('emdb_ids')}")
    print('Polymer Entities:')
    for eid, ent in data[pid]['entities'].items():
        desc = ent.get('rcsb_polymer_entity', {}).get('pdbx_description')
        poly_type = ent.get('entity_poly', {}).get('type')
        c_ids = ent.get('rcsb_polymer_entity_container_identifiers', {})
        asym_ids = c_ids.get('asym_ids', [])
        auth_asym_ids = c_ids.get('auth_asym_ids', [])
        uniprot = [x.get('database_accession') for x in c_ids.get('reference_sequence_identifiers', []) if x.get('database_name') == 'UniProt']
        genbank = [x.get('database_accession') for x in c_ids.get('reference_sequence_identifiers', []) if x.get('database_name') == 'GenBank']
        seq_len = ent.get('entity_poly', {}).get('rcsb_sample_sequence_length')
        print(f"  Entity {eid}:")
        print(f"    Description: {desc}")
        print(f"    Type: {poly_type}")
        print(f"    Length: {seq_len}")
        print(f"    label_asym_id (CIF chain): {asym_ids}")
        print(f"    auth_asym_id (PDB chain): {auth_asym_ids}")
        print(f"    UniProt: {uniprot} | GenBank: {genbank}")
