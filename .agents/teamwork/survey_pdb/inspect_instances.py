import requests
import json

targets = [
    ('7TRE', 'D', 'TPP1 in 7TRE'),
    ('7TRE', 'B', 'TERT in 7TRE'),
    ('7TRD', 'B', 'TERT in 7TRD'),
    ('7TRF', 'D', 'TERT in 7TRF'),
    ('2I46', 'A', 'TPP1 chain A in 2I46'),
    ('2I46', 'B', 'TPP1 chain B in 2I46')
]

for pid, asym, label in targets:
    url = f'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pid}/{asym}'
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        print('=============================================')
        auth_id = data.get('rcsb_polymer_entity_instance_container_identifiers', {}).get('auth_asym_id')
        print(f"{pid} asym_id='{asym}' ({label}): auth_asym_id='{auth_id}'")
        unmodeled = data.get('rcsb_unmodeled_polymer_instance_feature', [])
        print(f"Unmodeled features count: {len(unmodeled)}")
        for u in unmodeled:
            print(f"  unmodeled: type={u.get('type')} name={u.get('name')} positions={u.get('feature_positions')}")
