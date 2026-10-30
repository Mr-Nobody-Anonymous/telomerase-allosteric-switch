import os
from parse_2i46_seq import parse_simple_cif_loop

tert_entries = [
    ('7TRE', 'B', 'A'),
    ('7TRD', 'B', 'A'),
    ('7TRF', 'D', 'A')
]

for pid, asym, auth_asym in tert_entries:
    with open(f'cifs/{pid}.cif', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    tags, rows = parse_simple_cif_loop(lines, '_pdbx_poly_seq_scheme')
    tert_rows = [r for r in rows if r.get('_pdbx_poly_seq_scheme.asym_id') == asym]
    print(f'=== {pid} TERT (asym_id={asym}, auth_asym_id={auth_asym}) ===')
    print('Total sequence length in scheme:', len(tert_rows))
    
    # Modeled vs unmodeled
    modeled = [r for r in tert_rows if r.get('_pdbx_poly_seq_scheme.auth_seq_num') != '?']
    unmodeled = [r for r in tert_rows if r.get('_pdbx_poly_seq_scheme.auth_seq_num') == '?']
    print(f'Modeled residues: {len(modeled)}, Unmodeled: {len(unmodeled)}')
    
    # Check first and last modeled
    if modeled:
        print(f'First modeled: seq_id={modeled[0].get("_pdbx_poly_seq_scheme.seq_id")}, auth_seq_num={modeled[0].get("_pdbx_poly_seq_scheme.auth_seq_num")}, mon={modeled[0].get("_pdbx_poly_seq_scheme.mon_id")}')
        print(f'Last modeled: seq_id={modeled[-1].get("_pdbx_poly_seq_scheme.seq_id")}, auth_seq_num={modeled[-1].get("_pdbx_poly_seq_scheme.auth_seq_num")}, mon={modeled[-1].get("_pdbx_poly_seq_scheme.mon_id")}')
    
    # Check sequence correspondence: is seq_id == auth_seq_num?
    diffs = [r for r in modeled if r.get('_pdbx_poly_seq_scheme.seq_id') != r.get('_pdbx_poly_seq_scheme.auth_seq_num')]
    print(f'Residues where seq_id != auth_seq_num: {len(diffs)}')

    # Find unmodeled segments (loops/gaps)
    gaps = []
    in_gap = False
    gap_start = None
    for r in tert_rows:
        is_unobs = (r.get('_pdbx_poly_seq_scheme.auth_seq_num') == '?')
        seq_id = int(r.get('_pdbx_poly_seq_scheme.seq_id'))
        if is_unobs and not in_gap:
            in_gap = True
            gap_start = seq_id
        elif not is_unobs and in_gap:
            in_gap = False
            gaps.append((gap_start, seq_id - 1))
    if in_gap:
        gaps.append((gap_start, len(tert_rows)))
    print('Unmodeled gaps (seq_id ranges):', gaps)
