from parse_2i46_seq import parse_simple_cif_loop

rna_entries = [
    ('7TRE', 'A', 'B'),
    ('7TRD', 'A', 'B'),
    ('7TRF', 'C', 'B'),
    ('7TRC', 'D', 'B')
]

for pid, asym, auth_asym in rna_entries:
    with open(f'cifs/{pid}.cif', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    tags, rows = parse_simple_cif_loop(lines, '_pdbx_poly_seq_scheme')
    rna_rows = [r for r in rows if r.get('_pdbx_poly_seq_scheme.asym_id') == asym]
    modeled = [r for r in rna_rows if r.get('_pdbx_poly_seq_scheme.auth_seq_num') != '?']
    print(f"=== {pid} hTR RNA (asym_id={asym}, auth_asym_id={auth_asym}) ===")
    print(f"Total seq length: {len(rna_rows)}, Modeled nt count: {len(modeled)}")
    
    # find modeled continuous segments
    segments = []
    in_seg = False
    seg_start = None
    for r in rna_rows:
        is_mod = (r.get('_pdbx_poly_seq_scheme.auth_seq_num') != '?')
        auth_seq = r.get('_pdbx_poly_seq_scheme.auth_seq_num')
        if is_mod and not in_seg:
            in_seg = True
            seg_start = auth_seq
        elif not is_mod and in_seg:
            in_seg = False
            segments.append((seg_start, prev_auth))
        if is_mod:
            prev_auth = auth_seq
    if in_seg:
        segments.append((seg_start, prev_auth))
    print("Modeled segments (auth_seq_num ranges):", segments)
