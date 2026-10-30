import os

def parse_simple_cif_loop(lines, loop_prefix):
    # Find loop_
    in_loop = False
    collecting_tags = False
    tags = []
    rows = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'loop_':
            # check if following tags start with loop_prefix
            j = i + 1
            cur_tags = []
            while j < len(lines) and lines[j].strip().startswith('_'):
                cur_tags.append(lines[j].strip())
                j += 1
            if any(t.startswith(loop_prefix) for t in cur_tags):
                tags = cur_tags
                # now parse data rows until next loop_ or _tag or #
                data_tokens = []
                while j < len(lines):
                    dline = lines[j].strip()
                    if dline == 'loop_' or (dline.startswith('_') and not dline.startswith('_'+loop_prefix.split('.')[0])) or dline.startswith('#'):
                        break
                    # parse tokens handling quotes and semicolons
                    if dline.startswith(';'):
                        # multiline text
                        text = [dline[1:]]
                        j += 1
                        while j < len(lines) and not lines[j].strip().startswith(';'):
                            text.append(lines[j])
                            j += 1
                        data_tokens.append('\n'.join(text))
                    else:
                        # tokenize line
                        # regex for quoted strings or whitespace
                        import shlex
                        try:
                            tokens = shlex.split(dline)
                            data_tokens.extend(tokens)
                        except:
                            data_tokens.extend(dline.split())
                    j += 1
                
                # group tokens by len(tags)
                ncols = len(tags)
                for k in range(0, len(data_tokens), ncols):
                    chunk = data_tokens[k:k+ncols]
                    if len(chunk) == ncols:
                        rows.append(dict(zip(tags, chunk)))
                break
            i = j - 1
        i += 1
    return tags, rows

# Let's inspect 7TRE TPP1 sequence scheme
with open('cifs/7TRE.cif', 'r', encoding='utf-8') as f:
    lines = f.readlines()

tags, rows = parse_simple_cif_loop(lines, '_pdbx_poly_seq_scheme')
print('7TRE _pdbx_poly_seq_scheme total rows:', len(rows))
tpp1_rows = [r for r in rows if r.get('_pdbx_poly_seq_scheme.asym_id') == 'D']
print('7TRE TPP1 (asym_id D) rows:', len(tpp1_rows))

for r in tpp1_rows:
    seq_id = r.get('_pdbx_poly_seq_scheme.seq_id')
    mon_id = r.get('_pdbx_poly_seq_scheme.mon_id')
    auth_seq = r.get('_pdbx_poly_seq_scheme.auth_seq_num')
    auth_mon = r.get('_pdbx_poly_seq_scheme.auth_mon_id')
    pdb_strand = r.get('_pdbx_poly_seq_scheme.pdb_strand_id')
    pdb_seq = r.get('_pdbx_poly_seq_scheme.pdb_seq_num')
    if auth_seq in ['168', '169', '170', '171', '183', '212'] or seq_id in ['1', '2', '3', '4', '155', '163']:
        print(f"seq_id={seq_id:>3} mon_id={mon_id:>3} | auth_seq_num={auth_seq:>4} auth_mon_id={auth_mon:>3} pdb_strand={pdb_strand} pdb_seq={pdb_seq}")
