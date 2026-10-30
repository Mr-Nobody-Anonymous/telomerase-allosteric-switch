import gemmi

cif_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw"

st_e = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/7TRE.cif").sole_block())
st_d = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/7TRD.cif").sole_block())
st_f = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/7TRF.cif").sole_block())

tert_e = st_e[0]["A"]
rna_e = st_e[0]["B"]
tpp1_e = st_e[0]["C"]
dna_e = st_e[0]["D"]

print("=== 7TRE ATOM-PRECISE DISTANCES ===")
print(f"TERT residues: {len(tert_e)}")
print(f"RNA residues:  {len(rna_e)} (from {rna_e[0].name} {rna_e[0].seqid.num} to {rna_e[-1].name} {rna_e[-1].seqid.num})")
print(f"TPP1 residues: {len(tpp1_e)} (from {tpp1_e[0].name} {tpp1_e[0].seqid.num} to {tpp1_e[-1].name} {tpp1_e[-1].seqid.num})")
print(f"DNA residues:  {len(dna_e)} (from {dna_e[0].name} {dna_e[0].seqid.num} to {dna_e[-1].name} {dna_e[-1].seqid.num})")

key_res = [176, 756, 757, 759, 769, 774, 792, 794]

for resnum in key_res:
    r = next((r for r in tert_e if r.seqid.num == resnum), None)
    if not r:
        continue
    
    # Min distance to TPP1 (Chain C)
    min_tpp1 = 999.0
    tpp1_atom_pair = None
    for a1 in r:
        for tr in tpp1_e:
            for a2 in tr:
                d = a1.pos.dist(a2.pos)
                if d < min_tpp1:
                    min_tpp1 = d
                    tpp1_atom_pair = (a1.name, tr.name, tr.seqid.num, a2.name)
                    
    # Min distance to DNA (Chain D)
    min_dna = 999.0
    dna_atom_pair = None
    for a1 in r:
        for dr in dna_e:
            for a2 in dr:
                d = a1.pos.dist(a2.pos)
                if d < min_dna:
                    min_dna = d
                    dna_atom_pair = (a1.name, dr.name, dr.seqid.num, a2.name)
                    
    # Min distance to RNA (Chain B)
    min_rna = 999.0
    rna_atom_pair = None
    for a1 in r:
        for rr in rna_e:
            for a2 in rr:
                d = a1.pos.dist(a2.pos)
                if d < min_rna:
                    min_rna = d
                    rna_atom_pair = (a1.name, rr.name, rr.seqid.num, a2.name)
                    
    print(f"\nTERT {r.name:3s} {resnum:3d}:")
    print(f"  -> TPP1 (Chain C): {min_tpp1:6.3f} A  [{r.name}{resnum}:{tpp1_atom_pair[0]} <-> {tpp1_atom_pair[1]}{tpp1_atom_pair[2]}:{tpp1_atom_pair[3]}]")
    print(f"  -> DNA  (Chain D): {min_dna:6.3f} A  [{r.name}{resnum}:{dna_atom_pair[0]} <-> {dna_atom_pair[1]}{dna_atom_pair[2]}:{dna_atom_pair[3]}]")
    print(f"  -> RNA  (Chain B): {min_rna:6.3f} A  [{r.name}{resnum}:{rna_atom_pair[0]} <-> {rna_atom_pair[1]}{rna_atom_pair[2]}:{rna_atom_pair[3]}]")
