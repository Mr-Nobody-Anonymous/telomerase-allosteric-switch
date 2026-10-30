import gemmi

cif_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw"
out_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data"

def extract_subsystem(pdb_id, include_tpp1=False):
    st_src = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/{pdb_id}.cif").sole_block())
    
    st_new = gemmi.Structure()
    model = gemmi.Model("1")
    # Chain A: TEN (165-190) and TRAP (750-795)
    chain_a = gemmi.Chain("A")
    for r in st_src[0]["A"]:
        num = r.seqid.num
        if (165 <= num <= 190) or (750 <= num <= 795):
            chain_a.add_residue(r.clone())
    model.add_chain(chain_a)
    
    # Chain C: TPP1 (165-185, 210-225)
    if include_tpp1 and "C" in [c.name for c in st_src[0]]:
        chain_c = gemmi.Chain("C")
        for r in st_src[0]["C"]:
            num = r.seqid.num
            if (165 <= num <= 185) or (210 <= num <= 225):
                chain_c.add_residue(r.clone())
        model.add_chain(chain_c)
        
    st_new.add_model(model)
    return st_new

st_bound = extract_subsystem("7TRE", include_tpp1=True)
st_bound.write_pdb(f"{out_dir}/junction_bound_wt.pdb")

st_apo = extract_subsystem("7TRD", include_tpp1=False)
st_apo.write_pdb(f"{out_dir}/junction_apo_wt.pdb")

print("Successfully wrote junction PDBs:")
print(f"  Bound WT: {len(st_bound[0]['A'])} TERT residues, {len(st_bound[0]['C'])} TPP1 residues")
print(f"  Apo WT:   {len(st_apo[0]['A'])} TERT residues")
