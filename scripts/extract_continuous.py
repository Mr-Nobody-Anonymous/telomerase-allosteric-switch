import gemmi

cif_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw"
out_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data"

def extract_subsystem_continuous(pdb_id, include_tpp1=False):
    st_src = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/{pdb_id}.cif").sole_block())
    
    st_new = gemmi.Structure()
    model = gemmi.Model("1")
    
    # Chain A: TEN (165-190)
    chain_a = gemmi.Chain("A")
    for r in st_src[0]["A"]:
        if 165 <= r.seqid.num <= 190:
            chain_a.add_residue(r.clone())
    model.add_chain(chain_a)
    
    # Chain B: TRAP (750-795)
    chain_b = gemmi.Chain("B")
    for r in st_src[0]["A"]:
        if 750 <= r.seqid.num <= 795:
            chain_b.add_residue(r.clone())
    model.add_chain(chain_b)
    
    if include_tpp1 and "C" in [c.name for c in st_src[0]]:
        # Chain C: TPP1 (165-185)
        chain_c = gemmi.Chain("C")
        for r in st_src[0]["C"]:
            if 165 <= r.seqid.num <= 185:
                chain_c.add_residue(r.clone())
        model.add_chain(chain_c)
        
        # Chain D: TPP1 (210-225)
        chain_d = gemmi.Chain("D")
        for r in st_src[0]["C"]:
            if 210 <= r.seqid.num <= 225:
                chain_d.add_residue(r.clone())
        model.add_chain(chain_d)
        
    st_new.add_model(model)
    return st_new

st_bound = extract_subsystem_continuous("7TRE", include_tpp1=True)
st_bound.write_pdb(f"{out_dir}/junction_bound_continuous.pdb")

st_apo = extract_subsystem_continuous("7TRD", include_tpp1=False)
st_apo.write_pdb(f"{out_dir}/junction_apo_continuous.pdb")

print("Wrote continuous segment PDBs:")
print(f"  Bound: {[c.name + ': ' + str(len(c)) for c in st_bound[0]]}")
print(f"  Apo:   {[c.name + ': ' + str(len(c)) for c in st_apo[0]]}")
