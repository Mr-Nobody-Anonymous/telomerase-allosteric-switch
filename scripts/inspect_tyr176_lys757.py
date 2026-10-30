import gemmi

structures = ["7TRE", "7TRD", "7TRF"]

for pdb in structures:
    st = gemmi.make_structure_from_block(gemmi.cif.read_file(f"data/raw/{pdb}.cif").sole_block())
    tert = st[0]["A"]
    y176 = next(r for r in tert if r.seqid.num == 176)
    k757 = next(r for r in tert if r.seqid.num == 757)
    
    print(f"\n=================== {pdb} ===================")
    ca_y176 = y176["CA"][0]
    ca_k757 = k757["CA"][0]
    ca_dist = ca_y176.pos.dist(ca_k757.pos)
    print(f"Calpha - Calpha distance: {ca_dist:.3f} A")
    print(f"Tyr176 Calpha B-factor: {ca_y176.b_iso:.2f} A^2")
    print(f"Lys757 Calpha B-factor: {ca_k757.b_iso:.2f} A^2")
    
    # Tyr176 atoms present
    print(f"Tyr176 atoms: {[a.name for a in y176]}")
    print(f"Lys757 atoms: {[a.name for a in k757]}")
    
    # All pairs < 5.0 A
    pairs = []
    for a1 in y176:
        for a2 in k757:
            d = a1.pos.dist(a2.pos)
            if d < 5.0:
                pairs.append((a1.name, a2.name, d))
    pairs.sort(key=lambda x: x[2])
    print(f"Inter-residue atom contacts < 5.0 A (total {len(pairs)}):")
    for p in pairs:
        print(f"  Tyr176:{p[0]:4s} <---> Lys757:{p[1]:4s} = {p[2]:.3f} A")

    # Now let's examine the surrounding pocket of Tyr176 and Lys757
    print("\nContacts of Tyr176 sidechain (CB, CG, CD1, CD2, CE1, CE2, CZ, OH) with other residues (< 4.0 A):")
    sc_atoms = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
    y_contacts = []
    for a1 in y176:
        if a1.name in sc_atoms:
            for r in tert:
                if r.seqid.num != 176:
                    for a2 in r:
                        d = a1.pos.dist(a2.pos)
                        if d < 4.0:
                            y_contacts.append((a1.name, r.name, r.seqid.num, a2.name, d))
    y_contacts.sort(key=lambda x: x[4])
    for yc in y_contacts[:12]:
        print(f"  Tyr176:{yc[0]} -> {yc[1]} {yc[2]}:{yc[3]} = {yc[4]:.2f} A")
