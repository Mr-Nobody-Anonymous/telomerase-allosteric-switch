from pathlib import Path

import gemmi


def main():
    cif_path = Path("data/raw/2I46.cif")
    st = gemmi.read_structure(str(cif_path))
    st.setup_cell_images()
    cell = st.cell
    sg = st.find_spacegroup()
    print("Spacegroup:", sg.xhm())

    model = st[0]

    # Populate neighbor search
    ns = gemmi.NeighborSearch(model, cell, 5.0)
    ns.populate(include_h=False)

    # 1. Search inter-chain / symmetry contacts
    cs = gemmi.ContactSearch(5.0)
    # Ignore contacts within the same residue or same polymer if we want, or test everything
    # Let's inspect all contacts
    all_contacts = cs.find_contacts(ns)
    print(f"Total contacts <= 5.0 A across crystal: {len(all_contacts)}")

    target_resnums = set(range(165, 176)) # 165 to 175

    # Filter contacts where one partner is in Chain A residues 165-175
    de_contacts = []
    for c in all_contacts:
        p1 = c.partner1
        p2 = c.partner2
        
        # Check if p1 is target
        p1_is_target = (p1.chain.name == "A" and p1.residue.seqid.num in target_resnums)
        p2_is_target = (p2.chain.name == "A" and p2.residue.seqid.num in target_resnums)

        if not (p1_is_target or p2_is_target):
            continue

        # If both are in Chain A same ASU, check residue distance (intra-molecular)
        if p1.chain.name == "A" and p2.chain.name == "A" and c.image_idx == 0:
            # Intra-molecular contact within Chain A
            continue

        # Target is p1, partner is p2 (or vice versa)
        if p1_is_target:
            t = p1
            o = p2
        else:
            t = p2
            o = p1

        # Ignore water/solvent for lattice packing analysis unless noted
        is_water = o.residue.name in ["HOH", "SO4", "GOL"]

        de_contacts.append({
            "target_resnum": t.residue.seqid.num,
            "target_resname": t.residue.name,
            "target_atom": t.atom.name,
            "partner_chain": o.chain.name,
            "partner_resnum": o.residue.seqid.num,
            "partner_resname": o.residue.name,
            "partner_atom": o.atom.name,
            "dist": c.dist,
            "image_idx": c.image_idx,
            "is_water": is_water
        })

    print("\nContacts involving Chain A residues 165-175 with other molecules / symmetry mates:")
    print(f"Total such contacts: {len(de_contacts)} (non-water: {len([x for x in de_contacts if not x['is_water']] )})")

    de_contacts.sort(key=lambda x: x["dist"])
    
    # Non-water protein-protein lattice contacts
    protein_contacts = [x for x in de_contacts if not x["is_water"]]
    print("\n--- Protein-Protein Lattice Contacts (Non-Water) ---")
    if protein_contacts:
        for c in protein_contacts:
            print(f"Res {c['target_resnum']} {c['target_resname']}:{c['target_atom']} -> Chain {c['partner_chain']} Res {c['partner_resnum']} {c['partner_resname']}:{c['partner_atom']} dist={c['dist']:.3f} A (image {c['image_idx']})")
    else:
        print("NO protein-protein lattice contacts found <= 5.0 A!")

    # Water contacts
    water_contacts = [x for x in de_contacts if x["is_water"]]
    print(f"\n--- Water / Solvent Contacts ({len(water_contacts)} found) ---")
    for c in water_contacts[:15]:
        print(f"Res {c['target_resnum']} {c['target_resname']}:{c['target_atom']} -> {c['partner_resname']} {c['partner_resnum']}:{c['partner_atom']} dist={c['dist']:.3f} A")

    # Specifically for DE-loop residues 168-171:
    print("\n--- Summary specifically for DE-loop residues 168-171 ---")
    for r in [168, 169, 170, 171]:
        r_prot = [x for x in protein_contacts if x["target_resnum"] == r]
        r_wat = [x for x in water_contacts if x["target_resnum"] == r]
        print(f"Residue {r}: {len(r_prot)} protein lattice contacts, {len(r_wat)} solvent contacts.")
        for p in r_prot:
            print(f"   Protein contact: {p['target_atom']} -> Chain {p['partner_chain']} {p['partner_resname']}{p['partner_resnum']}:{p['partner_atom']} dist={p['dist']:.3f} A")

if __name__ == "__main__":
    main()
