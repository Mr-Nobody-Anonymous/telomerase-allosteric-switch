import urllib.request
import numpy as np

# Load 7TRE atoms and find all N, O pairs between TERT and TPP1 <= 4.0 A
# And check salt bridge residues:
# Acidic: ASP (OD1, OD2), GLU (OE1, OE2)
# Basic: LYS (NZ), ARG (NH1, NH2, NE), HIS (ND1, NE2)

from Bio.PDB import MMCIFParser

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('7TRE', r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\cifs\7TRE.cif')
model = structure[0]
tert = model['A']
tpp1 = model['C']

polar_types = {'N', 'O', 'S'}
acidic_atoms = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}
basic_atoms = {('LYS', 'NZ'), ('ARG', 'NH1'), ('ARG', 'NH2'), ('ARG', 'NE'), ('HIS', 'ND1'), ('HIS', 'NE2')}

polar_contacts = []
salt_bridges = []

for atom_t in tert.get_atoms():
    el_t = atom_t.element
    if el_t not in polar_types: continue
    res_t = atom_t.get_parent()
    resname_t = res_t.get_resname()
    auth_seq_t = res_t.id[1]
    name_t = atom_t.get_name()
    pos_t = atom_t.get_coord()
    
    for atom_p in tpp1.get_atoms():
        el_p = atom_p.element
        if el_p not in polar_types: continue
        res_p = atom_p.get_parent()
        resname_p = res_p.get_resname()
        auth_seq_p = res_p.id[1]
        name_p = atom_p.get_name()
        pos_p = atom_p.get_coord()
        
        dist = np.linalg.norm(pos_t - pos_p)
        if dist <= 4.0:
            polar_contacts.append({
                'tert_res': f"{resname_t}{auth_seq_t}",
                'tert_atom': name_t,
                'tpp1_res': f"{resname_p}{auth_seq_p}",
                'tpp1_atom': name_p,
                'dist': dist
            })
            
            # Check salt bridge
            # Case 1: TERT basic, TPP1 acidic
            is_sb1 = ((resname_t, name_t) in basic_atoms) and ((resname_p, name_p) in acidic_atoms)
            # Case 2: TERT acidic, TPP1 basic
            is_sb2 = ((resname_t, name_t) in acidic_atoms) and ((resname_p, name_p) in basic_atoms)
            if is_sb1 or is_sb2:
                salt_bridges.append({
                    'tert_res': f"{resname_t}{auth_seq_t}",
                    'tert_atom': name_t,
                    'tpp1_res': f"{resname_p}{auth_seq_p}",
                    'tpp1_atom': name_p,
                    'dist': dist,
                    'type': 'TERT(+):TPP1(-)' if is_sb1 else 'TERT(-):TPP1(+)'
                })

print(f"Total polar atom contacts (<= 4.0 A): {len(polar_contacts)}")
for pc in sorted(polar_contacts, key=lambda x: x['dist']):
    print(f"  TERT {pc['tert_res']}:{pc['tert_atom']} -- TPP1 {pc['tpp1_res']}:{pc['tpp1_atom']} = {pc['dist']:.2f} A")

print(f"\nSalt bridges (<= 4.0 A): {len(salt_bridges)}")
for sb in sorted(salt_bridges, key=lambda x: x['dist']):
    print(f"  TERT {sb['tert_res']}:{sb['tert_atom']} -- TPP1 {sb['tpp1_res']}:{sb['tpp1_atom']} = {sb['dist']:.2f} A ({sb['type']})")
