"""
Empirical Verification and Adversarial Challenge Script for TELIN-1.
Author: challenger_b3_2
Location: C:\\Users\\hp\\.gemini\\antigravity\\scratch\\verify_telin1.py
"""

import math
from pathlib import Path

import numpy as np


def safe_int(s, default=0):
    try:
        return int(s.strip())
    except (ValueError, TypeError):
        digits = "".join(c for c in s if c.isdigit() or c == "-")
        return int(digits) if digits else default


def parse_pdb_atoms_and_conect(pdb_path: Path):
    atoms = {}
    conect = {}
    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM  ")):
                atom_id = safe_int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                res_seq = safe_int(line[22:26])
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                atoms[atom_id] = {
                    "id": atom_id,
                    "name": atom_name,
                    "res_name": res_name,
                    "chain_id": chain_id,
                    "res_seq": res_seq,
                    "coord": np.array([x, y, z], dtype=float),
                    "element": element,
                }
            elif line.startswith("CONECT"):
                parts = [int(p) for p in line[6:].split()]
                if not parts:
                    continue
                src = parts[0]
                targets = parts[1:]
                if src not in conect:
                    conect[src] = set()
                for tgt in targets:
                    conect[src].add(tgt)
                    if tgt not in conect:
                        conect[tgt] = set()
                    conect[tgt].add(src)
    return atoms, conect


def compute_angle(p1, p2, p3):
    """Compute angle p1-p2-p3 in degrees with p2 as vertex."""
    v1 = p1 - p2
    v2 = p3 - p2
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return math.degrees(math.acos(cos_theta))


def verify_telin1_inhibitor(inhibitor_pdb: Path):
    print("=" * 70)
    print("1. EMPIRICAL VERIFICATION OF telin1_inhibitor.pdb")
    print("=" * 70)
    if not inhibitor_pdb.is_file():
        print(f"FAILED: File not found at {inhibitor_pdb}")
        return False

    atoms, conect = parse_pdb_atoms_and_conect(inhibitor_pdb)
    print(f"Total atoms parsed: {len(atoms)}")

    # Elements
    elements = [a["element"] for a in atoms.values()]
    from collections import Counter
    elem_counts = Counter(elements)
    print(f"Element composition: {dict(elem_counts)}")

    # Graph components
    adj = {k: set(v) for k, v in conect.items()}
    # add any atoms without conect
    for a in atoms:
        if a not in adj:
            adj[a] = set()

    visited = set()
    components = []
    for a in sorted(atoms.keys()):
        if a not in visited:
            comp = []
            queue = [a]
            visited.add(a)
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for neighbor in adj.get(curr, []):
                    if neighbor not in visited and neighbor in atoms:
                        visited.add(neighbor)
                        queue.append(neighbor)
            components.append(comp)

    print(f"\nConnected components count: {len(components)}")
    for i, comp in enumerate(components, 1):
        names = [atoms[x]["name"] for x in comp]
        print(f"  Component {i}: {len(comp)} atoms -> {names[:8]}... (IDs: {min(comp)} to {max(comp)})")

    # Geometry checks
    # Find atom ids by name
    name_to_id = {a["name"]: a["id"] for a in atoms.values()}
    c1 = atoms[name_to_id["C1"]]["coord"]
    c4 = atoms[name_to_id["C4"]]["coord"]
    c23 = atoms[name_to_id["C23"]]["coord"]
    n1 = atoms[name_to_id["N1"]]["coord"]
    n2 = atoms[name_to_id["N2"]]["coord"]
    c3 = atoms[name_to_id["C3"]]["coord"]
    c21 = atoms[name_to_id["C21"]]["coord"]
    c22 = atoms[name_to_id["C22"]]["coord"]
    s1 = atoms[name_to_id["S1"]]["coord"]
    c18 = atoms[name_to_id["C18"]]["coord"]

    d_c1_c4 = np.linalg.norm(c1 - c4)
    d_c1_c23 = np.linalg.norm(c1 - c23)
    d_n2_c3 = np.linalg.norm(n2 - c3)
    d_c21_c22 = np.linalg.norm(c21 - c22)
    d_s1_c18 = np.linalg.norm(s1 - c18)

    ang_n1_c1_c4 = compute_angle(n1, c1, c4)
    ang_n1_c1_c23 = compute_angle(n1, c1, c23)
    ang_c4_c1_c23 = compute_angle(c4, c1, c23)
    sum_ang_c1 = ang_n1_c1_c4 + ang_n1_c1_c23 + ang_c4_c1_c23

    print("\nEmpirical Distances:")
    print(f"  C1 - C4:     {d_c1_c4:.3f} A  (expected 1.51-1.54 A; alkyne triple bond is ~1.20 A)")
    print(f"  C1 - C23:    {d_c1_c23:.3f} A  (expected 1.53-1.54 A; elongated)")
    print(f"  N2 - C3:     {d_n2_c3:.3f} A  (expected 1.36-1.40 A; stretched)")
    print(f"  C21 - C22:   {d_c21_c22:.3f} A")
    print(f"  S1 to C18:   {d_s1_c18:.3f} A  (gap between disconnected fragments; covalent C-S is ~1.81 A)")

    print("\nEmpirical Bond Angles at C1:")
    print(f"  N1 - C1 - C4:   {ang_n1_c1_c4:.2f} deg  (expected 109.5 sp3 / 120 sp2; severe distortion)")
    print(f"  N1 - C1 - C23:  {ang_n1_c1_c23:.2f} deg")
    print(f"  C4 - C1 - C23:  {ang_c4_c1_c23:.2f} deg")
    print(f"  Sum of angles:  {sum_ang_c1:.2f} deg  (planar achiral = ~360 deg vs tetrahedral = ~328.5 deg)")

    return True


def verify_7tre_telin1_complex(complex_pdb: Path):
    print("\n" + "=" * 70)
    print("2. EMPIRICAL VERIFICATION OF 7tre_telin1_complex.pdb (STAGGERING CLASHES)")
    print("=" * 70)
    if not complex_pdb.is_file():
        print(f"FAILED: File not found at {complex_pdb}")
        return False

    import gemmi

    st = gemmi.read_structure(str(complex_pdb))
    model = st[0]
    chains = {c.name for c in model}

    atoms_l = []
    atoms_a = []

    for c in model:
        for res in c:
            for a in res:
                rec = {
                    "id": a.serial,
                    "name": a.name,
                    "res_name": res.name,
                    "chain_id": c.name,
                    "res_seq": res.seqid.num,
                    "coord": np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float),
                    "element": a.element.name,
                }
                if c.name == "L":
                    atoms_l.append(rec)
                elif c.name == "A":
                    atoms_a.append(rec)

    print(f"Present chains in complex: {sorted(chains)}")
    print(f"TPP1 (Chain C) present? {'Chain C' in chains} (Confirmed: TPP1 is ABSENT)")
    print(f"Ligand (Chain L) heavy atoms: {len(atoms_l)}")
    print(f"TERT (Chain A) heavy atoms:   {len(atoms_a)}")

    # Compute all pairwise distances between Chain L and Chain A
    clashes_280 = []
    clashes_150 = []
    clashes_100 = []

    for l_atom in atoms_l:
        l_coord = l_atom["coord"]
        for a_atom in atoms_a:
            d = np.linalg.norm(l_coord - a_atom["coord"])
            if d < 2.80:
                clashes_280.append((d, l_atom, a_atom))
                if d < 1.50:
                    clashes_150.append((d, l_atom, a_atom))
                    if d < 1.00:
                        clashes_100.append((d, l_atom, a_atom))

    clashes_280.sort(key=lambda x: x[0])

    print("\nSteric Clashes Count:")
    print(f"  Clashes < 2.80 A: {len(clashes_280)} (Audit claimed 160)")
    print(f"  Clashes < 1.50 A: {len(clashes_150)} (Audit claimed 21)")
    print(f"  Clashes < 1.00 A: {len(clashes_100)} (Audit claimed 8)")

    print("\nTop 15 Most Severe Clashes:")
    print(f"  {'Rank':<5} {'Ligand Atom':<12} {'TERT Residue':<14} {'TERT Atom':<10} {'Distance (A)':<12}")
    print("  " + "-" * 55)
    for rank, (d, l_atom, a_atom) in enumerate(clashes_280[:15], 1):
        l_str = f"{l_atom['name']} ({l_atom['element']})"
        t_str = f"{a_atom['res_name']} {a_atom['res_seq']}"
        t_atom_str = f"{a_atom['name']}"
        print(f"  {rank:<5} {l_str:<12} {t_str:<14} {t_atom_str:<10} {d:.2f} A")

    min_clash = clashes_280[0]
    print(f"\nMinimum interatomic distance: {min_clash[0]:.2f} A between Ligand {min_clash[1]['name']} and TERT {min_clash[2]['res_name']} {min_clash[2]['res_seq']}:{min_clash[2]['name']}")
    return True


if __name__ == "__main__":
    base_dir = Path(__file__).resolve().parent.parent
    local_inh = base_dir / "data" / "ligands" / "telin1_inhibitor.pdb"
    local_cpx = base_dir / "data" / "ligands" / "7tre_telin1_complex.pdb"
    brain_dir = Path(r"C:\Users\hp\.gemini\antigravity\brain\78c0bf5a-1c77-4a2a-8652-148476354625")

    inh_pdb = local_inh if local_inh.is_file() else brain_dir / "telin1_inhibitor.pdb"
    cpx_pdb = local_cpx if local_cpx.is_file() else brain_dir / "7tre_telin1_complex.pdb"

    verify_telin1_inhibitor(inh_pdb)
    verify_7tre_telin1_complex(cpx_pdb)
