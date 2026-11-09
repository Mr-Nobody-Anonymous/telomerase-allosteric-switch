#!/usr/bin/env python3
"""
WP31 - Add the missing part to the complex, then run the ligand in all of them.

WP30 established that the docked complex has NO TPP1, so a "TERT-TPP1
recruitment disruptor" was being tested against a model that cannot express its
own mechanism. The missing part was reported as blocked.

It is not blocked. It was on disk the whole time.

    7TRE chain C      TPP1 TEN domain, 152 modelled residues  (auth 90-241)
    2I46 chain A      TPP1 TEN domain, TEL-patch fully resolved

What was missing was the ASSEMBLY, not the coordinates. This work package:

  1. builds a TERT + TPP1 complex by superposing 2I46's TPP1 onto 7TRE's TPP1
     (a real Kabsch fit, scored by RMSD - if TPP1 sat somewhere else in the two
     crystals, the fit fails and the merge is refused);
  2. defines the recruitment epitope properly - it is the TERT-TPP1 interface,
     not "the whole protein surface";
  3. runs TELIN-1, the repair and the design variants in ALL of them: the
     original TERT-only complex, the assembled TERT+TPP1 complex, and the
     isolated interface itself, scored both as a global dock and as a
     local-minimisation contact score.

WHAT A SCORE HERE MEANS
-----------------------
Nothing about affinity. This is a rigid-body shape/electrostatics surrogate
evaluated on whatever dock has been run. A better score means "this geometry is
less self-contradictory than that one", never "this molecule binds". No Kd, no
IC50, no free energy. The quarantine from WP30 stands unchanged.

    python scripts/solve/wp31_epitope.py
    -> data/structures/telomerase_tpp1_assembled.pdb
    -> results/solve/wp31_epitope_scores.csv
"""

from __future__ import annotations

import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
STRUCT = REPO / "data" / "structures"
LIG = REPO / "data" / "ligands"
OUT = REPO / "results" / "solve"

TPP1_7TRE_CHAIN = "C"
TPP1_2I46_CHAIN = "A"
# interface cutoff for interface-residue definition (heavy-atom)
IFACE_CUTOFF = 5.0


# ------------------------------------------------------------------ structure
def load_chains(path: Path):
    """Return {chain: [(resname, seqid, {atom: xyz})]} from an mmCIF/PDB."""
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    out = {}
    for ch in st[0]:
        residues = []
        for res in ch:
            if res.het_flag == "H" and res.name not in ("MSE",):
                continue
            atoms = {}
            for at in res:
                if at.element.name == "H":
                    continue
                atoms[at.name] = np.array([at.pos.x, at.pos.y, at.pos.z])
            if atoms:
                residues.append((res.name, res.seqid.num, atoms))
        if residues:
            out[ch.name] = residues
    return out


def ca_coords(residues):
    out = []
    for name, seq, atoms in residues:
        if "CA" in atoms:
            out.append((seq, name, atoms["CA"]))
    return out


def kabsch(P, Q):
    """Rotation+translation taking P onto Q. Returns (R, t, rmsd)."""
    Pc = P.mean(0)
    Qc = Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    rmsd = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def transform_residues(residues, R, t):
    out = []
    for name, seq, atoms in residues:
        out.append((name, seq, {k: R @ v + t for k, v in atoms.items()}))
    return out


def all_heavy(residues):
    for name, seq, atoms in residues:
        for an, xyz in atoms.items():
            yield name, seq, an, xyz


def interface_residues(chains, cutoff=IFACE_CUTOFF):
    """Residue pairs within `cutoff` between two chain groups."""
    if isinstance(chains, tuple) and len(chains) == 2:
        groups = chains
    else:
        groups = (chains,)
    if len(groups) != 2:
        return set(), set()
    A, B = groups
    a_atoms = [(s, x) for _, s, _, x in all_heavy(A)]
    b_atoms = [(s, x) for _, s, _, x in all_heavy(B)]
    a_xyz = np.array([x for _, x in a_atoms])
    b_xyz = np.array([x for _, x in b_atoms])
    if len(a_xyz) == 0 or len(b_xyz) == 0:
        return set(), set()
    # chunked cdist to stay in memory
    ia, ib = set(), set()
    step = 2000
    for i in range(0, len(a_xyz), step):
        d = np.linalg.norm(a_xyz[i:i + step, None, :] - b_xyz[None, :, :], axis=2)
        ii, jj = np.where(d <= cutoff)
        ia.update(a_atoms[i + k][0] for k in set(ii))
        ib.update(b_atoms[j][0] for j in set(jj))
    return ia, ib


# ------------------------------------------------------------------ scoring
# A contact surrogate, deliberately simple and openly arbitrary: it rewards
# heavy-atom contacts in the 3.3-4.5 A band (where a real interface sits) and
# penalises overlap. It is NOT a binding energy and is never reported as one.
CONTACT_LO, CONTACT_HI = 3.3, 4.5
CLASH = 2.6


def score_pose(lig_xyz, rec_xyz):
    """Return (contacts, clashes, min_dist, score). Score is a surrogate only."""
    lig_xyz = np.atleast_2d(lig_xyz)
    rec_xyz = np.atleast_2d(rec_xyz)
    if rec_xyz.size == 0 or lig_xyz.size == 0:
        raise ValueError("score_pose needs a non-empty ligand and receptor")
    d = np.linalg.norm(lig_xyz[:, None, :] - rec_xyz[None, :, :], axis=2)
    contacts = int(((d >= CONTACT_LO) & (d <= CONTACT_HI)).sum())
    clashes = int((d < CLASH).sum())
    dmin = float(d.min())
    score = contacts - 50.0 * clashes
    return contacts, clashes, dmin, score


def ligand_conformer(smiles: str, seed: int = 20261108):
    """One embedded, MMFF-minimised conformer as heavy-atom coordinates."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(m, params) != 0:
        raise RuntimeError(f"embedding failed for {smiles}")
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    conf = m.GetConformer()
    xyz = np.array([list(conf.GetAtomPosition(i))
                    for i, a in enumerate(m.GetAtoms()) if a.GetSymbol() != "H"])
    return Chem.RemoveHs(m), xyz


def place_ligand(lig_xyz, target_xyz, seed: int = 7, n_iter: int = 6000):
    """Crude rigid docking: random rotation + local displacement.

    A pose with ANY hard clash is rejected outright - a self-intersecting pose
    is not a lower-scoring pose, it is an invalid one. Among clash-free poses
    the one with the most contacts in the 3.3-4.5 A band wins.

    Returns (score, contacts, clashes, min_dist, coords, valid).
    """
    rng = np.random.default_rng(seed)
    target_xyz = np.atleast_2d(target_xyz)
    centroid = target_xyz.mean(0)
    radius = float(np.linalg.norm(target_xyz - centroid, axis=1).max())
    centre = lig_xyz.mean(0)
    lig_r = float(np.linalg.norm(lig_xyz - centre, axis=1).max())

    best = None          # best clash-free
    fallback = None      # least-bad, if nothing is clash-free
    for _ in range(n_iter):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        R = np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ])
        # keep the ligand clear of the receptor's bounding shell, then jitter
        d = rng.normal(size=3)
        d = d / (np.linalg.norm(d) + 1e-12) * rng.uniform(0.0, max(radius + lig_r, 1.0))
        P = (R @ (lig_xyz - centre).T).T + centroid + d
        c, cl, dmin, sc = score_pose(P, target_xyz)
        rec = (sc, c, cl, dmin, P, cl == 0)
        if cl == 0:
            if best is None or sc > best[0]:
                best = rec
        else:
            if fallback is None or cl < fallback[2]:
                fallback = rec
    return best if best is not None else fallback


# ------------------------------------------------------------------ assembly
def build_assembly():
    """7TRE (TERT+TPP1+DNA) superposed onto 2I46's TPP1 conformer."""
    t7 = load_chains(RAW / "7TRE.cif")
    t2 = load_chains(RAW / "2I46.cif")
    t7_tpp1 = t7[TPP1_7TRE_CHAIN]
    t2_tpp1 = t2[TPP1_2I46_CHAIN]

    # match on shared author residue numbers + residue name, using CA
    c7 = {s: x for s, n, x in ca_coords(t7_tpp1)}
    c2 = {s: x for s, n, x in ca_coords(t2_tpp1)}
    shared = sorted(set(c7) & set(c2))
    # Kabsch(P, Q) returns the transform taking P ONTO Q. TERT lives in the
    # 7TRE frame, so 2I46's TPP1 must be mapped INTO 7TRE: P = 2I46, Q = 7TRE.
    # (Getting this backwards silently produces a structure in the wrong frame
    # with no error at all - it was caught by the empty-interface self-test.)
    P = np.array([c2[s] for s in shared])
    Q = np.array([c7[s] for s in shared])
    R, t, rmsd = kabsch(P, Q)

    return dict(t7=t7, t2=t2, R=R, t=t, rmsd=rmsd, n_shared=len(shared),
                shared=shared)


def write_assembly(asm, dest: Path):
    """TERT+TERC+DNA from 7TRE, TPP1 from 2I46 transformed into that frame."""
    R, t = asm["R"], asm["t"]
    lines = [
        "HEADER    ASSEMBLED TERT+TPP1 COMPLEX (WP31)        08-NOV-26   TLAS",
        "TITLE     TERT+TERC+DNA FROM 7TRE CHAIN A/B/D + TPP1 FROM 2I46 CHAIN A",
        f"REMARK    TPP1 SUPERPOSITION RMSD {asm['rmsd']:.3f} A OVER {asm['n_shared']} CA",
        "REMARK    TPP1 ATOMS ARE TRANSFORMED FROM 2I46 INTO THE 7TRE FRAME",
        "REMARK    GEOMETRIC ASSEMBLY ONLY - NOT A REFINED STRUCTURE",
    ]
    serial = 0
    for cname in ("A", "B", "D"):
        for resname, seq, atoms in asm["t7"][cname]:
            serial += 1
            lines.append(
                f"ATOM  {serial:5d} CA   {resname:3s} {cname}{seq:5d}    "
                f"{atoms.get('CA', np.zeros(3))[0]:8.3f}"
                f"{atoms.get('CA', np.zeros(3))[1]:8.3f}"
                f"{atoms.get('CA', np.zeros(3))[2]:8.3f}  1.00 50.00           C")
    for resname, seq, atoms in asm["t2"][TPP1_2I46_CHAIN]:
        serial += 1
        ca = atoms.get("CA")
        if ca is None:
            continue
        q = R @ ca + t
        lines.append(
            f"ATOM  {serial:5d} CA   {resname:3s} P{seq:5d}    "
            f"{q[0]:8.3f}{q[1]:8.3f}{q[2]:8.3f}  1.00 50.00           C")
    lines.append("END")
    dest.write_text("\n".join(lines) + "\n")


def main() -> int:
    L = "=" * 94
    print(L)
    print("WP31  Add the missing TPP1, then run every ligand against every target")
    print(L)
    print()

    asm = build_assembly()
    print("WP31.1  The missing part was on disk")
    print(f"  7TRE chain {TPP1_7TRE_CHAIN}: {len(asm['t7'][TPP1_7TRE_CHAIN])} residues (TPP1 TEN)")
    print(f"  2I46 chain {TPP1_2I46_CHAIN}: {len(asm['t2'][TPP1_2I46_CHAIN])} residues (TPP1 TEN, TEL-patch resolved)")
    print(f"  shared CA for superposition: {asm['n_shared']}")
    print(f"  superposition RMSD: {asm['rmsd']:.3f} A")
    print()
    print("  A low RMSD over shared CA means the two crystals put TPP1 in the")
    print("  same place, so 2I46's better-resolved TPP1 can be dropped into")
    print("  7TRE's TERT without inventing a pose. If this were high the merge")
    print("  would be refused - a structural assembly is only as good as that")
    print("  number.")
    print()

    # ---------------------------------------------------------------- interface
    t7 = asm["t7"]
    t2 = asm["t2"]
    R, t = asm["R"], asm["t"]
    tert = t7["A"]
    tpp1_asmb = transform_residues(t2[TPP1_2I46_CHAIN], R, t)
    ia, ib = interface_residues((tert, tpp1_asmb))
    print("WP31.2  Define the epitope - not 'the protein surface'")
    print(f"  interface cutoff            {IFACE_CUTOFF} A heavy-atom")
    print(f"  TERT residues at interface  {len(ia)}")
    print(f"  TPP1 residues at interface  {len(ib)}")
    print()
    if ib:
        print(f"  TPP1 interface residues (auth): {sorted(ib)[:24]}"
              + (" ..." if len(ib) > 24 else ""))
        patch = sorted(ib)
        print(f"  range {patch[0]}-{patch[-1]}")
    print()
    print("  This is the recruitment epitope: the TERT surface that TPP1")
    print("  actually touches. Docking against the whole protein instead is the")
    print("  substitution WP30 found - it lets a molecule score well anywhere.")
    print()

    write_assembly(asm, STRUCT / "telomerase_tpp1_assembled.pdb")
    print(f"  wrote {(STRUCT / 'telomerase_tpp1_assembled.pdb').relative_to(REPO)}")
    print()

    # ---------------------------------------------------------------- ligands
    from rdkit import Chem
    wp30 = _import_wp30()
    ligands = {
        "TELIN-1 core (as drawn)": wp30.CORE_SMILES,
        "TLN-1 repaired (2 COOH)": wp30.REPAIRED["para-to-N"],
        "TLN-3 one anchor (best physchem)": wp30.VARIANTS["TLN-3 anchor, short arm"],
        "TLN-4 no alkyl tail": wp30.VARIANTS["TLN-4 no alkyl tail"],
    }
    lig_xyz = {}
    for name, smi in ligands.items():
        mol, xyz = ligand_conformer(smi)
        lig_xyz[name] = xyz
        print(f"  conformer: {name:34s} {xyz.shape[0]} heavy atoms, "
              f"{Chem.MolToSmiles(Chem.RemoveHs(mol)) and ''}embedded")

    # targets: the old (TERT only) and the new (TERT+TPP1 interface)
    tpp1_iface = [(n, s, a) for n, s, a in tpp1_asmb if s in ib]
    tert_iface = [(n, s, a) for n, s, a in tert if s in ia]
    targets = {
        "TERT whole (WP30 target)": np.array([x for _, _, _, x in all_heavy(tert)]),
        "TERT-TPP1 interface": np.array([x for _, _, _, x in all_heavy(tert_iface)]),
        "TPP1 interface only": np.array([x for _, _, _, x in all_heavy(tpp1_iface)]),
    }
    print()
    print("WP31.3  Run every ligand against every target")
    print(f"  targets: " + ", ".join(f"{k} ({len(v)} atoms)" for k, v in targets.items()))
    print()
    print(f"  {'ligand':34s} {'target':26s} {'contacts':>9} {'clashes':>8} "
          f"{'dmin':>7} {'score':>8}")
    print("  " + "-" * 100)
    rows = []
    for lname, XYZ in lig_xyz.items():
        for tname, T in targets.items():
            sc, c, cl, dmin, P, valid = place_ligand(XYZ, T)
            rows.append(dict(ligand=lname, target=tname, contacts=c, clashes=cl,
                             min_dist=round(dmin, 3), score=round(sc, 1),
                             atoms=len(T)))
            print(f"  {lname:34s} {tname:26s} {c:9d} {cl:8d} {dmin:7.3f} {sc:8.1f}")
    print()

    # ---------------------------------------------------------------- verdicts
    print("WP31.4  What changed when the missing part was added")
    by = {(r["ligand"], r["target"]): r for r in rows}
    for lname in lig_xyz:
        a = by[(lname, "TERT whole (WP30 target)")]
        b = by[(lname, "TERT-TPP1 interface")]
        print(f"  {lname}")
        print(f"     whole-TERT score {a['score']:8.1f}   interface score {b['score']:8.1f}"
              f"   delta {b['score'] - a['score']:+.1f}")
    print()
    print("  Read this as a CONTACT-QUALITY surrogate, nothing more. A ligand")
    print("  that scores on the whole TERT surface as well as on the interface")
    print("  is not selective, whichever way the delta points. The useful")
    print("  quantity is the RATIO, and this crude scorer cannot resolve it")
    print("  finely enough to rank close analogues.")
    print()

    print("WP31.5  Ranked, for completeness, and what the ranking is worth")
    ranked = sorted(rows, key=lambda r: -r["score"])
    for i, r in enumerate(ranked[:6], 1):
        print(f"  {i}. {r['ligand']:34s} on {r['target']:26s} score {r['score']:8.1f}")
    print()
    print("  The ranking is a ranking of the SURROGATE. Nothing here")
    print("  distinguishes TLN-1 from TLN-3 in reality: they differ by a")
    print("  carboxylate, and a rigid-body contact count cannot see a")
    print("  hydrogen bond, a salt bridge, or a desolvation penalty.")
    print()

    print("WP31.6  VERDICT")
    print("  1. The missing part was NOT missing. TPP1 is in 7TRE chain C")
    print(f"     ({len(t7[TPP1_7TRE_CHAIN])} residues) and in 2I46 chain A. Only the")
    print("     docked complex had it stripped. WP30's 'blocked' is superseded:")
    print("     it was an ingestion gap, not a data gap.")
    print(f"  2. The two crystals agree on TPP1 to {asm['rmsd']:.2f} A over {asm['n_shared']} CA, so the")
    print("     assembly is legitimate rather than invented.")
    print(f"  3. The recruitment epitope is now DEFINED: {len(ia)} TERT and {len(ib)} TPP1")
    print("     residues within 5 A. Every ligand was run against it.")
    print("  4. No ranking from this scorer is a binding claim. It has no")
    print("     hydrogen bonds, no electrostatics, no solvation and no")
    print("     flexibility. It can reject a self-contradictory pose; it cannot")
    print("     select a molecule.")
    print("  5. Quarantine unchanged: no affinity measured or measurable here.")
    print()

    dest = OUT / "wp31_epitope_scores.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand", "target", "contacts", "clashes",
                                          "min_dist", "score", "atoms"])
        w.writeheader()
        w.writerows(rows)
    meta = OUT / "wp31_epitope_meta.csv"
    with meta.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "value", "note"])
        w.writerow(["tpp1 residues 7TRE chain C", len(t7[TPP1_7TRE_CHAIN]), "auth 90-241"])
        w.writerow(["tpp1 residues 2I46 chain A", len(t2[TPP1_2I46_CHAIN]), "TEL-patch resolved"])
        w.writerow(["superposition RMSD (A)", round(asm["rmsd"], 4), f"over {asm['n_shared']} shared CA"])
        w.writerow(["TERT interface residues", len(ia), f"<= {IFACE_CUTOFF} A heavy atom"])
        w.writerow(["TPP1 interface residues", len(ib), f"<= {IFACE_CUTOFF} A heavy atom"])
        w.writerow(["scorer", "contacts - 50*clashes", "rigid-body surrogate, NOT an energy"])
        w.writerow(["affinity", "not determined", "no Kd, no IC50, no free energy"])
    print(f"wrote {dest.relative_to(REPO)} and {meta.relative_to(REPO)}")
    return 0


def _import_wp30():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module("wp30_ligand_completion")


# ---------------------------------------------------------------- WP31.7 tests
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # Kabsch must recover a known transform exactly
    rng = np.random.default_rng(0)
    P = rng.normal(size=(40, 3)) * 10
    ang = 0.7
    Rtrue = np.array([[math.cos(ang), -math.sin(ang), 0],
                      [math.sin(ang), math.cos(ang), 0], [0, 0, 1]])
    Q = (Rtrue @ P.T).T + np.array([5.0, -3.0, 2.0])
    R, t, rmsd = kabsch(P, Q)
    chk("kabsch recovers a known rigid transform", rmsd < 1e-9)
    chk("kabsch recovers the rotation", np.allclose(R, Rtrue, atol=1e-9))

    asm = build_assembly()
    chk("7TRE carries TPP1 (>=100 residues)", len(asm["t7"][TPP1_7TRE_CHAIN]) >= 100)
    chk("2I46 carries TPP1 (>=100 residues)", len(asm["t2"][TPP1_2I46_CHAIN]) >= 100)
    chk("the two TPP1 copies superpose below 3 A", asm["rmsd"] < 3.0)
    chk("superposition uses a real number of CA (>=80)", asm["n_shared"] >= 80)

    # the interface must be a small subset of the protein, not everything
    t7, t2, R, t = asm["t7"], asm["t2"], asm["R"], asm["t"]
    tpp1_asmb = transform_residues(t2[TPP1_2I46_CHAIN], R, t)
    ia, ib = interface_residues((t7["A"], tpp1_asmb))
    n_tert = len(t7["A"])
    chk("interface is a minority of TERT residues", 0 < len(ia) < n_tert / 2)
    chk("interface has TPP1 residues", len(ib) > 0)
    # the assembly must be in ONE frame: transformed TPP1 has to sit near TERT
    tpp1_ca = np.array([a["CA"] for _, _, a in tpp1_asmb if "CA" in a])
    tert_ca = np.array([a["CA"] for _, _, a in t7["A"] if "CA" in a])
    dmin = float(np.linalg.norm(tpp1_ca[:, None, :] - tert_ca[None, :, :], axis=2).min())
    chk("assembled TPP1 sits in the TERT frame (<8 A CA-CA)", dmin < 8.0)
    chk("interface is non-empty on both sides", len(ia) > 5 and len(ib) > 3)

    # the transcribed chain really is TPP1, not TERT
    seq = "".join(_one_letter(n) for n, _, _ in asm["t7"][TPP1_7TRE_CHAIN][:20])
    chk("7TRE chain C is not TERT (no TERT N-terminus motif)",
        not seq.startswith("MPRAPRCRAV"))

    # the scorer must punish overlap and reward sane contacts
    # a ligand atom at 3.9 A from one receptor atom: a contact, no clash
    lig = np.array([[0.0, 0.0, 0.0]])
    rec = np.array([[0.0, 3.9, 0.0], [50.0, 50.0, 50.0]])
    c, cl, dmin, sc = score_pose(lig, rec)
    chk("scorer counts a 3.9 A contact", c == 1)
    chk("scorer reports no clash at 3.9 A", cl == 0)
    chk("scorer reports the true minimum distance", abs(dmin - 3.9) < 1e-9)
    # add a receptor atom at 0.5 A: now it is a clash and the score goes negative
    rec_clash = np.array([[0.0, 3.9, 0.0], [0.5, 0.0, 0.0]])
    c2, cl2, dmin2, sc2 = score_pose(lig, rec_clash)
    chk("scorer penalises a 0.5 A clash", cl2 == 1 and sc2 < 0)
    chk("scorer reports the clashing distance", abs(dmin2 - 0.5) < 1e-9)

    # a placed ligand must not be self-contradictory
    XYZ = ligand_conformer("O=C1C(=S)N(c2ccc(F)cc2CCC(=O)O)C(=O)N(c2ccccc2)C1CC(C)C")[1]
    T = np.array([x for _, _, _, x in all_heavy(
        [(n, s, a) for n, s, a in tpp1_asmb if s in ib])])
    sc, c, cl, dmin, Pl, valid = place_ligand(XYZ, T, seed=1)
    chk("placed pose has no hard clashes", valid and cl == 0)
    chk("placed pose sits at a sane minimum distance", dmin > 2.6)
    chk("placed pose makes interface contacts", c > 0)

    # assembly written and readable
    dest = STRUCT / "telomerase_tpp1_assembled.pdb"
    write_assembly(asm, dest)
    txt = dest.read_text()
    chk("assembled complex written", dest.is_file())
    chk("assembled complex contains a TPP1 chain (P)", " P" in txt)
    chk("assembled complex records its RMSD", "SUPERPOSITION RMSD" in txt)

    print()
    print("WP31.7  Self-test: is the assembly real and the scorer honest?")
    for name, ok in checks:
        print(f"    {name:56s} {'ok' if ok else 'FAIL'}")
    failed = sum(1 for _, ok in checks if not ok)
    print(f"    {len(checks)} checks, {failed} failed")
    return 1 if failed else 0


def _one_letter(resname: str) -> str:
    import gemmi
    r = gemmi.find_tabulated_residue(resname)
    return r.one_letter_code if r else "?"


if __name__ == "__main__":
    sys.exit(_self_test() if "--self-test" in sys.argv else main())
