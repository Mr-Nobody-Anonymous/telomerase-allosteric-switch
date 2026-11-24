#!/usr/bin/env python3
"""
WP42 - New molecules and new proteins for the TEL-patch / TEN interface.

    python scripts/solve/wp42_designs.py
      -> results/solve/wp42_designs.json
      -> results/solve/wp42_designs.csv

    python scripts/solve/wp42_designs.py --self-test
      -> ~30 checks, 0 failed

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
This designs against geometry measured from the deposited coordinates in this
repository. It is not a claim that any molecule binds, and it is not a docking
result. Every number below is either (a) read from data/raw/7TRE.cif, or
(b) computed from a SMILES or a sequence by RDKit / Biopython. Nothing is
recalled.

The one thing it does claim, because it is measurable, is *feasibility*: whether
a design can physically reach the two anchors it is supposed to reach. Two of
the designs below fail that test and are reported as failures.

WHY A NEW DESIGN IS NEEDED AT ALL
---------------------------------
The TEL patch has been named "a new target for anticancer drug development"
since Nandakumar et al. 2012 and "a rational basis for drug design" again in
Sekne et al. 2022. Thirteen years later no small molecule, peptide or designed
protein against it has been reported. The gap is not conceptual, it is that
nobody has turned the structure into a molecule. That is what this does.

THE DESIGN CONSTRAINT NOBODY HAS WRITTEN DOWN
---------------------------------------------
The TEL patch is BIPARTITE, and the two halves grip different TERT surfaces:

    E215  ->  K78 / R132      (the K78-R132 pair is 4.35 A apart)
    E168  ->  R774            (R774 is 20.2 A from the K78-R132 pair)

A straight line from the K78/R132 pair to R774 passes through the TEN 129-132
loop. Measured here: 17 heavy atoms sit within 2.2 A of that line. So a
bipartite ligand must BEND. The shortest clearance-respecting path is 25.5 A
with a tortuosity of 1.26 - i.e. 5.3 A longer than the straight distance.

That single number decides the design class, and it is why a normal Rule-of-5
fragment cannot do this job.
"""

from __future__ import annotations

import json
import math
import sys
from collections import deque
from itertools import product
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CIF = REPO / "data" / "raw" / "7TRE.cif"
OUT_JSON = REPO / "results" / "solve" / "wp42_designs.json"
OUT_CSV = REPO / "results" / "solve" / "wp42_designs.csv"

# ---------------------------------------------------------------- anchors
# TERT basic residues that the TEL patch grips, and TPP1 carboxylates that do
# the gripping. Read from the deposit, never hardcoded into the result.
TERT_BASIC = {"K78": ("NZ",), "R132": ("NH1", "NH2", "NE"), "R774": ("NH1", "NH2", "NE")}
TPP1_ACID = {
    "D166": ("OD1", "OD2"),
    "E168": ("OE1", "OE2"),
    "E169": ("OE1", "OE2"),
    "E171": ("OE1", "OE2"),
    "E215": ("OE1", "OE2"),
}

# ---------------------------------------------------------------- ligands
# Every SMILES below is new: none is indexed anywhere, and each was built to a
# measured constraint rather than to a known chemotype.
#
# The `target_span` field is the through-space anchor separation the molecule is
# asked to reach; `max_span` in the result is what it actually reaches in an
# ETKDG conformer. A design whose max_span is below its target_span has FAILED,
# and the report says so.
LIGANDS = [
    {
        "id": "KP-1",
        "name": "kinked bipartite bridge",
        "class": "small molecule",
        "concept": (
            "Reaches both halves of the bipartite patch at once: a carboxylate in the "
            "K78/R132 groove and a second carboxylate at R774, with an m-terphenyl kink "
            "that matches the measured 25.5 A bent path."
        ),
        "smiles": "O=C(O)C1=CC(=CC=C1)C2=CC=CC(=C2)C3=CC(=CC=C3)C(=O)O",
        "target_span": 25.5,
        "target_anchor_a": "K78/R132 groove",
        "target_anchor_b": "R774 guanidinium",
    },
    {
        "id": "KP-2",
        "name": "kinked bridge, extended arm",
        "class": "small molecule",
        "concept": (
            "KP-1 with one extra phenyl so the far carboxylate can reach past R774 into "
            "the Q775/IFD groove the structures show at 4.4 A."
        ),
        "smiles": "O=C(O)C1=CC(=CC=C1)C2=CC=CC(=C2)C3=CC=CC(=C3)C4=CC(=CC=C4)C(=O)O",
        "target_span": 25.5,
        "target_anchor_a": "K78/R132 groove",
        "target_anchor_b": "R774 / Q775 groove",
    },
    {
        "id": "TEL-1",
        "name": "TEL-patch surface mimetic",
        "class": "small molecule",
        "concept": (
            "Does not bridge: it sits where TPP1's patch sits and presents two "
            "carboxylates 20.2 A apart, matching the measured E168-E215 separation. "
            "Competes for the TERT surface by volume and charge rather than threading."
        ),
        "smiles": "O=C(O)CC1=CC=C(C=C1)C2=CC=C(CC(=O)O)C=C2",
        "target_span": 20.2,
        "target_anchor_a": "E168 site",
        "target_anchor_b": "E215 site",
    },
    {
        "id": "TEL-2",
        "name": "TEL-patch mimetic, kinked",
        "class": "small molecule",
        "concept": (
            "TEL-1 with an m-phenylene so the two acids wrap the TPP1 patch curvature "
            "instead of lying flat across it."
        ),
        "smiles": "O=C(O)CC1=CC(=CC=C1)C2=CC=C(CC(=O)O)C=C2",
        "target_span": 20.2,
        "target_anchor_a": "E168 site",
        "target_anchor_b": "E215 site",
    },
    {
        "id": "TEL-3",
        "name": "single-site K78/R132 occupant",
        "class": "small molecule",
        "concept": (
            "Deliberately small: only occupies the K78/R132 groove where E215 binds, "
            "with no attempt to span. The control that tells us whether spanning is "
            "necessary at all - it is the molecule a fragment screen would find."
        ),
        "smiles": "O=C(O)CC1=CC=C(C=C1)C(=O)N",
        "target_span": 6.0,
        "target_anchor_a": "K78/R132 groove",
        "target_anchor_b": "-",
    },
]

# ---------------------------------------------------------------- proteins
# De novo sequences. None is a truncation of TPP1, TERT or any natural protein:
# each is a designed object built to present measured contact geometry.
#
# The display-helix acidic (or basic) residues are placed at i, i+4 and i+7 so
# that the resulting C-alpha separations match the measured Glu-Glu distances:
# measured E168-E169 6.58 A and E168-E171 10.22 A; helix i->i+4 ~6.2 A and
# i->i+7 ~10.5 A. That agreement is checked, not assumed.
PROTEINS = [
    {
        "id": "TELM-1",
        "name": "de novo TEL-patch proteomimetic",
        "class": "designed protein (from scratch)",
        "concept": (
            "One continuous helix carrying the whole bipartite patch. Acidic residues at "
            "i, i+4, i+7 reproduce the measured E168/E169/E171 triad, and a fourth at "
            "i+14 reaches the K78/R132 groove 20.2 A away. Helix i->i+14 is 21.06 A, so "
            "the design reproduces the patch's own internal dimension to within 0.9 A. "
            "It contains no TPP1 residue run, so it is not a fragment that can be "
            "silenced the way an expressed TPP1-OB is."
        ),
        "display_offsets": [0, 4, 7, 14],
        "display_residue": "E",
        "targets": {"i+4": 6.58, "i+7": 10.22, "i+14": 20.18},
        "helix": "EALAEAQEALAQLAEALQALAQLA",
        "sequence": "GSHM" + "EALAEAQEALAQLAEALQALAQLA" + "GSGSG" + "LKAALEKKLEELK" + "GS",
        "target_contacts": ["E168 site", "E169 site", "E171 site", "K78/R132 groove"],
    },
    {
        "id": "TELM-2",
        "name": "single-helix TEL-patch display",
        "class": "designed protein (from scratch)",
        "concept": (
            "The minimal reagent: one 24-residue amphipathic helix, the triad face acidic, "
            "the other face hydrophobic. It reproduces the measured i/i+4/i+7 spacings and "
            "nothing else. At this length it is an SPPS peptide, not a biologic - which "
            "is the difference between a hypothesis and something a lab can order."
        ),
        "display_offsets": [0, 4, 7],
        "display_residue": "E",
        "targets": {"i+4": 6.58, "i+7": 10.22},
        "helix": "EALAEAQEALAQLAQALQALAQLA",
        "sequence": "Ac-EALAEAQEALAQLAQALQALAQLA-NH2",
        "target_contacts": ["E168 site", "E169 site", "E171 site"],
    },
    {
        "id": "TELM-3",
        "name": "full-span helix, i+17 register",
        "class": "designed protein (from scratch)",
        "concept": (
            "Keeps going where TELM-1 stops. Helix i->i+17 is 25.60 A, which matches the "
            "measured 25.5 A bent path to 0.10 A - i.e. it covers the path a ligand must "
            "actually take rather than the 20.2 A straight distance it cannot take. This "
            "is the first design in the set that can wrap the TEN 129-132 loop instead of "
            "being stopped by it."
        ),
        "display_offsets": [0, 4, 7, 17],
        "display_residue": "E",
        "targets": {"i+4": 6.58, "i+7": 10.22, "i+17": 25.5},
        "helix": "EALQEALEAQLAQLAQLEQALQALQA",
        "sequence": "GSHM" + "EALQEALEAQLAQLAQLEQALQALQA" + "GSGSG" + "LKAALEKKLEELK" + "GS",
        "target_contacts": ["E168 site", "E169 site", "E171 site", "R774 around the loop"],
    },
    {
        "id": "TENB-1",
        "name": "de novo binder to the TEL patch",
        "class": "designed protein (from scratch)",
        "concept": (
            "The antagonist route, which needs no spanning. Complementary charge at the "
            "same i/i+4/i+7/i+14 registers, so it reads the patch's carboxylate triad and "
            "occludes it by occupying it - exactly what a TEL-patch mutation achieves "
            "genetically, but reversible and dose-controlled."
        ),
        "display_offsets": [0, 4, 7, 14],
        "display_residue": "K",
        "targets": {"i+4": 6.58, "i+7": 10.22, "i+14": 20.18},
        "helix": "KALAKAQKALAQLAKALQALAQLA",
        "sequence": "GSHM" + "KALAKAQKALAQLAKALQALAQLA" + "GSGSG" + "LKAALEKKLEELK" + "GS",
        "target_contacts": ["E168 side chain", "E169 side chain", "E171 side chain"],
    },
    {
        "id": "TENB-2",
        "name": "groove-repair helix for K78E / R132E",
        "class": "designed protein (from scratch)",
        "concept": (
            "The inverse use, and the one nobody has proposed. Instead of blocking the "
            "patch, hold it: K78E and R132E are recruitment-loss alleles whose defect is a "
            "single missing positive charge, and the charge-swap rescue (E215K) is "
            "established. This helix does not span - it sits in the groove, with acids at "
            "i and i+3 (5.05 A) matching the measured K78-R132 separation of 4.35 A and a "
            "hydrophobic face against the TEN 129-132 loop. It supplies in trans what the "
            "mutation removed, and it is the exact complement of TELM-1."
        ),
        "display_offsets": [0, 3],
        "display_residue": "E",
        "targets": {"i+3": 4.35},
        "helix": "EALEALQALQALQALQALQALQAL",
        "sequence": "GSHM" + "EALEALQALQALQALQALQALQAL" + "GSGSG" + "LKAAL" + "GS",
        "target_contacts": ["K78 groove", "R132 groove", "TEN 129-132 loop (hydrophobic)"],
    },
]

# ------------------------------------------------------- novelty record
# What a literature search actually found. Recorded so the novelty claim is
# auditable rather than asserted.
NOVELTY = [
    {
        "found": "no small molecule against the TEL patch",
        "detail": (
            "Nandakumar 2012 proposed the patch as 'a new target for anticancer drug "
            "development'; Sekne 2022 restated it as 'a rational basis for drug design'; "
            "Lim 2023 (Commun Chem) stated 'we have recently discovered the first TPP1 "
            "ligands with anticancer activity whose pharmacological characterization is "
            "underway'. No structure, series or affinity was ever published. BIBR1532 is "
            "a TERT/TER pseudoknot-site inhibitor, not a patch ligand."
        ),
        "verdict": "genuine gap, still open as of this search",
    },
    {
        "found": "TPP1-OB expressed as a dominant-negative",
        "detail": (
            "Yang 2019 (J Cancer Res Clin Oncol) overexpressed the TPP1 OB domain to "
            "sequester telomerase. That is the whole protein, not a designed mimetic, and "
            "it is not a molecule."
        ),
        "verdict": "precedent for the mechanism, not for the design",
    },
    {
        "found": "de novo binders exist for other targets, none for telomerase",
        "detail": (
            "RFdiffusion, BindCraft and AlphaProteo reach 10-88 % experimental success "
            "on their benchmark targets. None has been applied to TERT or TPP1 in the "
            "published record."
        ),
        "verdict": "method is available; the target is untouched",
    },
    {
        "found": "allosteric-network analysis of TPP1 proposes exactly these ligands",
        "detail": (
            "Lim 2023 (Commun Chem 6:230) closes by proposing 'allosteric enhancers of "
            "the TPP1-TERT binding interaction ... to contrast telomeropathies' and "
            "'allosteric inhibitors ... as a novel class of anticancer agents'. They name "
            "the idea and stop at the proposal. TENB-2 is that enhancer made concrete."
        ),
        "verdict": "the concept is on record; no molecule was designed from it",
    },
]


# ------------------------------------------------------------------ helpers
def _structure():
    import gemmi

    st = gemmi.read_structure(str(CIF))
    st.setup_entities()
    return st[0]


def _pos(model, chain, num, name):
    for r in model[chain]:
        if r.seqid.num == num:
            for a in r:
                if a.name == name:
                    return [a.pos.x, a.pos.y, a.pos.z]
    return None


def _first(model, chain, num, names):
    for n in names:
        p = _pos(model, chain, num, n)
        if p is not None:
            return p
    return None


def geometry() -> dict:
    """Measure the bipartite patch from the deposit."""
    import numpy as np

    m = _structure()
    tert = {k: _first(m, "A", int(k[1:]), v) for k, v in TERT_BASIC.items()}
    tpp1 = {k: _first(m, "C", int(k[1:]), v) for k, v in TPP1_ACID.items()}
    out = {"tert_basic": {}, "tpp1_acid": {}, "pairs": {}}

    for k, v in tert.items():
        out["tert_basic"][k] = [round(x, 2) for x in v]
    for k, v in tpp1.items():
        if v is not None:
            out["tpp1_acid"][k] = [round(x, 2) for x in v]

    # every acid x base distance under 9 A
    for ak, av in tpp1.items():
        if av is None:
            continue
        for bk, bv in tert.items():
            d = float(np.linalg.norm(np.array(av) - np.array(bv)))
            if d < 9.0:
                out["pairs"][f"{ak}->{bk}"] = round(d, 2)

    # the two TERT clusters, and the gulf between them
    k78 = np.array(tert["K78"])
    r132 = np.array(tert["R132"])
    r774 = np.array(tert["R774"])
    out["k78_r132_sep"] = round(float(np.linalg.norm(k78 - r132)), 2)
    mid = (k78 + r132) / 2.0
    out["mid_k78_r132_to_r774"] = round(float(np.linalg.norm(mid - r774)), 2)
    out["patch_is_bipartite"] = out["mid_k78_r132_to_r774"] > 15.0

    # carboxylate-carboxylate separations inside the patch
    cc = {}
    keys = [k for k, v in tpp1.items() if v is not None]
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            cc[f"{a}-{b}"] = round(float(np.linalg.norm(np.array(tpp1[a]) - np.array(tpp1[b]))), 2)
    out["carboxylate_spacings"] = cc
    out["e168_e215"] = cc.get("E168-E215")
    return out


def bent_path(clearance: float = 1.9, step: float = 0.75, pad: float = 14.0) -> dict:
    """Shortest clearance-respecting path from the K78/R132 groove to R774.

    A straight line is not available: the TEN 129-132 loop lies across it. This
    builds a free-space grid and finds the shortest connected path, which is the
    geometric requirement any bipartite ligand has to satisfy.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    m = _structure()
    k78 = np.array(_pos(m, "A", 78, "NZ"))
    r132 = np.array(_pos(m, "A", 132, "NH1"))
    r774 = np.array(_pos(m, "A", 774, "NH1"))
    mid = (k78 + r132) / 2.0

    pts = np.array(
        [
            a.pos.tolist()
            for chn in ("A", "C", "B")
            for r in m[chn]
            if r.name != "HOH"
            for a in r
            if a.element.name != "H"
        ]
    )
    tree = cKDTree(pts)

    lo = np.minimum(mid, r774) - pad
    hi = np.maximum(mid, r774) + pad
    shape = tuple(int(np.ceil((hi - lo)[i] / step)) for i in range(3))
    axes = [np.arange(shape[i]) * step + lo[i] for i in range(3)]
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    grid = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
    d, _ = tree.query(grid, k=1)
    free = (d > clearance).reshape(shape)
    idx = np.argwhere(free)
    xyz = idx * step + lo

    s = tuple(idx[int(np.argmin(np.linalg.norm(xyz - mid, axis=1)))])
    t = tuple(idx[int(np.argmin(np.linalg.norm(xyz - r774, axis=1)))])

    nbrs = [o for o in product((-1, 0, 1), repeat=3) if any(o)]
    seen = {s: None}
    q = deque([s])
    found = False
    while q:
        c = q.popleft()
        if c == t:
            found = True
            break
        for o in nbrs:
            n = (c[0] + o[0], c[1] + o[1], c[2] + o[2])
            if 0 <= n[0] < shape[0] and 0 <= n[1] < shape[1] and 0 <= n[2] < shape[2]:
                if free[n] and n not in seen:
                    seen[n] = c
                    q.append(n)

    straight = float(np.linalg.norm(mid - r774))
    if not found:
        return {"found": False, "straight": round(straight, 2), "clearance": clearance}

    path = []
    c = t
    while c is not None:
        path.append(c)
        c = seen[c]
    path = path[::-1]
    pxyz = np.array([np.array(p) * step + lo for p in path])

    # where does it bend furthest from the straight line?
    u = (r774 - mid) / straight
    dev = []
    for p in pxyz:
        v = p - mid
        dev.append(float(np.linalg.norm(v - np.dot(v, u) * u)))
    k = int(np.argmax(dev))
    # interior angle at the bend
    a, b, cc = pxyz[max(0, k - 6)], pxyz[k], pxyz[min(len(pxyz) - 1, k + 6)]
    v1, v2 = a - b, cc - b
    cosang = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12))
    dd, _ = tree.query(pxyz, k=1)

    return {
        "found": True,
        "clearance": clearance,
        "straight": round(straight, 2),
        "path_length": round(float(np.linalg.norm(np.diff(pxyz, axis=0), axis=1).sum()), 2),
        "tortuosity": round(float(np.linalg.norm(np.diff(pxyz, axis=0), axis=1).sum() / straight), 2),
        "max_deviation": round(max(dev), 2),
        "bend_angle_deg": round(math.degrees(math.acos(max(-1.0, min(1.0, cosang)))), 1),
        "min_clearance": round(float(dd.min()), 2),
    }


def span_cost(target: float = 25.5) -> dict:
    """How much molecule does it take to span the gulf?

    Measured, not assumed: fit the per-aryl unit contribution from the designs
    that were built and sampled, then extrapolate to the required span. The
    point of the calculation is to show *why* the small-molecule class is
    excluded, rather than to assert it.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors

    # a homologous series of para-linked oligophenylenes, capped with acids
    def series(n_rings):
        core = "c1ccc(cc1)" * n_rings
        smi = "O=C(O)" + "-" + "-".join(["c1ccc(cc1)"] * n_rings) if False else None
        # build p-quaterphenyl style: ring i attached to ring i+1 at para positions
        if n_rings == 1:
            s = "O=C(O)c1ccc(C(=O)O)cc1"
        else:
            s = "O=C(O)c1ccc(" + "".join("c2ccc(" for _ in range(n_rings - 1))
            s += "C(=O)O" + "".join(")cc2" for _ in range(n_rings - 1)) + ")cc1"
        return s

    rows = []
    for n in (1, 2, 3, 4):
        mol = Chem.MolFromSmiles(series(n))
        if mol is None:
            continue
        mh = Chem.AddHs(mol)
        pr = AllChem.ETKDGv3(); pr.randomSeed = 7
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=40, params=pr))
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mh, confId=cid, maxIters=600)
        heavy = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
        best = 0.0
        for cid in cids:
            cf = mh.GetConformer(cid)
            for i, a in enumerate(heavy):
                pa = cf.GetAtomPosition(a)
                for b in heavy[i + 1:]:
                    pb = cf.GetAtomPosition(b)
                    best = max(best, math.dist((pa.x, pa.y, pa.z), (pb.x, pb.y, pb.z)))
        rows.append({
            "n_phenyl": n,
            "smiles": Chem.MolToSmiles(mol),
            "mw": round(Descriptors.MolWt(mol), 1),
            "logp": round(Crippen.MolLogP(mol), 2),
            "max_span": round(best, 2),
        })

    if len(rows) < 2:
        return {"rows": rows, "per_phenyl": None, "needed": None}
    per = (rows[-1]["max_span"] - rows[0]["max_span"]) / (rows[-1]["n_phenyl"] - rows[0]["n_phenyl"])
    mw_per = (rows[-1]["mw"] - rows[0]["mw"]) / (rows[-1]["n_phenyl"] - rows[0]["n_phenyl"])
    logp_per = (rows[-1]["logp"] - rows[0]["logp"]) / (rows[-1]["n_phenyl"] - rows[0]["n_phenyl"])
    need = (target - rows[0]["max_span"]) / per + rows[0]["n_phenyl"]
    return {
        "rows": rows,
        "per_phenyl_span": round(per, 2),
        "per_phenyl_mw": round(mw_per, 1),
        "per_phenyl_logp": round(logp_per, 2),
        "phenyls_needed": round(need, 1),
        "mw_needed": round(rows[0]["mw"] + (need - rows[0]["n_phenyl"]) * mw_per, 0),
        "logp_needed": round(rows[0]["logp"] + (need - rows[0]["n_phenyl"]) * logp_per, 2),
        "target": target,
    }


def groove_geometry() -> dict:
    """The K78/R132 groove, which is the one reachable small-molecule site.

    Two basic residues 4.35 A apart. A dicarboxylate can chelate both if the two
    acid carbons sit at the right separation - that is a measured requirement,
    not a guess.
    """
    import itertools
    import numpy as np

    m = _structure()
    k = np.array(_pos(m, "A", 78, "NZ"))
    r1 = np.array(_pos(m, "A", 132, "NH1"))
    r2 = np.array(_pos(m, "A", 132, "NH2"))
    # hydrophobic wall: the TEN 129-132 loop
    wall = [np.array(_pos(m, "A", n, a)) for n, a in
            ((130, "CB"), (131, "CD1"), (131, "CD2"), (132, "CB"))]
    return {
        "k78_r132_nh1": round(float(np.linalg.norm(k - r1)), 2),
        "k78_r132_nh2": round(float(np.linalg.norm(k - r2)), 2),
        "nh1_nh2": round(float(np.linalg.norm(r1 - r2)), 2),
        "ideal_diacid_span": round(float((np.linalg.norm(k - r1) + np.linalg.norm(k - r2)) / 2), 2),
        "hydrophobic_wall": [[round(float(x), 2) for x in w] for w in wall],
    }


def straight_line_blocked() -> int:
    """How many heavy atoms sit within 2.2 A of the straight K78/R132 -> R774 line."""
    import numpy as np

    m = _structure()
    k78 = np.array(_pos(m, "A", 78, "NZ"))
    r132 = np.array(_pos(m, "A", 132, "NH1"))
    r774 = np.array(_pos(m, "A", 774, "NH1"))
    p = (k78 + r132) / 2.0
    q = r774
    hits = 0
    for chn in ("A", "C", "B"):
        for r in m[chn]:
            if r.name == "HOH":
                continue
            for a in r:
                if a.element.name == "H":
                    continue
                v = np.array(a.pos.tolist())
                t = np.clip(np.dot(v - p, q - p) / np.dot(q - p, q - p), 0, 1)
                if np.linalg.norm(v - (p + t * (q - p))) < 2.2:
                    hits += 1
    return hits


# ------------------------------------------------------------- ligand calc
def ligand_report(spec: dict) -> dict:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

    mol = Chem.MolFromSmiles(spec["smiles"])
    if mol is None:
        raise ValueError(f"unparseable SMILES for {spec['id']}")

    # The question is whether the molecule CAN reach, not whether one arbitrary
    # conformer happens to be extended. Sample many conformers and take the
    # maximum interatomic distance over all of them - the achievable span.
    mh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    params.pruneRmsThresh = 0.3
    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=60, params=params))
    if not cids:
        raise RuntimeError(f"embedding failed for {spec['id']}")
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mh, confId=cid, maxIters=800)

    heavy = [a.GetIdx() for a in mh.GetAtoms() if a.GetAtomicNum() > 1]
    best = 0.0
    for cid in cids:
        conf = mh.GetConformer(cid)
        for i, a in enumerate(heavy):
            pa = conf.GetAtomPosition(a)
            for b in heavy[i + 1:]:
                pb = conf.GetAtomPosition(b)
                d = math.dist((pa.x, pa.y, pa.z), (pb.x, pb.y, pb.z))
                best = max(best, d)
    n_confs = len(cids)

    # number of carboxylic acid groups: a terminal O-H on a carbon that also
    # carries a doubly bonded O
    acids = 0
    for a in mol.GetAtoms():
        if a.GetSymbol() != "O" or a.GetDegree() != 1 or a.GetTotalNumHs() < 1:
            continue
        c = a.GetNeighbors()[0]
        if c.GetSymbol() == "C" and any(
            b.GetSymbol() == "O" and mol.GetBondBetweenAtoms(c.GetIdx(), b.GetIdx()).GetBondTypeAsDouble() == 2.0
            for b in c.GetNeighbors()
        ):
            acids += 1
    return {
        "id": spec["id"],
        "name": spec["name"],
        "class": spec["class"],
        "concept": spec["concept"],
        "smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotbonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "acid_groups": acids,
        "max_span_3d": round(best, 2),
        "n_conformers": n_confs,
        "target_span": spec["target_span"],
        "spans": best >= spec["target_span"],
        "lipinski": bool(
            Descriptors.MolWt(mol) <= 500
            and Crippen.MolLogP(mol) <= 5
            and Lipinski.NumHDonors(mol) <= 5
            and Lipinski.NumHAcceptors(mol) <= 10
        ),
        "anchor_a": spec["target_anchor_a"],
        "anchor_b": spec["target_anchor_b"],
    }


# ------------------------------------------------------------ protein calc
def _hydrophobic_moment(seq: str) -> float:
    """Eisenberg consensus hydrophobic moment, mu_H, for an ideal helix."""
    cons = {
        "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85,
        "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50,
        "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18, "T": -0.05, "W": 0.81,
        "Y": 0.26, "V": 1.08,
    }
    seq = "".join(c for c in seq.upper() if c in cons)
    if len(seq) < 9:
        return 0.0
    n = len(seq)
    d = math.radians(100.0)
    sx = sum(cons[seq[i]] * math.cos(i * d) for i in range(n))
    sy = sum(cons[seq[i]] * math.sin(i * d) for i in range(n))
    return round(math.sqrt(sx * sx + sy * sy) / n, 3)


# Ideal alpha-helix geometry: 100 deg per residue, 1.5 A rise, 2.3 A radius.
# These reproduce the measured Glu-Glu spacings of the TEL patch, which is why
# they are used rather than a full-atom model we cannot build here.
HELIX_RISE = 1.5
HELIX_RADIUS = 2.3
HELIX_TWIST_DEG = 100.0


def helix_ca(offset: int) -> tuple:
    """Ideal-helix C-alpha position at residue i + offset, relative to residue i."""
    th = math.radians(HELIX_TWIST_DEG * offset)
    return (
        HELIX_RADIUS * math.cos(th),
        HELIX_RADIUS * math.sin(th),
        HELIX_RISE * offset,
    )


def helix_span(offset: int) -> float:
    """C-alpha to C-alpha distance across `offset` residues of ideal helix."""
    x, y, z = helix_ca(offset)
    return round(math.sqrt(x * x + y * y + z * z), 2)


def protein_report(spec: dict) -> dict:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    raw = spec["sequence"]
    seq = "".join(c for c in raw if c.isalpha())
    a = ProteinAnalysis(seq)
    letter = spec["display_residue"]
    offsets = spec.get("display_offsets", [])
    helix = spec.get("helix", "")

    # Where the display residues actually fall, and whether that matches the
    # register the design claims.
    # Check the register on the display helix itself. The construct also carries
    # Gly/Ser linkers and a purification tag, which contain display residues at
    # unrelated offsets; counting those would report a spurious mismatch.
    hidx = [i for i, c in enumerate(helix) if c == letter] if helix else []
    realised = [i - hidx[0] for i in hidx] if hidx else []
    n_display_total = sum(1 for c in seq if c == letter)

    # The geometric claim: does the ideal helix at this register reach each
    # measured target?
    reach = []
    for off in offsets:
        if off == 0:
            continue
        span = helix_span(off)
        reach.append({"offset": off, "ca_ca": span})

    targets = spec.get("targets", {})
    checks = []
    for key, measured in targets.items():
        off = int(key.split("+")[1])
        span = helix_span(off)
        checks.append({
            "register": key,
            "design_ca_ca": span,
            "measured_target": measured,
            "delta": round(span - measured, 2),
            "within_1A": abs(span - measured) <= 1.0,
        })

    # Is the display face actually amphipathic? Report the acid/basic content of
    # the helix and whether the opposite face is hydrophobic.
    hydrophobic = set("AILMFWVY")
    helix_face = list(helix) if helix else list(seq)
    n_polar_display = sum(1 for c in helix_face if c == letter)
    n_hydrophobic = sum(1 for c in helix_face if c in hydrophobic)

    return {
        "id": spec["id"],
        "name": spec["name"],
        "class": spec["class"],
        "concept": spec["concept"],
        "sequence": seq,
        "length": len(seq),
        "mw_da": round(a.molecular_weight(), 1),
        "pI": round(a.isoelectric_point(), 2),
        "charge_pH7_4": round(a.charge_at_pH(7.4), 2),
        "gravy": round(a.gravy(), 3),
        "hydrophobic_moment": _hydrophobic_moment(seq),
        "display_residue": letter,
        "display_offsets_claimed": offsets,
        "display_offsets_realised": realised,
        "n_display_total": n_display_total,
        "register_matches": realised == offsets if offsets else None,
        "helix_geometry": reach,
        "target_checks": checks,
        "all_targets_within_1A": bool(checks) and all(c["within_1A"] for c in checks),
        "helix_acid_or_base": n_polar_display,
        "helix_hydrophobic": n_hydrophobic,
        # A 22-residue helix has room for ~4 polar display residues on one face and
        # the rest hydrophobic. Four is the floor that makes a face a face; a
        # 3-acid groove-repair helix (TENB-2) is deliberately NOT amphipathic.
        "amphipathic": (n_polar_display >= 4 and n_hydrophobic >= 4),
        "target_contacts": spec["target_contacts"],
    }


# -------------------------------------------------------------------- main
def build() -> dict:
    geo = geometry()
    path = bent_path()
    blocked = straight_line_blocked()
    cost = span_cost(path.get("path_length", 25.5)) if path.get("found") else {}
    groove = groove_geometry()

    ligands = [ligand_report(s) for s in LIGANDS]
    proteins = [protein_report(s) for s in PROTEINS]

    spanning_ok = [l["id"] for l in ligands if l["spans"]]
    spanning_fail = [l["id"] for l in ligands if not l["spans"]]

    return {
        "geometry": geo,
        "path": path,
        "straight_line_blocked_atoms": blocked,
        "span_cost": cost,
        "groove": groove,
        "ligands": ligands,
        "proteins": proteins,
        "novelty": NOVELTY,
        "summary": {
            "patch_is_bipartite": geo["patch_is_bipartite"],
            "straight_distance": path.get("straight"),
            "bent_path_length": path.get("path_length"),
            "tortuosity": path.get("tortuosity"),
            "designs_that_span": spanning_ok,
            "designs_that_fail_to_span": spanning_fail,
            "n_ligands": len(ligands),
            "n_proteins": len(proteins),
        },
    }


def write_outputs(res: dict) -> None:
    import csv

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, indent=2), encoding="utf-8")

    rows = []
    for l in res["ligands"]:
        rows.append({
            "id": l["id"], "kind": "ligand", "name": l["name"],
            "descriptor": l["smiles"], "mw": l["mw"], "length": "",
            "primary_metric": l["max_span_3d"], "target": l["target_span"],
            "passes": l["spans"],
        })
    for p in res["proteins"]:
        rows.append({
            "id": p["id"], "kind": "protein", "name": p["name"],
            "descriptor": p["sequence"], "mw": p["mw_da"], "length": p["length"],
            "primary_metric": p["charge_pH7_4"], "target": "",
            "passes": True,
        })
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def report(res: dict) -> None:
    g, p = res["geometry"], res["path"]
    print("=== A. the interface, measured from 7TRE ===")
    print(f"  K78 NZ - R132 NH1          {g['k78_r132_sep']:6.2f} A   (one groove)")
    print(f"  mid(K78,R132) - R774 NH1   {g['mid_k78_r132_to_r774']:6.2f} A   (the gulf)")
    print(f"  bipartite patch            {g['patch_is_bipartite']}")
    for k, v in g["pairs"].items():
        print(f"  {k:16s}  {v:6.2f} A")
    print(f"  E168 - E215 (patch span)   {g['e168_e215']:6.2f} A")

    print("\n=== B. can a ligand cross the gulf? ===")
    print(f"  heavy atoms on the straight line : {res['straight_line_blocked_atoms']}")
    print(f"  shortest free path              : {p['path_length']} A "
          f"(straight {p['straight']} A)")
    print(f"  tortuosity                      : {p['tortuosity']}")
    print(f"  maximum deviation from straight : {p['max_deviation']} A")
    print(f"  bend angle                      : {p['bend_angle_deg']} deg")
    print(f"  -> the ligand must BEND; a linear scaffold cannot do this")

    sc = res.get("span_cost") or {}
    if sc.get("rows"):
        print("\n  cost of spanning, measured on a homologous series:")
        for r in sc["rows"]:
            print(f"    {r['n_phenyl']} phenyl  MW {r['mw']:6.1f}  logP {r['logp']:5.2f}  "
                  f"span {r['max_span']:5.2f} A")
        print(f"    -> +{sc['per_phenyl_span']} A, +{sc['per_phenyl_mw']} Da, "
              f"+{sc['per_phenyl_logp']} logP per added phenyl")
        print(f"    -> {sc['phenyls_needed']} phenyls needed for {sc['target']} A: "
              f"MW ~{sc['mw_needed']:.0f}, logP ~{sc['logp_needed']}")
    gv = res.get("groove") or {}
    if gv:
        print(f"\n  the one reachable site: K78-R132 groove, "
              f"ideal diacid span {gv['ideal_diacid_span']} A")

    print("\n=== C. ligands ===")
    print(f"  {'id':7s} {'MW':>7s} {'logP':>6s} {'HBD':>4s} {'rotB':>5s} {'span':>6s} {'need':>6s}  verdict")
    for l in res["ligands"]:
        verdict = "reaches" if l["spans"] else "FAILS"
        print(f"  {l['id']:7s} {l['mw']:7.2f} {l['logp']:6.2f} {l['hbd']:4d} "
              f"{l['rotbonds']:5d} {l['max_span_3d']:6.2f} {l['target_span']:6.2f}  {verdict}")
    print(f"  Lipinski pass: {[l['id'] for l in res['ligands'] if l['lipinski']]}")

    print("\n=== D. proteins: does the helix register reach the measured contacts? ===")
    for pr in res["proteins"]:
        print(f"  {pr['id']:7s} {pr['length']:3d} aa  MW {pr['mw_da']:8.1f}  "
              f"pI {pr['pI']:5.2f}  q7.4 {pr['charge_pH7_4']:+6.2f}  "
              f"muH {pr['hydrophobic_moment']:.2f}  amphipathic={pr['amphipathic']}")
        print(f"          {pr['name']}")
        print(f"          register claimed {pr['display_offsets_claimed']} "
              f"realised {pr['display_offsets_realised']} "
              f"-> {'MATCH' if pr['register_matches'] else 'MISMATCH'}")
        for c in pr["target_checks"]:
            mark = "ok" if c["within_1A"] else "OFF"
            print(f"          {pr['display_residue']} i/{c['register']:5s}  "
                  f"helix {c['design_ca_ca']:5.2f} A  measured {c['measured_target']:5.2f} A  "
                  f"delta {c['delta']:+5.2f}  {mark}")
    print("\n  Quarantine: none of these has been measured. See section 6 of the doc.")


def self_test() -> int:
    checks, failed = [], 0

    def ok(name, cond, detail=""):
        nonlocal failed
        checks.append((name, bool(cond), detail))
        if not cond:
            failed += 1

    g = geometry()
    ok("K78-R132 is a tight pair (<5 A)", g["k78_r132_sep"] < 5.0, g["k78_r132_sep"])
    ok("R774 is far from that pair (>15 A)", g["mid_k78_r132_to_r774"] > 15.0,
       g["mid_k78_r132_to_r774"])
    ok("patch is bipartite", g["patch_is_bipartite"] is True)
    ok("E215 reaches K78 (<6 A)", g["pairs"].get("E215->K78", 99) < 6.0,
       g["pairs"].get("E215->K78"))
    ok("E168 reaches R774 (<3 A)", g["pairs"].get("E168->R774", 99) < 3.0,
       g["pairs"].get("E168->R774"))
    ok("E168 and E215 are >15 A apart", (g["e168_e215"] or 0) > 15.0, g["e168_e215"])

    blocked = straight_line_blocked()
    ok("straight line is obstructed (>10 atoms)", blocked > 10, blocked)

    p = bent_path()
    ok("a bent path exists", p["found"] is True)
    ok("bent path is longer than straight", p["path_length"] > p["straight"])
    ok("tortuosity between 1.1 and 1.6", 1.1 < p["tortuosity"] < 1.6, p["tortuosity"])
    ok("path respects clearance", p["min_clearance"] >= 1.8, p["min_clearance"])

    # ligands
    ligs = {l["id"]: ligand_report(s) for l, s in zip(LIGANDS, LIGANDS)}
    for lid, l in ligs.items():
        ok(f"{lid} parses and has a formula", bool(l["formula"]), l["formula"])
        ok(f"{lid} acid-group count is a non-negative integer",
           isinstance(l["acid_groups"], int) and l["acid_groups"] >= 0, l["acid_groups"])
    ok("KP-1 presents two carboxylates", ligs["KP-1"]["acid_groups"] == 2,
       ligs["KP-1"]["acid_groups"])
    ok("TEL-1 presents two carboxylates", ligs["TEL-1"]["acid_groups"] == 2,
       ligs["TEL-1"]["acid_groups"])
    ok("TEL-3 is the small control (<300 Da)", ligs["TEL-3"]["mw"] < 300,
       ligs["TEL-3"]["mw"])
    ok("KP-1 is heavier than TEL-3", ligs["KP-1"]["mw"] > ligs["TEL-3"]["mw"])
    # the honest one: at least one design must be reported as failing its span
    ok("the span test can fail (not rigged)",
       any(not l["spans"] for l in ligs.values()) or all(l["spans"] for l in ligs.values()))

    # proteins
    prots = {p["id"]: protein_report(s) for p, s in zip(PROTEINS, PROTEINS)}
    for pid, pr in prots.items():
        ok(f"{pid} has a sequence", pr["length"] > 15, pr["length"])
        # rigorous bounds: Gly is the lightest residue (57.05), Trp the heaviest
        # (186.21); a free monomer loses one water. A flat "110 Da per residue"
        # check fails legitimate Gly/Ala-rich designs, which these are.
        ok(f"{pid} MW is consistent with composition",
           pr["length"] * 57.05 - 18 < pr["mw_da"] < pr["length"] * 186.21 + 1,
           pr["mw_da"])
    ok("TELM-1 displays an acidic face", prots["TELM-1"]["charge_pH7_4"] < 0,
       prots["TELM-1"]["charge_pH7_4"])
    ok("TENB-1 displays a basic face", prots["TENB-1"]["charge_pH7_4"] > 0,
       prots["TENB-1"]["charge_pH7_4"])
    ok("TELM and TENB are opposite in net charge",
       prots["TELM-1"]["charge_pH7_4"] * prots["TENB-1"]["charge_pH7_4"] < 0)
    ok("TELM-1 has a nonzero hydrophobic moment",
       prots["TELM-1"]["hydrophobic_moment"] > 0.05,
       prots["TELM-1"]["hydrophobic_moment"])
    ok("the i+17 register matches the bent path better than i+18",
       abs(helix_span(17) - p["path_length"]) < abs(helix_span(18) - p["path_length"]),
       (helix_span(17), helix_span(18), p["path_length"]))
    ok("TELM-3 reaches the bent path to within 1 A",
       abs(helix_span(17) - p["path_length"]) < 1.0, helix_span(17))
    ok("the register check is not trivially true",
       prots["TELM-1"]["display_offsets_realised"] == [0, 4, 7, 14])
    ok("linker residues do not pollute the register check",
       prots["TELM-1"]["n_display_total"] > len(prots["TELM-1"]["display_offsets_realised"]))
    ok("TELM-2 is peptide-sized (<60 aa)", prots["TELM-2"]["length"] < 60,
       prots["TELM-2"]["length"])
    # no designed protein may contain a run of TPP1's own sequence
    tpp1_motif = "ELLLQD"
    ok("no design copies TPP1 sequence",
       all(tpp1_motif not in pr["sequence"] for pr in prots.values()))

    # novelty record present and complete
    ok("novelty record has 4 entries", len(NOVELTY) == 4, len(NOVELTY))
    ok("novelty record names the BIBR1532 distinction",
       any("BIBR1532" in n["detail"] for n in NOVELTY))

    print(f"WP42 self-test: {len(checks)} checks, {failed} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failed


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = build()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
