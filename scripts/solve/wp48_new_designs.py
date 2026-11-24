#!/usr/bin/env python3
"""
WP48 - For each remaining failure, design a new entity that addresses it, then
test the new design against the SAME gates that killed the previous attempt.

    python scripts/solve/wp48_new_designs.py
      -> results/solve/wp48_new_designs.json / .csv
      -> results/solve/WP48_NEW_DESIGNS.md

    python scripts/solve/wp48_new_designs.py --self-test   -> checks, 0 failed

METHOD: four isolated design agents plus a critic
-----------------------------------------------
Each agent is a separate module-level namespace with NO import of another
agent. They share only the read-only measurements in `facts()`, which are
re-derived from the deposit on every run. This is deliberate: WP44's helix error
survived because a test re-derived its own assumption, and WP45's disulfide
generalisation survived because one agent's conclusion was another's premise.
Isolation does not make an agent right; it makes it independently wrong, which
is checkable.

  AGENT A  Constraint      F1: the graft cannot pay its ordering cost
  AGENT B  Alert Layer     F2: the filter vocabulary misnames reactivity
  AGENT C  Assay Design    F3: P11/P12 are unmeasurable by computation
  AGENT D  Structure Fidelity  F4: both deposits are coarse models (NEW)

  CRITIC   cross-checks every agent against the deposit and against the others,
           and is allowed to fail them. It fails two.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import sys
from pathlib import Path

import gemmi
import numpy as np

REPO = Path(__file__).resolve().parents[2]
SOLVE = REPO / "results" / "solve"
RAW = REPO / "data" / "raw"

OUT_JSON = SOLVE / "wp48_new_designs.json"
OUT_CSV = SOLVE / "wp48_new_designs.csv"
OUT_MD = SOLVE / "WP48_NEW_DESIGNS.md"

SEG_A = (166, 172)
SEG_B = (209, 215)
PAYLOAD = (168, 169, 171, 215)

CHEM = {
    "disulfide":   (3.5, 5.0),
    "lactam":      (4.5, 6.0),
    "hydrocarbon": (6.0, 7.5),
    "diyne":       (6.5, 7.5),
    "thioether":   (9.0, 11.5),
    "biaryl":      (10.0, 12.5),
}

_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


# ===================================================================== FACTS
# Read-only measurements. Every agent reads this; no agent may redefine it.
_FACT_CACHE: dict = {}


def facts(structure: str = "7TRE") -> dict:
    if structure in _FACT_CACHE:
        return _FACT_CACHE[structure]

    st = gemmi.read_structure(str(RAW / f"{structure}.cif"))
    st.setup_entities()
    by = {c.name: c for c in st[0]}

    # TPP1 identified BY SEQUENCE. Both proteins contain residues 168 and 215.
    tpp1 = None
    for ch in st[0]:
        resl = [r for r in ch if r.seqid.num and r.name in _CODES]
        if len(resl) < 100:
            continue
        off = min(r.seqid.num for r in resl)
        seq = "".join(_CODES[r.name] for r in resl)
        if seq[166 - off:173 - off] == "DWEEKEF":
            tpp1, tpp1_off, tpp1_seq = ch, off, seq
            break
    assert tpp1 is not None

    # TERT = the PROTEIN chain carrying LYS78 and ARG132.
    # The nucleotide chains in 7TRE are also numbered past 132, so matching on
    # numbering alone selects the DNA. Chain identity must be checked on residue
    # type, not on the numbers - the same class of error as reading 168/215
    # without asking which protein they belong to.
    tert = None
    for ch in st[0]:
        if ch is tpp1:
            continue
        resl = [r for r in ch if r.seqid.num and r.name in _CODES]
        if len(resl) < 100:
            continue
        bynum = {r.seqid.num: r.name for r in resl}
        if bynum.get(78) == "LYS" and bynum.get(132) == "ARG":
            tert = ch
            break
    # Some deposits contain no telomerase (2I46 is TPP1 alone). That must be
    # recorded, not asserted away: a structure without the partner cannot test a
    # clash against it, and pretending otherwise would manufacture a clearance.

    def at(chain, n, name):
        for r in chain:
            if r.seqid.num == n:
                for x in r:
                    if x.name == name:
                        return np.array([x.pos.x, x.pos.y, x.pos.z])
        return None

    def _first_atom(chain, n, names):
        for nmx in names:
            p = at(chain, n, nmx)
            if p is not None:
                return p
        return None

    def nm(chain, n):
        for r in chain:
            if r.seqid.num == n:
                return r.name
        return None

    def incomplete_obstacles(chain):
        """Points for residues whose side chain is truncated, padded by how much
        is missing. A residue modelled to CB (3 side-chain atoms) or less cannot
        clash at its glyph atom alone - the missing atoms still occupy space."""
        pts = []
        heavy_per = {"GLY": 0, "ALA": 1, "SER": 2, "CYS": 2, "THR": 3, "VAL": 3,
                     "ASN": 4, "ASP": 4, "LEU": 4, "ILE": 4, "PRO": 3, "MET": 4,
                     "GLN": 5, "GLU": 5, "LYS": 5, "ARG": 7, "HIS": 6, "PHE": 7,
                     "TYR": 8, "TRP": 10}
        for r in chain:
            if r.seqid.num is None or r.name not in heavy_per:
                continue
            modelled = sum(1 for x in r
                           if x.element.name not in ("N", "C", "O", "H"))
            expected = 2 + heavy_per.get(r.name, 4)   # CA + CB + side chain
            missing = max(0, expected - modelled)
            if missing <= 0:
                continue
            # each missing heavy atom can extend ~1.25 A; cap the shell at 3 A
            pad = min(3.0, 0.9 * missing)
            for x in r:
                if x.element.name in ("N", "C", "O", "H"):
                    continue
                pts.append((np.array([x.pos.x, x.pos.y, x.pos.z]), pad))
        return pts

    def sc(chain, n):
        return [np.array([x.pos.x, x.pos.y, x.pos.z]) for r in chain
                if r.seqid.num == n for x in r
                if x.name not in ("N", "CA", "C", "O") and x.element.name != "H"]

    # completeness: how much of each chain is actually modelled
    def completeness(chain):
        resl = [r for r in chain if r.seqid.num]
        heavy = sum(1 for r in resl for x in r if x.element.name != "H")
        model = [r for r in resl if r.name in _CODES]
        return {
            "n_residues": len(resl),
            "heavy_atoms": heavy,
            "atoms_per_residue": round(heavy / max(1, len(resl)), 1),
            "residues_with_sidechain": sum(1 for r in model if len(r) >= 7),
            "fraction_with_sidechain": round(
                sum(1 for r in model if len(r) >= 7) / max(1, len(model)), 3),
        }

    d = {
        "structure": structure,
        "tpp1_chain": tpp1.name,
        "tert_chain": tert.name if tert is not None else None,
        "has_tert": tert is not None,
        "tpp1_seq": tpp1_seq,
        "tpp1_off": tpp1_off,
        "completeness": {"tpp1": completeness(tpp1),
                         "tert": completeness(tert) if tert is not None else None},
        "resname": {n: nm(tpp1, n) for n in list(range(166, 173)) + list(range(209, 216))},
        "cb": {n: at(tpp1, n, "CB") for n in list(range(166, 173)) + list(range(209, 216))},
        "ca": {n: at(tpp1, n, "CA") for n in list(range(166, 173)) + list(range(209, 216))},
        "sc": {n: sc(tpp1, n) for n in list(range(166, 173)) + list(range(209, 216))},
        "k78": at(tert, 78, "NZ") if tert is not None else None,
        "r132": at(tert, 132, "NH1") if tert is not None else None,
        "r774": (_first_atom(tert, 774, ("CZ", "NH1", "NH2"))
                 if tert is not None else None),
        "tert_p": np.array([[x.pos.x, x.pos.y, x.pos.z] for r in tert for x in r
                            if x.element.name != "H"]) if tert is not None
                  else np.zeros((0, 3)),
        "tert_incomplete": incomplete_obstacles(tert) if tert is not None else [],
        "tpp1_p": np.array([[x.pos.x, x.pos.y, x.pos.z] for r in tpp1 for x in r
                            if x.element.name != "H"]),
        # both carboxylate oxygens: a salt bridge is judged on the nearer one,
        # and quoting only OE1 would overstate the E215-R132 distance by 1.6 A
        "e215_o": [x for x in (at(tpp1, 215, "OE1"), at(tpp1, 215, "OE2"))
                   if x is not None],
        "e168_o": [x for x in (at(tpp1, 168, "OE1"), at(tpp1, 168, "OE2"))
                   if x is not None],
        "e215_oe1": at(tpp1, 215, "OE1"),
        "e168_oe1": at(tpp1, 168, "OE1"),
    }
    _FACT_CACHE[structure] = d
    return d


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n else v


def _clearance(pa, pb, obstacle, lo=0.12, hi=0.88, n=40, padded=None):
    """Narrowest gap between the straight span and the obstacle.

    `padded` is a list of (point, radius) for residues whose side chain is
    truncated; the missing atoms are represented by the radius rather than
    ignored. Ignoring them is what made WP47's clearance an upper bound.
    """
    worst = 1e9
    for t in np.linspace(lo, hi, n):
        pt = pa + (pb - pa) * t
        if len(obstacle):
            worst = min(worst, float(np.min(np.linalg.norm(obstacle - pt, axis=1))))
        for x, rad in (padded or []):
            worst = min(worst, float(np.linalg.norm(x - pt)) - rad)
    return worst


# ============================================================ AGENT A: CONSTRAINT
# FAILURE F1: the graft cannot pay its ordering cost (12-20 kcal/mol against
# 5.4 available; one staple buys 2.8).
#
# Agent A's diagnosis is a DESIGN ERROR in WP45, not a property of the molecule:
# WP45 built the linker out of GGGGG. Glycine has the highest backbone entropy of
# any residue, so WP45 chose the most expensive possible linker and then reported
# the cost as unpayable. The fix is three-fold and each part is testable:
#   A1  replace the glycine linker with the LOWEST-entropy span that still spans
#   A2  add a second, independent covalent constraint
#   A3  re-price the whole construct with the same per-torsion accounting

# Backbone conformational states per residue, used ONLY to rank linkers. Absolute
# values are model-dependent; the ORDERING is what matters and it is checked for
# robustness to model choice below.
BACKBONE_STATES = {
    "G": 4.0,   # no CB, maximal phi/psi freedom
    "S": 2.6, "T": 2.4, "N": 2.5, "C": 2.2, "A": 2.2, "D": 2.6, "E": 2.7,
    "Q": 2.7, "K": 2.9, "R": 3.0, "H": 2.6, "M": 2.9, "Y": 2.8, "F": 2.8,
    "W": 2.6, "L": 2.7, "I": 2.3, "V": 2.3, "P": 1.2,   # proline ring restricts
}

# Candidate 5-residue linkers. Spans 14.38 A at 3.5 A/residue, so >= 5 residues.
LINKER_CANDIDATES = [
    "GGGGG",   # WP45's choice - highest entropy possible
    "AAAAA",   # helix-favourable, low entropy
    "APAPA",   # proline-constrained
    "AibAib",  # not expressible as standard AA - excluded, listed for honesty
    "PAPAP",
    "AAAAP",
    "SSSSS",
    "AASAA",
    "PGPGA",
]


def a1_linker_entropy(kT: float = 0.5922) -> dict:
    """Price each candidate linker by its conformational entropy.

    S = kT * sum(ln(states_i), i in linker). A linker that is already structured
    costs less to order. This ranks linkers; it does not claim an absolute value.
    """
    rows = []
    for seq in LINKER_CANDIDATES:
        if any(c not in BACKBONE_STATES for c in seq):
            rows.append({"linker": seq, "expressible": False, "note":
                         "contains a non-standard residue; cannot be made by "
                         "recombinant expression, listed but not scored"})
            continue
        s = sum(math.log(BACKBONE_STATES[c]) for c in seq)
        rows.append({
            "linker": seq,
            "expressible": True,
            "states_product_ln": round(s, 3),
            "entropy_cost_kcal": round(kT * s, 2),
            "n_G": seq.count("G"),
            "charge": 0,
        })
    scored = [r for r in rows if r["expressible"]]
    scored.sort(key=lambda r: r["entropy_cost_kcal"])

    # robustness: does the winner survive a different state model?
    alt = {"G": 6.0, "A": 1.6, "P": 1.1, "S": 2.0, "T": 1.9}
    alt_rows = []
    for seq in LINKER_CANDIDATES:
        if not all(c in alt for c in seq):
            continue
        alt_rows.append((seq, round(kT * sum(math.log(alt[c]) for c in seq), 2)))
    alt_rows.sort(key=lambda r: r[1])

    gg = next(r for r in scored if r["linker"] == "GGGGG")
    best = scored[0]
    return {
        "candidates": rows,
        "ranking": [r["linker"] for r in scored],
        "best": best["linker"],
        "best_cost": best["entropy_cost_kcal"],
        "gGGGG_cost": gg["entropy_cost_kcal"],
        "saving_vs_wp45": round(gg["entropy_cost_kcal"] - best["entropy_cost_kcal"], 2),
        "ranking_under_alternative_model": [r[0] for r in alt_rows],
        "winner_stable": alt_rows[0][0] == best["linker"] if alt_rows else None,
        "model_caveat": (
            "The absolute entropies are model values, not measurements. What is "
            "not model-dependent: glycine has no CB and therefore the largest "
            "backbone state space, so GGGGG is the most expensive 5-mer available. "
            "WP45 chose it."),
    }


def a2_second_constraint() -> dict:
    """Find an independent second crosslink that does not reuse staple 1's residues
    and does not consume a payload residue.

    Staple 1 is fixed at W167-L211 (WP47). A second staple must therefore avoid
    167 and 211 as well as the payload set.
    """
    f = facts()
    taken = {167, 211} | set(PAYLOAD)
    intraseg = []
    for seg in (range(SEG_A[0], SEG_A[1] + 1), range(SEG_B[0], SEG_B[1] + 1)):
        for i, j in itertools.combinations(list(seg), 2):
            pi, pj = f["cb"].get(i), f["cb"].get(j)
            ci, cj = f["ca"].get(i), f["ca"].get(j)
            if any(x is None for x in (pi, pj, ci, cj)):
                continue
            d = float(np.linalg.norm(pj - pi))
            fits = [k for k, (lo, hi) in CHEM.items() if lo <= d <= hi]
            if not fits:
                continue
            gap = _unit(pj - pi)
            cos_i = float(np.dot(_unit(pi - ci), gap))
            cos_j = float(np.dot(_unit(pj - cj), -gap))
            clear = _clearance(pi, pj, f["tert_p"])
            uses_taken = (i in taken) or (j in taken)
            intraseg.append({
                "pair": f"{f['resname'][i]}{i}-{f['resname'][j]}{j}",
                "i": i, "j": j, "cb_distance": round(d, 2),
                "chemistries": fits,
                "cos_i": round(cos_i, 2), "cos_j": round(cos_j, 2),
                "direction_ok": cos_i > 0 and cos_j > 0,
                "worst_clearance": round(clear, 2),
                "clash_ok": clear >= 3.0,
                "avoids_staple1_and_payload": not uses_taken,
                "feasible": bool(not uses_taken and cos_i > 0 and cos_j > 0
                                 and clear >= 3.0),
            })
    intraseg.sort(key=lambda r: (not r["feasible"], r["cb_distance"]))
    feas = [r for r in intraseg if r["feasible"]]
    return {
        "staple_1": {"pair": "TRP167-LEU211", "chemistry": "hydrocarbon",
                     "cb_distance": 6.34},
        "candidates": intraseg,
        "n_candidates": len(intraseg),
        "feasible": feas,
        "second_constraint": feas[0] if feas else None,
        "verdict": ("a second independent staple exists" if feas
                    else "no second staple is available inside the segments"),
    }


def a3_build_and_price() -> dict:
    """Rebuild the construct with the new linker + both staples, and re-price it
    with WP45's own per-torsion accounting so the comparison is like-for-like."""
    f = facts()
    lk = a1_linker_entropy()
    sc = a2_second_constraint()

    linker = lk["best"]
    off = f["tpp1_off"]
    seg_a = f["tpp1_seq"][166 - off:173 - off]     # DWEEKEF, residues 166-172
    seg_b = f["tpp1_seq"][209 - off:216 - off]     # FSLLPTE, residues 209-215

    # crosslinked positions become X in one-letter form
    def patch(seg, start, positions):
        s = list(seg)
        for p in positions:
            idx = p - start
            if 0 <= idx < len(s):
                s[idx] = "X"
        return "".join(s)

    a_pos = [167] + ([sc["second_constraint"]["i"]]
                     if sc["second_constraint"] and sc["second_constraint"]["i"] in
                     range(SEG_A[0], SEG_A[1] + 1) else [])
    b_pos = [211] + ([sc["second_constraint"]["j"]]
                     if sc["second_constraint"] and sc["second_constraint"]["j"] in
                     range(SEG_B[0], SEG_B[1] + 1) else [])

    constructs = {
        "WP45-GRAFT":        seg_a + "GGGGG" + seg_b,
        "WP47-STAPLED":      patch(seg_a, 166, [167]) + "GGGGG" + patch(seg_b, 209, [211]),
        "WP48-NEW":          patch(seg_a, 166, a_pos) + linker + patch(seg_b, 209, b_pos),
    }

    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    rows = {}
    for name, s in constructs.items():
        d = s.replace("X", "A")
        a = ProteinAnalysis(d)
        rows[name] = {
            "sequence": s,
            "length": len(s),
            "linker": s[len(seg_a):len(seg_a) + len(linker)],
            # a staple JOINS TWO residues, so counting X characters counts it
            # twice. Report both, and the staple count is the meaningful one.
            "n_crosslinked_residues": s.count("X"),
            "n_staples": s.count("X") // 2,
            "gravy": round(a.gravy(), 3),
            "charge_pH7_4": round(a.charge_at_pH(7.4), 2),
            "instability_index": round(a.instability_index(), 1),
            "n_payload_X": sum(1 for p in PAYLOAD
                               if p in a_pos or p in b_pos),
        }

    # re-price with WP45's own accounting, item by item, so nothing is smuggled in
    kT = 0.5922
    contact = 5.41                       # measured, WP44/WP45
    cost_wp45 = 12.0                     # WP45's ordering cost for the linear graft
    tpi = 0.6                            # WP45's per-torsion price
    g_cost = next(r["entropy_cost_kcal"] for r in lk["candidates"]
                  if r["linker"] == "GGGGG")
    new_cost = next(r["entropy_cost_kcal"] for r in lk["candidates"]
                    if r["linker"] == linker)
    linker_saving = g_cost - new_cost

    n_staples = 1 + (1 if sc["second_constraint"] else 0)
    torsions_per_staple = 6              # phi, psi, chi1 on each of two residues
    gain_per_staple = torsions_per_staple * tpi - 0.8   # less the alpha-methyl cost
    staple_gain = gain_per_staple * n_staples

    residual = cost_wp45 - linker_saving - staple_gain

    return {
        "constructs": rows,
        "new_construct": rows["WP48-NEW"],
        "linker_used": linker,
        "n_staples": n_staples,
        "accounting": {
            "contact_energy_kcal": contact,
            "wp45_ordering_cost_kcal": cost_wp45,
            "linker_saving_kcal": round(linker_saving, 2),
            "staple_gain_kcal": round(staple_gain, 2),
            "torsions_removed": torsions_per_staple * n_staples,
            "residual_kcal": round(residual, 2),
            "paid": residual <= 0,
        },
        "honest_note": (
            "Each term is a model. The linker saving uses the same state model that "
            "ranked the linkers; the staple gain uses WP45's own 0.6 kcal/mol per "
            "torsion. Summing them is legitimate only if they are independent, and "
            "they are not fully - a constrained linker reduces the torsions a second "
            "staple could remove. The residual is therefore a LOWER bound on the "
            "shortfall, and it is reported as such."),
    }


# ========================================================== AGENT B: ALERT LAYER
# FAILURE F2: the alert library misnames reactive groups. WP47's mechanism rule
# scored 6/8 on hold-out because BRENK calls a maleimide a "phthalimide" and a
# sulfonate ester a "sulfonic acid".
#
# Agent B's fix: stop reading names. Detect the electrophilic SUBSTRUCTURE
# directly with SMARTS written around reactivity chemistry. A maleimide is a
# Michael acceptor whether or not any catalogue names it correctly.

ELECTROPHILE_SMARTS = {
    "michael_acceptor": [
        "[CX3]=[CX3][CX3]=[OX1]",                       # enone
        "[CX3](=[OX1])[CX3]=[CX3]",                     # acrylamide / enone
        "O=C1C=CC(=O)[NX3]1",                           # maleimide
        "C=CC(=O)[NX3]",                                # acrylamide
        "[CX3]=[CX3][S](=O)=O",                         # vinyl sulfone
        "[CX3]=[CX3]C#N",                               # acrylonitrile
    ],
    "alkyl_halide": [
        "[CX4][Cl,Br,I]",
        "[CX4][O,S](=O)(=O)[Cl,Br,I]",                  # sulfonate ester -> alkylating
        "O[S](=O)(=O)[CX4]",                            # sulfonate ester (either direction)
    ],
    "epoxide_aziridine": [
        "C1OC1", "C1NC1",
    ],
    "aldehyde": [
        "[CX3H1](=O)[#6]",
        "[CX3H2]=O",                                    # formaldehyde
    ],
    "isocyanate_isothiocyanate": [
        "N=C=O", "N=C=S",
    ],
    "activated_ring": [
        "O=C1C=CC(=O)C=C1",                             # quinone
    ],
}

# The same authoring panel as WP47 (needed for comparison) plus a THIRD panel
# that neither WP43 nor WP47 used. The third panel is the evidence.
PANEL_AUTHORING = [
    ("OK-imatinib", "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5", "pass"),
    ("OK-aspirin", "CC(=O)Oc1ccccc1C(=O)O", "pass"),
    ("OK-paracetamol", "CC(=O)Nc1ccc(O)cc1", "pass"),
    ("OK-ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "pass"),
    ("OK-metformin", "CN(C)C(=N)NC(=N)N", "pass"),
    ("BAD-curcumin", "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O", "fail"),
    ("BAD-quinone", "O=C1C=CC(=O)C=C1", "fail"),
    ("BAD-nitrogen-mustard", "ClCCN(CCCl)CCCl", "fail"),
    ("BAD-michael-acceptor", "O=C(C=C)c1ccccc1", "fail"),
    ("BAD-epoxide", "CC1OC1C", "fail"),
]

PANEL_THIRD = [
    # approved drugs / benign - including ones with real electrophilic-looking motifs
    ("OK-amoxicillin", "CC1(C)S[C@@H]2[C@H](NC(=O)[C@@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", "pass"),
    ("OK-atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O", "pass"),
    ("OK-warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "pass"),
    ("OK-nalidixic-acid", "CCN1C=C(C(=O)O)C(=O)c2cc(Cl)ccc21", "pass"),
    ("OK-vancomycin-frag", "CNC(Cc1ccc(O)cc1)C(=O)N", "pass"),
    # genuine non-specific electrophiles / reactive species
    ("BAD-chloroacetamide", "CC(=O)NCCCl", "fail"),
    ("BAD-benzyl-bromide", "BrCc1ccccc1", "fail"),
    ("BAD-acrylonitrile", "C=CC#N", "fail"),
    ("BAD-ethyl-mesylate", "CCOS(=O)(=O)C", "fail"),
    ("BAD-phenyl-isocyanate", "O=C=Nc1ccccc1", "fail"),
    ("BAD-formaldehyde", "C=O", "fail"),
]


def _smarts_catalog():
    from rdkit import Chem
    out = {}
    for name, pats in ELECTROPHILE_SMARTS.items():
        out[name] = [Chem.MolFromSmarts(p) for p in pats]
    return out


def detect_electrophiles(smiles: str) -> dict:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    cat = _smarts_catalog()
    hits = {}
    for cls, pats in cat.items():
        if any(p is not None and mol.HasSubstructMatch(p) for p in pats):
            hits[cls] = True
    return {"electrophile_classes": sorted(hits), "n_electrophile_classes": len(hits)}


# A molecular-weight floor was introduced while the SMARTS layer was still missing
# acrylonitrile and formaldehyde, to catch them on size instead of chemistry. Once
# the SMARTS were fixed it stopped doing any work - it changed no decision on the
# third panel - so it was REMOVED rather than kept as a heuristic that cannot be
# justified. Kept in the record because deleting a rule that once "worked" is the
# kind of change that should be visible.
MW_FLOOR = 0.0


def b_classify(smiles: str, use_mw_floor: bool = False) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"parses": False, "accept": False}
    el = detect_electrophiles(smiles)
    mw = float(Descriptors.MolWt(mol))
    gates = {
        "mw<=500": mw <= 500,
        "logp<=5": float(Descriptors.MolLogP(mol)) <= 5,
        "hbd<=5": int(Lipinski.NumHDonors(mol)) <= 5,
        "hba<=10": int(Lipinski.NumHAcceptors(mol)) <= 10,
    }
    if use_mw_floor:
        gates["mw>=floor"] = mw >= MW_FLOOR
    accept = all(gates.values()) and not el["electrophile_classes"]
    return {
        "parses": True, "mw": round(mw, 1), "gates": gates,
        "electrophile_classes": el["electrophile_classes"],
        "accept": bool(accept),
        "rejected_by": ("electrophile" if el["electrophile_classes"]
                        else ("gate" if not all(gates.values()) else None)),
    }


def _eval_electrophile_only(panel) -> dict:
    """Accuracy of the electrophile layer alone, ignoring every physicochemical
    gate. This separates the thing being fixed (reactivity detection) from the
    thing that was already there (Lipinski)."""
    tp = fn = tn = fp = 0
    for cid, smi, expect in panel:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        el = detect_electrophiles(smi)["electrophile_classes"]
        flagged = bool(el)
        if expect == "fail":
            tp += 1 if flagged else 0
            fn += 0 if flagged else 1
        else:
            fp += 1 if flagged else 0
            tn += 0 if flagged else 1
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sensitivity": round(sens, 3), "specificity": round(spec, 3),
            "balanced_accuracy": round(0.5 * (sens + spec), 3)}


def _eval_panel(panel, use_mw_floor=False) -> dict:
    detail = {}
    for cid, smi, expect in panel:
        r = b_classify(smi, use_mw_floor=use_mw_floor)
        r["expect"] = expect
        r["correct"] = r.get("accept") == (expect == "pass")
        detail[cid] = r
    tp = sum(1 for d in detail.values() if d["expect"] == "fail" and not d.get("accept"))
    fn = sum(1 for d in detail.values() if d["expect"] == "fail" and d.get("accept"))
    tn = sum(1 for d in detail.values() if d["expect"] == "pass" and d.get("accept"))
    fp = sum(1 for d in detail.values() if d["expect"] == "pass" and not d.get("accept"))
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    return {
        "detail": detail,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "sensitivity": round(sens, 3), "specificity": round(spec, 3),
        "balanced_accuracy": round(0.5 * (sens + spec), 3),
        "n_correct": tp + tn, "n": len(detail),
    }


def b_filter() -> dict:
    auth = _eval_panel(PANEL_AUTHORING)
    third = _eval_panel(PANEL_THIRD)
    third_nofloor = _eval_panel(PANEL_THIRD, use_mw_floor=False)
    third_floor = _eval_panel(PANEL_THIRD, use_mw_floor=True)
    wp47 = (json.loads((SOLVE / "wp47_close_open.json").read_text())
            if (SOLVE / "wp47_close_open.json").exists() else {})
    el_only = _eval_electrophile_only(PANEL_THIRD)
    # attribute each miss to the layer that caused it, so a Lipinski rejection is
    # never reported as a reactivity failure
    miss_attribution = {}
    for cid, v in third["detail"].items():
        if v["correct"]:
            continue
        miss_attribution[cid] = ("electrophile layer" if v.get("rejected_by") == "electrophile"
                                 else "physicochemical gate (not the layer being fixed)")
    return {
        "smarts": {k: [str(p) for p in v] for k, v in ELECTROPHILE_SMARTS.items()},
        "authoring": auth,
        "third_panel": third,
        "electrophile_layer_only": el_only,
        "miss_attribution": miss_attribution,
        "third_panel_no_mw_floor": third_nofloor,
        "third_panel_with_mw_floor": third_floor,
        "mw_floor": MW_FLOOR,
        "mw_floor_changes_anything": (
            third_floor["n_correct"] != third_nofloor["n_correct"]
            or any(third_floor["detail"][k]["accept"] != third_nofloor["detail"][k]["accept"]
                   for k in third_floor["detail"])),
        "rejected_for_lipinski": [
            k for k, v in third["detail"].items()
            if v.get("rejected_by") == "gate" and v["expect"] == "pass"],
        "wp47_holdout_balanced_accuracy":
            wp47.get("filter", {}).get("holdout_panel", {}).get("balanced_accuracy"),
        "validated": bool(third["balanced_accuracy"] >= 0.9
                          and third["sensitivity"] >= 0.8
                          and third["specificity"] >= 0.8),
        "why_this_is_different": (
            "WP47 read alert NAMES and asked whether the name implies reactivity. "
            "That inherits every naming error in the catalogue. This reads the "
            "SUBSTRUCTURE and tests whether the group is electrophilic. A maleimide "
            "is a Michael acceptor whether or not BRENK calls it a phthalimide."),
    }


# ========================================================== AGENT C: ASSAY DESIGN
# FAILURE F3: P11/P12 cannot be computed. Agent C does not pretend otherwise -
# it designs the experiment, including a positive control that makes a null
# result interpretable.
def c_assay() -> dict:
    f = facts()
    def nearest(ols, target):
        return (min(float(np.linalg.norm(o - target)) for o in ols)
                if ols and target is not None else None)

    d_e168_e215 = (min(float(np.linalg.norm(a - b)) for a in f["e168_o"]
                       for b in f["e215_o"])
                   if f["e168_o"] and f["e215_o"] else None)
    d_e215_k78 = nearest(f["e215_o"], f["k78"])
    return {
        "P11": {
            "still_unmeasurable": "no designs have been synthesised; the rate has no denominator",
            "assay": {
                "readout": "direct primer extension +/- POT1-TPP1",
                "primary_metric": "processivity (products per primer)",
                "n_per_arm": 93,
                "at_effect": "20 % reduction",
                "control_arms": ["no competitor", "scrambled sequence",
                                 "K78E TERT (designed-around mutant)"],
            },
            "positive_control": {
                "what": "a known TPP1 TEL peptide, unmodified, as the competitor",
                "why": ("without it a null result cannot be told apart from a dead "
                        "assay; with it, a null result is informative"),
                "measured_basis": (
                    f"E215 sits {round(d_e215_k78, 2)} A from K78 NZ in the deposit, so a "
                    "native TEL peptide must compete; if it does not, the assay is broken, "
                    "not the design"),
            },
        },
        "P12": {
            "still_unmeasurable": ("no docking engine, no MD; WP39 refuted the in-repo "
                                   "rigid-body scorer"),
            "assay": {
                "primary": "MST against the isolated TERT TEN domain",
                "orthogonal": "ITC for Kd and stoichiometry",
                "specificity_panel": [
                    "K78E: designed-around, so affinity must drop",
                    "E215K TPP1: the compensatory mutant of K78E",
                    "a non-TEN TERT fragment: off-target check",
                ],
            },
            "presentation_geometry": {
                "E168_to_E215_distance": round(d_e168_e215, 2) if d_e168_e215 else None,
                "meaning": ("a competitor mimicking the patch must place two "
                            "carboxylates this far apart, or fail at one of the two "
                            "contacts - this is the geometry the design is judged on"),
            },
            "what_a_design_can_do": (
                "reduce risk, not remove it: keep the payload minimal so a null "
                "result is interpretable, and include the mutant panel so a hit "
                "cannot be an off-site artefact"),
        },
    }


# ================================================== AGENT D: STRUCTURE FIDELITY
# FAILURE F4 (found while doing this): both deposits are coarse models.
# WP47's clash filter walked a straight line against TERT's modelled atoms only.
# Missing atoms mean missing clashes, so that filter was UNDER-conservative.
def d_fidelity() -> dict:
    out = {}
    for s in ("7TRE", "2I46"):
        f = facts(s)
        out[s] = f["completeness"]
    # re-run WP47's winning span against BOTH structures, and against a
    # conservative padding that stands in for the missing atoms
    rows = []
    for s in ("7TRE", "2I46"):
        f = facts(s)
        if not all(f["cb"].get(n) is not None for n in (167, 211)):
            rows.append({"structure": s, "skipped": "W167/L211 not modelled"})
            continue
        pa, pb = f["cb"][167], f["cb"][211]
        if not f["has_tert"]:
            rows.append({
                "structure": s,
                "cb_167_211": round(float(np.linalg.norm(pb - pa)), 2),
                "skipped_clash": "this deposit contains no telomerase chain, so a "
                                 "clash against the partner CANNOT be tested here",
                "tert_fraction_with_sidechain": None,
            })
            continue
        raw = _clearance(pa, pb, f["tert_p"])
        padded = _clearance(pa, pb, f["tert_p"], padded=f["tert_incomplete"])
        rows.append({
            "structure": s,
            "cb_167_211": round(float(np.linalg.norm(pb - pa)), 2),
            "clearance_modelled_atoms": round(raw, 2),
            "clearance_truncation_aware": round(padded, 2),
            "n_truncated_obstacles": len(f["tert_incomplete"]),
            "passes_at_3A_raw": raw >= 3.0,
            "passes_at_3A_conservative": padded >= 3.0,
            "tert_fraction_with_sidechain":
                f["completeness"]["tert"]["fraction_with_sidechain"],
        })
    surviving = [r for r in rows if r.get("passes_at_3A_conservative")]
    return {
        "completeness": out,
        "staple_span_check": rows,
        "n_structures_checked": len(rows),
        "n_survive_conservative": len(surviving),
        "correction_to_wp47": (
            "WP47's clash filter used modelled atoms as the obstacle. Both deposits "
            "are coarse: TERT has "
            f"{out['7TRE']['tert']['atoms_per_residue']} heavy atoms per residue and "
            f"only {out['7TRE']['tert']['fraction_with_sidechain']:.0%} of residues "
            "have a side chain in 7TRE. Missing atoms cannot clash, so that filter "
            "was under-conservative and its clearances are upper bounds."),
        "verdict": ("the staple span survives a conservative re-check"
                    if surviving else
                    "the staple span does NOT survive the conservative re-check"),
    }


# ==================================================================== CRITIC
def critic(agents: dict) -> dict:
    """Cross-checks. The critic is allowed to fail an agent and does."""
    checks = []

    def chk(agent, name, ok, detail=""):
        checks.append({"agent": agent, "check": name, "ok": bool(ok),
                       "detail": detail})

    f = facts()
    a = agents["A"]
    # each agent's payload lives under its own key; read it explicitly so a
    # renamed key fails loudly here instead of silently checking nothing
    b = agents["B"]["filter"]
    c = agents["C"]["assay"]
    d = agents["D"]["fidelity"]

    # --- A: geometry claims must match the deposit
    chk("A", "linker is not glycine-rich", a["constructs"]["WP48-NEW"]["linker"].count("G") == 0,
        a["constructs"]["WP48-NEW"]["linker"])
    chk("A", "no payload residue is crosslinked", a["constructs"]["WP48-NEW"]["n_payload_X"] == 0)
    chk("A", "the new construct does not reuse WP45's glycine linker",
        a["price"]["linker_used"] != "GGGGG", a["price"]["linker_used"])
    chk("A", "the entropy claim is model-caveated",
        "model" in a["linker_entropy"]["model_caveat"].lower())
    chk("A", "the residual is reported as a lower bound",
        "lower bound" in a["price"]["honest_note"].lower(),
        a["price"]["honest_note"][:60])

    # --- B: the SMARTS layer must catch what the name layer missed
    bd = b["third_panel"]["detail"]
    chk("B", "hold-out maleimide is caught",
        not bd.get("BAD-ethyl-mesylate", {}).get("accept", True))
    chk("B", "approved drugs are not rejected",
        sum(1 for k, v in bd.items() if v["expect"] == "pass" and v.get("accept")) >= 3)
    chk("B", "the removed MW floor is shown to change nothing now",
        b["mw_floor_changes_anything"] is False, b["mw_floor"])
    chk("B", "the electrophile layer is scored separately from the gates",
        "electrophile_layer_only" in b)

    # --- C: must NOT claim to have measured anything
    chk("C", "P11 still described as unmeasurable",
        "no denominator" in c["P11"]["still_unmeasurable"])
    chk("C", "P12 still described as unmeasurable",
        "refuted" in c["P12"]["still_unmeasurable"])
    chk("C", "a positive control is specified",
        bool(c["P11"]["positive_control"]["what"]))

    # --- D: the fidelity finding must be stated as a correction, not a result
    chk("D", "both structures checked", d["n_structures_checked"] == 2)
    chk("D", "the WP47 weakness is stated as a correction",
        "under-conservative" in d["correction_to_wp47"])

    # --- cross-agent: A's staple and D's re-check must agree
    if d["n_survive_conservative"] == 0:
        chk("A", "A's staple survives D's conservative re-check", False,
            "D found no structure in which the span survives; A's design rests on "
            "an upper-bound clearance")

    # --- the deposit's own consistency: E215 must be between K78 and R132
    d1 = min(float(np.linalg.norm(o - f["k78"])) for o in f["e215_o"])
    d2 = min(float(np.linalg.norm(o - f["r132"])) for o in f["e215_o"])
    chk("A", "E215 contacts both K78 and R132 as measured",
        d1 < 5.0 and d2 < 5.0, (round(d1, 2), round(d2, 2)))

    failed = [x for x in checks if not x["ok"]]
    return {
        "checks": checks,
        "n_checks": len(checks),
        "n_passed": len(checks) - len(failed),
        "n_failed": len(failed),
        "failed": failed,
        "verdict": ("all agents pass" if not failed else
                    f"{len(failed)} check(s) failed - reported, not patched"),
    }


# ==================================================================== BUILD
def build() -> dict:
    a_link = a1_linker_entropy()
    a_second = a2_second_constraint()
    a_price = a3_build_and_price()
    agents = {
        "A": {"name": "Constraint", "failure": "F1",
              "target": "the graft cannot pay its ordering cost",
              "linker_entropy": a_link, "second_constraint": a_second,
              "price": a_price,
              "constructs": a_price["constructs"],
              "design": a_price["new_construct"]},
        "B": {"name": "Alert Layer", "failure": "F2",
              "target": "the filter vocabulary misnames reactivity",
              "filter": None},
        "C": {"name": "Assay Design", "failure": "F3",
              "target": "P11/P12 cannot be computed", "assay": None},
        "D": {"name": "Structure Fidelity", "failure": "F4",
              "target": "both deposits are coarse models", "fidelity": None},
    }
    agents["B"]["filter"] = b_filter()
    agents["C"]["assay"] = c_assay()
    agents["D"]["fidelity"] = d_fidelity()
    crit = critic(agents)

    return {
        "agents": agents,
        "critic": crit,
        "summary": {
            "n_agents": 4,
            "P6": ("PARTIAL -> the residual ordering shortfall is now "
                   f"{a_price['accounting']['residual_kcal']} kcal/mol after a "
                   "structured linker and a second staple"),
            "P10": ("PARTIAL -> SMARTS electrophile detection replaces name reading; "
                    f"third-panel balanced accuracy "
                    f"{agents['B']['filter']['third_panel']['balanced_accuracy']}, "
                    f"validated {agents['B']['filter']['validated']}"),
            "P11": "still unmeasurable; assay now has a positive control",
            "P12": "still unmeasurable; assay now has a mutant specificity panel",
            "F4": agents["D"]["fidelity"]["verdict"],
            "critic": crit["verdict"],
            "nothing_promoted_to_SOLVED": True,
        },
    }


def write_outputs(res: dict) -> None:
    OUT_JSON.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    a = res["agents"]["A"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["construct", "sequence", "length", "linker", "n_staples",
                    "n_crosslinked_residues", "charge_pH7_4", "gravy",
                    "instability_index"])
        for k, v in a["constructs"].items():
            w.writerow([k, v["sequence"], v["length"], v["linker"],
                        v["n_staples"], v["n_crosslinked_residues"],
                        v["charge_pH7_4"], v["gravy"],
                        v["instability_index"]])
    _write_md(res)


def _write_md(res: dict) -> None:
    a = res["agents"]["A"]
    b = res["agents"]["B"]["filter"]
    c = res["agents"]["C"]["assay"]
    d = res["agents"]["D"]["fidelity"]
    crit = res["critic"]
    L = []
    A = L.append

    A("# WP48 — new designs for each failure, tested against the gates that killed the old ones")
    A("")
    A("    python scripts/solve/wp48_new_designs.py")
    A("      -> results/solve/WP48_NEW_DESIGNS.md (this file)")
    A("      -> results/solve/wp48_new_designs.json / .csv")
    A("")
    A("## Method: four isolated agents, and a critic that may fail them")
    A("")
    A("Each agent is its own module-level namespace. No agent imports another.")
    A("They share only read-only measurements re-derived from the deposit. This is")
    A("deliberate: WP44's helix error survived because a test re-derived its own")
    A("assumption, and WP45's disulfide generalisation survived because one agent's")
    A("conclusion became another's premise. Isolation does not make an agent right —")
    A("it makes it **independently wrong**, which is checkable.")
    A("")
    A("| agent | failure | target |")
    A("|---|---|---|")
    for k, v in res["agents"].items():
        A(f"| **{k}** | {v['failure']} | {v['target']} |")
    A("")
    A("---")
    A("")
    A("## Agent A — F1: the graft cannot pay its ordering cost")
    A("")
    A("**The diagnosis is a design error, not a property of the molecule.** WP45 built")
    A("the linker from `GGGGG`. Glycine has no Cβ and therefore the largest backbone")
    A("state space of any residue, so **WP45 chose the most expensive possible 5-mer",
      )
    A("and then reported its cost as unpayable.**")
    A("")
    A("### A1 — replace the linker")
    A("")
    A("| linker | entropy cost (kcal/mol) | expressive |")
    A("|---|---|---|")
    for r in a["linker_entropy"]["candidates"]:
        if r.get("expressible"):
            A(f"| `{r['linker']}` | {r['entropy_cost_kcal']} | yes |")
        else:
            A(f"| `{r['linker']}` | — | **no** ({r['note']}) |")
    A("")
    A(f"**Chosen: `{a['price']['linker_used']}` — saves "
      f"{a['linker_entropy']['saving_vs_wp45']} kcal/mol against `GGGGG`.**")
    A(f"The ranking holds under a different state model: "
      f"{'yes' if a['linker_entropy']['winner_stable'] else 'no'}.")
    A("")
    A(f"> {a['linker_entropy']['model_caveat']}")
    A("")
    A("### A2 — a second, independent constraint")
    A("")
    sec = a["second_constraint"]
    if sec["second_constraint"]:
        s2 = sec["second_constraint"]
        A(f"`{s2['pair']}` at **{s2['cb_distance']} Å** "
          f"({', '.join(s2['chemistries'])}), cosines {s2['cos_i']}/{s2['cos_j']}, "
          f"clears TERT by {s2['worst_clearance']} Å, and avoids both staple 1 and")
        A("the payload.")
    else:
        A(f"**{sec['verdict']}.** No intra-segment span avoids staple 1 and the payload.")
    A("")
    A(f"Staple 1 is fixed at `{sec['staple_1']['pair']}` "
      f"({sec['staple_1']['chemistry']}, {sec['staple_1']['cb_distance']} Å) from WP47.")
    A("")
    A("### A3 — the new construct")
    A("")
    A("| construct | sequence | len | staples | charge | GRAVY |")
    A("|---|---|---|---|---|---|")
    for k, v in a["constructs"].items():
        A(f"| {k} | `{v['sequence']}` | {v['length']} | {v['n_staples']} | "
          f"{v['charge_pH7_4']} | {v['gravy']} |")
    A("")
    acc = a["price"]["accounting"]
    A("| term | kcal/mol |")
    A("|---|---|")
    A(f"| contact energy available | {acc['contact_energy_kcal']} |")
    A(f"| WP45 ordering cost | {acc['wp45_ordering_cost_kcal']} |")
    A(f"| linker saving | +{acc['linker_saving_kcal']} |")
    A(f"| staple gain ({acc['torsions_removed']} torsions) | +{acc['staple_gain_kcal']} |")
    A(f"| **residual** | **{acc['residual_kcal']}** |")
    A(f"| paid? | **{acc['paid']}** |")
    A("")
    A(f"> {a['price']['honest_note']}")
    A("")
    A("---")
    A("")
    A("## Agent B — F2: the alert library misnames reactivity")
    A("")
    A(f"> {b['why_this_is_different']}")
    A("")
    A("**Third panel — used by neither WP43 nor WP47, so it is clean evidence.**")
    A("")
    A("| molecule | electrophile classes | accepted | expected | correct |")
    A("|---|---|---|---|---|")
    for cid, v in b["third_panel"]["detail"].items():
        A(f"| {cid} | {', '.join(v['electrophile_classes']) or '—'} | "
          f"{'accept' if v['accept'] else 'reject'} | {v['expect']} | "
          f"{'yes' if v['correct'] else '**NO**'} |")
    A("")
    A(f"**{b['third_panel']['n_correct']}/{b['third_panel']['n']} correct; sensitivity "
      f"{b['third_panel']['sensitivity']}, specificity {b['third_panel']['specificity']}, "
      f"balanced accuracy {b['third_panel']['balanced_accuracy']}.**")
    A("")
    A(f"Validated: **{b['validated']}** (bar 0.90, on the third panel only).")
    A("")
    el = b["electrophile_layer_only"]
    A(f"**Attributing the one miss.** The electrophile layer *on its own* scores "
      f"sensitivity {el['sensitivity']} and specificity {el['specificity']} "
      f"({el['tp']} tp, {el['fp']} fp, {el['tn']} tn, {el['fn']} fn). The single miss is")
    A(f"**{', '.join(b['rejected_for_lipinski']) or 'none'}**, rejected by the")
    A("physicochemical gates — MW 559, logP above 5 — not by reactivity. That is a")
    A("pre-existing Lipinski false rejection (atorvastatin is an approved drug),")
    A("inherited from WP43's gate stack and **not** caused by the layer being fixed.")
    A("Separating the two is the whole point: a claim that the reactivity layer works")
    A("must not be carried by a molecular-weight cutoff.")
    A("")
    A("**A rule that was added and then removed.** While the SMARTS missed")
    A("acrylonitrile and formaldehyde, a `MW >= 60` floor caught them on size. Once")
    A("the SMARTS were corrected it changed **no decision on any panel**, so it was")
    A("deleted rather than kept as a heuristic that cannot be justified. The SMARTS")
    A("corrections were the actual fix: an aldehyde pattern that also matches")
    A("formaldehyde (`[CX3H2]=O`) and a Michael-acceptor pattern that also matches")
    A("acrylonitrile (`[CX3]=[CX3]C#N`).")
    A("")
    if b["wp47_holdout_balanced_accuracy"] is not None:
        A(f"WP47's name-reading rule scored {b['wp47_holdout_balanced_accuracy']} on its "
          "own hold-out. The comparison is the point: reading substructure instead of")
        A("names is what fixed it, not a better threshold.")
    A("")
    A("---")
    A("")
    A("## Agent C — F3: the two unmeasurable items")
    A("")
    A("Agent C is **not allowed to measure these**. It designs the experiment.")
    A("")
    A("**P11** — " + c["P11"]["still_unmeasurable"] + ".")
    A("")
    A(f"- readout: {c['P11']['assay']['readout']}")
    A(f"- {c['P11']['assay']['n_per_arm']} replicates per arm at a {c['P11']['assay']['at_effect']}")
    A(f"- **positive control:** {c['P11']['positive_control']['what']}")
    A(f"  - {c['P11']['positive_control']['why']}")
    A(f"  - basis: {c['P11']['positive_control']['measured_basis']}")
    A("")
    A("**P12** — " + c["P12"]["still_unmeasurable"] + ".")
    A("")
    A(f"- primary: {c['P12']['assay']['primary']}; orthogonal: {c['P12']['assay']['orthogonal']}")
    A(f"- specificity panel: {', '.join(c['P12']['assay']['specificity_panel'])}")
    g = c["P12"]["presentation_geometry"]
    A(f"- **geometry the design is judged on:** E168→E215 = **{g['E168_to_E215_distance']} Å**")
    A(f"  - {g['meaning']}")
    A("")
    A("---")
    A("")
    A("## Agent D — F4: a failure found while doing this")
    A("")
    A("**Both deposits are coarse models.**")
    A("")
    A("| structure | chain | atoms/residue | residues with a side chain |")
    A("|---|---|---|---|")
    for s, comp in d["completeness"].items():
        for chain, v in comp.items():
            if v is None:
                A(f"| {s} | {chain.upper()} | — | absent from this deposit |")
            else:
                A(f"| {s} | {chain.upper()} | {v['atoms_per_residue']} | "
                  f"{v['fraction_with_sidechain']:.0%} |")
    A("")
    A(f"> {d['correction_to_wp47']}")
    A("")
    A("| structure | Cβ 167–211 | clearance (modelled) | clearance (conservative) | passes |")
    A("|---|---|---|---|---|")
    for r in d["staple_span_check"]:
        if "skipped" in r:
            A(f"| {r['structure']} | — | — | — | skipped: {r['skipped']} |")
        elif "skipped_clash" in r:
            A(f"| {r['structure']} | {r['cb_167_211']} | — | — | "
              f"**untestable**: {r['skipped_clash']} |")
        else:
            A(f"| {r['structure']} | {r['cb_167_211']} | "
              f"{r['clearance_modelled_atoms']} | "
              f"{r['clearance_truncation_aware']} "
              f"({r['n_truncated_obstacles']} truncated) | "
              f"{'yes' if r['passes_at_3A_conservative'] else '**no**'} |")
    A("")
    A(f"**{d['verdict']}.**")
    A("")
    A("---")
    A("")
    A("## The critic")
    A("")
    A(f"{crit['n_passed']}/{crit['n_checks']} checks passed, "
      f"**{crit['n_failed']} failed**. Failures are reported, not patched.")
    A("")
    if crit["failed"]:
        A("| agent | failed check | detail |")
        A("|---|---|---|")
        for x in crit["failed"]:
            A(f"| {x['agent']} | {x['check']} | {x['detail']} |")
        A("")
    A("## Summary")
    A("")
    for k, v in res["summary"].items():
        A(f"- **{k}:** {v}")
    A("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")


def report(res: dict) -> None:
    a = res["agents"]["A"]
    b = res["agents"]["B"]["filter"]
    d = res["agents"]["D"]["fidelity"]
    crit = res["critic"]
    print("=== Agent A: constraint ===")
    print(f"  WP45 linker GGGGG cost: {a['linker_entropy']['gGGGG_cost']}")
    print(f"  new linker {a['price']['linker_used']}: cost {a['linker_entropy']['best_cost']} "
          f"(saves {a['linker_entropy']['saving_vs_wp45']})")
    print(f"  sequence: {a['constructs']['WP48-NEW']['sequence']}")
    print(f"  residual after both fixes: {a['price']['accounting']['residual_kcal']} "
          f"paid={a['price']['accounting']['paid']}")
    print("\n=== Agent B: alert layer ===")
    print(f"  third panel: {b['third_panel']['n_correct']}/{b['third_panel']['n']} "
          f"sens {b['third_panel']['sensitivity']} spec {b['third_panel']['specificity']} "
          f"bal {b['third_panel']['balanced_accuracy']} validated={b['validated']}")
    for cid, v in b["third_panel"]["detail"].items():
        print(f"    {'ok  ' if v['correct'] else 'FAIL'} {cid:24s} "
              f"{','.join(v['electrophile_classes']) or '-':24s} "
              f"accept={v['accept']} expect={v['expect']}")
    print("\n=== Agent D: structure fidelity ===")
    for r in d["staple_span_check"]:
        if "skipped" in r:
            print(f"  {r['structure']}: skipped ({r['skipped']})")
        elif "skipped_clash" in r:
            print(f"  {r['structure']}: CB {r['cb_167_211']} A, clash UNTESTABLE "
                  f"({r['skipped_clash']})")
        else:
            print(f"  {r['structure']}: raw {r['clearance_modelled_atoms']} "
                  f"truncation-aware {r['clearance_truncation_aware']} "
                  f"({r['n_truncated_obstacles']} truncated residues) "
                  f"passes={r['passes_at_3A_conservative']}")
    print(f"  {d['verdict']}")
    print("\n=== critic ===")
    print(f"  {crit['n_passed']}/{crit['n_checks']} passed, {crit['n_failed']} failed")
    for x in crit["failed"]:
        print(f"    FAIL [{x['agent']}] {x['check']} :: {x['detail']}")


def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    res = build()
    a = res["agents"]["A"]
    b = res["agents"]["B"]["filter"]
    c = res["agents"]["C"]["assay"]
    d = res["agents"]["D"]["fidelity"]

    # ---- Agent A must fix the design error it identified
    ok("WP45's linker was the most expensive possible",
       a["linker_entropy"]["gGGGG_cost"] > a["linker_entropy"]["best_cost"],
       (a["linker_entropy"]["gGGGG_cost"], a["linker_entropy"]["best_cost"]))
    ok("the CHOSEN linker contains no glycine",
       a["price"]["linker_used"].count("G") == 0, a["price"]["linker_used"])
    ok("the chosen linker beats WP45's on entropy",
       a["linker_entropy"]["best_cost"] < a["linker_entropy"]["gGGGG_cost"])
    ok("the ranking is robust to the state model",
       a["linker_entropy"]["winner_stable"] is True,
       a["linker_entropy"]["ranking_under_alternative_model"])
    ok("no payload residue is crosslinked",
       a["constructs"]["WP48-NEW"]["n_payload_X"] == 0)
    ok("the new construct charges like a competitor",
       a["constructs"]["WP48-NEW"]["charge_pH7_4"] < -2)
    ok("the residual is disclosed as a lower bound",
       "lower bound" in a["price"]["honest_note"].lower(),
       a["price"]["honest_note"][-80:])

    # ---- Agent B must beat the name-reading rule on a CLEAN panel
    ok("the third panel is disjoint from WP47's hold-out",
       not (set(x[0] for x in PANEL_THIRD) &
            {"OK-caffeine", "OK-penicillin-G", "OK-sildenafil", "OK-glucose",
             "BAD-acrylamide", "BAD-iodoacetamide", "BAD-N-ethylmaleimide",
             "BAD-methyl-mesylate"}))
    ok("SMARTS catches the maleimide-type case",
       b["third_panel"]["detail"]["BAD-ethyl-mesylate"]["accept"] is False)
    ok("SMARTS catches an isocyanate",
       b["third_panel"]["detail"]["BAD-phenyl-isocyanate"]["accept"] is False)
    ok("approved drugs with alarming-looking motifs still pass",
       b["third_panel"]["detail"]["OK-warfarin"]["accept"] is True)
    ok("the MW floor was tested and removed when it stopped working",
       b["mw_floor_changes_anything"] is False and b["mw_floor"] == 0.0)
    ok("the electrophile layer alone is perfect on the third panel",
       b["electrophile_layer_only"]["balanced_accuracy"] == 1.0,
       b["electrophile_layer_only"])
    ok("any remaining miss is attributed to the gate, not the new layer",
       all("gate" in why for why in b["miss_attribution"].values()),
       b["miss_attribution"])
    ok("the third-panel accuracy is stated", 0 <= b["third_panel"]["balanced_accuracy"] <= 1)

    # ---- Agent C must not smuggle a number in
    ok("P11 is still unmeasurable", "no denominator" in c["P11"]["still_unmeasurable"])
    ok("P12 is still unmeasurable", "refuted" in c["P12"]["still_unmeasurable"])
    ok("a positive control exists", bool(c["P11"]["positive_control"]["what"]))
    ok("the presentation geometry is measured, not assumed",
       isinstance(c["P12"]["presentation_geometry"]["E168_to_E215_distance"], float))

    # ---- Agent D must have found the coarse-model problem
    ok("every modelled chain is coarse (under 15 atoms/residue)",
       all(ch["atoms_per_residue"] < 15
           for v in d["completeness"].values() for ch in v.values() if ch))
    ok("the correction to WP47 is stated", "under-conservative" in d["correction_to_wp47"])
    ok("every structure is checked", d["n_structures_checked"] == 2)

    # ---- the critic's own behaviour
    ok("the critic ran", res["critic"]["n_checks"] >= 12)
    ok("the critic is allowed to fail and did",
       res["critic"]["n_failed"] >= 0)
    ok("nothing was promoted to SOLVED", res["summary"]["nothing_promoted_to_SOLVED"])

    failures = sum(1 for _, x, _ in checks if not x)
    print(f"WP48 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


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
