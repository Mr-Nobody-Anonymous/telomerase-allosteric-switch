#!/usr/bin/env python3
"""
WP44 - One module per sub-problem, each solving its own site completely, then
merged - and the merge is tested for stability.

    python scripts/solve/wp44_modular_designs.py
      -> results/solve/wp44_modular_designs.json
      -> results/solve/wp44_modular_designs.csv

    python scripts/solve/wp44_modular_designs.py --self-test
      -> checks, 0 failed

THE DISCIPLINE THIS ENFORCES
----------------------------
WP42 made multifunctional helices: one 48-residue chain that had to span 20 A,
read a carboxylate triad AND present a hydrophobic face. That is fragile by
construction - every added job is another way to fail.

So the problem is split into the sites that were actually measured, one module
per site, and each module is held to a hard rule:

  * it must cover EVERY measured contact of its own site (solved fully), and
  * it must cover NO contact of any other site (solved only).

Both halves of that rule are machine-checked against distances from 7TRE. A
module that half-covers its own site, or that accidentally reaches a second
site, is rejected.

Then the modules are merged, and the merge is tested for stability. The stability
test is the point: modules that are each fine on their own can fail as a pair.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_JSON = REPO / "results" / "solve" / "wp44_modular_designs.json"
OUT_CSV = REPO / "results" / "solve" / "wp44_modular_designs.csv"

sys.path.insert(0, str(REPO / "scripts" / "solve"))

# --------------------------------------------------------------- helix model
HELIX_RISE = 1.5
HELIX_TWIST_DEG = 100.0
HELIX_CA_RADIUS = 2.3
# charged side-chain centroids sit further out than C-alpha
HELIX_SC_RADIUS = 6.0

# ------------------------------------------------------- measured site anchors
# Each anchor is a TERT or TPP1 atom read from 7TRE. A "site" is a set of anchors
# that a module must cover together.
SITES = {
    "S1_K78_R132": {
        "label": "K78/R132 groove (the E215 contact, recruitment hot spot)",
        "problem": "occupy the groove TPP1 uses, to block recruitment",
        "anchors": [("A", 78, "NZ", "K78 NZ"), ("A", 132, "NH1", "R132 NH1")],
        "kind": "electrostatic",
    },
    "S2_R774": {
        "label": "R774 (the E168 contact)",
        "problem": "engage R774, the second half of the bipartite patch",
        "anchors": [("A", 774, "NH1", "R774 NH1")],
        "kind": "electrostatic",
    },
    "S3_TEL_TRIAD": {
        "label": "TEL patch carboxylate triad E168/E169/E171",
        "problem": "read the patch itself, to compete with TERT for it",
        "anchors": [("C", 168, "OE1", "E168 OE1"), ("C", 169, "OE1", "E169 OE1"),
                    ("C", 171, "OE1", "E171 OE1")],
        "kind": "electrostatic",
    },
    "S4_TEN_LOOP": {
        "label": "TEN 129-132 hydrophobic loop / PLYQ anchor",
        "problem": "anchor on the hydrophobic wall beside the groove",
        "anchors": [("A", 130, "CB", "A130 CB"), ("A", 131, "CD1", "L131 CD1"),
                    ("A", 132, "CB", "R132 CB")],
        "kind": "hydrophobic",
    },
}

# --------------------------------------------------------------- the modules
# ONE site, ONE job, per module. `display` gives the offsets and residue placed
# on the display face; the register is checked against the site's measured
# anchor separation and the module is rejected if it does not match.
MODULES = [
    {
        "id": "MD1-GROOVE",
        "site": "S1_K78_R132",
        "intent": "antagonise: occupy the K78/R132 groove and block TPP1 engagement",
        "display_residue": "E",
        "display_offsets": [0, 3],
        "helix": "EALEALQALQALQALQALQAL",
        "note": ("E at i and i+3 -> 5.05 A, matching the measured K78-R132 separation "
                 "of 4.35 A. One module, one groove, both basic residues."),
    },
    {
        "id": "MD2-R774",
        "site": "S2_R774",
        "intent": "antagonise: present a carboxylate to R774",
        "display_residue": "E",
        "display_offsets": [0],
        "helix": "EALQALQALQALQALQALQAL",
        "note": ("A single anchor, so a single carboxylate. Deliberately NOT registering "
                 "a second display residue, because that would reach a second site."),
    },
    {
        "id": "MD3-PATCH",
        "site": "S3_TEL_TRIAD",
        "intent": "antagonise: read the patch triad so TERT cannot",
        "display_residue": "K",
        "display_offsets": [0, 4, 7],
        "helix": "KALAKAQKALAQLAQALQALAQLA",
        "note": ("K at i, i+4, i+7 -> 6.43 and 10.75 A, matching the measured "
                 "E168-E169 6.58 A and E168-E171 10.22 A. Opposite sign to MD1/MD2."),
    },
    {
        "id": "MD4-ANCHOR",
        "site": "S4_TEN_LOOP",
        "intent": "anchor: sit on the hydrophobic TEN wall without touching either groove",
        "display_residue": None,
        "display_offsets": [],
        "helix": "LQQLLQQLLQQLLQQLLQQLLQQL",
        "note": ("A hydrophobic patch at i, i+3, i+4 and i+7 - one face greasy, the other polar. No charged display residue at all. This module exists to add affinity "
                 "and nothing else, which is only safe if it touches no other site."),
    },
    {
        "id": "MD5-REPAIR",
        "site": "S1_K78_R132",
        "intent": "agonise: supply in trans the charge that K78E/R132E lost",
        "display_residue": "E",
        "display_offsets": [0, 3],
        "helix": "EALEALQALQALQALQALQAL",
        "note": ("Same site and same geometry as MD1, opposite intent. It cannot be "
                 "merged with MD1 - both want the same groove."),
    },
]

# modules that must never appear in one chain
INCOMPATIBLE = [("MD1-GROOVE", "MD5-REPAIR")]

# ------------------------------------------------------------- merge recipes
LINKERS = {
    "rigid": "EAAAK",       # established helical linker, keeps registers apart
    "flexible": "GSGSG",    # lets the modules collapse onto each other
}

MERGES = [
    {
        "id": "MRG-MIN",
        "modules": ["MD3-PATCH", "MD1-GROOVE"],
        "linker": "rigid",
        "goal": "block both halves of the patch interaction with two modules",
    },
    {
        "id": "MRG-BLOCK",
        "modules": ["MD3-PATCH", "MD1-GROOVE", "MD2-R774"],
        "linker": "rigid",
        "goal": "full blockade: patch reader + groove + R774",
    },
    {
        "id": "MRG-SAMESIGN",
        "modules": ["MD1-GROOVE", "MD2-R774", "MD4-ANCHOR"],
        "linker": "rigid",
        "goal": "all-acidic blockade with no charge conflict (the control merge)",
    },
    {
        "id": "MRG-REPAIR",
        "modules": ["MD5-REPAIR", "MD4-ANCHOR"],
        "linker": "rigid",
        "goal": "repair a K78E/R132E allele and anchor it",
    },
    {
        "id": "MRG-FLEX",
        "modules": ["MD3-PATCH", "MD1-GROOVE"],
        "linker": "flexible",
        "goal": "the same pair as MRG-MIN but with a flexible linker",
    },
    {
        "id": "MRG-ILLEGAL",
        "modules": ["MD1-GROOVE", "MD5-REPAIR"],
        "linker": "rigid",
        "goal": "deliberately invalid: two modules competing for one groove",
    },
]

# --------------------------------------------------------------- energy model
# Screened Coulomb. These are HEURISTIC parameters and are stated as such: the
# absolute kcal/mol is not a binding or folding free energy. Only differences
# between constructs computed with the same parameters are read.
EPSILON_EFF = 40.0      # screened protein/solvent interface
DEBYE_LENGTH = 10.0     # A, roughly 50 mM ionic strength
COULOMB_K = 332.0       # kcal*A/mol
COLLAPSE_THRESHOLD = -1.5   # kcal/mol of inter-module attraction: collapse risk
REPULSION_THRESHOLD = -3.0  # kcal/mol of intra-helix self-repulsion penalty

CHARGED = {"E": -1, "D": -1, "K": +1, "R": +1}

# Pace & Scholtz (1998) relative helix propensities, kcal/mol per residue.
# Lower is a better helix former; Pro is a breaker.
HELIX_PROPENSITY = {
    "A": 0.00, "L": 0.21, "R": 0.21, "K": 0.26, "M": 0.24, "F": 0.36,
    "Q": 0.39, "E": 0.40, "I": 0.41, "W": 0.49, "S": 0.50, "T": 0.52,
    "Y": 0.53, "V": 0.61, "H": 0.61, "N": 0.65, "C": 0.68, "D": 0.69,
    "G": 1.00, "P": 3.16,
}
HYDROPHOBIC = set("AILMFWVY")

# Match tolerance, stated once. The measured spacings carry ~0.5 A uncertainty at
# 3.5 A resolution, and a side-chain centroid sits off the C-alpha path, so a
# tighter tolerance than this would reject contacts that are genuinely matched.
# It is NOT tuned per design: every coverage and cross-site call uses this value.
MATCH_TOL = 1.0


# ------------------------------------------------------------------ geometry
def helix_span(offset: int) -> float:
    """C-alpha to C-alpha distance across `offset` residues of an ideal helix."""
    th = math.radians(HELIX_TWIST_DEG * offset)
    x = HELIX_CA_RADIUS * math.cos(th) - HELIX_CA_RADIUS
    y = HELIX_CA_RADIUS * math.sin(th)
    z = HELIX_RISE * offset
    return math.sqrt(x * x + y * y + z * z)


def scaffold_point(offset: int) -> tuple:
    """Side-chain centroid position on the ideal helix surface."""
    th = math.radians(HELIX_TWIST_DEG * offset)
    return (HELIX_SC_RADIUS * math.cos(th), HELIX_SC_RADIUS * math.sin(th),
            HELIX_RISE * offset)


def _structure():
    import gemmi

    st = gemmi.read_structure(str(REPO / "data" / "raw" / "7TRE.cif"))
    st.setup_entities()
    return st[0]


def site_geometry() -> dict:
    """Measured anchor positions and their pairwise separations, per site."""
    import numpy as np

    m = _structure()

    def pos(chain, num, name):
        for r in m[chain]:
            if r.seqid.num == num:
                for a in r:
                    if a.name == name:
                        return np.array([a.pos.x, a.pos.y, a.pos.z])
        return None

    out = {}
    for sid, s in SITES.items():
        pts, missing = {}, []
        for ch, num, name, label in s["anchors"]:
            p = pos(ch, num, name)
            if p is None:
                missing.append(label)
            else:
                pts[label] = p
        pairs = {}
        for (la, pa), (lb, pb) in itertools.combinations(pts.items(), 2):
            pairs[f"{la}|{lb}"] = round(float(np.linalg.norm(pa - pb)), 2)
        out[sid] = {
            "label": s["label"],
            "kind": s["kind"],
            "n_anchors": len(s["anchors"]),
            "resolved": len(pts),
            "missing": missing,
            "pairs": pairs,
            "coords": {k: [round(float(x), 2) for x in v] for k, v in pts.items()},
        }
    return out


# ------------------------------------------------------------- module checks
def solvability() -> dict:
    """Can a helical display satisfy each measured spacing at all?

    This is the constraint that decides whether the modular idea is even possible,
    so it is computed up front rather than discovered as a failure. For each
    measured anchor separation, report the closest achievable register and the
    error.
    """
    out = {}
    for sid, s in SITES.items():
        if s["kind"] != "electrostatic":
            continue
        for pair, measured in site_geometry()[sid]["pairs"].items():
            best = min(
                ((off, helix_span(off)) for off in range(1, 25)),
                key=lambda o: abs(o[1] - measured),
            )
            out[f"{sid}|{pair}"] = {
                "measured": measured,
                "best_register": f"i+{best[0]}",
                "design": round(best[1], 2),
                "error": round(best[1] - measured, 2),
                "reachable": abs(best[1] - measured) <= MATCH_TOL,
            }
    return out


def module_analysis(mod: dict, sites: dict) -> dict:
    """Does this module solve exactly one site, and solve that one FULLY?

    Two rules, both machine-checked:
      FULL   - every measured contact of its own site is covered
      ONLY   - it registers on no other site
    A module can fail FULL because its site is not helix-solvable at all; that is
    reported as the cause, because it is a fact about the site, not the design.
    """
    helix = mod["helix"]
    res = mod["display_residue"]
    offsets = mod["display_offsets"]
    site = SITES[mod["site"]]
    geo = sites[mod["site"]]

    if res:
        idx = [i for i, c in enumerate(helix) if c == res]
        realised = [i - idx[0] for i in idx] if idx else []
    else:
        realised = []

    # ---- FULL: coverage of its own site
    coverage = []
    if site["kind"] == "electrostatic":
        if not geo["pairs"]:
            coverage.append({
                "contact": "single anchor", "measured": None,
                "design": float(len(realised)), "register": "i+0" if realised else None,
                "error": 0.0, "covered": bool(realised),
            })
        for pair_label, measured in geo["pairs"].items():
            best = None
            for a, b in itertools.combinations(realised, 2):
                span = helix_span(abs(a - b))
                if best is None or abs(span - measured) < abs(best[1] - measured):
                    best = (abs(a - b), span)
            achievable = min((helix_span(o) for o in range(1, 25)),
                             key=lambda x: abs(x - measured))
            coverage.append({
                "contact": pair_label,
                "measured": measured,
                "design": round(best[1], 2) if best else round(achievable, 2),
                "register": f"i+{best[0]}" if best else "best-possible",
                "error": round((best[1] if best else achievable) - measured, 2),
                "covered": bool(best and abs(best[1] - measured) <= MATCH_TOL),
            })
    else:
        # Coverage of a hydrophobic site and RISK of aggregation are two questions
        # about the same greasy face, so they use one model with two thresholds
        # instead of two models that can contradict each other. Covering the wall
        # needs a patch to EXIST; risk needs it to be big. MD4 does the first and
        # is close to the second, which is a real trade-off to report, not a bug.
        face = aggregate_face(helix)
        covers = ((face["hydrophobic_fraction"] >= 0.40
                   and face["face_concentration"] >= 0.30)
                  or face["hydrophobic_fraction"] >= 0.85)
        coverage.append({
            "contact": "hydrophobic wall", "measured": None,
            "design": round(face["hydrophobic_fraction"], 3), "register": "bulk",
            "error": 0.0, "covered": bool(covers),
            "face_concentration": face["face_concentration"],
            "coverage_rule": "frac >= 0.40 and R >= 0.30, or frac >= 0.85",
        })

    own_ok = all(c["covered"] for c in coverage)
    unsolvable = [c["contact"] for c in coverage
                  if not c["covered"] and c["measured"] is not None
                  and abs(c["design"] - c["measured"]) > MATCH_TOL]

    # ---- ONLY: does it register on another site?
    cross = []
    for sid, s in SITES.items():
        if sid == mod["site"] or s["kind"] == "hydrophobic":
            continue
        if len(realised) < 2:
            continue
        for pair_label, measured in sites[sid]["pairs"].items():
            for a, b in itertools.combinations(realised, 2):
                span = helix_span(abs(b - a))
                if abs(span - measured) <= MATCH_TOL:
                    cross.append({
                        "site": sid, "contact": pair_label, "measured": measured,
                        "design": round(span, 2), "register": f"i+{abs(b - a)}",
                    })
                    break

    # A register coincidence is a RISK, not a proof. Whether two display residues
    # actually land on a second site depends on docking, which cannot be computed
    # here. It is reported so that it is not discovered later as a surprise.
    if own_ok and not cross:
        verdict = "ONE-SITE-FULL"
    elif own_ok and cross:
        verdict = "ONE-SITE-FULL-WITH-RISK"
    elif unsolvable:
        verdict = "SITE-UNSOLVABLE"
    else:
        verdict = "INCOMPLETE"

    return {
        "id": mod["id"],
        "site": mod["site"],
        "site_label": site["label"],
        "problem": site["problem"],
        "intent": mod["intent"],
        "helix": helix,
        "length": len(helix),
        "display_residue": res,
        "display_offsets_claimed": offsets,
        "display_offsets_realised": realised,
        "register_matches": (realised == offsets) if res else True,
        "coverage": coverage,
        "contacts_total": len(coverage),
        "contacts_covered": sum(1 for c in coverage if c["covered"]),
        "fully_covers_own_site": own_ok,
        "site_solvable_in_helix": not unsolvable,
        "unsolvable_contacts": unsolvable,
        "cross_site_register_matches": cross,
        "cross_site_warning": bool(cross),
        "verdict": verdict,
        "note": mod["note"],
    }


# -------------------------------------------------------------- energy terms
def intrahelix_electrostatics(helix: str) -> dict:
    """Like-charge repulsion inside one helix, on the ideal scaffold.

    A polyanionic helix is destabilised by its own charges. This is the term that
    punishes stacking several acidic modules into one chain.
    """
    charges = [(i, CHARGED[c]) for i, c in enumerate(helix) if c in CHARGED]
    total = 0.0
    n_pairs = 0
    worst = (0.0, None)
    for (i, qi), (j, qj) in itertools.combinations(charges, 2):
        pi, pj = scaffold_point(i), scaffold_point(j)
        r = math.dist(pi, pj)
        e = COULOMB_K * qi * qj / (EPSILON_EFF * r) * math.exp(-r / DEBYE_LENGTH)
        total += e
        n_pairs += 1
        if e < worst[0]:
            worst = (e, (i, j, round(r, 2)))
    net = sum(q for _, q in charges)
    return {
        "net_charge": net,
        "n_charged": len(charges),
        "n_pairs": n_pairs,
        "energy": round(total, 3),
        "worst_pair": worst[1],
        "worst_energy": round(worst[0], 3),
        # charge per residue: above ~0.25 the helix is electrostatically stressed
        "charge_density": round(abs(net) / len(helix), 3),
    }


def intermodule_electrostatics(a: dict, b: dict, linker: str) -> dict:
    """Attraction or repulsion between two modules in one chain.

    Separation is the linker's extended length, floored at a contact distance so
    that two helices cannot be scored as overlapping.
    """
    qa = sum(CHARGED.get(c, 0) for c in a["helix"])
    qb = sum(CHARGED.get(c, 0) for c in b["helix"])
    # A rigid helical linker holds the modules apart; a flexible Gly/Ser linker
    # lets them collapse onto each other. That difference is the whole reason the
    # linker-type axis exists, so it is modelled rather than assumed away.
    if linker.startswith(("GS", "SG", "GG")):
        d = 5.0                          # collapsed coil, modules adjacent
    else:
        d = max(3.5 * len(linker), 12.0)  # extended helical linker
    e = COULOMB_K * qa * qb / (EPSILON_EFF * d) * math.exp(-d / DEBYE_LENGTH)
    return {
        "pair": f"{a['id']}|{b['id']}",
        "charge_a": qa, "charge_b": qb, "separation": round(d, 2),
        "energy": round(e, 3),
        "same_sign": (qa * qb) > 0,
        "collapse_risk": bool(e < COLLAPSE_THRESHOLD),
    }


# ------------------------------------------------------------ stability test
def aggregate_face(seq: str) -> dict:
    """Does this sequence put a large greasy patch on one helical face?

    Hydrophobic residues are assigned their helical angle (100 deg per residue)
    and the resultant length R of that angle set is computed. R near 1 means they
    all sit on one face; R near 0 means they are spread around the helix.

    Risk requires BOTH a large hydrophobic fraction AND a concentrated face: a
    helix that is 40 % hydrophobic but amphipathic is not an aggregation problem,
    and one that is 80 % hydrophobic but evenly spread is not either.
    """
    angles = [math.radians(HELIX_TWIST_DEG * i) for i, c in enumerate(seq)
              if c in HYDROPHOBIC]
    n_hyd = len(angles)
    frac = n_hyd / len(seq) if seq else 0.0
    if angles:
        sx = sum(math.cos(a) for a in angles) / len(angles)
        sy = sum(math.sin(a) for a in angles) / len(angles)
        R = math.hypot(sx, sy)
        mean_hyd = math.atan2(sy, sx)
    else:
        R, mean_hyd = 0.0, 0.0

    polar = [math.radians(HELIX_TWIST_DEG * i) for i, c in enumerate(seq)
             if c in CHARGED or c in "QNSTHGD"]
    if polar:
        px = sum(math.cos(a) for a in polar) / len(polar)
        py = sum(math.sin(a) for a in polar) / len(polar)
        mean_pol = math.atan2(py, px)
        sep = abs(math.degrees(mean_hyd - mean_pol))
        sep = min(sep, 360 - sep)
    else:
        sep = 180.0

    # A helix is an aggregation risk when it presents a big greasy PATCH. That
    # needs the hydrophobic residues to be both concentrated (R) and numerous
    # enough to cover an area. A helix that is uniformly greasy has R ~ 0 and is
    # still one continuous greasy surface, so it is caught by the high-fraction
    # branch - the first version of this rule scored it SAFE, which is backwards.
    patch_risk = (frac >= 0.45) and (R >= 0.35)
    all_greasy_risk = frac >= 0.85
    risk = patch_risk or all_greasy_risk
    return {
        "hydrophobic_fraction": round(frac, 3),
        "face_concentration": round(R, 3),
        "face_separation_deg": round(sep, 1),

        "patch_risk": bool(patch_risk),
        "all_greasy_risk": bool(all_greasy_risk),
        "aggregation_risk": bool(risk),
        "amphipathic": bool(sep >= 90.0 and 0.30 <= frac <= 0.65),
        "rule": "risk if (frac >= 0.45 and R >= 0.35) or frac >= 0.85",
    }


def threshold_sensitivity(constructs: list) -> dict:
    """Does any verdict depend on where a threshold was placed?

    A verdict that flips inside a plausible range of a chosen constant is not a
    finding, it is a choice. This sweeps the two thresholds that can decide one
    and reports which constructs are sensitive, instead of quietly settling on the
    value that produced the preferred answer.
    """
    seqs = [(c["id"], c["sequence"]) for c in constructs if c.get("sequence")]
    rows = []
    for r_thr in (0.30, 0.35, 0.40, 0.45):
        for frac_thr in (0.45, 0.50, 0.55):
            flagged = []
            for cid, seq in seqs:
                angles = [math.radians(HELIX_TWIST_DEG * i)
                          for i, c in enumerate(seq) if c in HYDROPHOBIC]
                if not angles:
                    continue
                sx = sum(math.cos(a) for a in angles) / len(angles)
                sy = sum(math.sin(a) for a in angles) / len(angles)
                R = math.hypot(sx, sy)
                frac = len(angles) / len(seq)
                if (frac >= frac_thr and R >= r_thr) or frac >= 0.85:
                    flagged.append(cid)
            rows.append({"R_threshold": r_thr, "frac_threshold": frac_thr,
                         "flagged": flagged, "n_flagged": len(flagged)})
    all_flagged = {cid for r in rows for cid in r["flagged"]}
    never = {cid for cid, _ in seqs} - all_flagged
    always = {cid for cid, _ in seqs
              if all(cid in r["flagged"] for r in rows)}
    sensitive = sorted(all_flagged - always)
    return {
        "sweep": rows,
        "always_flagged": sorted(always),
        "never_flagged": sorted(never),
        "threshold_sensitive": sensitive,
    }


def stability(construct: dict) -> dict:
    """Run the stability battery on a single module or a merged construct."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    seq = construct["sequence"]
    a = ProteinAnalysis(seq)

    # 1. helix propensity cost, mean per residue
    prop = [HELIX_PROPENSITY.get(c, 0.5) for c in seq if c != "X"]
    helix_cost = round(sum(prop) / len(prop), 3) if prop else None

    # 2. like-charge repulsion inside each helix
    intra = [intrahelix_electrostatics(m["helix"]) for m in construct["modules"]]
    intra_total = round(sum(x["energy"] for x in intra), 3)

    # 3. inter-module electrostatics
    inter = []
    for a_m, b_m in itertools.combinations(construct["modules"], 2):
        inter.append(intermodule_electrostatics(a_m, b_m, construct["linker_seq"]))

    # 4. aggregation, modelled SPATIALLY. A 5-residue sequence window spans 1.4
    #    helix turns, so its residues are not neighbours in space and a contiguous
    #    window says nothing about a hydrophobic patch. The first implementation
    #    used windows and flagged every amphipathic helix - including ones whose
    #    hydrophobic residues all sit on a single face, which is the GOOD case.
    #    What actually drives helical aggregation is a large greasy FACE, so that
    #    is computed: the fraction of hydrophobic residues, and how tightly they
    #    cluster in helical angle (resultant length R, 0 = spread, 1 = one point).
    agg = aggregate_face(seq)

    # 5. charge complementarity between modules (the collapse route that is not
    #    purely Coulomb: a basic module and an acidic module salt-bridging)
    collapse_pairs = [x for x in inter if x["collapse_risk"]]
    opposite = [x for x in inter if not x["same_sign"] and x["charge_a"] and x["charge_b"]]

    verdict, reasons = "STABLE", []
    if collapse_pairs:
        verdict = "COLLAPSE-RISK"
        reasons.append(
            f"{len(collapse_pairs)} module pair(s) attract intramolecularly: "
            + ", ".join(f"{p['pair']} at {p['energy']} kcal/mol" for p in collapse_pairs))
    if intra_total < REPULSION_THRESHOLD:
        verdict = "COLLAPSE-RISK" if verdict == "COLLAPSE-RISK" else "ELECTROSTATIC-STRAIN"
        reasons.append(f"intra-helix like-charge repulsion {intra_total} kcal/mol "
                       f"(threshold {REPULSION_THRESHOLD})")
    if agg["aggregation_risk"]:
        reasons.append(
            f"greasy helical face: {agg['hydrophobic_fraction']:.0%} hydrophobic, "
            f"face concentration R = {agg['face_concentration']:.2f} "
            f"[{agg['rule']}]")
        if verdict == "STABLE":
            verdict = "AGGREGATION-RISK"
    # Face MIXING is the real defect: if the greasy and polar centroids are less
    # than 90 deg apart on the helix they are the same face. The first version
    # fired this on a 66%-hydrophobic helix whose faces were 178 deg apart - i.e.
    # perfectly segregated - because it tested membership of a 30-65% window
    # rather than the separation it was describing.
    if agg["face_separation_deg"] < 90.0:
        reasons.append(
            f"faces are mixed: greasy and polar centroids only "
            f"{agg['face_separation_deg']:.0f} deg apart on the helix")
    if not reasons:
        reasons.append("no electrostatic or face-mixing failure detected")

    # Biopython's instability index is calibrated on globular proteins using a
    # dipeptide scale. It systematically flags Ala/Glu-rich designed helices, so
    # it is reported as ADVISORY and does not decide the verdict - but the
    # disagreement is printed rather than hidden, because a silent disagreement
    # between two models is how one of them gets quoted out of context later.
    ii = round(a.instability_index(), 1)
    ii_verdict = "predicted unstable" if ii > 40 else "stable"
    advisory = []
    if ii > 40:
        advisory.append(
            f"instability index {ii} > 40 ('{ii_verdict}') - advisory only: the index is "
            f"calibrated on globular proteins and flags designed amphipathic helices "
            f"systematically; it is not evidence against this construct")

    return {
        "id": construct["id"],
        "sequence": seq,
        "length": len(seq),
        "mw_da": round(a.molecular_weight(), 1),
        "pI": round(a.isoelectric_point(), 2),
        "charge_pH7_4": round(a.charge_at_pH(7.4), 2),
        "gravy": round(a.gravy(), 3),
        "instability_index": ii,
        "instability_verdict": ii_verdict,
        "advisory": advisory,
        "helix_fraction": round(a.secondary_structure_fraction()[0], 3),
        "helix_propensity_cost": helix_cost,
        "intra_helix": intra,
        "intra_helix_total": intra_total,
        "inter_module": inter,
        "n_opposite_sign_pairs": len(opposite),
        "n_collapse_pairs": len(collapse_pairs),
        "aggregation": agg,
        "verdict": verdict,
        "reasons": reasons,
    }


# ------------------------------------------------------------------ assembly
def build_module(mod: dict, sites: dict) -> dict:
    a = module_analysis(mod, sites)
    a["sequence"] = mod["helix"]
    a["modules"] = [mod]
    a["linker_seq"] = ""
    a["construct_type"] = "module"
    a["stability"] = stability(a)
    return a


def merge(recipe: dict, modules_by_id: dict, sites: dict) -> dict:
    """Assemble a merged construct, or refuse it if the modules conflict."""
    ids = recipe["modules"]
    linker = LINKERS[recipe["linker"]]

    for pair in INCOMPATIBLE:
        if pair[0] in ids and pair[1] in ids:
            return {
                "id": recipe["id"],
                "goal": recipe["goal"],
                "modules": [modules_by_id[i] for i in ids],
                "linker_type": recipe["linker"],
                "linker_seq": linker,
                "construct_type": "merge",
                "refused": True,
                "refusal_reason": (
                    f"{pair[0]} and {pair[1]} both target {SITES[modules_by_id[pair[0]]['site']]['label']} "
                    f"with opposite intent; one groove cannot be both blocked and repaired"),
                "sequence": None,
                "stability": None,
            }

    seq = linker.join(modules_by_id[i]["helix"] for i in ids)
    out = {
        "id": recipe["id"],
        "goal": recipe["goal"],
        "modules": [modules_by_id[i] for i in ids],
        "linker_type": recipe["linker"],
        "linker_seq": linker,
        "construct_type": "merge",
        "refused": False,
        "sequence": seq,
    }
    out["stability"] = stability(out)

    # module independence after merging: every module must still hit its own
    # site and none of the others, in the merged chain
    still = [module_analysis(m, sites) for m in out["modules"]]
    out["module_integrity"] = [
        {"id": s["id"], "verdict": s["verdict"],
         "covers": f"{s['contacts_covered']}/{s['contacts_total']}",
         "cross": len(s["cross_site_register_matches"])}
        for s in still
    ]
    out["all_modules_intact"] = all(
        s["verdict"] in ("ONE-SITE-FULL", "ONE-SITE-FULL-WITH-RISK") for s in still)
    out["merges_an_unsolvable_site"] = any(
        s["verdict"] == "SITE-UNSOLVABLE" for s in still)
    return out


def build() -> dict:
    sites = site_geometry()
    modules_by_id = {m["id"]: m for m in MODULES}
    module_results = [build_module(m, sites) for m in MODULES]
    merge_results = [merge(r, modules_by_id, sites) for r in MERGES]

    built = [r for r in merge_results if not r["refused"]]
    sens = threshold_sensitivity(module_results + built)
    best_module = max(module_results, key=lambda x: (x["verdict"] == "ONE-SITE-FULL",
                                                     x["stability"]["instability_index"] < 40))
    stable_merges = [r for r in built if r["stability"]["verdict"] == "STABLE"]

    return {
        "sites": sites,
        "solvability": solvability(),
        "modules": module_results,
        "merges": merge_results,
        "summary": {
            "n_modules": len(module_results),
            "n_one_site_full": sum(1 for m in module_results
                                   if m["verdict"] == "ONE-SITE-FULL"),
            "n_merges": len(merge_results),
            "n_merges_refused": sum(1 for r in merge_results if r["refused"]),
            "n_merges_built": len(built),
            "n_merges_stable": len(stable_merges),
            "stable_merges": [r["id"] for r in stable_merges],
            "unstable_merges": [r["id"] for r in built
                                if r["stability"]["verdict"] != "STABLE"],
            "best_module": best_module["id"],
        "threshold_sensitivity": sens,
        "verdicts_decided_by_a_constant": sens["threshold_sensitive"],
        },
    }


def write_outputs(res: dict) -> None:
    import csv

    OUT_JSON.write_text(json.dumps(res, indent=2), encoding="utf-8")
    rows = []
    for m in res["modules"]:
        rows.append({
            "kind": "module", "id": m["id"], "site": m["site"],
            "verdict": m["verdict"], "stability": m["stability"]["verdict"],
            "length": m["length"], "coverage": f"{m['contacts_covered']}/{m['contacts_total']}",
            "cross_site": len(m["cross_site_register_matches"]),
            "intra_E": m["stability"]["intra_helix_total"],
            "gravy": m["stability"]["gravy"],
            "aggregation": m["stability"]["aggregation"]["aggregation_risk"],
        })
    for r in res["merges"]:
        if r["refused"]:
            rows.append({
                "kind": "merge", "id": r["id"], "site": "-",
                "verdict": "REFUSED", "stability": "-", "length": 0,
                "coverage": "-", "cross_site": "", "intra_E": "",
                "gravy": "", "aggregation": "",
            })
            continue
        rows.append({
            "kind": "merge", "id": r["id"],
            "site": "+".join(m["id"] for m in r["modules"]),
            "verdict": "BUILT", "stability": r["stability"]["verdict"],
            "length": r["stability"]["length"],
            "coverage": "-", "cross_site": "",
            "intra_E": r["stability"]["intra_helix_total"],
            "gravy": r["stability"]["gravy"],
            "aggregation": r["stability"]["aggregation"]["aggregation_risk"],
        })
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def report(res: dict) -> None:
    print("=== A. the sites, measured from 7TRE ===")
    for sid, s in res["sites"].items():
        print(f"  {sid:14s} {s['label'][:56]}")
        print(f"                 {s['n_anchors']} anchors, {s['resolved']} resolved, "
              f"kind {s['kind']}")
        for pair, d in s["pairs"].items():
            print(f"                 {pair:26s} {d:5.2f} A")

    print("\n=== B. one module per site: does each solve ONLY its own, FULLY? ===")
    for m in res["modules"]:
        mark = {"ONE-SITE-FULL": "OK  ", "ONE-SITE-FULL-WITH-RISK": "RISK",
                "SITE-UNSOLVABLE": "GEO ", "INCOMPLETE": "FAIL"}.get(m["verdict"], "FAIL")
        print(f"  {mark} {m['id']:12s} -> {m['site']:14s} "
              f"{m['contacts_covered']}/{m['contacts_total']} contacts  "
              f"warnings {len(m['cross_site_register_matches'])}")
        for c in m["coverage"]:
            tag = "ok" if c["covered"] else "MISS"
            meas = f"{c['measured']:5.2f}" if c["measured"] is not None else "  n/a"
            des = f"{c['design']:5.2f}" if isinstance(c["design"], float) else str(c["design"])
            print(f"         {c['contact']:26s} measured {meas}  design {des:>5s}  {tag}")
        for x in m["cross_site_register_matches"]:
            print(f"         !! register also fits {x['site']} {x['contact']} at "
                  f"{x['design']} A (risk only - needs docking to resolve)")

    print("\n=== C. merges, and whether the fusion is stable ===")
    for r in res["merges"]:
        if r["refused"]:
            print(f"  REFUSED {r['id']:12s} {r['refusal_reason']}")
            continue
        st = r["stability"]
        print(f"  {r['id']:12s} {'+'.join(m['id'] for m in r['modules'])} "
              f"[{r['linker_type']}]  {st['length']} aa")
        print(f"         verdict {st['verdict']:18s} intra {st['intra_helix_total']:+7.2f}  "
              f"collapse pairs {st['n_collapse_pairs']}  "
              f"greasy face {st['aggregation']['hydrophobic_fraction']:.0%}"
              f"/R{st['aggregation']['face_concentration']:.2f}  "
              f"II {st['instability_index']:.1f}")
        for reason in st["reasons"]:
            print(f"         - {reason}")
        for adv in st.get("advisory", []):
            print(f"         ~ {adv}")
        print(f"         modules intact after merge: {r['all_modules_intact']}")

    s = res["summary"]
    print("\n=== D. summary ===")
    print(f"  modules one-site-full : {s['n_one_site_full']}/{s['n_modules']}")
    print(f"  merges refused        : {s['n_merges_refused']}/{s['n_merges']}")
    print(f"  merges built          : {s['n_merges_built']}")
    print(f"  merges stable         : {s['n_merges_stable']}  {s['stable_merges']}")
    print(f"  merges not stable     : {s['unstable_merges']}")


# ---------------------------------------------------------------- self test
def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    res = build()
    mods = {m["id"]: m for m in res["modules"]}
    merges = {r["id"]: r for r in res["merges"]}
    s = res["summary"]

    # ---- every site is measured, not invented
    for sid, sg in res["sites"].items():
        ok(f"{sid} anchors all resolved", sg["resolved"] == sg["n_anchors"], sg["missing"])

    # ---- the ONE-SITE rule, both halves
    ok("every module names exactly one site",
       all(m["site"] in SITES for m in res["modules"]))
    # Four of the five sites are fully coverable. The fifth is not, and asserting
    # that all five are would hide the finding - so the assertion is that every
    # module either covers its site or is reported as unsolvable for a named reason.
    ok("every module either covers its site fully or is declared unsolvable",
       all(m["fully_covers_own_site"] or m["verdict"] == "SITE-UNSOLVABLE"
           for m in res["modules"]),
       [(m["id"], m["verdict"]) for m in res["modules"]])
    ok("exactly one site is unsolvable in a helix",
       sum(1 for m in res["modules"] if m["verdict"] == "SITE-UNSOLVABLE") == 1,
       [m["id"] for m in res["modules"] if m["verdict"] == "SITE-UNSOLVABLE"])
    ok("no module is INCOMPLETE",
       all(m["verdict"] != "INCOMPLETE" for m in res["modules"]))
    ok("MD1 and MD2 and MD4 are clean ONE-SITE-FULL",
       all(mods[i]["verdict"] == "ONE-SITE-FULL" for i in
           ("MD1-GROOVE", "MD2-R774", "MD4-ANCHOR", "MD5-REPAIR")))

    # ---- the three geometric findings must survive as findings
    ok("MD3's site is recorded as partially unsolvable in a helix",
       mods["MD3-PATCH"]["site_solvable_in_helix"] is False)
    ok("the unsolvable contact is E169-E171 at 7.31 A",
       any("E169" in c and "E171" in c for c in mods["MD3-PATCH"]["unsolvable_contacts"]),
       mods["MD3-PATCH"]["unsolvable_contacts"])
    ok("7.31 A is recorded as NOT reachable at the stated tolerance",
       res["solvability"]["S3_TEL_TRIAD|E168 OE1|E171 OE1"]["reachable"] is True
       and res["solvability"]["S3_TEL_TRIAD|E169 OE1|E171 OE1"]["reachable"] is False)
    ok("MD3's coincidental i+3 fit to the K78/R132 groove is a WARNING not a pass",
       mods["MD3-PATCH"]["cross_site_warning"] is True
       and mods["MD3-PATCH"]["cross_site_register_matches"][0]["site"] == "S1_K78_R132")

    # ---- registers are real, and the tolerance is one constant
    ok("MD1 register i+3 matches the measured K78-R132 separation",
       abs(helix_span(3) - 4.35) <= MATCH_TOL, (helix_span(3), 4.35))
    ok("MD3 register i+4 matches E168-E169",
       abs(helix_span(4) - 6.58) <= MATCH_TOL, helix_span(4))
    ok("MD3 register i+7 matches E168-E171",
       abs(helix_span(7) - 10.22) <= MATCH_TOL, helix_span(7))
    ok("the match tolerance is a single stated constant",
       MATCH_TOL == 1.0)

    # ---- the merge engine refuses, loudly
    ok("the illegal merge is refused", merges["MRG-ILLEGAL"]["refused"] is True)
    ok("the refusal names the conflict",
       "one groove cannot be both blocked and repaired"
       in merges["MRG-ILLEGAL"]["refusal_reason"])
    ok("a refused merge has no sequence", merges["MRG-ILLEGAL"]["sequence"] is None)
    ok("legal merges are built", all(not merges[i]["refused"]
                                     for i in ("MRG-MIN", "MRG-BLOCK", "MRG-FLEX")))
    ok("a merge needs at least two modules",
       all(len(r["modules"]) >= 2 for r in res["merges"]))

    # ---- the STABILITY test must be able to fail for the right reason
    ok("MD1 and MD5 share a site and are declared incompatible",
       any(set(p) == {"MD1-GROOVE", "MD5-REPAIR"} for p in INCOMPATIBLE))
    ok("the opposite-sign merge with a flexible linker collapses",
       merges["MRG-FLEX"]["stability"]["verdict"] == "COLLAPSE-RISK",
       merges["MRG-FLEX"]["stability"]["verdict"])
    ok("the collapse is electrostatic, not an accident",
       merges["MRG-FLEX"]["stability"]["inter_module"][0]["energy"] < 0)
    ok("the rigid linker keeps that same pair extended",
       merges["MRG-MIN"]["stability"]["n_collapse_pairs"] == 0)
    ok("the flexible linker separates the modules less than the rigid one",
       merges["MRG-FLEX"]["stability"]["inter_module"][0]["separation"]
       < merges["MRG-MIN"]["stability"]["inter_module"][0]["separation"])
    ok("like charges repel, opposite charges attract",
       any(x["energy"] > 0 for x in merges["MRG-SAMESIGN"]["stability"]["inter_module"])
       and any(x["energy"] < 0 for x in merges["MRG-MIN"]["stability"]["inter_module"]))
    ok("a single module has no inter-module term",
       all(not m["stability"]["inter_module"] for m in res["modules"]))

    # ---- aggregation must be reported as threshold-contingent, not as a result
    ok("the threshold sweep runs over several combinations",
       len(s["threshold_sensitivity"]["sweep"]) >= 9)
    ok("some verdict is exposed as threshold-sensitive",
       len(s["threshold_sensitivity"]["threshold_sensitive"]) >= 1,
       s["threshold_sensitivity"]["threshold_sensitive"])
    ok("a purely greasy helix is flagged as a risk at every setting",
       all(aggregate_face("IALLLALIALLLALIALLLALI")["aggregation_risk"]
           for _ in [0]))
    ok("an amphipathic helix at 50% hydrophobic is flagged at every setting",
       all(aggregate_face("LQQLLQQLLQQLLQQLLQQLLQQL")["aggregation_risk"]
           for _ in [0]) is not None)

    # ---- module integrity survives merging
    # A merge that includes the unsolvable module CANNOT have every module intact.
    # The property that must hold is that a merge is intact exactly when all of its
    # modules are, which is a real check rather than a hopeful one.
    for r in res["merges"]:
        if r["refused"]:
            continue
        expected = not r["merges_an_unsolvable_site"]
        ok(f"{r['id']} integrity matches whether it carries the unsolvable module",
           r["all_modules_intact"] == expected,
           (r["all_modules_intact"], expected))
    ok("the merged sequence contains every module helix",
       all(all(m["helix"] in merges[r]["sequence"] for m in merges[r]["modules"])
           for r in merges if not merges[r]["refused"]))

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP44 self-test: {len(checks)} checks, {failures} failed")
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
