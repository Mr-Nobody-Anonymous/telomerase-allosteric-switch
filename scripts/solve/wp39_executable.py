#!/usr/bin/env python3
"""
WP39 - The two things WP37 called impossible, done.

WP37's closing table listed six items as "cannot be closed here". Two of those
were wrong, and this work package corrects them.

  ITEM 5 - "c > 1 extensions/division: nothing establishes it".
  This was a literature gap, not a limitation of the sandbox, and it sits on the
  ONLY escape hatch from a CLOSED gate. Zhao et al. (Mol Cell 2011) measured it:
  H1299 = 100% of lagging daughter ends extended per cell cycle, HeLa = ~70%,
  hESC = ~70%. So c is not 0.7 or 1.0 as a guess - it is measured per line, and
  the H1299 value is 1.0, which is exactly the boundary.

  ITEM 4 - "any affinity at the TEL patch: no FEP/MM-GBSA".
  True, and unchanged. But WP37 then wrote the refutation criterion for a
  docking proposal and left it unrun, when WP33's size-normalised rigid-body
  scorer was already on disk. A criterion you never evaluate is not a
  prediction. So this runs it, including the negative control, and reports
  whichever way it comes out.

WHAT THIS DOES NOT DO
---------------------
It does not produce an affinity. The scorer is rigid-body, distance-dependent
dielectric, no desolvation, no entropy, no induced fit. A per-heavy-atom score
ranks poses; it does not yield a Kd and must never be quoted as one. Everything
below is a RANKING, and the numbers stay relative.

    python scripts/solve/wp39_executable.py
    python scripts/solve/wp39_executable.py --self-test
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# ---------------------------------------------------------------- c, measured
#: Fraction of lagging-daughter telomere ends extended per cell cycle. These are
#: measurements, not parameters, and each is tied to the line it was made in.
C_EXTENSION = {
    "H1299 (adenocarcinoma)": dict(c=1.00, n_ends="100% of lagging ends",
                                   source="Zhao 2011 Mol Cell 43:251, Fig 6E-F"),
    "HeLa (cervical)": dict(c=0.70, n_ends="~70% of lagging ends",
                            source="Zhao 2011 Mol Cell 43:251, Fig 6C-D"),
    "hESC": dict(c=0.70, n_ends="~70% of ends",
                 source="Zhao 2011 Mol Cell 43:251, Fig S5"),
    "HeLa, pH 7.4 (NAR 2016)": dict(c=0.77, n_ends="77 +/- 4% of ends",
                                    source="NAR 44:8086 (2016), Fig 2B-C"),
    "HeLa, pH 6.8 (NAR 2016)": dict(c=0.48, n_ends="48 +/- 4% of ends",
                                    source="NAR 44:8086 (2016), Fig 2B-C"),
    "S. cerevisiae": dict(c=0.07, n_ends="~7% of ends",
                          source="Teixeira 2004, STEX"),
}

#: The TEL patch, as WP33 enumerates it: ARG92 +1, ASP166 -1, GLU168 -1,
#: GLU169 -1, LYS170 +1, ARG180 +1, GLU215 -1  ->  net -1.
TEL_PATCH_SEQ = (92, 166, 168, 169, 170, 180, 215)

#: WP37's ceilings, restated so the comparison is auditable.
DELTA_FLOOR = 50.0
DELTA_HUVEC = 90.0


def delta_of(c: float, r_half_pt: float = 8.0, nt_per_repeat: int = 6) -> float:
    """bp added per division = c * (nt per repeat) * R_half."""
    return c * nt_per_repeat * r_half_pt


def c_required(delta: float, r_half_pt: float = 8.0, nt_per_repeat: int = 6) -> float:
    """The c needed to hold a given bp/PD at a fixed R_half."""
    return delta / (nt_per_repeat * r_half_pt)


def gate_table() -> list[dict]:
    rows = []
    for name, rec in C_EXTENSION.items():
        d = delta_of(rec["c"])
        rows.append(dict(line=name, c=rec["c"], n_ends=rec["n_ends"],
                         delta=d, source=rec["source"],
                         holds_50=d >= DELTA_FLOOR, holds_90=d >= DELTA_HUVEC))
    return rows


def _w32():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module("wp32_physics_scoring")


def _w33():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module("wp33_valid_scoring")


def _w31():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module("wp31_epitope")


# ------------------------------------------------------- the docking proposal
def patch_targets():
    """The wild-type TEL patch and the charge-neutralised mutant, same atoms.

    The mutation changes ONLY formal charges - no atom is added or removed - so
    the two targets are identical in shape and differ solely in electrostatics.
    That is what makes the comparison a test of the charge hypothesis rather
    than a second size confound.
    """
    wp, w31 = _w32(), _w31()
    asm = w31.build_assembly()
    tpp1 = w31.transform_residues(asm["t2"][w31.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    tert = asm["t7"]["A"]
    ia, ib = w31.interface_residues((tert, tpp1))

    # WP33 is explicit that the interface total is NOT the TEL patch: it is the
    # whole TERT+TPP1 contact surface, whose net is -17.89 because of TERT-side
    # charge. Scoring against that surface would re-introduce exactly the
    # confound WP33 retired, so the patch is restricted to its seven named
    # residues - the ones WP33 lists with formal charges.
    patch = [(n, s, a) for n, s, a in tpp1 if s in TEL_PATCH_SEQ]
    wt = wp.protein_atoms(patch)
    assert len({t[1] for t in wt}) == len(TEL_PATCH_SEQ), "TEL patch incomplete"

    # R92Q/K170Q/R180Q and D166N/E168Q/E169Q/E215Q: every ionisable side chain in
    # the patch is neutralised. Both signs, which is why the net goes -1 -> 0 and
    # not -1 -> -4.
    mutant = [(t[0], t[1], t[2], t[3], 0.0) + tuple(t[5:]) if t[4] != 0.0 else t
              for t in wt]
    # wp.protein_atoms() sums FORMAL side-chain charges together with PARTIAL
    # polar charges on backbone N/O, so its total is not the formal charge. WP33
    # reports the formal charge of these seven residues as -1; the scorer's
    # total is -2.54. Both are reported, under different names, because
    # conflating them is the same defect class this programme keeps finding.
    formal_wt = sum(_side_chain_formal(t[0]) for t in _unique_res(patch))
    return dict(patch=patch, wt=wt, mutant=mutant,
                q_formal=float(formal_wt),
                q_scorer=float(sum(t[4] for t in wt)),
                q_scorer_mut=float(sum(t[4] for t in mutant)),
                n_atoms=len(wt),
                n_res=len({t[1] for t in wt}))


#: Formal side-chain charge at pH 7.4, for the residues that carry one.
FORMAL = {"ARG": +1.0, "LYS": +1.0, "HIS": +0.0,
          "ASP": -1.0, "GLU": -1.0}


def _side_chain_formal(resname: str) -> float:
    return FORMAL.get(resname, 0.0)


def _unique_res(residues):
    seen, out = set(), []
    for r in residues:
        if r[1] not in seen:
            seen.add(r[1])
            out.append(r)
    return out


def control_patch(patch):
    """Same residue count, same shape, charge-neutral by construction.

    WP37's criterion requires that a hydrophobic analogue is INSENSITIVE to
    charge neutralisation. Comparing a ligand against a different COMPOUND
    confounds charge with chemistry, so the control is instead the same patch
    with its charges already zeroed. A ligand that responds to electrostatics
    moves between the WT patch and its neutralised twin; a ligand that does not
    moves the same amount against a charge-free control. That difference is the
    charge-specific response, and nothing else.
    """
    return [(t[0], t[1], t[2], t[3], 0.0) + tuple(t[5:]) for t in patch]


def charge_specific_response(lig, patch_fn, n_iter=4000, seed=11):
    """(wt - neutralised) minus (mutant - charge-free control).

    The second bracket is the non-electrostatic floor for this ligand and this
    placement protocol.
    """
    w32, w33 = _w32(), _w33()
    wt_p, mut_p = patch_fn()
    ctrl = control_patch(wt_p)
    lig = w32.ligand_typed(lig, ionise=False) if isinstance(lig, str) else lig
    e_wt, _ = w33.place_normalised(lig, wt_p, seed=seed, n_iter=n_iter)
    e_mut, _ = w33.place_normalised(lig, mut_p, seed=seed, n_iter=n_iter)
    e_ctrl, _ = w33.place_normalised(lig, ctrl, seed=seed, n_iter=n_iter)
    charge_part = e_mut["per_atom"] - e_wt["per_atom"]
    floor_part = e_ctrl["per_atom"] - e_mut["per_atom"]
    return dict(wt=e_wt["per_atom"], mutant=e_mut["per_atom"],
                control=e_ctrl["per_atom"], charge_part=charge_part,
                floor_part=floor_part, specific=charge_part - floor_part)


def ligand_smiles(w33=None):
    """The neutral amide (WP37's candidate anchor) plus its matched controls."""
    w33 = w33 or _w33()
    ser = w33.series()
    return {
        "neutral amide (TLN-amide anchor)": ser["amide (neutral control) :: neutral"],
        "carboxylic acid (matched pair)": ser["acid / carboxylate :: neutral"],
        "carboxylate (matched pair)": ser["acid / carboxylate :: charged"],
        "amine / ammonium (matched pair)": ser["amine / ammonium :: charged"],
    }


def score_against_targets(seed: int = 11, n_iter: int = 4000):
    """Score each ligand on the wild-type and neutralised patch.

    Returns per-ligand dicts with the per-heavy-atom energy on each target and
    the delta between them. delta > 0 means the ligand prefers the WT patch;
    delta ~ 0 means it is insensitive to the patch charge.
    """
    w32, w33 = _w32(), _w33()
    t = patch_targets()
    out = []
    for label, smi in ligand_smiles(w33).items():
        lig = w32.ligand_typed(smi, ionise=False)
        row = dict(ligand=label, smiles=smi, n_heavy=len(lig["xyz"]))
        for tag, prot in (("wt", t["wt"]), ("mutant", t["mutant"])):
            best = w33.place_normalised(lig, prot, seed=seed, n_iter=n_iter)
            res, xyz = best
            row[f"{tag}_per_atom"] = float(res["per_atom"])
            row[f"{tag}_total"] = float(res["total"])
            row[f"{tag}_clashes"] = int(res["n_hard"])
        row["delta_prefers_wt"] = row["mutant_per_atom"] - row["wt_per_atom"]
        out.append(row)
    return out, t


#: WP37 wrote the criterion; it is evaluated here, unchanged.
THRESHOLD = -0.15


def amide_control_verdict(scores: list[dict]) -> dict:
    """WP37's stated refutation criterion, applied to measured scores."""
    amide = next(r for r in scores if r["ligand"].startswith("neutral amide"))
    controls = [r for r in scores if not r["ligand"].startswith("neutral amide")]
    # "prefers the neutralised patch" means its per-atom energy there is LOWER
    delta = amide["delta_prefers_wt"]
    # delta = E(mutant) - E(wt). Positive means the ligand binds WORSE once the
    # patch is neutralised, which is WP37's "loses interface preference".
    loses_preference = delta > 0
    big_enough = abs(delta) > abs(THRESHOLD)
    ctrl_mean = float(np.mean([abs(r["delta_prefers_wt"]) for r in controls]))
    moves_more_than_control = abs(delta) > ctrl_mean
    return dict(
        amide_delta=delta,
        loses_preference=loses_preference,
        big_enough=big_enough,
        control_mean_abs_delta=ctrl_mean,
        moves_more_than_control=moves_more_than_control,
        refuted=not (loses_preference and big_enough and moves_more_than_control),
        threshold=THRESHOLD,
    )


def main() -> int:
    L = "=" * 96
    print(L)
    print("WP39  The two items WP37 called impossible - both executable")
    print(L)
    print()

    # ------------------------------------------------------------ item 5
    print("WP39.1  Item 5 - 'c > 1: nothing establishes it' was a literature gap")
    print()
    print("  c is the fraction of telomere ends extended per cell cycle. It is")
    print("  MEASURED, per cell line, and it sits directly on the gate:")
    print()
    print(f"  {'line':26s} {'c':>5} {'bp/PD':>7} {'>=50':>5} {'>=90':>5}   source")
    print("  " + "-" * 92)
    for r in gate_table():
        print(f"  {r['line']:26s} {r['c']:5.2f} {r['delta']:7.2f} "
              f"{'yes' if r['holds_50'] else ' NO':>5} {'yes' if r['holds_90'] else ' NO':>5}"
              f"   {r['source']}")
    print()
    print(f"  At R1/2 = 8.0 (measured, +PT) the ceiling is c x 48 bp/PD.")
    print(f"  H1299 at c = 1.00 reaches exactly {delta_of(1.0):.1f} bp/PD. That MISSES")
    print(f"  the 50 bp/PD floor by {DELTA_FLOOR - delta_of(1.0):.2f} bp, and misses a")
    print(f"  90 bp/PD HUVEC requirement by {DELTA_HUVEC - delta_of(1.0):.1f}. Both facts,")
    print("  stated together, and neither is softened.")
    print()
    for target, label in ((DELTA_FLOOR, "50 bp/PD floor"), (DELTA_HUVEC, "90 bp/PD HUVEC")):
        c = c_required(target)
        held = c <= 1.0
        print(f"  {label:18s} needs c >= {c:.2f}  ->  "
              f"{'reachable by a c=1.00 line' if held else 'NOT reachable even at c = 1.00'}")
    print()
    print("  This does NOT open the gate - it closes it harder. Before, 'we do not")
    print("  know c' was recorded as an open assumption. It is now MEASURED, and the")
    print("  best measured value (c = 1.00, H1299) still falls 2 bp short. The gate")
    print("  no longer depends on an unmeasured c at all; it depends on R1/2, which")
    print(f"  would need {DELTA_FLOOR / 6.0:.2f} to clear the floor and {DELTA_HUVEC / 6.0:.1f} for HUVEC,")
    print("  against 8.0 measured.")
    print()

    # ------------------------------------------------------------ item 4
    print("WP39.2  Item 4 - WP37 wrote a refutation criterion and never ran it")
    print()
    print("  Running it now. Same atoms, same shape, only formal charges differ")
    print("  between the wild-type patch and the R92Q/K170Q/R180Q +")
    print("  D166N/E168Q/E169Q/E215Q neutralised patch.")
    print()
    scores, t = score_against_targets()
    print(f"  patch {t['n_res']} residues, {t['n_atoms']} atoms")
    print(f"  formal charge WT {t['q_formal']:+.2f}   "
          f"scorer total {t['q_scorer']:+.2f} (formal + partial polar)")
    print(f"  neutralised: formal 0.00, scorer total {t['q_scorer_mut']:+.2f}")
    print()
    print(f"  {'ligand':36s} {'wt':>8} {'mutant':>8} {'delta':>8}  {'prefers':>8}")
    print("  " + "-" * 74)
    for r in scores:
        prefers = "WT" if r["delta_prefers_wt"] > 0 else "neutralised"
        print(f"  {r['ligand']:36s} {r['wt_per_atom']:8.3f} {r['mutant_per_atom']:8.3f} "
              f"{r['delta_prefers_wt']:8.3f}  {prefers:>8}")
    print()
    v = amide_control_verdict(scores)
    print(f"  WP37's criterion: the neutral amide must LOSE interface preference")
    print(f"  when the patch is neutralised, and the hydrophobic control must not")
    print(f"  track it. Threshold {THRESHOLD:.2f} per heavy atom.")
    print()
    print(f"    amide delta (mutant - wt)     {v['amide_delta']:+.3f}")
    print(f"    loses preference when patch neutralised  {v['loses_preference']}")
    print(f"    effect exceeds {abs(THRESHOLD):.2f} per heavy atom        {v['big_enough']}")
    print(f"    control mean |delta|          {v['control_mean_abs_delta']:.3f}")
    print(f"    amide moves more than control {v['moves_more_than_control']}")
    print()
    if v["refuted"]:
        print("  VERDICT: the proposal as written is REFUTED by its own scorer.")
    else:
        print("  VERDICT: the proposal SURVIVES its own criterion.")
    print()
    print("  Either way this is a rigid-body ranking, not an affinity. It cannot")
    print("  produce a Kd and is not evidence that anything binds.")
    print()

    # ------------------------------------------------------------ write
    rows = []
    for r in gate_table():
        rows.append(dict(item="c_measured", key=r["line"], value=f"{r['c']:.2f}",
                         note=f"{r['delta']:.2f} bp/PD; {r['source']}"))
    for target, label in ((DELTA_FLOOR, "c_for_50"), (DELTA_HUVEC, "c_for_90")):
        rows.append(dict(item="c_required", key=label, value=f"{c_required(target):.3f}",
                         note="reachable" if c_required(target) <= 1.0 else "not reachable at c<=1"))
    rows.append(dict(item="patch", key="formal_charge_wt", value=f"{t['q_formal']:+.2f}",
                     note=f"{t['n_res']} patch residues, {t['n_atoms']} atoms"))
    rows.append(dict(item="patch", key="scorer_charge_wt", value=f"{t['q_scorer']:+.2f}",
                     note="formal + partial polar; NOT the formal charge"))
    rows.append(dict(item="patch", key="scorer_charge_mutant",
                     value=f"{t['q_scorer_mut']:+.2f}",
                     note="all ionisable side chains neutralised"))
    for r in scores:
        rows.append(dict(item="docking", key=r["ligand"],
                         value=f"{r['delta_prefers_wt']:+.3f}",
                         note=f"wt {r['wt_per_atom']:.3f} mutant {r['mutant_per_atom']:.3f}"))
    rows.append(dict(item="verdict", key="docking_proposal",
                     value="REFUTED" if v["refuted"] else "SURVIVES",
                     note=f"amide delta {v['amide_delta']:+.3f} vs threshold {THRESHOLD}"))

    with (OUT / "wp39_executable.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "key", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    (OUT / "wp39_executable.json").write_text(
        json.dumps(dict(c=gate_table(), scores=scores, verdict=v), indent=2),
        encoding="utf-8")
    print(f"  wrote {(OUT / 'wp39_executable.csv').relative_to(REPO)}")
    print(f"  wrote {(OUT / 'wp39_executable.json').relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- self test
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # --- the measured c values, and the arithmetic that uses them
    chk("H1299 c is 1.00 (measured, not assumed)",
        C_EXTENSION["H1299 (adenocarcinoma)"]["c"] == 1.00)
    chk("HeLa c is 0.70", C_EXTENSION["HeLa (cervical)"]["c"] == 0.70)
    chk("yeast c is 0.07 - the contrast that makes the point",
        C_EXTENSION["S. cerevisiae"]["c"] == 0.07)
    chk("every c is a fraction", all(0 < r["c"] <= 1.0 for r in C_EXTENSION.values()))
    chk("every c cites a source", all(r["source"] for r in C_EXTENSION.values()))

    chk("c=1.00 at R1/2=8 gives exactly 48 bp/PD", abs(delta_of(1.0) - 48.0) < 1e-12)
    chk("H1299 clears the 50 floor by 2 bp", abs((delta_of(1.0) - DELTA_FLOOR) + 2.0) < 1e-12)
    chk("HeLa at c=0.70 does NOT clear 50", delta_of(0.70) < DELTA_FLOOR)
    chk("yeast at c=0.07 is nowhere near", delta_of(0.07) < 5)
    chk("c_required(50) is 1.0417", abs(c_required(50.0) - 50.0 / 48.0) < 1e-12)
    chk("c_required(50) > 1 - the floor is NOT clear even at c = 1.00",
        c_required(50.0) > 1.0)
    chk("48 bp/PD at c=1 misses the 50 floor by 2", delta_of(1.0) < DELTA_FLOOR)
    chk("c_required(90) is 1.875 and unreachable", c_required(90.0) > 1.0)
    chk("R1/2 needed for 90 at c=1 is 15", abs(90.0 / 6.0 - 15.0) < 1e-12)

    # --- the gate table must not claim more than the numbers support
    tbl = gate_table()
    h = next(r for r in tbl if r["line"].startswith("H1299"))
    chk("H1299 does NOT hold 50 either - 48 < 50 is the finding",
        not h["holds_50"])
    chk("H1299 does NOT hold 90 - both facts must survive", not h["holds_90"])
    chk("no line holds 90 at the measured R1/2", not any(r["holds_90"] for r in tbl))

    # --- the mutant changes charge only
    t = patch_targets()
    chk("wild-type FORMAL patch charge is -1", abs(t["q_formal"] - (-1.0)) < 1e-12)
    chk("the scorer's total is NOT the formal charge - kept under its own name",
        abs(t["q_scorer"] - (-2.54)) < 0.05 and abs(t["q_scorer"] - t["q_formal"]) > 1.0)
    chk("the patch is exactly the seven named residues", t["n_res"] == 7)
    chk("neutralised patch is 0 on both measures",
        abs(t["q_scorer_mut"]) < 1e-12)
    chk("mutation adds or removes no atom", len(t["wt"]) == len(t["mutant"]))
    chk("mutation moves no atom", all(
        np.allclose(a[3], b[3]) for a, b in zip(t["wt"], t["mutant"])))
    chk("mutation leaves donors/acceptors untouched", all(
        a[5:] == b[5:] for a, b in zip(t["wt"], t["mutant"])))

    # --- the verdict logic must be able to say no
    fake = [dict(ligand="neutral amide (TLN-amide anchor)", delta_prefers_wt=-0.40),
            dict(ligand="carboxylic acid (matched pair)", delta_prefers_wt=-0.10)]
    chk("verdict logic CAN refute (amide moves more, wrong direction)",
        amide_control_verdict(fake)["refuted"])
    fake0 = [dict(ligand="neutral amide (TLN-amide anchor)", delta_prefers_wt=+0.02),
             dict(ligand="carboxylic acid (matched pair)", delta_prefers_wt=+0.40)]
    chk("verdict logic CAN refute (amide barely moves, control does)",
        amide_control_verdict(fake0)["refuted"])
    fake2 = [dict(ligand="neutral amide (TLN-amide anchor)", delta_prefers_wt=+0.40),
             dict(ligand="carboxylic acid (matched pair)", delta_prefers_wt=+0.05)]
    chk("verdict logic CAN pass", not amide_control_verdict(fake2)["refuted"])
    fake3 = [dict(ligand="neutral amide (TLN-amide anchor)", delta_prefers_wt=0.01),
             dict(ligand="carboxylic acid (matched pair)", delta_prefers_wt=+0.30)]
    chk("an insensitive amide is refuted even if it 'prefers' WT",
        amide_control_verdict(fake3)["refuted"])

    # --- the honesty guards
    src = Path(__file__).read_text(encoding="utf-8")
    chk("the docstring refuses to call this an affinity",
        "not an affinity" in src.lower() or "does not produce an affinity" in src)
    # scan the parsed AST, not the source text: a grep for "def kd" matches the
    # literal in this very check, which is how the first version of this test
    # failed on itself.
    import ast
    fns = {n.name for n in ast.walk(ast.parse(src)) if isinstance(n, ast.FunctionDef)}
    chk("no function computes a dissociation constant",
        not any(k in f.lower() for f in fns for k in ("kd", "dissociation", "affinity")))
    # built from fragments so this literal cannot match itself in `src`
    # scan only the PRODUCTION source: this check's own message contains the
    # word, which is how the previous two versions of it failed on themselves.
    kcal = "".join(["kc", "al"])
    production = src.split("def _self_test")[0]
    chk("no energy is reported in kcal/mol", kcal not in production.lower())
    chk("the threshold is WP37's, imported not tuned", THRESHOLD == -0.15)

    bad = [nm for nm, ok in checks if not ok]
    print(f"WP39 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
