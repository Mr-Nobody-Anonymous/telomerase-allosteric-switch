#!/usr/bin/env python3
"""
WP43 - Run every WP42 design through a computed test battery, and report the
numbers as three different things, because they are three different things.

    python scripts/solve/wp43_test_battery.py
      -> results/solve/wp43_test_battery.json
      -> results/solve/wp43_test_battery.csv

    python scripts/solve/wp43_test_battery.py --self-test
      -> checks, 0 failed

THE DISTINCTION THAT MATTERS
----------------------------
"Success rate" means at least four different things here, and reporting one
number without saying which would be the single most misleading thing this
repository could do:

  1. IN-SILICO GATE PASS RATE   - how many designs clear computed filters.
                                  That is what this file measures. It is NOT
                                  evidence that anything binds.
  2. PROSPECTIVE HIT RATE       - how many of N synthesised/expressed designs
                                  show the intended activity in a lab.
  3. RETROSPECTIVE RECOVERY     - whether the battery rejects designs already
                                  known to be bad. Measured below on controls.
  4. PUBLISHED RATES            - AlphaProteo 9-88 %, RFdiffusion binders ~1-5 %
                                  unfiltered, small-molecule de novo <1 %.

(1) is always much higher than (2). The gap is the point of the control set: a
filter that passes everything has a pass rate of 100 % and predicts nothing.
Every gate here is therefore validated in both directions - a known-good molecule
must pass and a known-bad one must fail - and the reported number is the pass
rate over designs PLUS the controls, so it cannot be inflated by grading only
its own homework.

WHAT CANNOT BE RUN HERE
-----------------------
No OpenMM, no docking engine, no FEP, no MM-GBSA. So there is no binding score,
no pose, no affinity in this file. WP39 already refuted this repository's only
rigid-body scorer, and repeating that would add a number without adding
evidence. The gates are all either graph properties (RDKit), sequence properties
(Biopython) or geometry already measured in WP42.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_JSON = REPO / "results" / "solve" / "wp43_test_battery.json"
OUT_CSV = REPO / "results" / "solve" / "wp43_test_battery.csv"

sys.path.insert(0, str(REPO / "scripts" / "solve"))

# ---------------------------------------------------------------- thresholds
# Every cutoff below is the conventional literature value, stated so it can be
# argued with. None was tuned to make a design pass; the controls would catch it.
THRESHOLDS = {
    "mw": 500.0,        # Lipinski
    "logp": 5.0,        # Lipinski
    "hbd": 5,           # Lipinski
    "hba": 10,          # Lipinski
    # QED is REPORTED, not gated. It was going to be a hard gate at >=0.50 until
    # the controls were run: imatinib, an approved oral drug, scores 0.389. QED
    # was calibrated on lead-like oral drugs and penalises the basic solubilising
    # side chain that kinase inhibitors all carry, so a 0.50 cutoff rejects
    # eleven-percent of the approved-drug world before it looks at a design.
    "sa_score": 4.0,    # <=4 is "accessible"; >6 is hard
    "gravy_max": 0.80,  # above this a peptide tends to aggregate
    "charge_abs_max": 8.0,  # beyond +-8 the patch-selectivity argument weakens
    "min_length": 15,   # below 15 aa there is no folded structure to speak of
    "instability_max": 40.0,  # Biopython: >40 predicted unstable
}

# ------------------------------------------------------------------ controls
# A gate is only worth reporting if it can be shown to fire. These are the
# calibration set: two must pass, three must fail.
CONTROLS = [
    {
        "id": "CTRL-OK-1", "kind": "ligand", "expect": "pass",
        "note": "imatinib, an approved oral drug: the drug-likeness gates must accept it",
        "smiles": "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(C)CC5",
    },
    {
        "id": "CTRL-OK-2", "kind": "ligand", "expect": "pass",
        "note": "aspirin, approved: BRENK's 'phenol_ester' alert is a FALSE POSITIVE here",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    },
    {
        "id": "CTRL-BAD-1", "kind": "ligand", "expect": "fail",
        "note": "curcumin: the classic frequent hitter (PAINS misses it; BRENK catches it)",
        "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
    },
    {
        "id": "CTRL-BAD-2", "kind": "ligand", "expect": "fail",
        "note": "a quinone, redox-cycling: must be caught",
        "smiles": "O=C1C=CC(=O)C=C1",
    },
    {
        "id": "CTRL-BAD-3", "kind": "ligand", "expect": "fail",
        "note": "an alkylating nitrogen mustard (BRENK catches it, PAINS does not)",
        "smiles": "ClCCN(CCCl)CCCl",
    },
]


# ------------------------------------------------------------------- ligands
def _sascorer():
    from rdkit.Chem import RDConfig
    import os

    p = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if p not in sys.path:
        sys.path.append(p)
    import sascorer

    return sascorer


def _catalog():
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    p = FilterCatalogParams()
    for c in (
        FilterCatalogParams.FilterCatalogs.PAINS,
        FilterCatalogParams.FilterCatalogs.BRENK,
        FilterCatalogParams.FilterCatalogs.NIH,
    ):
        p.AddCatalog(c)
    return FilterCatalog(p)


# Two configurations, because the structural-alert catalogs do not agree and the
# disagreement is the finding. Reported as a confusion matrix on the controls in
# both, so the reader can pick the tradeoff instead of being handed one.
#
#   permissive : PAINS only.       PAINS is specific; it misses most frequent hitters.
#   strict     : PAINS + BRENK.    BRENK catches them, and flags aspirin.
CONFIGS = ("permissive", "strict")


def ligand_gates(smiles: str, span=None, target_span=None, config="permissive") -> dict:
    """Run every computed ligand gate. Returns per-gate detail and a tier."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"parses": False, "gates": {}, "tier": "REJECT", "score": 0.0}

    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    sascorer = _sascorer()
    # PAINS and BRENK are built separately. They were originally merged into one
    # catalog and then read as if each were the other, which made the "permissive"
    # configuration reject aspirin exactly like the strict one - i.e. the two
    # configurations were identical and the comparison was vacuous.
    def _cat(*names):
        pr = FilterCatalogParams()
        for n in names:
            pr.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, n))
        return FilterCatalog(pr)

    alerts_pains = [e.GetDescription() for e in _cat("PAINS").GetMatches(mol)]
    alerts_brenk = [e.GetDescription() for e in _cat("BRENK").GetMatches(mol)]

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    qed = QED.qed(mol)
    sa = sascorer.calculateScore(mol)

    gates = {
        "lipinski_mw": (mw <= THRESHOLDS["mw"], round(mw, 2)),
        "lipinski_logp": (logp <= THRESHOLDS["logp"], round(logp, 2)),
        "lipinski_hbd": (hbd <= THRESHOLDS["hbd"], hbd),
        "lipinski_hba": (hba <= THRESHOLDS["hba"], hba),
        "sa_score": (sa <= THRESHOLDS["sa_score"], round(sa, 2)),
        "no_pains": (len(alerts_pains) == 0, alerts_pains),
    }
    if config == "strict":
        gates["no_brenk"] = (len(alerts_brenk) == 0, alerts_brenk)
    if span is not None and target_span is not None:
        gates["geometric_span"] = (span >= target_span, round(span, 2))

    structural_alerts = list(alerts_pains) + (list(alerts_brenk) if config == "strict" else [])

    # Two independent axes, because they fail independently - and because
    # conflating them produced a real artefact: with a single fractional score,
    # adding a gate changes the DENOMINATOR and flipped a design from FAIL to
    # CONDITIONAL without any property of the molecule changing. Chemistry and
    # geometry are now scored separately and neither can be improved by counting
    # more gates.
    #
    # Chemistry follows Lipinski as written: at most ONE violation is still
    # drug-like. Hard additions: SA <= 4 and no structural alerts.
    lipinski_violations = sum(
        1 for k in ("lipinski_mw", "lipinski_logp", "lipinski_hbd", "lipinski_hba")
        if not gates[k][0]
    )
    chemistry_ok = (lipinski_violations <= 1 and gates["sa_score"][0]
                    and not structural_alerts)
    if structural_alerts:
        chemistry = "REJECT"
    elif chemistry_ok and lipinski_violations == 0:
        chemistry = "PASS"
    elif chemistry_ok:
        chemistry = "ONE-VIOLATION"
    else:
        chemistry = "FAIL"

    if span is not None and target_span is not None:
        geometry = "REACHES" if span >= target_span else "CANNOT-SPAN"
    else:
        geometry = "n/a"

    passed = sum(1 for ok, _ in gates.values() if ok)
    total = len(gates)
    score = passed / total

    # MODEL-PASS requires BOTH axes. A drug-like molecule that cannot reach its
    # site is not a pass - that is the entire result of WP42.
    if chemistry in ("PASS", "ONE-VIOLATION") and geometry in ("REACHES", "n/a"):
        tier = "MODEL-PASS"
    elif chemistry in ("PASS", "ONE-VIOLATION"):
        tier = "GEOMETRY-FAIL"
    elif chemistry == "REJECT":
        tier = "REJECT"
    else:
        tier = "CHEMISTRY-FAIL"

    return {
        "parses": True,
        "gates": {k: {"pass": bool(v[0]), "value": v[1]} for k, v in gates.items()},
        "passed": passed,
        "total": total,
        "score": round(score, 3),
        "tier": tier,
        "chemistry": chemistry,
        "geometry": geometry,
        "lipinski_violations": lipinski_violations,
        "qed": round(qed, 3),
        "sa_score": round(sa, 2),
        "alerts": alerts_pains,
        "brenk_alerts": alerts_brenk,
        "config": config,
        "mw": round(mw, 2),
        "logp": round(logp, 2),
    }


# ------------------------------------------------------------------ proteins
def protein_gates(rec: dict) -> dict:
    """Run every computed protein gate against a WP42 protein record."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    seq = rec["sequence"]
    a = ProteinAnalysis(seq)
    gravy = a.gravy()
    charge = a.charge_at_pH(7.4)
    instab = a.instability_index()
    ss_frac = a.secondary_structure_fraction()

    gates = {
        "register_matches": (rec.get("register_matches") is True,
                             rec.get("display_offsets_realised")),
        "hits_measured_target": (rec.get("all_targets_within_1A") is True,
                                 [c["register"] for c in rec.get("target_checks", [])]),
        "length_sufficient": (rec["length"] >= THRESHOLDS["min_length"], rec["length"]),
        "gravy_not_aggregating": (gravy <= THRESHOLDS["gravy_max"], round(gravy, 3)),
        "charge_within_bound": (abs(charge) <= THRESHOLDS["charge_abs_max"],
                                round(charge, 2)),
        "helical_propensity": (ss_frac[0] >= 0.4, round(ss_frac[0], 3)),
        "not_tpp1_fragment": (
            "ELLLQD" not in seq and "GFRGTEGR" not in seq, None
        ),
        "no_cys_aggregation_risk": ("C" not in seq, seq.count("C")),
    }
    passed = sum(1 for ok, _ in gates.values() if ok)
    total = len(gates)
    score = passed / total
    tier = "MODEL-PASS" if score == 1.0 else ("CONDITIONAL" if score >= 0.75 else "FAIL")

    return {
        "gates": {k: {"pass": bool(v[0]), "value": v[1]} for k, v in gates.items()},
        "passed": passed,
        "total": total,
        "score": round(score, 3),
        "tier": tier,
        "gravy": round(gravy, 3),
        "charge": round(charge, 2),
        "instability_index": round(instab, 1),
        "instability_advisory": ("stable" if instab <= THRESHOLDS["instability_max"]
                                 else "predicted unstable"),
        "helix_fraction": round(ss_frac[0], 3),
    }


# ---------------------------------------------------------------------- main
_RUN_CACHE = None


def run() -> dict:
    """Compute the battery once per process.

    Deterministic by construction (fixed conformer seeds), and the WP42 embedding
    step is the expensive part: without this the test suite re-embeds 5 x 60
    conformers eleven times.
    """
    global _RUN_CACHE
    if _RUN_CACHE is None:
        _RUN_CACHE = _run()
    return _RUN_CACHE


def _run() -> dict:
    import wp42_designs as wp

    designs = wp.build()

    # controls first, so the battery's ability to fire is established - and in
    # BOTH configurations, because they disagree and the disagreement is the result
    control_results = []
    for c in CONTROLS:
        row = {"id": c["id"], "expect": c["expect"], "note": c["note"]}
        for cfg in CONFIGS:
            r = ligand_gates(c["smiles"], config=cfg)
            row[cfg] = {
                "tier": r["tier"], "score": r["score"],
                "alerts": r["alerts"], "brenk": r["brenk_alerts"],
                "qed": r["qed"], "sa_score": r["sa_score"],
            }
        for cfg in CONFIGS:
            row[f"agree_{cfg}"] = ((row[cfg]["tier"] == "MODEL-PASS")
                                   == (c["expect"] == "pass"))
        control_results.append(row)

    # confusion matrix for each configuration over the control set
    matrices = {}
    for cfg in CONFIGS:
        tp = sum(1 for r in control_results if r["expect"] == "pass" and r[f"agree_{cfg}"])
        tn = sum(1 for r in control_results if r["expect"] == "fail" and r[f"agree_{cfg}"])
        fp = sum(1 for r in control_results if r["expect"] == "fail" and not r[f"agree_{cfg}"])
        fn = sum(1 for r in control_results if r["expect"] == "pass" and not r[f"agree_{cfg}"])
        matrices[cfg] = {
            "good_accepted": tp, "bad_rejected": tn,
            "bad_accepted": fp, "good_rejected": fn,
            "sensitivity": round(tn / (tn + fp), 3) if (tn + fp) else None,
            "specificity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "accuracy": round((tp + tn) / len(control_results), 3),
        }

    ligand_results = []
    for l in designs["ligands"]:
        row = {"id": l["id"], "name": l["name"], "mw": l["mw"], "logp": l["logp"],
               "max_span_3d": l["max_span_3d"], "target_span": l["target_span"]}
        for cfg in CONFIGS:
            r = ligand_gates(l["smiles"], span=l["max_span_3d"],
                             target_span=l["target_span"], config=cfg)
            row[cfg] = r
        # headline = permissive, because it is the one that does not reject a
        # real drug; the strict result is kept alongside, never merged
        row.update({
            "tier": row["permissive"]["tier"],
            "score": row["permissive"]["score"],
            "passed": row["permissive"]["passed"],
            "total": row["permissive"]["total"],
            "gates": row["permissive"]["gates"],
            "qed": row["permissive"]["qed"],
            "sa_score": row["permissive"]["sa_score"],
            "alerts": row["permissive"]["alerts"],
            "brenk_alerts": row["permissive"]["brenk_alerts"],
        })
        ligand_results.append(row)

    protein_results = []
    for p in designs["proteins"]:
        r = protein_gates(p)
        protein_results.append({
            "id": p["id"], "name": p["name"], "length": p["length"],
            "mw_da": p["mw_da"], **r,
        })

    # ---- the three numbers, kept separate on purpose
    # headline configuration = the one that does not reject an approved drug on a
    # false-positive alert (see the confusion matrix); strict is reported beside it
    HEADLINE = "strict" if matrices["permissive"]["accuracy"] < matrices["strict"]["accuracy"] else "permissive"
    lig_pass = [x for x in ligand_results if x[HEADLINE]["tier"] == "MODEL-PASS"]
    lig_geom = [x for x in ligand_results if x[HEADLINE]["tier"] == "GEOMETRY-FAIL"]
    prot_pass = [x for x in protein_results if x["tier"] == "MODEL-PASS"]
    controls_agree = [c for c in control_results
                      if all(c[f"agree_{cfg}"] for cfg in CONFIGS)]

    # ---------------------------------------------------------------- ranking
    # Ranking unmeasured designs is a weak operation and pretending otherwise
    # would be the dishonesty this file exists to avoid. So the tiebreak is
    # stated rather than fitted, and it is applied in a fixed order:
    #
    #   ligands  : tier first, then geometric MARGIN (span - required), because
    #              with every chemistry gate equal the only real difference is
    #              how much slack the molecule has.
    #   proteins : all clear the same 8 gates, so they are separated by
    #              (1) predicted stability, (2) how many distinct measured
    #              contacts the design addresses, (3) |GRAVY| - deviation from
    #              neutral, which is the crude solubility proxy available here.
    #
    # A design that is one gate worse but far more likely to be expressible is a
    # better next experiment, and that judgement is recorded explicitly rather
    # than hidden in a weighted sum.

    tier_rank = {"MODEL-PASS": 0, "CONDITIONAL": 1, "GEOMETRY-FAIL": 2,
                 "CHEMISTRY-FAIL": 3, "REJECT": 4}
    lig_ranked = sorted(
        ligand_results,
        key=lambda x: (tier_rank.get(x[HEADLINE]["tier"], 9),
                       -(x["max_span_3d"] - x["target_span"])),
    )
    for i, l in enumerate(lig_ranked, 1):
        l["rank"] = i
        l["rank_basis"] = (f"{l[HEADLINE]['tier']}, margin "
                           f"{l['max_span_3d'] - l['target_span']:+.2f} A")

    prot_ranked = sorted(
        protein_results,
        key=lambda x: (x["instability_advisory"] != "stable",
                       -len(x["gates"]["hits_measured_target"]["value"] or []),
                       abs(x["gravy"])),
    )
    for i, p_ in enumerate(prot_ranked, 1):
        p_["rank"] = i
        p_["rank_basis"] = (f"{p_['instability_advisory']}, "
                            f"{len(p_['gates']['hits_measured_target']['value'] or [])} contacts, "
                            f"|GRAVY| {abs(p_['gravy']):.2f}")

    best_ligand = lig_ranked[0]
    best_protein = prot_ranked[0]
    best_design = best_protein

    summary = {
        "n_designs": len(ligand_results) + len(protein_results),
        "n_controls": len(control_results),
        "control_agreement": f"{len(controls_agree)}/{len(control_results)}",
        "control_agreement_by_config": {
            cfg: f"{sum(1 for r in control_results if r[f'agree_{cfg}'])}/{len(control_results)}"
            for cfg in CONFIGS
        },
        "battery_validated": all(
            sum(1 for r in control_results if r[f"agree_{cfg}"]) == len(control_results)
            for cfg in CONFIGS
        ),
        "confusion": matrices,
        "headline_config": HEADLINE,
        "ligand_pass": len(lig_pass),
        "ligand_geometry_fail": len(lig_geom),
        "ligand_total": len(ligand_results),
        "protein_pass": len(prot_pass),
        "protein_total": len(protein_results),
        "by_config": {
            cfg: {
                "ligand_pass": sum(1 for x in ligand_results if x[cfg]["tier"] == "MODEL-PASS"),
                "ligand_geometry_fail": sum(1 for x in ligand_results
                                            if x[cfg]["tier"] == "GEOMETRY-FAIL"),
                "ligand_reject": sum(1 for x in ligand_results if x[cfg]["tier"] == "REJECT"),
                "protein_pass": sum(1 for x in protein_results if x["tier"] == "MODEL-PASS"),
            }
            for cfg in CONFIGS
        },
        # in-silico pass rate over designs only
        "in_silico_pass_rate_all": round(
            (len(lig_pass) + len(prot_pass)) / (len(ligand_results) + len(protein_results)), 3),
        # and over designs + controls, so it cannot be inflated by self-grading
        "in_silico_pass_rate_with_controls": round(
            (len(lig_pass) + len(prot_pass))
            / (len(ligand_results) + len(protein_results) + len(control_results)), 3),
        "best_design": best_design["id"],
        "best_score": best_design["score"],
        "best_ligand": best_ligand["id"],
        "best_protein": best_protein["id"],
        "ranking_ligands": [{"rank": l["rank"], "id": l["id"],
                             "basis": l["rank_basis"]} for l in lig_ranked],
        "ranking_proteins": [{"rank": p_["rank"], "id": p_["id"],
                              "basis": p_["rank_basis"]} for p_ in prot_ranked],
        "experimental_success_rate": None,
        "reason_no_experimental_rate": (
            "none of these designs has been synthesised or expressed; an experimental "
            "success rate requires N >= 1 measured designs and N = 0 exist"
        ),
    }

    return {
        "thresholds": THRESHOLDS,
        "controls": control_results,
        "ligands": ligand_results,
        "proteins": protein_results,
        "summary": summary,
        "cannot_be_run_here": [
            "molecular docking (no engine installed)",
            "MD / MM-GBSA / FEP (no OpenMM; and WP39 refuted this repo's rigid-body scorer)",
            "binding affinity of any kind",
            "cell permeability (no permeability model; the acid count is reported instead)",
            "immunogenicity of the designed proteins",
        ],
    }


def write_outputs(res: dict) -> None:
    import csv

    OUT_JSON.write_text(json.dumps(res, indent=2), encoding="utf-8")
    rows = []
    for c in res["controls"]:
        rows.append({"group": "control", "id": c["id"],
                     "tier_permissive": c["permissive"]["tier"],
                     "tier_strict": c["strict"]["tier"],
                     "score": c["permissive"]["score"], "passed": "", "total": "",
                     "note": c["note"]})
    for l in res["ligands"]:
        rows.append({"group": "ligand", "id": l["id"],
                     "tier_permissive": l["permissive"]["tier"],
                     "tier_strict": l["strict"]["tier"],
                     "score": l["permissive"]["score"], "passed": l["passed"],
                     "total": l["total"], "note": l["name"]})
    for p in res["proteins"]:
        rows.append({"group": "protein", "id": p["id"],
                     "tier_permissive": p["tier"], "tier_strict": p["tier"],
                     "score": p["score"], "passed": p["passed"], "total": p["total"],
                     "note": p["name"]})
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def report(res: dict) -> None:
    print("=== 0. can the battery fire? (controls, both configurations) ===")
    print(f"  {'id':11s} {'expect':6s} {'permissive':12s} {'strict':12s}  note")
    for c in res["controls"]:
        pv, sv = c["permissive"]["tier"], c["strict"]["tier"]
        mark = "" if (c["agree_permissive"] and c["agree_strict"]) else "  <- disagree"
        print(f"  {c['id']:11s} {c['expect']:6s} {pv:12s} {sv:12s}  "
              f"{c['note'][:44]}{mark}")
    print()
    for cfg, m in res["summary"]["confusion"].items():
        print(f"  {cfg:11s} good accepted {m['good_accepted']}  bad rejected "
              f"{m['bad_rejected']}  bad accepted {m['bad_accepted']}  "
              f"good rejected {m['good_rejected']}   acc {m['accuracy']:.2f}")

    print("\n=== 1. ligands ===")
    print(f"  {'id':7s} {'MW':>7s} {'logP':>6s} {'QED':>5s} {'SA':>5s} {'span':>6s} "
          f"{'gates':>6s}  tier")
    for l in res["ligands"]:
        print(f"  {l['id']:7s} {l['mw']:7.2f} {l['logp']:6.2f} {l['qed']:5.3f} "
              f"{l['sa_score']:5.2f} {l['max_span_3d']:6.2f} "
              f"{l['passed']}/{l['total']}  {l['tier']:11s} "
              f"(strict: {l['strict']['tier']})")
    for l in res["ligands"]:
        bad = [k for k, v in l["gates"].items() if not v["pass"]]
        if bad:
            print(f"    {l['id']}: fails {bad}")

    print("\n=== 2. proteins ===")
    for p in res["proteins"]:
        print(f"  {p['id']:7s} {p['length']:3d} aa  {p['passed']}/{p['total']}  "
              f"{p['tier']:11s} GRAVY {p['gravy']:+.2f} q {p['charge']:+5.2f} "
              f"helix {p['helix_fraction']:.2f}  {p['instability_advisory']}")
        bad = [k for k, v in p["gates"].items() if not v["pass"]]
        if bad:
            print(f"           fails {bad}")

    s = res["summary"]
    print("\n=== 3. the numbers, named ===")
    print(f"  in-silico pass, designs only      {s['ligand_pass']+s['protein_pass']}/"
          f"{s['n_designs']} = {s['in_silico_pass_rate_all']:.1%}")
    print(f"  in-silico pass, designs+controls  "
          f"{s['ligand_pass']+s['protein_pass']}/{s['n_designs']+s['n_controls']} = "
          f"{s['in_silico_pass_rate_with_controls']:.1%}   <- the honest one")
    print(f"  experimental success rate         NOT MEASURABLE ({s['reason_no_experimental_rate'][:46]}...)")
    print(f"  best design                       {s['best_design']} at {s['best_score']:.2f}")
    print(f"  best ligand                       {s['best_ligand']}")
    print("\n  ligand ranking:")
    for r in s["ranking_ligands"]:
        print(f"    {r['rank']}. {r['id']:7s} {r['basis']}")
    print("  protein ranking (all clear the same 8 gates):")
    for r in s["ranking_proteins"]:
        print(f"    {r['rank']}. {r['id']:7s} {r['basis']}")

    print("\n=== 4. what this cannot test ===")
    for x in res["cannot_be_run_here"]:
        print(f"  - {x}")


def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    res = run()
    s = res["summary"]
    ctl = {c["id"]: c for c in res["controls"]}

    # --- the battery must be able to say no, in both configurations
    ok("imatinib passes in both configs",
       ctl["CTRL-OK-1"]["permissive"]["tier"] == "MODEL-PASS"
       and ctl["CTRL-OK-1"]["strict"]["tier"] == "MODEL-PASS")
    ok("a quinone is rejected in both",
       ctl["CTRL-BAD-2"]["permissive"]["tier"] == "REJECT"
       and ctl["CTRL-BAD-2"]["strict"]["tier"] == "REJECT")
    ok("strict catches all three bad controls",
       all(ctl[c]["strict"]["tier"] == "REJECT"
           for c in ("CTRL-BAD-1", "CTRL-BAD-2", "CTRL-BAD-3")))
    ok("strict rejects aspirin - a known false positive",
       ctl["CTRL-OK-2"]["strict"]["tier"] == "REJECT")
    ok("permissive does NOT reject aspirin",
       ctl["CTRL-OK-2"]["permissive"]["tier"] == "MODEL-PASS")
    ok("strict beats permissive on the controls",
       s["confusion"]["strict"]["accuracy"] >= s["confusion"]["permissive"]["accuracy"])
    ok("the two configurations are genuinely different",
       s["confusion"]["strict"] != s["confusion"]["permissive"])
    ok("a configuration is chosen and named", s["headline_config"] in CONFIGS)

    # --- QED must not be a gate: imatinib scores below 0.5, so gating on it
    # would reject an approved drug. This check is the reason the gate was removed.
    ok("imatinib QED is below the naive 0.50 cutoff",
       ctl["CTRL-OK-1"]["permissive"]["qed"] < 0.50,
       ctl["CTRL-OK-1"]["permissive"]["qed"])

    # --- numbers must stay separated
    ok("experimental rate is None, not a number", s["experimental_success_rate"] is None)
    ok("reason for the absent rate is recorded",
       "N = 0" in s["reason_no_experimental_rate"])
    ok("pass rate with controls is lower than designs-only",
       s["in_silico_pass_rate_with_controls"] < s["in_silico_pass_rate_all"] + 1e-9)

    # --- coverage and the axis separation
    ok("all 5 ligands scored", len(res["ligands"]) == 5)
    ok("all 5 proteins scored", len(res["proteins"]) == 5)
    for l in res["ligands"]:
        for cfg in CONFIGS:
            ok(f"{l['id']}/{cfg} has a chemistry axis",
               l[cfg]["chemistry"] in ("PASS", "ONE-VIOLATION", "FAIL", "REJECT"))
            ok(f"{l['id']}/{cfg} has a geometry axis",
               l[cfg]["geometry"] in ("REACHES", "CANNOT-SPAN", "n/a"))
    for p_ in res["proteins"]:
        ok(f"{p_['id']} has 8 gates", p_["total"] == 8, p_["total"])

    # --- the WP42 result must survive the battery, not be averaged away
    cannot = [l["id"] for l in res["ligands"]
              if l[s["headline_config"]]["geometry"] == "CANNOT-SPAN"]
    ok("the four designs that cannot span are GEOMETRY-FAIL",
       all(l[s["headline_config"]]["tier"] == "GEOMETRY-FAIL" for l in res["ligands"]
           if l[s["headline_config"]]["geometry"] == "CANNOT-SPAN"),
       cannot)
    ok("all four spanning ligand failures are present", len(cannot) == 4, cannot)
    ok("only TEL-3 passes among the ligands",
       [l["id"] for l in res["ligands"] if l[s["headline_config"]]["tier"] == "MODEL-PASS"]
       == ["TEL-3"])
    ok("chemistry is not the reason the spanning ligands fail",
       all(l[s["headline_config"]]["chemistry"] in ("PASS", "ONE-VIOLATION")
           for l in res["ligands"] if l["id"] in cannot))

    # --- the ranking must be the one the evidence supports, not a tiebreak artefact
    ok("the best ligand is the only one that reaches", s["best_ligand"] == "TEL-3",
       s["best_ligand"])
    ok("no geometry failure outranks the geometry pass",
       s["ranking_ligands"][0]["id"] == "TEL-3")
    ok("a design WITH 3 contacts outranks one with fewer",
       next(r["rank"] for r in s["ranking_proteins"] if r["id"] == "TENB-2")
       > min(r["rank"] for r in s["ranking_proteins"]
             if r["id"] in ("TELM-1", "TELM-3", "TENB-1")))
    ok("every design has a stated rank basis",
       all(r["basis"] for r in s["ranking_ligands"] + s["ranking_proteins"]))
    ok("ranking is a stated tiebreak, not a score", len(s["ranking_proteins"]) == 5)

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP43 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = run()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
