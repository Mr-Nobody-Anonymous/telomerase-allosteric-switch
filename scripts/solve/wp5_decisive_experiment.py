#!/usr/bin/env python3
"""
WP5 - The decisive experiment: design and power.

WP2-WP3 leave exactly one open distinction.  Put every published measurement in one table and
ask which model can produce it:

    Retention-only (R)   TPP1 changes the DNA/output affinity of the enzyme, with an effect that
                         is expressed only while the enzyme is cycling.
    Allostery (A)        TPP1 additionally changes a rate at a site it does not contact -
                         the active site is the cleanest one.

Every published quantity so far is reproduced by R with no free parameters (WP3).  A is not
excluded by any of them.  The experiment that separates them must therefore perturb TPP1 while
holding the DNA *output* fixed, and measure a rate at the active site.

This script:
  5.1  enumerates candidate experiments and scores them against the two models;
  5.2  computes the power of the surviving design;
  5.3  states the pre-declared decision rule.

Run:  python scripts/solve/wp5_decisive_experiment.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"


# --------------------------------------------------------------------------------------
# 5.1  candidate experiments scored against the two surviving models
# --------------------------------------------------------------------------------------
CANDIDATES = [
    dict(
        experiment="Primer titration +/- TPP1 (the repository's proposed decisive test)",
        retention="saturating [P] plateau unchanged; apparent K_app decreases",
        allostery="plateau increases",
        separates="C1 vs B, within the retention family",
        note="Cannot separate R from A: both are retention-family predictions. WP2 shows the "
             "assay runs at [P]/K_d = 100, where the affinity term is already silent.",
        verdict="not decisive for the open question",
    ),
    dict(
        experiment="Resting-state primer dissociation +/- TPP1 (already done, L&C 2010)",
        retention="unchanged",
        allostery="unchanged",
        separates="nothing",
        note="Both models predict no change, and none is observed. This is why the assay cannot "
             "be used as evidence about mechanism in either direction.",
        verdict="uninformative for the open question",
    ),
    dict(
        experiment="Steady-state processivity +/- TPP1 (already done, many papers)",
        retention="2.5-6x increase, reproduced to 0.4 SD in WP3",
        allostery="same prediction",
        separates="nothing",
        note="A scalar observable; both models can hit any value by parameter choice.",
        verdict="uninformative",
    ),
    dict(
        experiment="Single-round translocation efficiency +/- TPP1 (already done, L&C 2010)",
        retention="57% -> 99% (a retention prediction: fewer failures to commit)",
        allostery="same prediction",
        separates="nothing on its own",
        note="Strong constraint on the retention model's parameters, but A predicts it too.",
        verdict="constrains R, does not separate",
    ),
    dict(
        experiment="Single-turnover nucleotide incorporation +/- TPP1, no primer",
        retention="k_chem unchanged (TPP1 acts on output affinity, which has no substrate)",
        allostery="k_chem changes",
        separates="R vs A directly",
        note="THE decisive measurement. TPP1 contacts TERT at TEN/IFD-TRAP, 38.9 A from the "
             "catalytic site (WP4). With no primer there is no retention route to the active "
             "site, so any change in k_chem must be propagated through the protein.",
        verdict="DECISIVE",
    ),
    dict(
        experiment="Single-turnover incorporation +/- TPP1 at saturating primer",
        retention="k_chem unchanged; only the probability of retaining the primer changes",
        allostery="k_chem changes",
        separates="R vs A",
        note="A useful positive control for the no-primer condition: if TPP1 changes k_chem here "
             "but not without primer, the effect requires a bound substrate and is not pure "
             "allostery.",
        verdict="decisive (control)",
    ),
    dict(
        experiment="Hydrogen-deuterium exchange +/- TPP1, no DNA",
        retention="no change outside the TPP1 footprint",
        allostery="protection/deprotection at the active site or RT/CTE",
        separates="R vs A",
        note="Complementary; sees conformational ensembles rather than rates, so it can be "
             "positive while k_chem is unchanged (a real possibility, and the reason to pair "
             "them).",
        verdict="decisive (complementary)",
    ),
]


def power_for_kchem_effect(sigma_log: float, n: int, alpha: float = 0.05) -> float:
    """
    Two-sample t-test on log k_chem, two-sided alpha, power 0.8 target.
    Returns the minimum detectable fold-change.
    """
    from scipy import stats
    # solve for delta such that power = 0.8
    d = 0.1
    while d < 10:
        ncp = d / (sigma_log * math.sqrt(2.0 / n))
        crit = stats.t.ppf(1 - alpha / 2, df=2 * n - 2)
        power = 1 - stats.nct.cdf(crit, df=2 * n - 2, nc=ncp) + stats.nct.cdf(-crit, df=2 * n - 2, nc=ncp)
        if power >= 0.8:
            return math.exp(d)
        d *= 1.02
    return float("inf")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    print("!" * 92)
    print("!! SUPERSEDED: the no-primer single-turnover experiment this script ranks first is")
    print("!! NOT POSSIBLE.  Telomerase nucleotide incorporation requires a primer base-paired to")
    print("!! the RNA template, so with no primer there is no k_chem to measure.  That candidate")
    print("!! is withdrawn; see results/solve/CORRECTIONS.md (C1) for the corrected design")
    print("!! (translocation at matched DNA residence time).  The power table below remains valid")
    print("!! for any matched-occupancy k_chem readout.")
    print("!" * 92)
    print()
    print("=" * 92)
    print("WP5.1  Candidate experiments against the two surviving models")
    print("=" * 92)
    for c in CANDIDATES:
        print(f"  {c['verdict'].upper():32s} {c['experiment']}")
        print(f"       retention predicts : {c['retention']}")
        print(f"       allostery predicts : {c['allostery']}")
        print(f"       separates          : {c['separates']}")
        print(f"       note               : {c['note']}")
        print()

    print("=" * 92)
    print("WP5.2  Power of the decisive design")
    print("=" * 92)
    print("  Readout: single-turnover k_chem (min^-1), TPP1-saturated vs mock.")
    print("  A single-turnover incorporation measurement on reconstituted telomerase reaches a")
    print("  within-experiment CV of 10-20% when counts are high; 15% is assumed below.")
    print()
    print(f"  {'CV':>6} {'n/group':>8} {'min detectable fold-change':>28}")
    rows = []
    for cv in (0.10, 0.15, 0.20, 0.30):
        sigma_log = math.sqrt(math.log(1 + cv ** 2))
        for n in (3, 4, 5, 6, 8):
            fc = power_for_kchem_effect(sigma_log, n)
            print(f"  {cv:>6.2f} {n:>8d} {fc:>26.2f}x")
            rows.append(dict(cv=cv, n_per_group=n, min_detectable_fold=round(fc, 3)))
        print()

    print("  Design consequence: at CV = 15% and n = 5 per condition the design detects a")
    print("  1.35x change in k_chem at 80% power (1.6x at CV=20%).  A TPP1-dependent allosteric effect on the")
    print("  chemical step would have to be smaller than that to escape detection - in which")
    print("  case it is also too small to explain the 2.7x processivity effect, which WP3")
    print("  accounts for without it.")
    print()

    print("=" * 92)
    print("WP5.3  Pre-declared decision rule")
    print("=" * 92)
    rule = [
        ("k_chem(TPT) / k_chem(mock) >= 1.6x or <= 1/1.6x, 95% CI excluding 1, "
         "with no primer", "allostery DEMONSTRATED: TPP1 changes a rate at a site it does "
                           "not contact"),
        ("CI excludes 1 only when primer is present", "retention-dependent, not allosteric: "
                                                      "the effect requires substrate, so it "
                                                      "is an affinity/first-order term"),
        ("CI includes 1 in both conditions", "retention-only CONFIRMED to the design's "
                                             "resolution; no mechanical component at the "
                                             "catalytic step"),
        ("no measurable single-turnover rate under either condition", "assay failed - the "
                                                                      "result is uninformative, "
                                                                      "not negative"),
    ]
    for cond, out in rule:
        print(f"  if {cond}\n     -> {out}\n")

    with open(OUT / "wp5_candidates.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(CANDIDATES[0]))
        w.writeheader()
        w.writerows(CANDIDATES)
    with open(OUT / "wp5_power.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT/'wp5_candidates.csv'}")
    print(f"wrote {OUT/'wp5_power.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
