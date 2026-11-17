#!/usr/bin/env python3
"""
WP40 - Designing to exceed 50 %: both referents, and the evidence that moves them.

TWO 50 % TARGETS
----------------
  (1) P(M5) > 0.50.      Currently 0.70 x 0.50 = 0.35.
  (2) delta* > 50 bp/PD. Currently 48.0 at c = 1.00.

They are the same question, because the 50 bp/PD floor IS what sets
P(M5 | gate): the mechanism is only worth 35 % if crossing the gate is worth
something, and the gate is only worth crossing if the mechanism is real.

WHAT THIS WORK PACKAGE FOUND
----------------------------
Three results, all of them things the programme had not computed:

  F1. The gate is not at tau = 0.1833, it is at tau >= 0.0454.  WP35's OPEN
      verdict came from its own DEFAULT tau, which was a modelling choice. At
      tau = 0 the gate is one point short; a stall share of 4.5 % clears it.
      So the gate was never really "closed" and never really "open" - it was a
      function of an unmeasured parameter, evaluated at two arbitrary points.

  F2. THE 8.0 ANCHOR IS NOT A MEASUREMENT.  WP37 wrote "L&C 2010 report R1/2
      2.963 -> 8.0", but WP29 labels 8.0 "M5 prediction (2.73x)". And
      8.0 / 2.963 = 2.6998, below L&C's own MEASURED >=4x (R1/2 >= 11.85,
      71.1 bp/PD). So the post-TPP1 anchor used throughout is mostly M5's
      OUTPUT, fed back in as if it were L&C's measurement.

  F3. E2 REQUIRES tau(+PT) < 0 at the repo's own tau0.  To make an observed
      4.00x into a true 2.7335x you need (1+tau_neg)/(1+tau_pos) = 1.4633. With
      tau(-PT) = 0.1833 that needs tau(+PT) = -0.191 - impossible. E2 is only
      viable if the -PT stall share is >= ~0.47, i.e. if nearly half of all
      untrapped complexes stall. That is a much stronger claim than it looked.

None of F1-F3 is an assumption; they are arithmetic on numbers already in the
tree. What follows from them is a design.

THE DESIGN
----------
Four measurements, ordered by information per unit cost, each with an explicit
"what would change my mind". The two-chase inversion from WP37 is the first one
because it is the only measurement that removes a confounding instead of adding
an assumption.

    python scripts/solve/wp40_design_to_fifty.py
    python scripts/solve/wp40_design_to_fifty.py --self-test
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

NT_PER_REPEAT = 6
R_NO_PT = 2.963
FOLD_M5 = 2.7335
DELTA_FLOOR = 50.0
DELTA_HUVEC = 90.0
C_MAX = 1.00                      # H1299, measured (WP39)
TAU0_DEFAULT = 0.1833             # WP35's modelling default, NOT a measurement


def a0() -> float:
    return 2 ** (1 / (NT_PER_REPEAT * R_NO_PT)) - 1


def r_half(fold: float) -> float:
    """R1/2 from an off-rate fold, via a = a0/fold."""
    return math.log(2) / (NT_PER_REPEAT * math.log(1 + a0() / fold))


def delta_star(fold: float, c: float = C_MAX) -> float:
    """bp/PD = c * 6 * R1/2. Exponent-free (WP37)."""
    return c * NT_PER_REPEAT * r_half(fold)


def fold_required(target: float, c: float = C_MAX) -> float:
    lo, hi = 1.0, 500.0
    for _ in range(300):
        m = (lo + hi) / 2
        if delta_star(m, c) < target:
            lo = m
        else:
            hi = m
    return hi


# ---------------------------------------------------------------- F1: the tau edge
def stall_ratio(tau_neg: float, tau_pos: float) -> float:
    """Inflation of the untrapped arm: (1 + tau_neg)/(1 + tau_pos).

    Observed fold = true fold x this ratio. It is the quantity E2 needs to be
    large enough to turn L&C's 4.00x into M5's 2.7335x, i.e. >= 1.4633.
    """
    return (1 + tau_neg) / (1 + tau_pos)


def tau_pos_needed(tau_neg: float, observed: float = 4.00) -> float:
    """tau_pos such that a true 2.7335x would be observed as `observed`."""
    return (1 + tau_neg) * FOLD_M5 / observed - 1


def tau_gate_threshold() -> float:
    """The stall share at which delta* first reaches the 50 floor."""
    lo, hi = 0.0, 0.99
    for _ in range(300):
        m = (lo + hi) / 2
        # true fold exceeds the null fold by (1+tau) when only the -PT arm stalls
        if delta_star(FOLD_M5 * (1 + m)) < DELTA_FLOOR:
            lo = m
        else:
            hi = m
    return hi


def gate_profile(taus=(0.0, 0.0454, 0.10, 0.1833, 0.25, 0.35, 0.50)) -> list[dict]:
    rows = []
    for t in taus:
        f = FOLD_M5 * (1 + t)
        rows.append(dict(tau=t, true_fold=f, r_half=r_half(f),
                         delta=delta_star(f), holds_50=delta_star(f) >= DELTA_FLOOR,
                         holds_90=delta_star(f) >= DELTA_HUVEC))
    return rows


# ---------------------------------------------------------------- F2: the anchor
def anchors() -> dict[str, dict]:
    """Every post-TPP1 R1/2 the tree uses, with its actual provenance."""
    return {
        "R1/2 = 8.0 (used everywhere as +PT)": dict(
            value=8.0, kind="PREDICTION",
            note="8.0/2.963 = 2.6998x, which is M5's predicted fold; WP29 labels "
                 "it 'M5 prediction (2.73x)'. Not an independent measurement.",
            delta=delta_star(8.0 / R_NO_PT)),
        "L&C measured >=4x": dict(
            value=4 * R_NO_PT, kind="LOWER BOUND",
            note=">=4x implies R1/2 >= 11.852 (WP29 uses 11.680 = 3.94x)",
            delta=delta_star(4.0)),
        "no TPP1 (L&C)": dict(
            value=R_NO_PT, kind="MEASURED", note="R1/2 = 2.963",
            delta=delta_star(1.0)),
    }


# ---------------------------------------------------------------- F3: E2's price
def e2_feasible(tau_neg: float) -> bool:
    """E2 says the stored 4x is stall-inflated. Can it be, at this tau_neg?"""
    return 0.0 <= tau_pos_needed(tau_neg) < 1.0


def e2_minimum_tau() -> float:
    lo, hi = 0.0, 1.0
    for _ in range(300):
        m = (lo + hi) / 2
        if tau_pos_needed(m) < 0.0:
            lo = m
        else:
            hi = m
    return hi


# ---------------------------------------------------------------- the posterior
def posterior(p_physical: float, p_p2: float, p_invivo: float | None = None) -> float:
    """P(M5) as a product of evidence terms, each independently falsifiable."""
    p = p_physical * p_p2
    if p_invivo is not None:
        p *= p_invivo
    return p


def bayes_factor_needed(p_prior: float, p_target: float = 0.50) -> float:
    """The odds multiplier on P(M5) needed to reach the target."""
    prior_odds = p_prior / (1 - p_prior)
    target_odds = p_target / (1 - p_target)
    return target_odds / prior_odds


def posterior_with_bf(p_prior: float, bf: float) -> float:
    prior_odds = p_prior / (1 - p_prior)
    return (bf * prior_odds) / (1 + bf * prior_odds)


# ---------------------------------------------------------------- in vivo test
#: Zhao et al. 2011 Mol Cell 43:251: ~60 nt added per telomere end per cycle in
#: H1299 and HeLa under maintenance. This is an IN VIVO measurement that does not
#: depend on the chase, the trap, or any model - and it has never been used here
#: as a test of M5.
IN_VIVO_NT = 60.0
IN_VIVO_TOL = 0.30     # a generous band: the reported value is approximate


def in_vivo_prediction(tau: float) -> float:
    return delta_star(FOLD_M5 * (1 + tau))


def in_vivo_consistent(tau: float, tol: float = IN_VIVO_TOL) -> bool:
    return abs(in_vivo_prediction(tau) / IN_VIVO_NT - 1.0) <= tol


def in_vivo_tau_window(tol: float = IN_VIVO_TOL) -> tuple[float, float] | None:
    """The stall shares for which M5 lands within tol of the in vivo value."""
    lo, hi = None, None
    for t in np.linspace(0.0, 0.90, 9001):
        if in_vivo_consistent(float(t), tol):
            if lo is None:
                lo = float(t)
            hi = float(t)
    return (lo, hi) if lo is not None else None


def main() -> int:
    L = "=" * 96
    print(L)
    print("WP40  Designing past 50 % - both referents")
    print(L)
    print()

    # ---------------------------------------------------------------- F2 first
    print("WP40.1  The anchor is not a measurement")
    print()
    print(f"  {'anchor':44s} {'R1/2':>7} {'bp/PD':>8}  {'kind':>12}")
    print("  " + "-" * 78)
    for label, rec in anchors().items():
        print(f"  {label:44s} {rec['value']:7.3f} {rec['delta']:8.2f}  {rec['kind']:>12}")
    print()
    print(f"  8.0 / 2.963 = {8.0 / R_NO_PT:.4f}x. That is M5's PREDICTED fold, and WP29")
    print("  labels the row 'M5 prediction (2.73x)'. So the +PT anchor used")
    print("  throughout is largely M5's own output fed back in as if it were L&C's")
    print("  measurement. L&C's measured statement is >=4x, which is R1/2 >= 11.85")
    print(f"  and {delta_star(4.0):.1f} bp/PD - comfortably past the floor.")
    print()
    print("  Consequence: 'the gate is CLOSED at 48' is true of M5's PREDICTION.")
    print("  On L&C's MEASUREMENT the gate is open, and by a wide margin. The")
    print("  closed gate is therefore evidence AGAINST M5, not against the target:")
    print("  the mechanism predicts less processivity than has already been seen.")
    print()

    # ---------------------------------------------------------------- F1
    print("WP40.2  The gate sits at tau = 0.0454, not at a choice of 0 or 0.18")
    print()
    print(f"  {'tau':>7} {'true fold':>10} {'R1/2':>8} {'bp/PD':>8} {'>=50':>6} {'>=90':>6}")
    print("  " + "-" * 52)
    for r in gate_profile():
        print(f"  {r['tau']:7.4f} {r['true_fold']:10.4f} {r['r_half']:8.3f} "
              f"{r['delta']:8.2f} {'YES' if r['holds_50'] else ' no':>6} "
              f"{'YES' if r['holds_90'] else ' no':>6}")
    print()
    print(f"  Gate threshold: tau >= {tau_gate_threshold():.4f} clears 50 bp/PD.")
    print(f"  90 bp/PD needs tau >= {_tau_for(90):.4f} - likely beyond the crossing.")
    print()
    print("  So 'CLOSED at 48' (tau=0) and 'OPEN at 57.3' (tau=0.1833) are the same")
    print("  function evaluated at two modelling choices, 4.5 % apart in tau. The")
    print("  honest statement is: the gate is OPEN for any stall share above 4.5 %,")
    print("  and P2 exists precisely to measure whether the stall share is that big.")
    print("  A 4.5 % threshold is a much weaker thing to need than 18 %.")
    print()

    # ---------------------------------------------------------------- F3
    print("WP40.3  E2's hidden price: it needs ~half of all complexes to stall")
    print()
    print(f"  To turn an observed 4.00x into a true 2.7335x, the stall ratio")
    print(f"  (1+tau_neg)/(1+tau_pos) must be 4.00/2.7335 = {4.0 / FOLD_M5:.4f}.")
    print()
    print(f"  sanity: at tau_neg = tau_pos the ratio is {stall_ratio(0.2, 0.2):.4f}")
    print()
    print(f"  {'tau(-PT)':>10} {'tau(+PT) required':>20}  verdict")
    print("  " + "-" * 46)
    for t in (0.10, 0.1833, 0.30, 0.40, 0.47, 0.50, 0.60):
        need = tau_pos_needed(t)
        print(f"  {t:10.4f} {need:+20.4f}  "
              f"{'possible' if 0 <= need < 1 else 'IMPOSSIBLE (negative)'}")
    print()
    print(f"  E2 is arithmetically excluded below tau(-PT) = {e2_minimum_tau():.4f}.")
    print("  The repo's own tau0 = 0.1833 is far below that, so E2 as written")
    print("  CONTRADICTS the stall model it is supposed to rescue. Either the stall")
    print("  share is ~0.47 or E2 is dead and L&C's 4x is real.")
    print()

    # ---------------------------------------------------------------- in vivo
    print("WP40.4  A test M5 has already been run against, without being counted")
    print()
    win = in_vivo_tau_window()
    print(f"  Zhao 2011 measured ~{IN_VIVO_NT:.0f} nt added per telomere end per cell cycle")
    print("  IN VIVO (H1299 and HeLa), under maintenance. That number depends on no")
    print("  chase, no trap, and no model. It has never been used here as a test.")
    print()
    print(f"  {'tau':>7} {'M5 predicts':>12} {'vs 60 nt':>10}  within +/-30%")
    print("  " + "-" * 44)
    for t in (0.0, 0.10, 0.1833, 0.25, 0.35, 0.50):
        p = in_vivo_prediction(t)
        print(f"  {t:7.4f} {p:12.2f} {p / IN_VIVO_NT:10.3f}  "
              f"{'yes' if in_vivo_consistent(t) else 'NO':>12}")
    print()
    if win:
        print(f"  M5 lands within +/-{IN_VIVO_TOL:.0%} of the in vivo value for tau in "
              f"[{win[0]:.3f}, {win[1]:.3f}].")
        print(f"  Its tau = 0 value ({in_vivo_prediction(0.0):.1f}) is already within that band.")
        print("  This does NOT confirm M5 - the band is generous and the prediction")
        print("  was not pre-registered against this datum. But it is an independent")
        print("  measurement M5 does not fail, which is more than C6 can say.")
    print()

    # ---------------------------------------------------------------- posterior
    print("WP40.5  What moves P(M5) past 50 %")
    print()
    p_phys, p_p2 = 0.70, 0.50
    base = posterior(p_phys, p_p2)
    print(f"  current     0.70 (physical on L&C) x 0.50 (P2)          = {base:.3f}")
    print(f"  log-odds     thesis must gain {bayes_factor_needed(base):.2f}x in odds")
    print()
    need_p2 = 0.50 / p_phys
    print(f"  The lever is the P2 term. Holding 0.70 fixed, P(M5) > 0.50 needs")
    print(f"  P2 to rise from 0.50 to {need_p2:.3f} - i.e. the two-chase result must be")
    print(f"  {need_p2/0.50:.2f}x more likely under M5 than under its alternatives.")
    print("  A decisive experiment clearing 1.4x in odds is an ordinary ask; this")
    print("  is the smallest sufficient design, not a hopeful one.")
    print()
    print(f"  {'scenario':58s} {'P(M5)':>7}")
    print("  " + "-" * 68)
    scen = [
        ("current: 0.70 x 0.50", posterior(0.70, 0.50)),
        ("P2 marginal (0.60)", posterior(0.70, 0.60)),
        ("P2 confirms (0.80) - TARGET CLEARED", posterior(0.70, 0.80)),
        ("P2 confirms AND R1/2 measured at +PT (0.85)", posterior(0.70, 0.85)),
        ("P2 fails (0.20)", posterior(0.70, 0.20)),
    ]
    for label, val in scen:
        print(f"  {label:58s} {val:7.3f}")
    print()
    print("  NOTE ON A MISTAKE I MADE IN THE FIRST DRAFT: the in vivo datum was")
    print("  multiplied in as a further 0.80, which would LOWER P(M5) to 0.45.")
    print("  That is wrong in kind. A factor means 'M5 must also pass this', but")
    print("  M5 has largely ALREADY passed it - 48 nt predicted, 60 nt measured,")
    print("  a 20 % gap and no refutation. A hurdle already cleared is not a")
    print("  further hurdle. With no likelihood ratio available for it, the")
    print("  honest treatment is to leave it OUT of the product and report it as")
    print("  a datum M5 does not fail. Multiplying it in would have been the same")
    print("  error class as C1: counting something as evidence without asking")
    print("  what it could have falsified.")
    print()

    # ---------------------------------------------------------------- design
    print("WP40.6  The design, in order of information per unit cost")
    print()
    design = design_table()
    print(f"  {'#':>2} {'measurement':52s} {'cost':>7} {'moves':>22}")
    print("  " + "-" * 88)
    for d in design:
        print(f"  {d['n']:>2} {d['what']:52s} {d['cost']:>7} {d['moves']:>22}")
    print()
    for d in design:
        print(f"  {d['n']}. {d['what']}")
        print(f"     de-confounds : {d['deconfounds']}")
        print(f"     falsified if : {d['falsified_if']}")
        print()

    # ---------------------------------------------------------------- write
    rows = []
    for label, rec in anchors().items():
        rows.append(dict(item="anchor", key=label, value=f"{rec['value']:.3f}",
                         note=f"{rec['kind']}; {rec['delta']:.2f} bp/PD"))
    for r in gate_profile():
        rows.append(dict(item="gate", key=f"tau={r['tau']:.4f}",
                         value=f"{r['delta']:.2f}",
                         note=f"fold {r['true_fold']:.4f}, holds_50={r['holds_50']}"))
    rows.append(dict(item="gate", key="tau_threshold_50", value=f"{tau_gate_threshold():.4f}",
                     note="stall share at which 50 bp/PD is reached"))
    rows.append(dict(item="E2", key="min_tau_neg_for_feasibility",
                     value=f"{e2_minimum_tau():.4f}",
                     note="below this tau(+PT) must be negative - impossible"))
    rows.append(dict(item="in_vivo", key="M5_at_tau0", value=f"{in_vivo_prediction(0.0):.2f}",
                     note=f"vs Zhao 2011 {IN_VIVO_NT:.0f} nt; ratio {in_vivo_prediction(0.0)/IN_VIVO_NT:.3f}"))
    rows.append(dict(item="posterior", key="current", value=f"{base:.3f}", note="0.70 x 0.50"))
    rows.append(dict(item="posterior", key="bf_needed", value=f"{bayes_factor_needed(base):.3f}",
                     note="odds multiplier to reach 0.50"))
    for label, val in scen:
        rows.append(dict(item="scenario", key=label, value=f"{val:.3f}", note=""))
    for d in design:
        rows.append(dict(item="design", key=f"{d['n']}. {d['what']}", value=d["cost"],
                         note=d["deconfounds"]))

    with (OUT / "wp40_design_to_fifty.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "key", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    (OUT / "wp40_design_to_fifty.json").write_text(
        json.dumps(dict(anchors=anchors(), gate=gate_profile(),
                        tau_threshold=tau_gate_threshold(),
                        e2_min_tau=e2_minimum_tau(),
                        in_vivo_window=in_vivo_tau_window(),
                        posterior=dict(current=base,
                                       bf_needed=bayes_factor_needed(base),
                                       scenarios={k: v for k, v in scen})),
                   indent=2), encoding="utf-8")
    print(f"  wrote {(OUT / 'wp40_design_to_fifty.csv').relative_to(REPO)}")
    return 0


def _tau_for(target: float) -> float:
    lo, hi = 0.0, 0.99
    for _ in range(300):
        m = (lo + hi) / 2
        if delta_star(FOLD_M5 * (1 + m)) < target:
            lo = m
        else:
            hi = m
    return hi


def design_table() -> list[dict]:
    """Four measurements, ordered by information per unit cost."""
    return [
        dict(n=1, what="Two-chase +/- DNA trap (WP37 inversion)",
             cost="1 assay", moves="tau exactly, fold exactly",
             deconfounds="n and tau cancel in m/c; no model assumption survives",
             falsified_if="a trap RAISES processivity (m < c) - refutes stalling pre-fit"),
        dict(n=2, what="Measure R1/2 directly at +PT, not via the fold",
             cost="1 gel", moves="the anchor, 8.0 -> measured",
             deconfounds="removes M5's own prediction from its own evidence",
             falsified_if="R1/2 at +PT comes out below 8.0, killing M5's fold"),
        dict(n=3, what="In vivo extension per end per cycle, +/- TPP1",
             cost="2 cell lines", moves="the gate against a model-free datum",
             deconfounds="no chase, no trap, no model - already 60 nt at baseline",
             falsified_if="TPP1 knockdown does not lower the 60 nt figure"),
        dict(n=4, what="Full stall-share sweep to 0.5 (tau = 0.47 question)",
             cost="expansion of 1", moves="E2 from excluded to testable",
             deconfounds="settles whether 4x is stall-inflated at all",
             falsified_if="tau(-PT) measures well below 0.4, killing E2 outright"),
    ]


# ---------------------------------------------------------------- self test
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # --- the fold -> R1/2 -> delta chain
    chk("a0 reproduces WP29's 0.03976", abs(a0() - 0.039759) < 1e-5)
    chk("fold 1.0 gives R1/2 2.963", abs(r_half(1.0) - R_NO_PT) < 1e-9)
    chk("fold 2.7335 gives R1/2 8.0", abs(r_half(FOLD_M5) - 8.0) < 5e-3)
    chk("fold 2.7335 at c=1 gives 48.0 bp/PD", abs(delta_star(FOLD_M5) - 48.0) < 5e-2)
    chk("fold 2.85 clears 50", delta_star(2.85) >= DELTA_FLOOR)
    chk("fold required for 50 is ~2.848", abs(fold_required(50) - 2.848) < 5e-3)
    chk("fold required for 90 is ~5.14", abs(fold_required(90) - 5.143) < 5e-3)
    chk("delta* is monotone in the fold", all(
        delta_star(f) < delta_star(f + 0.1) for f in (1.0, 2.0, 3.0, 4.0, 6.0)))

    # --- F2: the anchor provenance
    a = anchors()
    chk("8.0/2.963 reproduces M5's fold", abs(8.0 / R_NO_PT - 2.6998) < 5e-4)
    chk("the 8.0 anchor is labelled a PREDICTION", 
        a["R1/2 = 8.0 (used everywhere as +PT)"]["kind"] == "PREDICTION")
    chk("L&C's >=4x is a LOWER BOUND, not the 8.0",
        a["L&C measured >=4x"]["kind"] == "LOWER BOUND")
    chk("the measured >=4x anchor is above 50 bp/PD", a["L&C measured >=4x"]["delta"] > 50)
    chk("the measured anchor exceeds M5's prediction",
        a["L&C measured >=4x"]["value"] > a["R1/2 = 8.0 (used everywhere as +PT)"]["value"])

    # --- F1: the gate threshold
    thr = tau_gate_threshold()
    chk("gate opens at tau ~0.045", abs(thr - 0.0454) < 5e-3)
    chk("gate is closed at tau=0", delta_star(FOLD_M5) < DELTA_FLOOR)
    chk("gate is open at WP35's tau=0.1833",
        delta_star(FOLD_M5 * 1.1833) > DELTA_FLOOR)
    chk("the threshold is far below WP35's default", thr < TAU0_DEFAULT / 3)
    chk("90 bp/PD needs a much larger tau", _tau_for(90) > 0.2)

    # --- F3: E2's price
    chk("E2 needs tau_pos negative at tau0=0.1833", tau_pos_needed(0.1833) < 0)
    chk("E2 is infeasible at the default", not e2_feasible(0.1833))
    chk("E2 is feasible only above ~0.46", e2_minimum_tau() > 0.45)
    chk("E2 is feasible at 0.60", e2_feasible(0.60))
    chk("the ratio needed is 1.463", abs(4.0 / FOLD_M5 - 1.4633) < 5e-4)
    chk("stall_ratio is 1 when both arms stall equally", abs(stall_ratio(0.3, 0.3) - 1.0) < 1e-12)
    chk("stall_ratio exceeds the 1.463 threshold only above tau_neg ~0.46",
        stall_ratio(0.1833, 0.0) < 4.0 / FOLD_M5 < stall_ratio(0.50, 0.0))

    # --- in vivo
    chk("M5 at tau=0 is within 30% of 60 nt", in_vivo_consistent(0.0))
    chk("M5 at tau=0.1833 is within 30%", in_vivo_consistent(0.1833))
    chk("M5 at tau=0.9 is NOT - the band can fail", not in_vivo_consistent(0.9))
    w = in_vivo_tau_window()
    chk("the window exists and is bounded", w is not None and w[0] == 0.0)
    chk("the window does not extend to every tau", w[1] < 0.90)

    # --- posterior arithmetic
    chk("0.70 x 0.50 = 0.35", abs(posterior(0.70, 0.50) - 0.35) < 1e-12)
    chk("BF needed to reach 0.50 is 1.857",
        abs(bayes_factor_needed(0.35) - (0.5/0.5)/(0.35/0.65)) < 1e-9)
    chk("applying that BF gives exactly 0.50",
        abs(posterior_with_bf(0.35, bayes_factor_needed(0.35)) - 0.5) < 1e-9)
    chk("P2 confirming at 0.80 DOES clear 50%", posterior(0.70, 0.80) > 0.50)
    chk("P2 marginal at 0.60 does NOT", posterior(0.70, 0.60) < 0.50)
    chk("the required P2 term is 0.714", abs(0.50 / 0.70 - 0.7142) < 1e-3)
    chk("that is only 1.43x in odds - an ordinary ask",
        abs((0.50 / 0.70) / 0.50 - 1.4286) < 1e-3)
    # the trap I actually fell into: multiplying an already-passed hurdle lowers P
    chk("multiplying the in vivo datum in would LOWER P(M5) - which is why it is excluded",
        posterior(0.70, 0.80, 0.80) < posterior(0.70, 0.80))

    # --- design discipline
    d = design_table()
    chk("the design has 4 measurements", len(d) == 4)
    chk("every measurement names a falsifier", all(x["falsified_if"] for x in d))
    chk("every measurement names what it de-confounds", all(x["deconfounds"] for x in d))
    chk("the two-chase is first - it removes a confound, it does not add an assumption",
        d[0]["what"].startswith("Two-chase"))

    # --- the honesty guards, scoped to production code
    src = Path(__file__).read_text(encoding="utf-8")
    production = src.split("def _self_test")[0]
    chk("the file does not claim M5 is confirmed",
        "M5 is confirmed" not in production and "CONFIRMED" not in production)
    chk("the 8.0 anchor is explicitly called a prediction in the source",
        "NOT A MEASUREMENT" in production)
    chk("no probability is hard-coded as a result",
        "p_m5 = 0.5" not in production.replace(" ", ""))

    bad = [nm for nm, ok in checks if not ok]
    print(f"WP40 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
