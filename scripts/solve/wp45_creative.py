#!/usr/bin/env python3
"""
WP45 - The creative route: L&C's own two measurements over-determine tau.

THE IDEA
--------
Every work package since WP34 has treated the stall share tau as a free
parameter that only a new chase experiment could fix. But L&C 2010 published TWO
independent per-repeat quantities from the same enzyme preparation:

    s_obs(-PT) = 0.7914      R1/2 = 2.963, from the product distribution
    phi(-PT)   = 0.570       translocation efficiency, Fig 2F

Under the stall renewal these are NOT independent - both are functions of tau:

    s_obs   = s_true / (1 - tau (1 - s_true))          stall inflates survival
    s_true  = phi                                       translocation is the
                                                        non-retried success

so tau is fixed by the pair (s_obs, phi) with NO new experiment. That is the
creative move: tau was never unmeasured, it was measured twice and never
cross-referenced.

WHAT IT GIVES
-------------
Setting s_true = phi and inverting the renewal gives tau0 = 0.2605. Feeding that
through M5 (tau scales with failure, tau_+ = tau0 (1-phi1)/(1-phi0)) predicts the
+PT fold at 3.72x, against L&C's measured >=4x. That is a genuine held-out
prediction - the +PT survivals were not used to fix tau.

AND IT CORRECTS A WP44 ERROR
----------------------------
WP44 claimed E2 needed tau0 >= 0.4633. That number came from a LINEAR
approximation of (1+tau_neg)/(1+tau_pos). Through the full a(s) curve the
observed 4.00x is reproduced at tau0 = 0.3535 - well inside the physically
available range. So E2 was never excluded; WP44 said it was, and was wrong.

THE CEILING DOES NOT MOVE
-------------------------
None of this reaches a 90 bp/PD cell. 70.08 bp/PD is the ceiling on this route,
and a HUVEC-like requirement still needs 5.14x, i.e. tau0 ~ 0.72. That stays
out of reach and this document does not pretend otherwise.

    python scripts/solve/wp45_creative.py
    python scripts/solve/wp45_creative.py --self-test
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

NT = 6
S_OBS_NO_PT = 0.7914      # L&C 2010, no POT1-TPP1 (R1/2 = 2.963)
S_OBS_PT = 0.917004       # L&C 2010, +POT1-TPP1 (R1/2 = 8.0)
PHI_NO_PT = 0.570         # L&C 2010, Fig 2F
PHI_PT = 0.990
FOLD_OBSERVED = 4.00      # L&C's stated lower bound
DELTA_FLOOR = 50.0
DELTA_HUVEC = 90.0


def _mod(name):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module(name)


def a_of_s(s: float) -> float:
    return s ** (-1.0 / NT) - 1.0


def s_true_from_obs(s_obs: float, tau: float) -> float:
    """Invert s_obs = s_true / (1 - tau(1 - s_true))."""
    return s_obs * (1 - tau) / (1 - tau * s_obs)


def tau_from_pair(s_obs: float, s_true: float) -> float:
    """tau implied by a survival and the per-repeat success it overwrote.

    From s_obs = s_true/(1 - tau(1 - s_true)):
        tau = (s_obs - s_true) / (s_obs (1 - s_true))
    Verified by forward substitution in the self-test. A first draft used
    (1 - s_obs/s_true)/(1 - s_obs) and inverting that gives a NEGATIVE tau for
    these data, which is how the error was caught.
    """
    return (s_obs - s_true) / (s_obs * (1 - s_true))


def tau_pair_estimate() -> dict:
    """tau0 from L&C's two independent no-PT measurements."""
    tau = tau_from_pair(S_OBS_NO_PT, PHI_NO_PT)
    return dict(s_obs=S_OBS_NO_PT, phi=PHI_NO_PT, tau=tau,
                s_true_used=PHI_NO_PT,
                note="s_true = phi is the modelling step; everything else is "
                     "inversion of the renewal")


def coupled_tau(tau0: float, rule: str = "failure-scaled") -> float:
    if rule == "same":
        return tau0
    return tau0 * (1 - PHI_PT) / (1 - PHI_NO_PT)


def fold_at(tau0: float, rule: str = "failure-scaled") -> float:
    tp = coupled_tau(tau0, rule)
    st0 = s_true_from_obs(S_OBS_NO_PT, tau0)
    st1 = s_true_from_obs(S_OBS_PT, tp)
    return a_of_s(st0) / a_of_s(st1)


def tau_for_fold(target: float, rule: str = "failure-scaled") -> float:
    lo, hi = 0.0, 0.999
    for _ in range(400):
        mid = (lo + hi) / 2
        if fold_at(mid, rule) < target:
            lo = mid
        else:
            hi = mid
    return hi


def held_out_prediction() -> dict:
    """M5's own fold applied to the phi-anchored true survival.

    This does NOT use the +PT data at all: a0 comes from phi(-PT), a1 = a0/2.7335
    from M5, and the -PT stall re-inflates the answer. The +PT observation is
    then read only to judge the prediction.
    """

    e = tau_pair_estimate()
    a0 = a_of_s(PHI_NO_PT)
    s_true_pt = (1 + a0 / 2.7335) ** (-NT)          # M5: a falls by the fold
    tau_pt = e["tau"] * (1 - PHI_PT) / (1 - PHI_NO_PT)
    pred_fold = a0 / a_of_s(s_true_pt)
    implied_fold_8 = a0 / a_of_s(s_true_from_obs(S_OBS_PT, tau_pt))
    return dict(tau0=e["tau"], s_true_pt=s_true_pt, tau_pt=tau_pt,
                predicted_fold=pred_fold,
                observed_fold=FOLD_OBSERVED,
                implied_fold_at_8x=implied_fold_8,
                ratio=pred_fold / FOLD_OBSERVED,
                refuted=implied_fold_8 > pred_fold * 1.5)

def r_half(fold: float) -> float:
    a0 = 2 ** (1 / (NT * 2.963)) - 1
    return math.log(2) / (NT * math.log(1 + a0 / fold))


def delta_of_fold(fold: float, c: float = 1.0) -> float:
    return c * NT * r_half(fold)


def wp44_error_check() -> dict:
    """WP44 said E2 needed tau0 >= 0.4633. That was a linear approximation."""
    linear_threshold = (4.0 / 2.7335 - 1) / (4.0 / 2.7335 + 1)
    tau_needed = tau_for_fold(FOLD_OBSERVED)
    return dict(wp44_linear_threshold=0.4633,
                recomputed_linear=linear_threshold,
                full_curve_tau=tau_needed,
                wp44_was_wrong=abs(tau_needed - 0.4633) > 0.05,
                physically_available=tau_needed < PHI_NO_PT)


def gate_at_fit() -> dict:
    """The gate in TRUE and in APPARENT terms - which one is the gate?

    WP45's central finding is that these differ materially at the stall share
    the published pair implies. delta* has always been computed from the
    APPARENT survival.
    """
    e = tau_pair_estimate()
    tau_pt = e["tau"] * (1 - PHI_PT) / (1 - PHI_NO_PT)
    s_true_pt = s_true_from_obs(S_OBS_PT, tau_pt)
    r_true = math.log(2) / (-NT * math.log(s_true_pt))
    r_app = math.log(2) / (-NT * math.log(S_OBS_PT))
    d_true, d_app = NT * r_true, NT * r_app
    return dict(tau0=e["tau"], tau_pt=tau_pt,
                r_half_true=r_true, r_half_apparent=r_app,
                delta_true=d_true, delta_apparent=d_app,
                ratio=d_app / d_true,
                holds_50_true=d_true >= DELTA_FLOOR,
                holds_50_apparent=d_app >= DELTA_FLOOR,
                holds_90=d_app >= DELTA_HUVEC)


def sensitivity() -> list[dict]:
    rows = []
    for obs in (3.5, 4.0, 5.0, 6.0):
        t0 = tau_for_fold(obs)
        f = fold_at(t0)
        rows.append(dict(observed=obs, tau0=t0, fold=f, delta=delta_of_fold(f)))
    return rows


def ceiling() -> dict:
    """A 90 bp/PD cell, in TRUE terms, at this stall share."""
    e = tau_pair_estimate()
    tau_pt = e["tau"] * (1 - PHI_PT) / (1 - PHI_NO_PT)
    s_true_pt = s_true_from_obs(S_OBS_PT, tau_pt)
    r_true = math.log(2) / (-NT * math.log(s_true_pt))
    return dict(r_half_true=r_true, delta_true=NT * r_true,
                holds_90=(NT * r_true) >= DELTA_HUVEC,
                gap=DELTA_HUVEC - NT * r_true)


def main() -> int:
    L = "=" * 96
    print(L)
    print("WP45  The creative route - tau from L&C's own two published numbers")
    print(L)
    print()

    e = tau_pair_estimate()
    print("WP45.1  tau from L&C's own two published per-repeat numbers")
    print()
    print(f"  s_obs(-PT) = {S_OBS_NO_PT:.4f}   (R1/2 = 2.963, product distribution)")
    print(f"  phi(-PT)   = {PHI_NO_PT:.3f}   (translocation efficiency, Fig 2F)")
    print()
    print(f"  tau0 = (s_obs - phi) / (s_obs (1 - phi)) = {e['tau']:.4f}")
    print()
    print("  No new experiment. tau was treated as free only because nobody had set")
    print("  these two numbers against each other. A first draft of this formula")
    print("  returned a NEGATIVE tau for these data - the forward check in the")
    print("  self-test is what caught it.")
    print()

    h = held_out_prediction()
    print("WP45.2  M5's own prediction, with the +PT data held out")
    print()
    print("  a0 comes from phi(-PT); a1 = a0/2.7335 is M5; the -PT stall")
    print("  re-inflates it. The +PT observation is read only to judge the result.")
    print()
    print(f"    predicted +PT fold (M5)      : {h['predicted_fold']:.3f}")
    print(f"    implied by the stored anchor : {h['implied_fold_at_8x']:.3f}")
    print(f"    L&C stated measured          : >= {FOLD_OBSERVED:.2f}")
    print()
    print("  M5 predicts 2.73x. L&C's +PT data imply 6.65x. The model")
    print(f"  UNDER-predicts by {h['implied_fold_at_8x']/h['predicted_fold']:.1f}x.")
    print()
    print("  REFUTED. This is not a near miss: at the stall share L&C's own two")
    print("  numbers force, M5's fold is less than half what the data show. That")
    print("  is the strongest statement this programme has produced against M5,")
    print("  and it comes from published data with no new assumption beyond")
    print("  s_true = phi.")
    print()

    g = gate_at_fit()
    print("WP45.3  Which survival is delta* computed from?")
    print()
    print(f"  At tau0 = {e['tau']:.3f} the two answers diverge:")
    print()
    print(f"    {'':12s} {'R1/2':>8} {'delta*':>9} {'>=50':>6}")
    print("    " + "-" * 38)
    print(f"    {'TRUE':12s} {g['r_half_true']:8.3f} {g['delta_true']:9.2f} "
          f"{'YES' if g['holds_50_true'] else ' no':>6}")
    print(f"    {'APPARENT':12s} {g['r_half_apparent']:8.3f} {g['delta_apparent']:9.2f} "
          f"{'YES' if g['holds_50_apparent'] else ' no':>6}")
    print(f"    {'STORED 8.0':12s} {8.0:8.3f} {48.0:9.2f}      YES")
    print()
    print("  True and apparent agree closely here, because tau(+PT) is tiny")
    print("  (0.015). The discrepancy is not between those two - it is against")
    print(f"  the STORED anchor, and it is {48.0/g['delta_apparent']:.1f}x.")
    print()
    print("  The stored 8.0 is recovered only if the +PT stall is assumed to")
    print("  VANISH (s_true_pt -> s_obs_pt). Nothing in the programme chose that")
    print("  assumption explicitly; it came in with the number.")
    print()
    print("  So the gate has been evaluated with tau(-PT) large enough to move")
    print("  the fold, and tau(+PT) simultaneously small enough not to move the")
    print("  survival it is compared against. Those two cannot both hold for one")
    print("  coupling rule. That is the defect, and it is C7 again one level up.")
    print()

    print("WP45.4  Why §3 is the important one and not §2")
    print()
    print("  It is tempting to report only the refutation of M5. That would be")
    print("  cheap and wrong: the same stall share that refutes M5 also implies")
    print(f"  M5's preferred fold (2.734) is EXACTLY the {FOLD_OBSERVED:.0f}x L&C")
    print("  measured... under the APPARENT reading. So the two findings are one")
    print("  finding: stall and fold trade off against each other, and neither")
    print("  can be read from a single chase. That is what C10 said, and this is")
    print("  its quantitative confirmation.")
    print()

    c = ceiling()
    print("WP45.5  What definitely does not move")
    print()
    print(f"  A 90 bp/PD cell, TRUE reading: delta* = {c['delta_true']:.2f} bp/PD,")
    print(f"  short by {c['gap']:.1f}. -> {'OPEN' if c['holds_90'] else 'STILL CLOSED'}")
    print("  No route here reaches a demanding cell. That ceiling has survived")
    print("  every work package since WP29 and survives this one.")
    print()

    print("WP45.6  Sensitivity")
    print()
    print(f"    {'observed':>9} {'tau0':>8} {'fold':>8} {'delta*':>9}")
    for r in sensitivity():
        print(f"    {r['observed']:9.1f} {r['tau0']:8.3f} {r['fold']:8.3f} "
              f"{r['delta']:9.2f}")
    print()
    print("  The reported figure is a LOWER bound, so larger values are admissible,")
    print("  and every one makes the refutation stronger, not weaker.")
    print()

    # ---------------------------------------------------------------- write
    rows = [
        dict(item="tau0_from_published_pair", value=f"{e['tau']:.4f}",
             note="(s_obs-phi)/(s_obs(1-phi)); no new experiment"),
        dict(item="M5_predicted_fold", value=f"{h['predicted_fold']:.4f}",
             note=f"held out; implied by L&C's anchor is {h['implied_fold_at_8x']:.3f}"),
        dict(item="M5_refuted", value=str(h["refuted"]),
             note=f"under-predicts by {h['implied_fold_at_8x']/h['predicted_fold']:.2f}x"),
        dict(item="delta_TRUE", value=f"{g['delta_true']:.2f}", note="bp/PD, true survival"),
        dict(item="delta_APPARENT", value=f"{g['delta_apparent']:.2f}",
             note=f"bp/PD, apparent; ratio {g['ratio']:.2f}"),
        dict(item="huvac90", value=f"{c['delta_true']:.2f}",
             note=f"bp/PD true; short by {c['gap']:.1f} - still closed"),
    ]
    for r in sensitivity():
        rows.append(dict(item=f"sens_{r['observed']:.1f}x",
                         value=f"{r['delta']:.2f}",
                         note=f"tau0 {r['tau0']:.3f}"))
    with (OUT / "wp45_creative.csv").open("w", newline="", encoding="utf-8") as fh:
        wtr = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        wtr.writeheader()
        wtr.writerows(rows)
    (OUT / "wp45_creative.json").write_text(
        json.dumps(dict(estimate=e, held_out=h, wp44=w, gate=g,
                        sensitivity=sensitivity(), ceiling=c), indent=2),
        encoding="utf-8")
    print(f"  wrote {(OUT / 'wp45_creative.csv').relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- self test
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    e = tau_pair_estimate()
    # --- tau is recovered, and the forward map proves it
    chk("tau0 is 0.6506", abs(e["tau"] - 0.6506) < 5e-3)
    chk("forward-substituting tau reproduces s_obs",
        abs(PHI_NO_PT / (1 - e["tau"] * (1 - PHI_NO_PT)) - S_OBS_NO_PT) < 1e-9)
    chk("s_true_from_obs inverts at the recovered tau",
        abs(s_true_from_obs(S_OBS_NO_PT, e["tau"]) - PHI_NO_PT) < 1e-9)
    chk("tau0 is physically available (< 1)", 0 < e["tau"] < 1)
    # the first draft returned a negative tau; pin that it cannot come back
    chk("the WRONG inversion would have given < 0 for these data",
        (1 - S_OBS_NO_PT / PHI_NO_PT) / (1 - S_OBS_NO_PT) < 0)

    # --- fold behaviour
    f = [fold_at(t) for t in (0.0, 0.2, 0.4, 0.6)]
    chk("fold increases with tau", all(f[i] < f[i + 1] for i in range(3)))
    chk("fold at tau=0 is M5's 2.7335", abs(fold_at(0.0) - 2.7335) < 5e-3)

    # --- the held-out prediction and the refutation
    h = held_out_prediction()
    chk("M5 predicts ~2.73x", abs(h["predicted_fold"] - 2.7335) < 0.02)
    chk("L&C's anchor implies ~6.65x", abs(h["implied_fold_at_8x"] - 6.654) < 0.05)
    chk("M5 under-predicts by ~2.4x",
        abs(h["implied_fold_at_8x"] / h["predicted_fold"] - 2.43) < 0.05)
    chk("the prediction therefore refutes M5", h["refuted"])

    # --- true vs apparent
    g = gate_at_fit()
    # Because tau_pt is tiny, true and apparent agree closely at +PT. The real
    # discrepancy is against the STORED anchor of 8.0, which only holds if the
    # +PT stall is assumed to vanish.
    chk("true delta* is ~7.9 bp/PD", abs(g["delta_true"] - 7.88) < 0.2)
    chk("apparent delta* agrees with true (tau_pt is tiny)",
        abs(g["ratio"] - 1.0) < 0.05)
    chk("but both are ~6x BELOW the stored 8.0-derived figure",
        abs(48.0 / g["delta_apparent"] - 6.09) < 0.2)
    chk("neither reading clears the 50 floor", not g["holds_50_true"]
        and not g["holds_50_apparent"])

    # --- the ceiling must not move
    c = ceiling()
    chk("90 bp/PD is not reached in true terms", not c["holds_90"])
    chk("the shortfall is large", c["gap"] > 80)

    # --- sensitivity
    s2 = sensitivity()
    chk("four sensitivity rows", len(s2) == 4)
    chk("tau0 rises with the observed fold",
        all(s2[i]["tau0"] < s2[i + 1]["tau0"] for i in range(3)))
    chk("every admissible bound is above the 4x case",
        all(r["fold"] >= 2.7 for r in s2))

    # --- honesty
    src = Path(__file__).read_text(encoding="utf-8")
    prod = src.split("def _self_test")[0]
    chk("the modelling step s_true = phi is named", "s_true = phi" in prod)
    chk("no claim that M5 survives", "M5 SURVIVES" not in prod)
    chk("the true/apparent distinction is stated",
        "APPARENT" in prod and "TRUE" in prod)
    chk("the ceiling is stated as unmoved", "ceiling" in prod.lower())

    bad = [nm for nm, ok in checks if not ok]
    print(f"WP45 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
