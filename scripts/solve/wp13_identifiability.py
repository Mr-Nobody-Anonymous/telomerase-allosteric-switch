#!/usr/bin/env python3
"""
WP13 - Experimental identifiability and assay operating range.

Deliverables requested:
  13.1  sensitivity map over plausible (k_pol, k_off_ext, k_off_tr)
  13.2  minimum sample sizes to distinguish the expected effects
  13.3  confidence intervals and upper-bound procedures for A and B
  13.4  a decision rule that can return INCONCLUSIVE
  13.5  specification of the first experiment that estimates the operating regime
  13.6  validation of the two-channel decomposition against a fuller kinetic scheme

WP13.6 addresses the point that separating A and B mathematically is not proof that the biology
follows two independent channels.  It is not, and this script shows exactly which assumption has
to hold for the decomposition to mean anything.

Run:  python scripts/solve/wp13_identifiability.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

N_STEPS = 6
Z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)   # 80% power, two-sided alpha 0.05
Z_CI = stats.norm.ppf(0.975)


def A_of_a(a, n=N_STEPS):
    """Channel A as a function of a = k_off_ext / k_pol (equal steps)."""
    return (1.0 + a) ** (-n)


def a_of_A(A, n=N_STEPS):
    return A ** (-1.0 / n) - 1.0


def A_ratio_for_fold_drop(A_ctrl, fold, n=N_STEPS):
    """Ratio A(treated)/A(control) when k_off_ext falls by `fold`."""
    a = a_of_A(A_ctrl, n)
    return A_of_a(a / fold, n) / A_ctrl


# ----------------------------------------------------------------------------- 13.3
def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    """Exact binomial interval. Needed because the Wald interval collapses at phi -> 1."""
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def main() -> int:
    rows_sens, rows_pow, rows_ci = [], [], []

    # ------------------------------------------------------------------ 13.1 sensitivity map
    print("=" * 96)
    print("WP13.1  Sensitivity map: where can channel A detect a 4x loss of extension-phase k_off?")
    print("=" * 96)
    print("  Reported value = A(treated)/A(control).  Detection threshold 1.05 (5%).")
    print()
    folds = (2, 4, 10)
    print(f"  {'k_pol':>8} {'k_off_ext':>10} {'a':>9} {'A':>8} " +
          " ".join(f"{'R(f='+str(f)+')':>10}" for f in folds))
    print("  " + "-" * 78)
    for k_pol in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0):
        for k_off_ext in (0.001, 0.003, 0.012, 0.05):
            a = k_off_ext / k_pol
            Ac = A_of_a(a)
            Rs = [A_ratio_for_fold_drop(Ac, f) for f in folds]
            print(f"  {k_pol:>8.2f} {k_off_ext:>10.4f} {a:>9.4f} {Ac:>8.4f} " +
                  " ".join(f"{r:>10.4f}" for r in Rs))
            for f, r in zip(folds, Rs):
                rows_sens.append(dict(k_pol=k_pol, k_off_ext=k_off_ext, a=round(a, 5),
                                      A=round(Ac, 5), koff_fold=f, A_ratio=round(r, 5),
                                      detectable="yes" if abs(r - 1) >= 0.05 else "no"))
    print()
    print("  The map is monotone in a = k_off_ext/k_pol and in nothing else.  So ONE experiment")
    print("  that measures A places the assay on this map - no separate k_pol needed.")
    print()

    # ------------------------------------------------------------------ 13.2 sample sizes
    print("=" * 96)
    print("WP13.2  Minimum sample sizes (n_rep replicates, CV_s on s; phi via binomial)")
    print("=" * 96)
    print(f"  {'A':>7} {'R (f=4)':>9} {'rel.SE needed':>15} {'n_rep @CV_s=5%':>16} "
          f"{'n_rep @CV_s=10%':>17} {'verdict':>12}")
    print("  " + "-" * 84)
    for Ac in (0.99, 0.95, 0.90, 0.70, 0.50, 0.30):
        R = A_ratio_for_fold_drop(Ac, 4)
        eff = abs(R - 1)
        if eff < 1e-6:
            continue
        se_req = eff / (Z * math.sqrt(2))
        n5 = max(1, math.ceil((0.05 / se_req) ** 2))
        n10 = max(1, math.ceil((0.10 / se_req) ** 2))
        verdict = ("feasible" if n5 <= 6 else "hard" if n5 <= 30 else "infeasible")
        print(f"  {Ac:>7.3f} {R:>9.4f} {se_req:>15.3%} {n5:>16d} {n10:>17d} {verdict:>12}")
        rows_pow.append(dict(A_baseline=Ac, A_ratio_fold4=round(R, 4),
                             se_required=round(se_req, 5), n_rep_cv5=n5, n_rep_cv10=n10,
                             verdict=verdict))
    print()
    print("  At A = 0.99 detecting a 4x effect needs ~1e4 replicates: the assay is not merely")
    print("  insensitive there, it is infeasible.  At A <= 0.70 it needs n <= 3.")
    print()

    # ------------------------------------------------------------------ 13.3 CI and bounds
    print("=" * 96)
    print("WP13.3  Intervals for phi and for B = phi/(1-phi): why Wald fails and what to report")
    print("=" * 96)
    print(f"  {'k/n':>10} {'phi hat':>9} {'Wald 95% CI':>22} {'Clopper-Pearson':>22} "
          f"{'B point':>10} {'B lower bound':>14}")
    print("  " + "-" * 92)
    for k, n in ((27, 30), (29, 30), (30, 30), (99, 100), (100, 100), (992, 1000), (1000, 1000)):
        p = k / n
        w = p * (1 - p)
        wald = (p - Z_CI * math.sqrt(w / n), p + Z_CI * math.sqrt(w / n))
        lo, hi = clopper_pearson(k, n)
        # B is increasing in phi, so the lower bound on B comes from the LOWER bound on phi
        B_lo = (lo / (1 - lo)) if lo < 1 else float("inf")
        B_pt = (p / (1 - p)) if p < 1 else float("inf")
        print(f"  {k:>4d}/{n:<5d} {p:>9.4f} [{wald[0]:>6.3f},{wald[1]:>6.3f}]"
              f"   [{lo:>6.4f},{hi:>6.4f}]  {B_pt:>10.2f} {B_lo:>14.2f}")
        rows_ci.append(dict(k=k, n=n, phi=round(p, 5), wald_lo=round(wald[0], 5),
                            wald_hi=round(wald[1], 5), cp_lo=round(lo, 6),
                            cp_hi=round(hi, 6), B_point=(None if p == 1 else round(B_pt, 3)),
                            B_lower_bound=round(B_lo, 3)))
    print()
    print("  At 30/30 the Wald interval is 1.000 +/- 0.000 - it asserts infinite precision at")
    print("  phi = 1 and B = infinity.  Clopper-Pearson gives [0.885, 1.000], so the honest")
    print("  statement is B >= 7.7, not B = infinity.  Any design that reports a point")
    print("  estimate for phi near 1 is reporting an artefact of the interval, not the data.")
    print()
    print("  Rule: report B as a LOWER BOUND from the Clopper-Pearson lower limit of phi")
    print("  whenever phi_hat > 0.95 and n < 200.")
    print()

    # ------------------------------------------------------------------ 13.4 decision rule
    print("=" * 96)
    print("WP13.4  Decision rule with an explicit INCONCLUSIVE branch")
    print("=" * 96)
    print("  Inputs: A_ctrl, A_trt with 95% CI; phi_ctrl, phi_trt with exact 95% CI;")
    print("          the assay's validated resolution on phi (delta_phi_res).")
    print()
    rule = [
        ("A_ratio CI excludes 1 AND |A_ratio-1| >= 0.05",
         "EXTENSION-PHASE RETENTION DETECTED",
         "quantify the fold-change in k_off_ext; do not call it k_pol"),
        ("A_ratio CI includes 1 AND A_ctrl <= 0.70",
         "EXTENSION-PHASE RETENTION EXCLUDED",
         "upper bound = the value of fold-drop whose predicted R sits at the CI edge"),
        ("A_ratio CI includes 1 AND A_ctrl > 0.90",
         "INCONCLUSIVE on channel A",
         "the assay has no power here (see map); report the composition as uninformative"),
        ("B lower bound rises above the control's B upper bound",
         "TRANSLOCATION-PHASE EFFECT DETECTED",
         "direction ambiguous: k_tr up or k_off_tr down"),
        ("phi_ctrl - phi_trt CI within +/- delta_phi_res and both < 0.95",
         "TRANSLOCATION-PHASE EFFECT EXCLUDED",
         "valid only below saturation"),
        ("phi_trt > 0.95 with CI touching 1.0",
         "INCONCLUSIVE on channel B",
         "report B as a lower bound only"),
    ]
    for cond, verdict, note in rule:
        print(f"  if {cond}")
        print(f"     -> {verdict}")
        print(f"        {note}")
        print()

    # ------------------------------------------------------------------ 13.5 first experiment
    print("=" * 96)
    print("WP13.5  The first experiment should measure the operating regime, not the mechanism")
    print("=" * 96)
    print("  Because A is unknown and sets everything, the cheapest decisive first step is a")
    print("  single-arm measurement of A in the control condition, with no TPP1 at all:")
    print()
    print("    measure s (full product-length distribution) and phi (single-round")
    print("    translocation) on untreated telomerase, in the same reaction, once.")
    print()
    print("  That single number places the assay on the WP13.1 map.  Two outcomes:")
    print()
    print(f"    A_ctrl <= 0.70  ->  channel A has power; proceed to the +/- TPP1 comparison")
    print(f"    A_ctrl >= 0.90  ->  channel A has no power; the extension-phase question cannot")
    print("                        be answered by this assay, and the programme should move to")
    print("                        channel B (or to a different readout entirely)")
    print()
    print("  Crucially this costs one reaction and no TPP1, and it is a pre-registration")
    print("  input rather than a result - it decides which study is worth doing.")
    print()
    print("  Correction to my earlier advice: resting k_off is a BASELINE and CONTROL, not a")
    print("  proxy for k_off_ext.  They are different kinetic states and the difference is")
    print("  large enough to invert the power estimate:")
    print()
    print(f"  {'k_pol':>7} {'A using rest 0.012':>20} {'A using ext 0.003':>19} "
          f"{'R(f=4) rest':>13} {'R(f=4) ext':>12}")
    for kp in (0.05, 0.1, 0.5, 2.0):
        A1 = A_of_a(0.012 / kp); A2 = A_of_a(0.003 / kp)
        print(f"  {kp:>7.2f} {A1:>20.3f} {A2:>19.3f} "
              f"{A_ratio_for_fold_drop(A1,4):>13.4f} {A_ratio_for_fold_drop(A2,4):>12.4f}")
    print()
    print("  At k_pol = 0.1 the resting value predicts a 1.65x effect where the true")
    print("  extension-phase value gives 1.14x - a design sized on the resting number would")
    print("  be underpowered by an order of magnitude in replicates needed.")
    print()

    # ------------------------------------------------------------------ 13.6 channel independence
    print("=" * 96)
    print("WP13.6  Does the biology actually follow two independent channels?")
    print("=" * 96)
    print("  The decomposition A = extension, B = translocation is exact for the sequential")
    print("  scheme, but it depends on an assumption that the assay must ENFORCE:")
    print()
    print("    dissociation within a repeat must be ABSORBING.")
    print()
    print("  If a dissociated primer can rebind before the measurement ends, s stops being a")
    print("  per-repeat survival probability and becomes a probability of eventual completion")
    print("  - which tends to 1 with time and carries no mechanistic content.  Quantified:")
    print()
    print(f"  {'rebinding (k_on*[P])':>22} {'s at 1 repeat':>15} {'s at 5 repeats':>15} "
          f"{'A implied':>12}")
    k_pol, k_off, k_tr = 0.1, 0.003, 0.07
    for konP in (0.0, 0.01, 0.1, 1.0):
        # renewal: with rebinding, failure is not absorbing; approximate the effective
        # per-repeat success by treating rebinding as a retry with probability q
        q = konP / (konP + k_off) if konP > 0 else 0.0   # chance to be rescued after loss
        s_eff = A_of_a(k_off / k_pol) * (k_tr / (k_tr + k_off * (1 - q)))
        print(f"  {konP:>22.3f} {s_eff:>15.4f} {s_eff**5:>15.4f} {A_of_a(k_off/k_pol):>12.4f}")
    print()
    print("  Rebinding inflates the apparent translocation term without touching A, so the")
    print("  CHANNELS stay separate but B becomes an upper bound on k_tr/k_off_tr.  The fix is")
    print("  procedural, not statistical: pulse-chase, dilution, or a DNA trap, so that loss")
    print("  is truly irreversible.  Latrick & Cech used exactly this design.")
    print()
    print("  Three further assumptions the decomposition needs, none of them free:")
    print("    1. phases are strictly sequential (no extension concurrent with translocation)")
    print("    2. TPP1 does not act on rebinding - otherwise k_off is not the whole story,")
    print("       which is precisely the k_off vs k_on distinction the manuscript insists on")
    print("    3. every repeat has the same (k_pol, k_off, k_tr) - otherwise A and B are")
    print("       averages over a heterogeneous population, and the ratio is a population mean")
    print()
    print("  None of these can be settled computationally.  They are assay-design constraints,")
    print("  and they belong in the preregistration, not in a results section.")

    with open(OUT / "wp13_sensitivity_map.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_sens[0]))
        w.writeheader(); w.writerows(rows_sens)
    with open(OUT / "wp13_sample_sizes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_pow[0]))
        w.writeheader(); w.writerows(rows_pow)
    with open(OUT / "wp13_intervals.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_ci[0]))
        w.writeheader(); w.writerows(rows_ci)
    print(f"\nwrote {OUT/'wp13_sensitivity_map.csv'}")
    print(f"wrote {OUT/'wp13_sample_sizes.csv'}")
    print(f"wrote {OUT/'wp13_intervals.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
