#!/usr/bin/env python3
"""
WP15 - Experimental design validation package.

Three priorities:
  15.0  Replace "survives the tested profiles" with a PROOF over all profiles   (new)
  15.1  Active-fraction calibration w: estimand, identifiability, uncertainty
  15.2  Trap validation: multiple methods, five readouts, method-dependence
  15.3  Hierarchical power analysis with purification as the experimental unit

Run:  python scripts/solve/wp15_design_validation.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats, optimize

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
RNG = np.random.default_rng(1515)
Z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)


# ------------------------------------------------------------------ 15.0 proven bounds
def h(u, f):
    return np.log(1.0 + (f - 1.0) * np.exp(-u))


def log_R(u, f, n):
    """log of A(treated)/A(control) for per-step ratios given by u_i = log(1+a_i)."""
    return n * np.log(f) - h(u, f).sum()


def proven_bounds(A0, f, n=6):
    """Exact min and max of R over ALL vectors u >= 0 with sum(u) = -log A0."""
    L = -math.log(A0)
    R_max = f ** n * math.exp(-n * h(L / n, f))      # uniform
    R_min = f * math.exp(-h(L, f))                    # concentrated on one step
    return R_min, R_max


def main() -> int:
    rows_bounds, rows_act, rows_trap, rows_hier = [], [], [], []

    # ------------------------------------------------------------------ 15.0
    print("=" * 96)
    print("WP15.0  From 'tested profiles' to a PROOF over all profiles")
    print("=" * 96)
    print("  With u_i = log(1 + a_i) and sum_i u_i = L = -log(A0):")
    print()
    print("      log R(f) = n*log(f) - sum_i h(u_i),     h(u) = log(1 + (f-1)e^-u)")
    print()
    print("  h is convex in u, so -sum h attains its maximum at u uniform and its MINIMUM at")
    print("  the most concentrated feasible vector.  Bounds, exactly:")
    print()
    print("      R_max = f^n * exp(-n*h(L/n))      [all steps equal]")
    print("      R_min = f   * exp(-h(L))          [all survival in ONE step]")
    print()
    print(f"  {'A0':>6} {'R_min':>9} {'R_max':>9} {'spread':>8}  interpretation")
    for A0 in (0.99, 0.95, 0.90, 0.80, 0.70, 0.50):
        Rmin, Rmax = proven_bounds(A0, 4.0)
        print(f"  {A0:>6.2f} {Rmin:>9.4f} {Rmax:>9.4f} {100*(Rmax/Rmin-1):>7.2f}%")
        rows_bounds.append(dict(A0=A0, fold=4.0, R_min=round(Rmin, 6), R_max=round(Rmax, 6),
                                spread_pct=round(100 * (Rmax / Rmin - 1), 3)))
    print()
    # numerical confirmation over a wide family
    n, f, A0 = 6, 4.0, 0.70
    L = -math.log(A0)
    Rmin, Rmax = proven_bounds(A0, f)
    worst_lo, worst_hi = 1e9, -1e9
    for scale in (0.2, 0.3, 1.0, 5.0, 20.0):
        for _ in range(40000):
            w = RNG.dirichlet(np.ones(n) * scale)
            r = math.exp(log_R(w * L, f, n))
            worst_lo, worst_hi = min(worst_lo, r), max(worst_hi, r)
    for _ in range(40000):
        w = np.exp(RNG.normal(0, 4, n)); w /= w.sum()
        r = math.exp(log_R(w * L, f, n))
        worst_lo, worst_hi = min(worst_lo, r), max(worst_hi, r)
    print(f"  240,000 random profiles across 5 Dirichlet concentrations + heavy log-normal:")
    print(f"    observed range [{worst_lo:.6f}, {worst_hi:.6f}]")
    print(f"    proven bounds   [{Rmin:.6f}, {Rmax:.6f}]   -> contained: "
          f"{worst_lo >= Rmin - 1e-9 and worst_hi <= Rmax + 1e-9}")
    print()
    print("  CONCLUSION: at A0 = 0.70 a 4-fold fall in extension-phase k_off gives R in")
    print(f"  [{Rmin:.4f}, {Rmax:.4f}] for EVERY within-population heterogeneity — a spread of")
    print(f"  {100*(Rmax/Rmin-1):.2f}%.  This is no longer conditional on the tested profiles.")
    print("  The earlier caveat ('tested profiles do not establish universality') is discharged")
    print("  for within-population heterogeneity; it still stands for POPULATION MIXTURES,")
    print("  which are a different model and remain covered by 15.1.")
    print()

    # ------------------------------------------------------------------ 15.1 active fraction
    print("=" * 96)
    print("WP15.1  Active fraction w: what is estimable, and what is assumed")
    print("=" * 96)
    print("  Observed A is a survival-weighted mean of subpopulation survivals:")
    print("      A_obs = sum_j w_j * A_j")
    print("  This is one equation in {w_j, A_j}.  It is NOT identifiable from A alone —")
    print("  confirming the audit: w cannot be inferred from A without circular reasoning.")
    print()
    def n_for_mixture(R):
        eff = abs(R - 1)
        if eff <= 1e-9:
            return None
        return max(3, math.ceil((0.05 * Z / (eff / math.sqrt(2))) ** 2))

    print("  What IS estimable, given the mixture model:")
    print()
    # two-population model: A_obs = w*A_r + (1-w)*1
    print("    A_obs = w*A_r + (1-w)      =>      w = (1 - A_obs)/(1 - A_r)")
    print("  so w is identified ONLY up to A_r, the responsive subpopulation's survival.")
    print()
    print("  The dangerous case is not a low A_obs.  It is a LOW A_obs produced by a mixture")
    print("  that contains a subpopulation in the no-power regime — the gate passes and the")
    print("  study proceeds while part of the enzyme cannot be measured at all.")
    print()

    def n_for(A0, fold=4.0):
        Rmin, _ = proven_bounds(A0, fold)
        eff = Rmin - 1
        if eff <= 1e-9:
            return None
        return max(3, math.ceil((0.05 * Z / (eff / math.sqrt(2))) ** 2))

    print(f"  {'subpopulation A values':>26} {'w':>6} {'A_obs':>7} {'gate':>6} "
          f"{'R_obs':>7} {'n_obs':>6} {'n for worst subpop':>19}")
    print("  " + "-" * 92)
    cases = [
        (0.30, 0.95, 0.5, "50/50 of A=0.30 and A=0.95"),
        (0.30, 0.99, 0.5, "50/50 of A=0.30 and A=0.99"),
        (0.20, 0.95, 0.6, "60/40 of A=0.20 and A=0.95"),
        (0.30, 0.70, 0.5, "50/50 of A=0.30 and A=0.70 (both powered)"),
        (0.70, 1.00, 0.5, "50/50 of A=0.70 and 100% inactive"),
    ]
    for a1, a2, w, lab in cases:
        A_obs = w * a1 + (1 - w) * a2
        R = (w * (proven_bounds(a1, 4.0)[0]) + (1 - w) * (proven_bounds(a2, 4.0)[0])) / A_obs
        pass_gate = A_obs <= 0.70
        n_obs = n_for_mixture(R)
        cand = [x for x in (n_for(a1), n_for(a2), n_for_mixture(R)) if x]
        n_worst = max(cand) if cand else None
        print(f"  {lab:>26} {w:>6.2f} {A_obs:>7.3f} {'PASS' if pass_gate else 'fail':>6} "
              f"{R:>7.4f} {str(n_obs):>6} {str(n_worst):>19}")
        rows_act.append(dict(mixture=lab, w_a1=w, A_1=a1, A_2=a2, A_obs=round(A_obs, 4),
                             gate="pass" if pass_gate else "fail", R_obs=round(R, 4),
                             n_obs=n_obs, n_worst_subpopulation=n_worst))
    print()
    print("  Rows 1-3 are the hazard: A_obs passes the A <= 0.70 gate, yet the mixture")
    print("  contains a subpopulation at A >= 0.95 that would need n = 26 to 688 to see ITS")
    print("  effect.  The aggregate looks well powered (n = 3) while part of the enzyme is")
    print("  unmeasurable.  Row 5 is the mirror case: the gate rejects a preparation that is")
    print("  half fully responsive.")
    print()
    print()
    print("  Two consequences:")
    print("    1. the aggregate A can pass the gate while a subpopulation is unmeasurable, so")
    print("       w must be known BEFORE the study is sized, not inferred from A after;")
    print("    2. w must come from an independent measurement. Three options, ranked:")
    print("       (a) active-site titration / activity standard curve against a purified")
    print("           standard -> w with a CI.  Direct, but calibration-dependent.")
    print("       (b) TPP1 binding stoichiometry (bound TPP1 per TERT) -> an upper bound on w.")
    print("       (c) product-length curvature across a dilution series: a mixture shows")
    print("           non-geometric decay that a homogeneous population cannot.  Useful as a")
    print("           consistency check, NOT as the primary estimate — it shares information")
    print("           with A and is therefore partly circular.")
    print()
    print("  Required reporting: w with a 95% CI, and the method.  If the CI includes values")
    print("  below 0.8, the feasibility gate is NOT applied.")
    print()

    # small numerical demonstration of (c) curvature
    print("  Demonstration of (c): per-repeat survival measured across an input-concentration")
    print("  series is flat for a homogeneous population and rises toward 1 for a mixture")
    print("  (more enzyme molecules available to rescue a given primer).")
    print()
    print(f"  {'input conc (rel.)':>18} {'A homogeneous':>15} {'A 60/40 mixture':>17}")
    A_h, A_mix = 0.70, 0.60
    for c in (0.25, 0.5, 1.0, 2.0, 4.0):
        # crude occupancy model: fraction of primers with an enzyme engaged
        occ = c / (c + 1.0)
        print(f"  {c:>18.2f} {A_h:>15.3f} {occ*A_mix + (1-occ)*1.0:>17.3f}")
    print("  (indicative only; the exact curvature depends on the assembly model)")
    print()

    # ------------------------------------------------------------------ 15.2 trap validation
    print("=" * 96)
    print("WP15.2  Trap validation: five readouts, two methods, and method-dependence")
    print("=" * 96)
    print("  Acceptance: each readout within +/- 10% of the no-trap control, for BOTH methods.")
    print("  Any readout failing in either method blocks the B-channel conclusion.")
    print()
    readouts = ["basal nucleotide addition", "TPP1 binding (Kd)", "active concentration",
                "primer dissociation", "product-length distribution (s)"]
    methods = ["dilution (100x)", "excess DNA trap", "heparin trap"]
    print(f"  {'readout':36s} " + " ".join(f"{m:>18}" for m in methods))
    print("  " + "-" * 96)
    # illustrative: expected disturbance by mechanism of action
    expected = {
        "basal nucleotide addition": (0.02, 0.05, 0.25),
        "TPP1 binding (Kd)": (0.03, 0.06, 0.15),
        "active concentration": (0.05, 0.03, 0.10),
        "primer dissociation": (0.90, 0.95, 0.40),
        "product-length distribution (s)": (0.04, 0.07, 0.20),
    }
    for r in readouts:
        vals = expected[r]
        cells = " ".join(f"{(1-v)*100:>17.1f}%" for v in vals)
        print(f"  {r:36s} {cells}")
        rows_trap.append(dict(readout=r, dilution_pct_retained=round((1 - vals[0]) * 100, 1),
                              dna_trap_pct_retained=round((1 - vals[1]) * 100, 1),
                              heparin_pct_retained=round((1 - vals[2]) * 100, 1)))
    print()
    print("  Read as: percent of the no-trap value retained.  Values are ILLUSTRATIVE —")
    print("  the point is the structure of the requirement, not the numbers. Two properties")
    print("  matter and they are different:")
    print()
    print("    dilution and excess-DNA trap remove rebinding without binding the enzyme;")
    print("    heparin binds the enzyme and is the most likely to alter k_pol or Kd.")
    print()
    print("  Because a trap is applied to BOTH arms, its systematic effect does not cancel.")
    print("  Preference order: dilution > excess DNA >> heparin.  If two methods disagree")
    print("  beyond their CIs, the inferred dissociation is method-dependent and neither is")
    print("  reportable as k_tr/k_off_tr; report agreement/disagreement as the result.")
    print()

    # ------------------------------------------------------------------ 15.3 hierarchical power
    print("=" * 96)
    print("WP15.3  Hierarchical power: purification as the experimental unit")
    print("=" * 96)
    print("  Model: paired within purification (each prep gives control and treated), with a")
    print("  purification-level random slope on the effect:")
    print()
    print("      log A_trt[p] - log A_ctrl[p] = mu + b_p + e_p,")
    print("      b_p ~ N(0, tau^2)   [prep-to-prep variation in the EFFECT]")
    print("      e_p ~ N(0, 2*sigma^2)  [within-prep measurement error, both arms]")
    print()
    print("  The paired estimate has variance  (tau^2 + 2*sigma^2) / P  over P purifications.")
    print("  If tau^2 > 0, extra REACTIONS per preparation buy nothing: only more")
    print("  purifications reduce the error.  This is the point the audit raised.")
    print()
    print("  Shown at two operating points: A=0.70 (large effect, heterogeneity irrelevant)")
    print("  and A=0.90 (small effect, where prep-level variance dominates).")
    print()
    sigma = 0.05
    for A0 in (0.70, 0.90):
        logR_eff = math.log(proven_bounds(A0, 4.0)[0])
        print(f"  A0 = {A0:.2f}   (log effect = {logR_eff:.4f}, sigma = {sigma})")
        print(f"    {'tau (prep SD on effect)':>24} {'P needed':>10} {'P if tau=0':>12} "
              f"{'inflation':>11}")
        for tau in (0.0, 0.02, 0.05, 0.10, 0.20):
            var = tau ** 2 + 2 * sigma ** 2
            n_units = max(2, math.ceil((Z ** 2) * var / (logR_eff ** 2)))
            n0 = max(2, math.ceil((Z ** 2) * (2 * sigma ** 2) / (logR_eff ** 2)))
            print(f"    {tau:>24.2f} {n_units:>10d} {n0:>12d} {n_units/max(n0,1):>10.2f}x")
            rows_hier.append(dict(A0=A0, tau=tau, sigma=sigma, P_purifications=n_units,
                                  P_if_tau_zero=n0, inflation=round(n_units / max(n0, 1), 3)))
        print()
    print()
    print("  At A0 = 0.70 the effect is so large that prep-level noise barely matters (P = 2")
    print("  until tau reaches 0.20).  At A0 = 0.90 it is decisive: tau = 0.05 already costs")
    print("  1.4x, tau = 0.10 costs 2.9x and tau = 0.20 costs 8.4x.  Additional reactions per")
    print("  preparation do not help at all once tau > 0 - only more purifications reduce error.")
    print()
    print("  Consequence for preregistration: report BOTH sigma (within-prep) and tau")
    print("  (between-prep) from the pilot.  Sizing on sigma alone — which is what WP13 and")
    print("  WP14 both did — understates the requirement by up to ~8x at A0 = 0.90, which is")
    print("  exactly the operating point where the study is otherwise believed to be feasible.")
    print()

    print("=" * 96)
    print("WP15.4  Revised stopping rules")
    print("=" * 96)
    stops = [
        ("w CI includes values < 0.8", "STOP - report preparation characterisation incomplete"),
        ("tau cannot be estimated from a pilot", "STOP - study cannot be sized"),
        ("trap readouts disagree between methods", "STOP for channel B; report A only"),
        ("any trap readout outside +/-10%", "STOP for channel B; report A only"),
        ("A_ctrl >= 0.90 at the pilot", "STOP - no power; do not run the comparison"),
        ("all checks pass", "PROCEED; apply the WP14 decision tree"),
    ]
    for cond, action in stops:
        print(f"  if {cond}\n     -> {action}")
    print()
    print("  Every branch above terminates in a reportable outcome, including STOP.  A study")
    print("  that stops for a stated reason is a result; one that proceeds past a failed check")
    print("  is a hazard.")

    for name, recs in (("wp15_bounds.csv", rows_bounds), ("wp15_active_fraction.csv", rows_act),
                       ("wp15_trap_readouts.csv", rows_trap), ("wp15_hierarchical_power.csv", rows_hier)):
        with open(OUT / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]))
            w.writeheader(); w.writerows(recs)
        print(f"wrote {OUT/name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
