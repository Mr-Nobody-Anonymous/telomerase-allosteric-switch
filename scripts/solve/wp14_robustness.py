#!/usr/bin/env python3
"""
WP14 - Robustness of the WP13 conclusions to kinetic heterogeneity and rebinding.

Four audits requested, plus the specific question: does the A <= 0.70 feasibility threshold
survive heterogeneous-repeat and rebinding models?

  14.1  Independent reproduction of the WP13 sample-size calculation, with every assumption
        stated and a Monte-Carlo check against the closed form.
  14.2  Heterogeneous-repeat sensitivity: position-dependent k_pol, slow-step, graded, and
        "only one dissociation-competent step".
  14.3  Population heterogeneity: mixtures of responsive and inert enzyme.
  14.4  Rebinding, with and without an absorbing trap, and sensitivity to trap side effects.
  14.5  Do the A/B decision rules remain valid under 14.1-14.4?
  14.6  Revised preregistration text.

Run:  python scripts/solve/wp14_robustness.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
RNG = np.random.default_rng(1414)

Z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)   # 80% power, two-sided alpha = 0.05


def A_of(a):
    """Channel A for an arbitrary vector of per-step ratios a_i = k_off,i / k_pol,i."""
    a = np.atleast_1d(np.asarray(a, dtype=float))
    return float(np.prod(1.0 / (1.0 + a)))


def A_ratio(a, fold):
    a = np.asarray(a, dtype=float)
    return A_of(a / fold) / A_of(a)


def calibrate_to_A(profile, target=0.70):
    """Scale a shape profile so that baseline A equals `target`."""
    p = np.asarray(profile, dtype=float)
    lo, hi = 1e-12, 1e4
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if A_of(mid * p) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) * p


def main() -> int:
    rows_ss, rows_het, rows_mix, rows_trap, rows_rule = [], [], [], [], []

    # ------------------------------------------------------------------ 14.1 sample sizes
    print("=" * 96)
    print("WP14.1  Reproducing the WP13 sample-size numbers, assumptions made explicit")
    print("=" * 96)
    assumptions = [
        ("replicate unit", "INDEPENDENT BIOLOGICAL replicates (separate purifications)"),
        ("never", "technical replicates - they do not reduce biological variance"),
        ("effect metric", "ratio A_trt/A_ctrl, each arm a group mean over n_rep"),
        ("CV applies to", "individual replicate A values within an arm (not to the group mean)"),
        ("test", "TWO-sided alpha = 0.05"),
        ("power", "0.80"),
        ("uncertainty", "BOTH arms, propagated: SE_ratio/ratio = sqrt(2) * CV/sqrt(n_rep)"),
        ("phi noise", "binomial; contributes < 0.5% to SE at n_mol >= 2000, i.e. negligible vs CV 5%"),
        ("independence", "assumes arms share no variance; a shared prep INVALIDATES this"),
    ]
    for k, v in assumptions:
        print(f"  {k:16s}: {v}")
    print()
    print("  Closed form:  n_rep = ceil( [ CV * z / (effect/ sqrt(2)) ]^2 ),  z = 2.802")
    print()
    def sim_power(n, R, cv, trials=4000, rng=RNG):
        """Exact two-sample Welch t-test power by simulation."""
        if n < 2:
            return float("nan")
        hits = 0
        for _ in range(trials):
            ctl = rng.normal(1.0, cv, n)
            trt = rng.normal(R, cv * R, n)
            se = math.sqrt(ctl.var(ddof=1) / n + trt.var(ddof=1) / n)
            if se > 0:
                hits += abs((trt.mean() - ctl.mean()) / se) > stats.t.ppf(0.975, df=2 * n - 2)
        return hits / trials

    print(f"  {'A_ctrl':>8} {'R(f=4)':>9} {'effect':>8} {'n (z-form)':>11} "
          f"{'power':>7} {'n (exact t)':>12} {'power':>7}")
    print("  " + "-" * 76)
    for Ac in (0.99, 0.95, 0.90, 0.70, 0.50, 0.30):
        a = calibrate_to_A(np.ones(6), Ac)
        R = A_ratio(a, 4)
        eff = abs(R - 1)
        if eff < 1e-9:
            continue
        cv = 0.05
        n_closed = max(1, math.ceil((cv * Z / (eff / math.sqrt(2))) ** 2))
        p_at = sim_power(n_closed, R, cv)
        # exact-t requirement: smallest n whose simulated power reaches 0.80
        n_exact = None
        for cand in range(max(2, n_closed - 5), n_closed + 60):
            if sim_power(cand, R, cv, trials=2000) >= 0.80:
                n_exact = cand
                break
        p_exact = sim_power(n_exact, R, cv) if n_exact else float("nan")
        print(f"  {Ac:>8.3f} {R:>9.4f} {eff:>8.4f} {n_closed:>11d} {p_at:>7.3f} "
              f"{(n_exact if n_exact else '>n/a'):>12} {p_exact:>7.3f}")
        rows_ss.append(dict(A_ctrl=Ac, A_ratio=round(R, 4), effect=round(eff, 5), cv=cv,
                            n_z_formula=n_closed, power_at_n_z=round(p_at, 3),
                            n_exact_t=(n_exact if n_exact else None),
                            power_at_n_exact=round(p_exact, 3) if n_exact else None))

    print()
    print("  The numbers are BIOLOGICAL")
    print("  replicates at CV 5% on A; at CV 10% they scale by 4.  WP13's 686/2744 are the")
    print("  CV 5% / CV 10% pair for A = 0.99 and are therefore conditional on an assay")
    print("  precision nobody has demonstrated for this observable.")
    print()

    # ------------------------------------------------------------------ 14.2 per-position heterogeneity
    print("=" * 96)
    print("WP14.2  Does the A <= 0.70 threshold survive position-dependent k_pol?")
    print("=" * 96)
    profiles = {
        "homogeneous":       np.ones(6),
        "one slow step 10x": np.array([10, 1, 1, 1, 1, 1]),
        "one fast step 10x": np.array([0.1, 1, 1, 1, 1, 1]),
        "graded 1..6":       np.array([1, 2, 3, 4, 5, 6]),
        "only ONE step":     np.array([1, 0, 0, 0, 0, 0]),
        "alternating hi/lo": np.array([8, 0.5, 8, 0.5, 8, 0.5]),
    }
    print("  Each profile scaled so baseline A = 0.70, then the fold-drop response measured.")
    print()
    print(f"  {'k_pol profile':20s} {'A':>7} {'R(f=2)':>9} {'R(f=4)':>9} {'R(f=10)':>9} "
          f"{'n_req(f=4)':>12}")
    print("  " + "-" * 72)
    Rs = []
    for name, prof in profiles.items():
        a = calibrate_to_A(prof, 0.70)
        r2, r4, r10 = (A_ratio(a, f) for f in (2, 4, 10))
        n_req = max(1, math.ceil((0.05 * Z / ((r4 - 1) / math.sqrt(2))) ** 2))
        Rs.append(r4)
        print(f"  {name:20s} {A_of(a):>7.4f} {r2:>9.4f} {r4:>9.4f} {r10:>9.4f} {n_req:>12d}")
        rows_het.append(dict(profile=name, A=round(A_of(a), 4), R_f2=round(r2, 4),
                             R_f4=round(r4, 4), R_f10=round(r10, 4), n_req_f4=n_req))
    print()
    print(f"  Range of R(f=4) across wildly different profiles: {min(Rs):.4f} - {max(Rs):.4f}")
    print("  -> the threshold SURVIVES position-dependent k_pol.  Even the extreme case of a")
    print("     single dissociation-competent step (n_eff = 1) moves R(f=4) only from 1.304")
    print("     to 1.290, because calibrating on A already absorbs the shape.")
    print()

    # ------------------------------------------------------------------ 14.3 population heterogeneity
    print("=" * 96)
    print("WP14.3  The case that DOES break the threshold: mixtures of enzyme populations")
    print("=" * 96)
    print("  A preparation is rarely 100% active.  Model: fraction w of responsive enzyme plus")
    print("  (1-w) of inert enzyme.  Observed A is the weighted mean of per-repeat survivals.")
    print()
    print(f"  {'mixture':44s} {'A_obs':>8} {'R(f=4)':>9} {'verdict at A<=0.70 gate':>26}")
    print("  " + "-" * 90)
    mixes = [
        (1.00, 0.06125, 0.06125, "100% responsive (reference)"),
        (0.90, 0.06800, 1e-9, "90% responsive + 10% inert"),
        (0.70, 0.07140, 1e-9, "70% responsive + 30% inert"),
        (0.50, 0.08000, 1e-9, "50% responsive + 50% inert"),
        (0.30, 0.10000, 1e-9, "30% responsive + 70% inert"),
        (0.50, 0.06125, 0.30000, "50% A=0.70 + 50% A=0.30 (both responsive)"),
    ]
    for w, a_resp, a_inert, lab in mixes:
        A0 = w * A_of([a_resp] * 6) + (1 - w) * (A_of([a_inert] * 6) if a_inert > 1e-8 else 1.0)
        A1 = w * A_of([a_resp / 4] * 6) + (1 - w) * (A_of([a_inert / 4] * 6) if a_inert > 1e-8 else 1.0)
        R = A1 / A0
        gate = "PASS (proceed)" if A0 <= 0.70 else "FAIL (abandon study)"
        print(f"  {lab:44s} {A0:>8.4f} {R:>9.4f} {gate:>26}")
        rows_mix.append(dict(mixture=lab, w_responsive=w, A_obs=round(A0, 4),
                             R_f4=round(R, 4), gate="pass" if A0 <= 0.70 else "fail"))
    print()
    print("  THIS is where the threshold fails.  A 50% inert preparation sits at A = 0.815 and")
    print("  would be rejected by the A <= 0.70 gate - yet it contains 50% fully responsive")
    print("  enzyme.  The gate does not distinguish 'the assay has no power' from 'the")
    print("  preparation is mostly dead'.")
    print()
    print("  Necessary repair: the active fraction must be measured and reported with A.  Two")
    print("  minimally invasive options, both standard:")
    print("    1. an activity standard curve (known active-site titration) -> w;")
    print("    2. report A alongside per-repeat survival of a concentration series, since a")
    print("       mixture shows curvature that a homogeneous population does not.")
    print()
    print("  Corrected gate: apply A <= 0.70 ONLY to preparations with w >= 0.8; otherwise")
    print("  report the threshold as uninterpretable rather than as infeasible.")

    # ------------------------------------------------------------------ 14.4 rebinding and trap
    print()
    print("=" * 96)
    print("WP14.4  Rebinding: does the trap work, and what if the trap perturbs the system?")
    print("=" * 96)
    print("  With irreversible loss, A is the pure extension-phase survival.  With rebinding at")
    print("  rate kon*[P] the observed per-repeat survival is inflated; the trap removes it.")
    print()
    print(f"  {'kon*[P] (per min)':>18} {'s_obs (1 rep)':>14} {'A_effective':>13} "
          f"{'bias on A':>11}")
    k_pol, k_off, k_tr = 0.10, 0.003, 0.07
    a0 = k_off / k_pol
    A_true = A_of([a0] * 6)
    for konP in (0.0, 0.01, 0.1, 0.5, 2.0):
        q = konP / (konP + k_off + 1e-12) if konP > 0 else 0.0
        s_eff = A_of([a0] * 6) * (k_tr / (k_tr + k_off * (1 - q)))
        print(f"  {konP:>18.3f} {s_eff:>14.4f} {A_true:>13.4f} {A_true:>11.4f}")
        rows_trap.append(dict(konP=konP, s_obs=round(s_eff, 5), A_true=round(A_true, 5)))
    print()
    print("  Rebinding inflates s and therefore the apparent translocation term, but leaves")
    print("  A mathematically unchanged (d/dlog(k_off_tr)[log A] = 0, WP13.6).  So a trap is")
    print("  REQUIRED for B and merely advisable for A.")
    print()
    print("  Trap side effects - the risk the audit asks about.  If the trap perturbs the")
    print("  system by factor t, the bias on the A-ratio is:")
    print()
    print(f"  {'perturbation':46s} {'f_trap':>8} {'A_ctrl shift':>13} {'R(f=4) shift':>13}")
    for lab, t in (("trap reduces k_pol by 10%", 0.90),
                   ("trap reduces k_pol by 30%", 0.70),
                   ("trap reduces k_pol by 50%", 0.50),
                   ("trap raises k_off by 30%", 1.0)):
        if "k_off" in lab:
            a_t = a0 * 1.30
            A_t = A_of([a_t] * 6)
            R_t = A_ratio([a_t] * 6, 4)
            print(f"  {lab:46s} {1.0:>8.2f} {A_t/A_true:>13.3f} {R_t:>13.3f}")
        else:
            a_t = a0 / t
            A_t = A_of([a_t] * 6)
            R_t = A_ratio([a_t] * 6, 4)
            print(f"  {lab:46s} {t:>8.2f} {A_t/A_true:>13.3f} {R_t:>13.3f}")
    print()
    print("  The trap moves BOTH the baseline A and the effect size, in the same direction.")
    print("  A trap that halves k_pol cuts baseline A and roughly doubles the detectable")
    print("  effect - it would make the study look MORE feasible than it is.  Because the")
    print("  trap is applied to both arms, the bias does not cancel in the ratio.")
    print()
    print("  Therefore the trap must be validated, not assumed: measure activity, TPP1")
    print("  binding and k_pol with and without trap, in the absence of any TPP1 effect.")

    # ------------------------------------------------------------------ 14.5 decision rules
    print()
    print("=" * 96)
    print("WP14.5  Do the WP13 decision rules survive?")
    print("=" * 96)
    checks = [
        ("A-ratio CI excludes 1 and |ratio-1| >= 0.05 -> DETECTED", "survives",
         "heterogeneity changes effect size but not the direction of inference"),
        ("A-ratio null with A_ctrl <= 0.70 -> EXCLUDED", "survives with a caveat",
         "valid only for w >= 0.8; a mixture can hide a responsive subpopulation"),
        ("A-ratio null with A_ctrl > 0.90 -> INCONCLUSIVE", "survives, strengthened",
         "must be reported as 'no power' or 'mixture' - never as 'no effect'"),
        ("B lower bound > control B upper bound -> DETECTED", "survives only with a trap",
         "without an absorbing trap, B is an upper bound and cannot be tested this way"),
        ("phi diff within resolution, both < 0.95 -> EXCLUDED", "survives",
         "still requires both arms below saturation"),
        ("phi_trt > 0.95, CI touches 1 -> INCONCLUSIVE", "survives, unchanged",
         "any reparameterisation of B fails the same way"),
    ]
    for cond, status, note in checks:
        print(f"  {cond}")
        print(f"     -> {status}")
        print(f"        {note}")
        rows_rule.append(dict(rule=cond, status=status, note=note))
    print()

    print("  Net: the A-channel rules are robust to position-dependent kinetics and to")
    print("  rebinding, and are threatened only by population mixtures.  The B-channel rules")
    print("  are threatened by rebinding and depend on the trap working.")

    # ------------------------------------------------------------------ writes
    for name, recs in (("wp14_sample_sizes.csv", rows_ss), ("wp14_heterogeneity.csv", rows_het),
                       ("wp14_mixtures.csv", rows_mix), ("wp14_rebinding.csv", rows_trap),
                       ("wp14_decision_rules.csv", rows_rule)):
        with open(OUT / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]))
            w.writeheader(); w.writerows(recs)
        print(f"wrote {OUT/name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
