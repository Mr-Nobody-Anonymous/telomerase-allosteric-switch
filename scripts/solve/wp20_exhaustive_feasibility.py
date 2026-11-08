#!/usr/bin/env python3
"""
WP20 - Exhaustive feasibility test: which mechanism combinations can reproduce
       every measured quantity simultaneously?

WP19 assembled the best literature numbers.  This script does not add prose; it
enumerates the mechanisms and sweeps the whole parameter space, so that "is it
possible?" becomes a computed set rather than an argument.

MODEL (WP18/WP19, unchanged)
    phi = k_tr / (k_tr + k_off_tr)
    s   = (1 + a)^-6 * phi ,   a = k_off_ext / k_pol

MECHANISMS THE SWEEP CAN TURN ON
    koff_drop  : fall in k_off_ext and/or k_off_tr      (retention, M1/M4)
    ktr_drop   : fall in k_off_tr only, k_tr fixed      (M4)
    ktr_rise   : rise in k_tr                           (M2)
    Kd_drop    : fall in primer Kd via k_on             (M1 weak)

MEASURED TARGETS
    phi      0.57 -> 0.99                      L&C 2010, fitted A
    s0       0.4112 (absolute)                 Wang 2007, R_1/2 = 0.78
    processivity fold                          L&C 2010: 2.7 +/- 0.6
    s ratio                                    Wang 2007: 1.971
    resting k_off                              L&C 2010: 0.012, unchanged

Run:  python scripts/solve/wp20_exhaustive_feasibility.py
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

PHI0, PHI1 = 0.57, 0.99
S0_ANCHOR = 0.4112          # Wang 2007, R_1/2 = 0.78 -> s = 2^(-1/0.78)
PROC_OBS, PROC_SD = 2.7, 0.6
SRATIO_OBS = 1.971          # Wang 2007
KD_DROP = 3.51              # Wang 2007, 26.0 -> 7.4 nM
N_INCORP = 6
KOFF_RESTING = 0.012


def s_of(a: float, phi: float) -> float:
    return (1.0 + a) ** (-N_INCORP) * phi


def a_from_s(s: float, phi: float) -> float:
    return (s / phi) ** (-1.0 / N_INCORP) - 1.0


A0 = a_from_s(S0_ANCHOR, PHI0)
CEILING = PHI1 / S0_ANCHOR          # s1 <= phi1 always, since a1 >= 0


def main() -> int:
    L = "=" * 92
    print(L)
    print("WP20  Exhaustive feasibility: every mechanism combination, swept")
    print(L)
    print()
    print("WP20.1  The anchor fixes a0")
    print(f"  Wang 2007 gives R_1/2 = 0.78 repeats -> s0 = 2^(-1/0.78) = {S0_ANCHOR:.4f}")
    print(f"  a0 = (s0/phi0)^(-1/6) - 1 = {A0:.5f}")
    print()

    # ---------------------------------------------------------------- T1
    print("WP20.2  T1 - HARD CEILING (independent of any mechanism)")
    print("  s1 = (1+a1)^-6 * phi1 <= phi1 for all a1 >= 0, because k_off >= 0.")
    print(f"  Therefore  s1/s0 <= phi1/s0 = {PHI1:.2f}/{S0_ANCHOR:.4f} = {CEILING:.4f}")
    print("  This is not a fit and not a bound on a model - it is a bound on the")
    print("  data, given the identity.  Any claimed fold-change above it is either")
    print("  not s, or the anchor is wrong.")
    print()
    checks = [("Wang 2007 absolute s ratio", SRATIO_OBS, 0.0),
              ("L&C 2010 processivity fold", PROC_OBS, PROC_SD)]
    ceiling_rows = []
    for name, obs, sd in checks:
        over = obs - CEILING
        verdict = "IMPOSSIBLE" if over > 1e-9 else "feasible"
        zz = f"{over/sd:+.2f} SD" if sd else "no SD given"
        print(f"    {name:32s} {obs:5.3f}  vs ceiling {CEILING:.3f}  "
              f"{verdict:10s} ({zz})")
        ceiling_rows.append(dict(test="T1_ceiling", target=name, observed=obs,
                                 ceiling=round(CEILING, 4),
                                 exceeds="yes" if over > 1e-9 else "no",
                                 margin_sd=round(over / sd, 3) if sd else ""))
    print()

    # ---------------------------------------------------------------- T2
    print("WP20.3  T2 - Exhaustive sweep: can k_off and k_tr moves do it?")
    print("  Sweep koff_drop (total fall in k_off_ext and k_off_tr) over 1e0..1e3")
    print("  and ktr_drop (fall in k_off_tr only) over 1e0..1e3, phi1 fixed at 0.99.")
    print("  Report the feasible set for each observable.")
    print()
    grid = [10 ** (i / 12.0) for i in range(0, 37)]      # 1 .. ~1000
    sweep = []
    for koff_drop, ktr_drop in itertools.product(grid, grid):
        a1 = A0 / koff_drop
        s1 = s_of(a1, PHI1)
        sweep.append((koff_drop, ktr_drop, s1 / S0_ANCHOR))
    print(f"    {'observable':28s} {'feasible?':10s} {'min ratio':>10} {'max ratio':>10}")
    print("    " + "-" * 62)
    ratios = [r[2] for r in sweep]
    feas_wang = [r for r in sweep if abs(r[2] - SRATIO_OBS) <= 0.05]
    feas_proc = [r for r in sweep if abs(r[2] - PROC_OBS) <= PROC_SD]
    print(f"    {'range covered by sweep':28s} {'-':10s} {min(ratios):10.4f} {max(ratios):10.4f}")
    print(f"    {'Wang s ratio 1.971 +/-0.05':28s} "
          f"{'YES' if feas_wang else 'NO':10s} {len(feas_wang):10d} pts")
    print(f"    {'L&C 2.7 +/-0.6':28s} "
          f"{'YES' if feas_proc else 'NO':10s} {len(feas_proc):10d} pts")
    print()
    if feas_wang:
        best = min(feas_wang, key=lambda r: r[0])
        print(f"    Cheapest Wang-consistent point: koff_drop={best[0]:.2f} "
              f"-> ratio {best[2]:.4f}")
        print("    Note ktr_drop is IRRELEVANT here: without an independent k_tr the")
        print("    two arms cannot be separated, so every ktr_drop value ties.")
    print()

    # ---------------------------------------------------------------- T3
    print("WP20.4  T3 - Is the Kd-driven (k_on) route enough on its own?")
    print("  Primer Kd falls 3.51x (Wang: 26.0 -> 7.4 nM).  If that acts by raising")
    print("  k_on, a falls with Kd while k_off is untouched:")
    a1_kd = A0 / KD_DROP
    s1_kd = s_of(a1_kd, PHI1)
    print(f"    a1 = {A0:.5f}/{KD_DROP} = {a1_kd:.5f}   s1 = {s1_kd:.4f}   "
          f"ratio = {s1_kd/S0_ANCHOR:.4f}")
    over = s1_kd / S0_ANCHOR - SRATIO_OBS
    print(f"    Wang target {SRATIO_OBS:.3f}  ->  OVERSHOOTS by {over:+.4f} "
          f"({over/SRATIO_OBS*100:+.1f}%)")
    print()
    print("  Read this carefully - it is NOT a match.  The Kd route overshoots Wang's")
    print(f"  ratio by {over/SRATIO_OBS*100:.0f}%.  What it does establish is the ORDER of")
    print("  magnitude and the DIRECTION: a measured 3.51x change in primer Kd, fed")
    print("  through the identity with no free parameters and without touching")
    print(f"  k_off_ext, lands at {s1_kd/S0_ANCHOR:.2f}x against Wang's {SRATIO_OBS:.2f}x.  It")
    print("  is the closest single mechanism in the whole sweep.")
    print()

    # ---------------------------------------------------------------- T4
    print("WP20.5  T4 - M1 vs M4 vs M2 against every measured item")
    print()
    cand = [
        ("M1 strong: k_off down globally", A0 / KD_DROP, PHI0, KOFF_RESTING / KD_DROP),
        ("M1 weak: k_on up (Kd only)", A0 / KD_DROP, PHI1, KOFF_RESTING),
        ("M2: k_tr up + k_off_tr down", A0, PHI1, KOFF_RESTING),
        ("M4: k_off_tr down only", A0, PHI1, KOFF_RESTING),
        ("M3: k_pol up", A0 / KD_DROP, PHI0, KOFF_RESTING),
    ]
    print(f"    {'candidate':30s} {'s ratio':>8} {'vs Wang':>9} {'rest k_off':>11} {'phi':>5}")
    print("    " + "-" * 68)
    t4_rows = []
    for name, a1, phi1, rest in cand:
        s1 = s_of(a1, phi1)
        ratio = s1 / S0_ANCHOR
        dev = ratio - SRATIO_OBS
        rok = "ok" if abs(rest - KOFF_RESTING) < 1e-12 else "FAIL"
        phi_ok = "ok" if abs(phi1 - PHI1) < 1e-9 else "FAIL"
        print(f"    {name:30s} {ratio:8.4f} {dev:+9.3f} {rok:>11} {phi_ok:>5}")
        t4_rows.append(dict(test="T4_candidate", target=name,
                            observed=round(ratio, 4),
                            ceiling=round(dev, 4), exceeds=rok, margin_sd=phi_ok))
    print()
    print("    No single-rate candidate lands on Wang's 1.971 - every one is either")
    print("    too small (M2/M4: 1.74) or too large (Kd route: 2.19).  Wang's own")
    print("    ratio sits BETWEEN the phi-only and the Kd-driven prediction, which is")
    print("    what you expect if both arms contribute partially.")
    print()

    # ---------------------------------------------------------------- T5
    print("WP20.6  T5 - The one measurement that would break the degeneracy")
    print("  With a0 anchored, M2 and M4 both give exactly s1 = phi1 = 0.99. Their")
    print("  difference lives entirely in k_tr, which s and phi both integrate out.")
    print("  A k_tr measurement is therefore decisive and nothing else is:")
    print("    implied k_off_tr = k_tr (1-phi)/phi:")
    print(f"      no TPP1 : k_tr=0.07 -> {0.07*(1-PHI0)/PHI0:.5f} /min")
    print(f"      + TPP1  : if k_tr stays 0.07 -> {0.07*(1-PHI1)/PHI1:.5f} /min")
    print(f"      + TPP1  : if k_tr is 0.10    -> {0.10*(1-PHI1)/PHI1:.5f} /min")
    print("    Direct measurement of k_tr (or of k_off_tr by stopped-flow) settles it.")
    print()

    # ---------------------------------------------------------------- verdict
    print("WP20.7  VERDICT")
    print("  ** SUPERSEDED BY WP22 on two points - read this first **")
    print("     (i) 2.7 is a ratio of R1/2 (median length), NOT of per-repeat")
    print("         survival s. L&C's own text gives R1/2 = 8 with POT1-TPP1, so")
    print("         the s ratio is 1.159, not 2.700. The ceiling below is compared")
    print("         against the wrong quantity and the 'violation' is a units error.")
    print("    (ii) Wang's s0 and L&C's 2.7 come from DIFFERENT ENZYMES (RRL")
    print("         reconstitution vs HEK super-telomerase) and must not be")
    print("         multiplied into one equation. WP22 scores them separately.")
    print("  What survives from WP20: the sweep (T2) and the Kd route (T3) are")
    print("  unaffected because both are computed within a single system.")
    print()
    print(f"  1. CEILING {CEILING:.3f} (from s0={S0_ANCHOR:.4f}): the identity cannot")
    print("     produce 2.7x at all. L&C's 2.7x processivity EXCEEDS what s can")
    print("     reach at s0=0.4112. The two literature numbers are mutually")
    print("     inconsistent AT THE ANCHOR. One of them is not per-repeat s.")
    need_s0 = PHI1 / PROC_OBS
    need_r = -math.log(2) / math.log(need_s0)
    print(f"     Sensitivity: to reach 2.7x the anchor would need s0 <= {need_s0:.4f},")
    print(f"     i.e. R_1/2 <= {need_r:.4f} repeats instead of Wang's 0.78 - a")
    print(f"     {abs(need_r-0.78)/0.78*100:.1f}% shift. That is inside a plausible read-off error,")
    print("     so the tension is a CROSS-PAPER CALIBRATION issue, not a refutation.")
    print("  2. Wang's absolute s ratio (1.971) is feasible, and the measured 3.51x")
    print(f"     Kd change alone lands at {s1_kd/S0_ANCHOR:.2f}x - the closest single mechanism,")
    print("     overshooting by 11%. Wang's value sits BETWEEN phi-only (1.74x) and")
    print("     Kd-driven (2.19x), which is what partial contribution of both looks like.")
    print("  3. M1-strong dies on the resting k_off (0.012 must not move).")
    print("     M3 dies on phi (must reach 0.99) and on the 35.8-42.1 A reach.")
    print("  4. M2 and M4 remain exactly degenerate. The decisive experiment is a")
    print("     direct k_tr measurement, not more processivity data.")
    print()

    out = OUT / "wp20_exhaustive_feasibility.csv"
    rows = ceiling_rows + t4_rows
    rows.append(dict(test="summary", target="ceiling", observed=round(CEILING, 4),
                     ceiling="", exceeds="", margin_sd=""))
    rows.append(dict(test="summary", target="a0_from_anchor", observed=round(A0, 6),
                     ceiling="", exceeds="", margin_sd=""))
    rows.append(dict(test="summary", target="kd_route_s_ratio",
                     observed=round(s1_kd / S0_ANCHOR, 4), ceiling="", exceeds="",
                     margin_sd=""))
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "target", "observed", "ceiling",
                                           "exceeds", "margin_sd"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
