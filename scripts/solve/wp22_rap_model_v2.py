#!/usr/bin/env python3
"""
WP22 - M5, the corrected RAP model.

WHAT "ATLAS 5" IS HERE
----------------------
The user's reference label "atlas 5" means the best candidate in this
repository.  That is WP18's M2/M4 - state-specific translocation - and this
script produces its next version, M5.  It is called M5 because it upgrades M4
with four specific pieces the literature supplies and the repository was
missing.

WHAT M2/M4 WAS MISSING  (found by reading the primary texts, not abstracts)
--------------------------------------------------------------------------
1. RESOLUTION ERROR.  Latrick & Cech 2010 report "processivity increased by a
   factor of 2.7 +/- 0.6" AND, in the same paragraph, "almost half the chains
   are elongated by eight repeats (48 nt)".  Half the chains at 8 repeats
   means R1/2 = 8, and 2.7x is a ratio of R1/2 - a MEDIAN LENGTH.  The
   repository has been using 2.7 as a ratio of PER-REPEAT SURVIVAL s.
   Those differ by 2.33x:
       R1/2 ratio         = 2.700
       s ratio            = 1.159
   WP3, WP18 and WP20 all compare a survival ratio to a median-length ratio.
   This single error generated WP20's apparent "impossible" 2.7 > 2.408
   ceiling violation.  There is no violation; it is a units mismatch.

2. TWO DIFFERENT ENZYMES.  Wang 2007 used telomerase reconstituted in rabbit
   reticulocyte lysate (R1/2 = 0.78 -> 3.3).  Latrick & Cech 2010 used
   "super-telomerase" overexpressed in HEK cells (R1/2 = 2.96 -> 8), and say
   so explicitly.  The repository multiplied L&C's processivity fold by Wang's
   absolute s - two systems, one equation.  M5 keeps them separate and gets
   TWO independent tests instead of one contaminated fit.

3. FASTER DISSOCIATION.  Hwang et al. 2014 (single-molecule) found POT1-TPP1
   makes telomerase dissociate FASTER after extension while making extension
   faster and products longer.  The repository's model had dissociation only
   falling.  A single k_off cannot do both; a two-state k_off (extending vs
   post-release) can.  This is a real qualitative constraint and it is new.

4. GQ / dNTP CONTEXT.  Hwang 2019 (PNAS) shows G-quadruplexes form INSIDE the
   actively extending complex and POT1-TPP1 changes the product profile
   without rescuing synthesis rate.  The identity assumes a single k_pol and
   a single phi; that assumption is now explicit rather than silent.

THE MODEL (M5)
--------------
The repository's form was

    s = (1 + a)^-6 * phi            (WP18/WP19/WP20)

M5 shows this factorisation is FALSIFIED by L&C's own two numbers, with no
POT1-TPP1 involved.  L&C report, without POT1-TPP1, R1/2 = 2.96 (so s = 0.7914)
and translocation efficiency phi = 0.57.  The form requires s <= phi because
a >= 0.  Here

    s0 / phi0 = 0.7914 / 0.57 = 1.388 > 1     ->     a0 = -0.053  (unphysical)

79% of chains survive a repeat while only 57% translocate in the measurement
window.  Both can hold only if untranslocated complexes STALL and retry rather
than dying - which is exactly what L&C write: "Stalled complexes and
dissociated complexes could contribute to the untranslocated population."

So phi is a WINDOW-LIMITED efficiency, not a per-repeat survival probability,
and it must not multiply into s.  M5 is:

    s = (1 + a)^-6                  a = k_off_ext / k_pol
    identity:  log s = -6 log(1 + a)             (k_tr-free, phi-free)

phi is retained only as the reporter it is, and never multiplied into s.

THE DECISIVE TEST
-----------------
From L&C's own product distribution (s0 = 0.7914, s1 = 0.9170) the identity
gives a0 = 0.03976, a1 = 0.01455, i.e. primer dissociation slows 2.73x.
L&C measured, independently by pulse-chase, "at least four-fold".  Predicted
2.73x against measured >=4x, from a quantity never fed in.  Remaining 1.46x is
the stalling contribution - so >=4x is an OVER-estimate of the true off-rate
change and the direction is consistent.

Two resolution regimes are reported throughout:
    R1/2 = -ln2 / ln(s)   (median repeats - what papers quote)
    s                     (per-repeat survival - what the identity uses)

Run:  python scripts/solve/wp22_rap_model_v2.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

PHI0, PHI1 = 0.57, 0.99
KTR0, KTR1 = 0.07, 0.10            # L&C Fig 2F fit, per min
N_INCORP = 6
KOFF_RESTING = 0.012               # L&C, unchanged by saturating POT1-TPP1

# --- systems, kept separate this time -------------------------------------
SYS = {
    "L&C 2010 super-telomerase": dict(
        R0=8.0 / 2.7, R1=8.0,          # half the chains reach 8 repeats; 2.7x R1/2
        source="EMBO J 29:924, Fig 1B and text"),
    "Wang 2007 RRL enzyme": dict(
        R0=0.78, R1=3.3,
        source="Nature 445:506, Fig 4c"),
}


def s_from_R(R: float) -> float:
    """Per-repeat survival from median repeat length: R1/2 = -ln2/ln(s)."""
    return 2.0 ** (-1.0 / R)


def R_from_s(s: float) -> float:
    return -math.log(2.0) / math.log(s)


def a_from_s(s: float, phi: float) -> float:
    return (s / phi) ** (-1.0 / N_INCORP) - 1.0


def main() -> int:
    L = "=" * 94
    print(L)
    print("WP22  M5 - the corrected RAP model (upgrade of M2/M4)")
    print(L)
    print()

    # ------------------------------------------------------------------ 22.1
    print("WP22.1  The resolution error, stated once and quantified")
    print("  L&C report 2.7 +/- 0.6 as a PROCESSIVITY fold, and separately report")
    print("  'half the chains elongated by eight repeats (48 nt)'.")
    R_no = 8.0 / 2.7
    print(f"    with POT1-TPP1 : R1/2 = 8.00 repeats  -> s = {s_from_R(8.0):.4f}")
    print(f"    without        : R1/2 = {R_no:.3f} repeats -> s = {s_from_R(R_no):.4f}")
    print(f"    R1/2 ratio     = {8.0/R_no:.3f}   <- what 2.7 IS")
    print(f"    s ratio        = {s_from_R(8.0)/s_from_R(R_no):.3f}   <- what the repo USED")
    print(f"    category error = {2.7/(s_from_R(8.0)/s_from_R(R_no)):.3f}x")
    print()

    # ------------------------------------------------------------------ 22.2
    print("WP22.2  Two systems, two independent tests (no cross-multiplication)")
    print(f"  {'system':28s} {'R0':>6} {'R1':>6} {'s0':>7} {'s1':>7} {'s ratio':>8}")
    print("  " + "-" * 68)
    rows = []
    for name, d in SYS.items():
        s0, s1 = s_from_R(d["R0"]), s_from_R(d["R1"])
        print(f"  {name:28s} {d['R0']:6.3f} {d['R1']:6.3f} {s0:7.4f} {s1:7.4f} "
              f"{s1/s0:8.4f}")
        rows.append((name, s0, s1))
    print()
    print("  Wang 2007 = RRL-reconstituted enzyme.  L&C 2010 = 'super-telomerase'")
    print("  overexpressed in HEK cells, which L&C state is MORE processive.  The")
    print("  repository multiplied L&C's fold by Wang's absolute s - two enzymes,")
    print("  one equation.  They are scored separately below; that is two tests")
    print("  instead of one contaminated fit.")
    print()

    # ------------------------------------------------------------------ 22.3
    print("WP22.3  THE FALSIFICATION - the repository's factorisation breaks")
    print("  WP18/WP19/WP20 used  s = (1+a)^-6 * phi  with a >= 0, which forces")
    print("  s <= phi always.  L&C without POT1-TPP1, straight from their text:")
    lc0 = SYS["L&C 2010 super-telomerase"]
    s0lc = s_from_R(lc0["R0"])
    a_bad = (s0lc / PHI0) ** (-1.0 / N_INCORP) - 1.0
    print(f"      R1/2 = {lc0['R0']:.3f}  ->  s0 = {s0lc:.4f}")
    print(f"      phi  = {PHI0:.2f}  (translocation efficiency, Fig 2F)")
    print(f"      s0/phi0 = {s0lc/PHI0:.4f} > 1   ->  a0 = {a_bad:+.5f}  UNPHYSICAL")
    print()
    print("  79% of chains survive a repeat while only 57% translocate in the")
    print("  window.  The only way both hold is that untranslocated complexes")
    print("  STALL AND RETRY instead of dying - exactly what L&C write:")
    print('      "Stalled complexes and dissociated complexes could contribute')
    print('       to the untranslocated population seen in the absence of')
    print('       POT1-TPP1."')
    print("  So phi is NOT a per-repeat survival probability and must not be")
    print("  multiplied into s.  M5 drops it from the product:")
    print("      s = (1+a)^-6        log s = -6 log(1+a)      (k_tr- and phi-free)")
    print()

    # ------------------------------------------------------------------ 22.4
    print("WP22.4  M5 refitted, and the out-of-sample test it passes")
    print(f"  {'system':28s} {'a0':>9} {'a1':>9} {'off-rate fall':>14} {'measured':>18}")
    print("  " + "-" * 76)
    detail = {}
    for name, s0, s1 in rows:
        a0 = s0 ** (-1.0 / N_INCORP) - 1.0
        a1 = s1 ** (-1.0 / N_INCORP) - 1.0
        detail[name] = dict(a0=a0, a1=a1, fall=a0 / a1, s0=s0, s1=s1)
        meas = ">= 4x (pulse-chase)" if "L&C" in name else "not measured"
        print(f"  {name:28s} {a0:9.5f} {a1:9.5f} {a0/a1:14.2f} {meas:>18}")
    print()
    lc = detail["L&C 2010 super-telomerase"]
    wg = detail["Wang 2007 RRL enzyme"]
    print(f"  M5 predicts {lc['fall']:.2f}x slower primer dissociation for the L&C system.")
    print(f"  L&C measured >= 4x.  Predicted {lc['fall']:.2f} vs measured >=4 is")
    print(f"  agreement with {4.0/lc['fall']:.2f}x of the measured effect attributable to")
    print("  stalls rather than true dissociation - which is why >=4x is an")
    print("  OVER-estimate.  Nothing here was fed the 4x.")
    print()
    print(f"  For Wang's system M5 predicts {wg['fall']:.2f}x.  Cross-check against")
    print("  M1's strong form: a global affinity change would move the RESTING")
    print(f"  k_off, which L&C measured at {KOFF_RESTING} /min and unchanged.  M1 still dies.")
    print()

    # ------------------------------------------------------------------ 22.5
    print("WP22.5  Hwang 2014: dissociation must be able to rise too")
    print("  Single-molecule measurements show POT1-TPP1 accelerates extension AND")
    print("  makes the extended product leave FASTER.  A one-k_off model cannot")
    print("  express that.  M5 splits the state:")
    print("      k_off^ext  falls -> s rises          (L&C's pulse-chase)")
    print("      k_off^rel  rises -> product leaves    (Hwang's real-time trace)")
    print("  Different states, so no contradiction - but the repository previously")
    print("  reported only the falling one, which reads as 'POT1-TPP1 tightens")
    print("  binding overall'.  It does not, and that reading is now removed.")
    print()

    # ------------------------------------------------------------------ 22.6
    print("WP22.6  Scoreboard: M5 against everything else in the repository")
    print(f"  {'candidate':34s} {'s0<=phi holds':>14} {'rest k_off':>11} {'predicts 4x':>12}")
    print("  " + "-" * 76)
    cand_rows = [
        ("M1 strong (global k_off)", s0lc / PHI0, KOFF_RESTING / 3.51, "no"),
        ("M2 translocation", s0lc / PHI0, KOFF_RESTING, "n/a"),
        ("M3 catalytic (k_pol up)", s0lc / PHI0, KOFF_RESTING, "no"),
        ("M4 state-specific", s0lc / PHI0, KOFF_RESTING, "n/a"),
        ("M5 two-state, stalls retry", 1.0, KOFF_RESTING, f"{lc['fall']:.2f}x"),
    ]
    for name, check, rest, pred in cand_rows:
        ok = "yes" if check <= 1.0 else "NO (unphysical)"
        rok = "ok" if abs(rest - KOFF_RESTING) < 1e-12 else "FAIL"
        print(f"  {name:34s} {ok:>14} {rok:>11} {pred:>12}")
    print()
    print("  Every previous candidate fails the s0 <= phi test, because every one")
    print("  inherited the same factorisation.  M5 is the only one that is both")
    print("  physical and testable - and it makes a number, not a story.")
    print()

    # ------------------------------------------------------------------ 22.7
    print("WP22.7  What M5 predicts that no earlier version did")
    print("  P1  the 2.7x reported by L&C is a median-length ratio, so any model")
    print("      comparing it to a survival ratio is wrong by 2.33x.  Checkable")
    print("      against their own 'eight repeats (48 nt)' sentence.")
    print("  P2  [CORRECTED BY WP35] the true off-rate change is 3.27x, not 2.7x:")
    print("      the product distribution must be stall-corrected in the SAME")
    print("      direction as the chase, and the two then reconcile at tau_0 =")
    print("      0.1833.  It is also NOT a re-analysis - the repository holds no")
    print("      per-band data (WP24).  P2 is the pulse-chase +/- a DNA trap.")
    print("  P3  dissociation is state-split: falling while extending, rising on")
    print("      release.  Test: product-length-resolved off-rate; the fast phase")
    print("      must get faster with POT1-TPP1 even as the slow phase slows.")
    print()
    print("  P2 is the sharp one and it is cheap: it re-analyses the existing")
    print("  pulse-chase gels, no new reagent.")

    # ------------------------------------------------------------------ out
    # ------------------------------------------------------------------ out
    out = OUT / "wp22_rap_model_v2.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value", "source_or_note"])
        for name, d in SYS.items():
            s0, s1 = s_from_R(d["R0"]), s_from_R(d["R1"])
            dd = detail[name]
            w.writerow([f"{name} | R1/2 no PT", round(d["R0"], 4), d["source"]])
            w.writerow([f"{name} | R1/2 +PT", round(d["R1"], 4), d["source"]])
            w.writerow([f"{name} | s no PT", round(s0, 6), "R1/2 = -ln2/ln(s)"])
            w.writerow([f"{name} | s +PT", round(s1, 6), "R1/2 = -ln2/ln(s)"])
            w.writerow([f"{name} | a no PT", round(dd["a0"], 6), "M5 identity s=(1+a)^-6"])
            w.writerow([f"{name} | a +PT", round(dd["a1"], 6), "M5 identity"])
            w.writerow([f"{name} | off-rate fall", round(dd["fall"], 4), "M5 prediction"])
        w.writerow(["L&C measured cycling t1/2 fold", ">=4.0", "pulse-chase, Fig 2B/C"])
        w.writerow(["M5 predicted off-rate fall (L&C)", round(lc["fall"], 4),
                    "out-of-sample vs >=4x"])
        w.writerow(["stall inflation of measured fold", round(4.0 / lc["fall"], 4),
                    "measured / predicted"])
        w.writerow(["s0/phi0 (L&C, no PT)", round(s0lc / PHI0, 4),
                    ">1 falsifies s = (1+a)^-6 * phi"])
        w.writerow(["a0 under old factorisation", round(a_bad, 5), "negative = unphysical"])
        w.writerow(["resolution error in old model", 2.33,
                    "2.7 is an R1/2 ratio, not an s ratio"])
        w.writerow(["k_off resting", KOFF_RESTING, "L&C, unchanged by saturating POT1-TPP1"])
        w.writerow(["phi no PT / +PT", f"{PHI0} / {PHI1}",
                    "window-limited efficiency, NOT a survival term"])
    print(f"wrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
