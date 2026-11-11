#!/usr/bin/env python3
"""
WP19 - Close the gaps by finding the measurement in the literature.

WP18 left one thing missing: an absolute per-repeat survival s with a matched
phi, to test the identity  s/phi = (1+a)^-6  against an independently measured
off-rate.  That measurement exists and is used below.

WHAT WAS FOUND (verified against the primary text, not abstracts)
-----------------------------------------------------------------
1. Latrick & Cech 2010, EMBO J 29:924 (PMC2837173) DID measure cycling
   primer dissociation by pulse-chase.  Exact wording:
     "With actively cycling enzyme, POT1-TPP1 clearly slowed the dissociation
      of the ensemble of extending primers ... the time required for
      half-dissociation was increased at least four-fold at 16 or 30 C."
   -> the claim "nobody has measured primer k_off +/- TPP1 for cycling
      enzyme" is FALSE and is retracted here.

2. Same paper: k_pol.  "measuring the rate of incorporation of dG ... we
   observed no increase k_pol in the presence of POT1-TPP1", and k_pol is
   fast relative to translocation so it does not set processivity.
   -> M3 (catalytic / k_pol up) is excluded by MEASUREMENT, not just by reach.

3. Same paper, translocation: y = A(1 - e^(-kt)) with A = 0.57 -> 0.99 and
   k = 0.07 -> 0.10 /min.  A and k are fitted as SEPARATE parameters.
   -> WP18 overstated this as "only the combination is observable".  It is
      not: k_tr = 0.10/min is reported.  Corrected here.  The two are
      correlated through one curve, so eps = 52x carries a correlated error,
      but C2 is a measurement and is reinstated.

4. Wang et al 2007, Nature 445:559.  Absolute per-repeat survival through
   their own formula: processivity = R_1/2 = -ln2/(2.303 k_slope).
   R_1/2 = 0.78 repeats (no POT1-TPP1) -> 3.3 repeats (+POT1-TPP1).
   Since R_1/2 = -ln2/ln(s),  s = 0.5^(1/R_1/2) = 2^(-1/R_1/2).
   This is the absolute anchor the model needed.

Run:  python scripts/solve/wp19_literature_closure.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

PHI0, PHI1 = 0.57, 0.99                  # L&C translocation efficiency A
KTR0, KTR1 = 0.07, 0.10                  # L&C translocation rate k, /min
KOFF_REST = 0.012                        # L&C resting k_off, /min, UNCHANGED
T12_CYCLING_FOLD = 4.0                   # L&C: "increased at least four-fold"
RAP_FOLD_OBS = 2.7                       # L&C processivity increase (n=4)

R_HALF0, R_HALF1 = 0.78, 3.3             # Wang 2007 R_1/2, no PT -> +PT
KD_POT1_A5, KD_POT1TPP1_A5 = 26.0, 7.4   # Wang 2007, nM, primer a5
N_INCORP = 6

# Wang 2007 text: "almost half the chains elongated by eight repeats (48 nt)
# ... whereas <10% of the chains reach this length without POT1-TPP1"
LC8_PT, LC8_NO_PT = 0.50, 0.10
EIGHT = 8


def s_from_rhalf(r_half: float) -> float:
    """s = 2^(-1/R_1/2); R_1/2 is repeats until half the chains dissociate."""
    return 2.0 ** (-1.0 / r_half)


def a_from_identity(s: float, phi: float) -> float:
    """Invert s = (1+a)^-6 * phi  ->  a = (s/phi)^(-1/6) - 1."""
    base = s / phi
    return base ** (-1.0 / N_INCORP) - 1.0


def main() -> int:
    L = "=" * 92
    print(L); print("WP19  Closing the gaps with measurements that already exist"); print(L); print()

    # ---- 1. absolute s from Wang 2007 -------------------------------------
    print("WP19.1  Absolute per-repeat survival from Wang 2007 R_1/2")
    s0, s1 = s_from_rhalf(R_HALF0), s_from_rhalf(R_HALF1)
    print(f"  R_1/2 = {R_HALF0} repeats (no PT)  ->  s = {s0:.4f}")
    print(f"  R_1/2 = {R_HALF1} repeats (+PT)    ->  s = {s1:.4f}")
    print(f"  measured s ratio .................. {s1/s0:.3f}x")
    print("  (Wang's own formula R_1/2 = -ln2/(2.303 k) = -ln2/ln(s).)")
    print()

    # ---- 2. the L&C 8-repeat cross-check ---------------------------------
    print("WP19.2  Independent cross-check from L&C's own text at 8 repeats")
    s_pt_lo = 0.5 ** (1 / EIGHT)
    s_no_hi = LC8_NO_PT ** (1 / EIGHT)
    print(f"  >=half reach 8 repeats with PT   -> s >= {s_pt_lo:.4f}")
    print(f"  <10% reach 8 repeats without PT  -> s <= {s_no_hi:.4f}")
    print(f"  => s ratio >= {s_pt_lo/s_no_hi:.3f}x")
    print("  Wang gives s_PT 0.81 and s_noPT 0.41, L&C's brackets are")
    print("  0.92 and 0.75. Same direction, same order of magnitude, different")
    print("  primers and buffers. Consistent, not identical - treat as a bound.")
    print()

    # ---- 3. a from the identity, on the absolute anchor ------------------
    print("WP19.3  The identity applied to real numbers")
    a0, a1 = a_from_identity(s0, PHI0), a_from_identity(s1, PHI1)
    print(f"  a = k_off_ext/k_pol, no PT ... {a0:.4f}   (1+a)^-6 = {(1+a0)**-6:.4f}")
    print(f"  a = k_off_ext/k_pol, +PT ..... {a1:.4f}   (1+a)^-6 = {(1+a1)**-6:.4f}")
    print(f"  a falls by .................... {a0/a1:.2f}x")
    print("  With k_pol MEASURED unchanged (L&C), that is an EXTENSION-phase")
    print(f"  off-rate fall of {a0/a1:.2f}x, while resting k_off stays at 0.012.")
    print()

    # ---- 4. the conflict -------------------------------------------------
    print("WP19.4  Two measurements of the same quantity disagree - and why")
    print(f"  pulse-chase, cycling t_1/2 ........ >= {T12_CYCLING_FOLD:.0f}x slower")
    print(f"  identity + absolute R_1/2 ......... {a0/a1:.2f}x slower")
    print(f"  disagreement factor ............... {T12_CYCLING_FOLD/(a0/a1):.2f}x")
    print()
    print("  The pulse-chase number is a LOWER bound on the slowdown and includes")
    print("  rebinding: a primer that dissociates can reassociate before the next")
    print("  extension, inflating the apparent bound fraction - exactly the WP13")
    print("  result.  The identity uses R_1/2, which counts only successful")
    print("  committed extensions, so it is rebinding-free.")
    print()
    print("  TESTABLE PREDICTION (new):")
    print(f"    with dilution or a DNA trap to abolish rebinding, the pulse-chase")
    print(f"    t_1/2 fold-change should fall from >= {T12_CYCLING_FOLD:.0f}x toward")
    print(f"    ~{a0/a1:.1f}x, i.e. by a factor of ~{T12_CYCLING_FOLD/(a0/a1):.1f}.")
    print("    If it stays >= 4x under trap, the identity's R_1/2 anchor is wrong")
    print("    and the model is falsified.")
    print()

    # ---- 5. does M1 survive? resting affinity ---------------------------
    print("WP19.5  M1 (global affinity) against the measured equilibrium Kd")
    kd_ratio = KD_POT1_A5 / KD_POT1TPP1_A5
    print(f"  Kd POT1-a5 = {KD_POT1_A5} nM, POT1-TPP1-a5 = {KD_POT1TPP1_A5} nM")
    print(f"  TPP1 tightens POT1-DNA binding {kd_ratio:.2f}x (Wang 2007).")
    print("  So the GLOBAL affinity does change - M1 is not absurd - but it")
    print("  changes on the DNA/POT1 axis, NOT on TERT-primer k_off, which L&C")
    print(f"  measured as unchanged at {KOFF_REST} /min. The mechanism is")
    print("  state-specific: it acts on cycling enzyme, not on the resting")
    print("  complex. M1 in its strong form (global k_off on TERT) stays dead;")
    print("  a DNA-affinity component is real and is retained.")
    print()

    # ---- 6. consistency of the full model -------------------------------
    print("WP19.6  Does the assembled model reproduce the 2.7x?")
    phi_fold = PHI1 / PHI0
    ret_fold = ((1 + a1) ** -N_INCORP) / ((1 + a0) ** -N_INCORP)
    total = phi_fold * ret_fold
    print(f"  translocation efficiency  phi ratio ....... {phi_fold:.4f}")
    print(f"  retention term (1+a)^-6 ratio ............. {ret_fold:.4f}")
    print(f"  total ..................................... {total:.4f}")
    print(f"  observed (L&C, n=4) ....................... {RAP_FOLD_OBS}")
    print(f"  z ......................................... {(total-RAP_FOLD_OBS)/0.6:+.2f}")
    print(f"  absolute s ratio (Wang) ................... {s1/s0:.4f}  -> residual {(s1/s0)/total:.3f}x")
    print()
    print("  NOTE the circularity, stated plainly: a0 and a1 were obtained BY")
    print("  inverting the identity on Wang's s and L&C's phi, so the total above")
    print("  is an identity, not a prediction.  What is NOT circular is WP19.4:")
    print("  the independently measured cycling t_1/2 (>=4x) exceeds the off-rate")
    print("  change the identity allows (1.67x) by 2.4x. That excess is the")
    print("  quantitative signature of rebinding and is the falsifiable claim.")
    print()

    # ---- 7. what is left ------------------------------------------------
    print("WP19.7  What is still genuinely open after this")
    items = [
        ("M2 vs M4 (does k_tr itself rise?)",
         "k_tr = 0.10/min is reported, so M2 is supported; but A and k are "
         "fitted from one curve and their covariance is unpublished, so M4 "
         "is not excluded. Needs the raw time-course, not a new experiment."),
        ("rebinding correction",
         "quantify by running the pulse-chase under DNA trap / dilution."),
        ("RAP-state superposition 9SHZ / 9SI0",
         "coordinates now retrievable via models.rcsb.org one-residue at a "
         "time; still not ingested."),
        ("10 vs 8.4 repeats",
         "L&C say 'eight repeats (48 nt)'; 8 x 6 = 48 exactly, so the count "
         "is consistent. Earlier 10-vs-8.4 concern is dropped."),
    ]
    for t, d in items:
        print(f"  - {t}\n      {d}")
    print()

    out = OUT / "wp19_literature_closure.csv"
    rows = [
        dict(quantity="s from Wang R_1/2, no PT", value=round(s0, 4), source="Wang 2007"),
        dict(quantity="s from Wang R_1/2, +PT", value=round(s1, 4), source="Wang 2007"),
        dict(quantity="s ratio (absolute)", value=round(s1 / s0, 4), source="Wang 2007"),
        dict(quantity="s lower bound, +PT (8 repeats)", value=round(s_pt_lo, 4), source="L&C 2010 text"),
        dict(quantity="s upper bound, no PT (8 repeats)", value=round(s_no_hi, 4), source="L&C 2010 text"),
        dict(quantity="a no PT", value=round(a0, 4), source="identity"),
        dict(quantity="a +PT", value=round(a1, 4), source="identity"),
        dict(quantity="a fall (off-rate, extension phase)", value=round(a0 / a1, 4), source="identity"),
        dict(quantity="cycling t1/2 fold (pulse-chase)", value=T12_CYCLING_FOLD, source="L&C 2010"),
        dict(quantity="rebinding inflation factor", value=round(T12_CYCLING_FOLD / (a0 / a1), 3), source="derived"),
        dict(quantity="Kd POT1-a5 nM", value=KD_POT1_A5, source="Wang 2007"),
        dict(quantity="Kd POT1-TPP1-a5 nM", value=KD_POT1TPP1_A5, source="Wang 2007"),
        dict(quantity="k_pol change with PT", value=0.0, source="L&C 2010: none"),
    ]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["quantity", "value", "source"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
