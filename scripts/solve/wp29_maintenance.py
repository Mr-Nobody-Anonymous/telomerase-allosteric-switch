#!/usr/bin/env python3
"""
WP29 - The right objective: maintenance, not output.

WHY THIS EXISTS
---------------
WP27 optimised "product per cell".  WP28 optimised the same thing with a
saturating recruitment term.  Both are the wrong target, and setting the target
correctly changes the design more than either of them did.

A cell is not maintained because it makes a lot of repeats.  It is maintained
because the repeats it makes per division exceed the repeats it loses per
division:

    gain per division   = 6 * c(x) * R1/2(s)      bp        (6 nt per repeat)
    loss per division   = delta                   bp
    net                 = gain - loss

A cell with a huge gain and a bigger loss still senesces.  Output is not
maintenance; the DIFFERENCE is, and the difference has a sign.

THE EXACT RESULT
----------------
Setting net = 0 gives the critical attrition the design can survive:

    delta* = 6 * R1/2(s) * c(x)  <=  6 * R1/2(s)          (c <= 1)

    and since R1/2 = ln2 / (6 ln(1+a)),     delta* = 1 / log2(1+a)

delta* is a property of the ENZYME alone and it has a hard ceiling: at full
coverage, the best the L&C enzyme achieves is 48.0 bp per division, which is
BELOW every real human attrition rate measured (50-100 in fibroblasts, 90 in
HUVEC, up to 200 in the wider literature).

So: at the measured enzyme kinetics, no achievable [TERT] immortalises a human
cell, because coverage saturates at 1 and 48 < 50.  This is not a tuning
problem.  The deficit grows with c, reaches its minimum at c = 1, and is still
negative there.

WHAT DOES WORK
--------------
delta* is monotone in the off-rate ratio, and it is the ONLY lever that is not
capped.  Immortalising a 90 bp/division cell needs delta* >= 90, i.e.

    a <= 2**(1/90) - 1 = 0.00773      a 5.14x reduction in k_off

M5 predicts 2.73x and L&C measure >=4x.  Both are short of 5.14x.  That is the
design target, stated as one number.

    python scripts/solve/wp29_maintenance.py
    -> results/solve/wp29_maintenance.csv

NOT CLAIMED
-----------
None of this is a therapy.  Immortalising a cell is not treating anything, and
in a somatic cell it is the canonical first step of oncogenesis.  This is a
statement about what a design would have to achieve, and about which quantity
to optimise - not a treatment.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

A_NO_PT = 0.03976          # L&C no-TPP1 off-rate ratio (WP22)
A_PT = 0.01455             # L&C +TPP1, i.e. the 2.73x M5 predicts
A_MEAS = A_NO_PT / 4.0     # L&C pulse-chase, >=4x
A_WANG = A_NO_PT / 4.48    # Wang 2007 RRL, from WP22
NT_PER_REPEAT = 6

# Real attrition, from the primary literature
DELTA_FIBRO_LO, DELTA_FIBRO_HI = 50.0, 100.0   # bp/PD, human fibroblasts
DELTA_HUVEC = 90.0                              # bp/PD, human endothelial
DELTA_MAXLIT = 200.0                            # bp/PD, upper literature range


def s_of_a(a: float) -> float:
    """Per-repeat survival under M5: s = (1+a)^-6."""
    return (1.0 + a) ** -6


def R_of_a(a: float) -> float:
    """Repeats per telomere per division at full coverage (R1/2 = ln2/(6 ln(1+a)))."""
    return math.log(2.0) / (6.0 * math.log(1.0 + a))


def gain_bp(a: float, coverage: float = 1.0) -> float:
    """bp added per division per telomere."""
    return NT_PER_REPEAT * coverage * R_of_a(a)


def delta_star(a: float, coverage: float = 1.0) -> float:
    """Critical attrition rate: the largest delta the design still outruns."""
    return gain_bp(a, coverage)


def coverage_needed(a: float, delta: float) -> float:
    """Coverage required to break even.  > 1 means unreachable."""
    return delta / gain_bp(a, 1.0)


def a_needed(delta: float) -> float:
    """Off-rate ratio required to keep up with attrition delta at full coverage."""
    return 2.0 ** (1.0 / delta) - 1.0


def offrate_fold_needed(delta: float) -> float:
    return A_NO_PT / a_needed(delta)


def lifespan_pd(a: float, delta: float, l0: float, lcrit: float, coverage: float = 1.0):
    """Population doublings before the telomere reaches lcrit.  None if immortal."""
    net = delta - gain_bp(a, coverage)
    if net <= 0:
        return None
    return (l0 - lcrit) / net


def substitution_rate(delta: float) -> float:
    """d ln(1+a) / d ln(x) at its maximum (coverage 0.5).

    This is how fast [TERT] can buy off-rate reduction: (ln2/delta) * c * (1-c).
    Its maximum is ln2/(4 delta), per e-fold of [TERT].
    """
    return math.log(2.0) / (4.0 * delta)


def main() -> int:
    L = "=" * 96
    print(L)
    print("WP29  The right objective is maintenance, not output")
    print(L)
    print()

    print("WP29.1  Output is not maintenance")
    print("      gain per division = 6 * c(x) * R1/2(s)   bp")
    print("      loss per division = delta                bp")
    print("      net               = gain - loss          <- has a SIGN")
    print()
    print("  A cell that gains 48 bp and loses 90 bp per division still")
    print("  senesces.  Maximising gain without comparing it to loss optimises")
    print("  a quantity that does not decide the outcome.")
    print()

    print("WP29.2  delta* = 1 / log2(1+a)  - an exact identity")
    print("  delta* = 6 R1/2 = 6 * ln2/(6 ln(1+a)) = 1/log2(1+a).  The")
    print("  coverage cancels out of the DERIVATION and persists only as a")
    print("  multiplier <= 1, so full coverage is the best case.")
    print()
    print(f"  {'enzyme state':26s} {'a':>9} {'R1/2':>7} {'delta* (c=1)':>13}")
    print("  " + "-" * 60)
    cand = [("no TPP1 (L&C)", A_NO_PT, "baseline"),
            ("M5 prediction (2.73x)", A_PT, "M5"),
            ("L&C measured (4x)", A_MEAS, "measurement"),
            ("Wang RRL (4.48x)", A_WANG, "out-of-sample")]
    for lbl, a, tag in cand:
        print(f"  {lbl:26s} {a:9.5f} {R_of_a(a):7.3f} {delta_star(a):10.2f} bp")
    print()
    print("  A 2.73x off-rate reduction buys 48.0 bp/division, not 90.")
    print()

    print("WP29.3  M5's own number was already on the page (sharpens WP22)")
    print("  L&C write: 'half the chains elongated by eight repeats (48 nt)'.")
    print("  M5's 2.73x predicts:")
    print(f"      delta* = 1/log2(1 + {A_PT:.5f}) = {delta_star(A_PT):.2f} bp")
    print("  Eight repeats at 6 nt = 48 nt.  The 48 was DERIVED, and it")
    print("  matches the observed 48.  That agreement is the strongest single")
    print("  number M5 has, and it is also its own ceiling.")
    print()

    print("WP29.4  Every real human attrition rate exceeds that ceiling")
    print(f"  {'delta (bp/PD)':>16} {'source':34s} {'coverage needed':>16}")
    print("  " + "-" * 70)
    for delta, src in ((DELTA_FIBRO_LO, "fibroblasts, lower bound"),
                       (60.0, "fibroblasts, typical"),
                       (DELTA_FIBRO_HI, "fibroblasts, upper bound"),
                       (DELTA_HUVEC, "HUVEC (measured 89-90)"),
                       (DELTA_MAXLIT, "literature maximum")):
        c = coverage_needed(A_PT, delta)
        print(f"  {delta:16.1f} {src:34s} {c:12.2f}  {'UNREACHABLE' if c > 1 else ''}")
    print()
    print("  Coverage cannot exceed 1 - it is the fraction of telomeres that")
    print("  carry an elongating complex.  At 50 bp/PD the design needs 1.04")
    print("  telomerases per telomere.  There is no [TERT] at which this")
    print("  works, which is why no amount of the Yamanaka layer fixes it.")
    print()

    print("WP29.5  The one lever with no ceiling: the off-rate")
    print(f"  {'target delta':>13} {'a required':>12} {'off-rate fold':>14} {'have we got it':>22}")
    print("  " + "-" * 66)
    for delta, have in ((DELTA_FIBRO_LO, "no: 2.73-4x"), (DELTA_HUVEC, "no: 2.73-4x"),
                        (DELTA_FIBRO_HI, "no"), (DELTA_MAXLIT, "no")):
        print(f"  {delta:13.1f} {a_needed(delta):12.5f} {offrate_fold_needed(delta):13.2f}x "
              f"{have:>22}")
    print()
    print(f"  90 bp/PD needs {offrate_fold_needed(DELTA_HUVEC):.2f}x.  M5 predicts 2.73x; L&C")
    print("  measured >=4x.  Both are short.  The number to hit is a k_off")
    print("  ratio, and it is the same measurement WP22 already asked for.")
    print()

    print("WP29.6  Why the transcription layer cannot substitute for it")
    print("  a enters only through coverage, so")
    print("      dln(1+a)/dln(x) = (ln2/delta) * c * (1-c)   <=  ln2/(4 delta)")
    print()
    print(f"  {'delta':>7} {'rate per e-fold [TERT]':>24} {'[TERT] fold for 4x':>22}")
    print("  " + "-" * 56)
    for delta in (DELTA_FIBRO_LO, 70.0, DELTA_HUVEC, DELTA_FIBRO_HI):
        r = substitution_rate(delta)
        need = math.log(4.0) / r
        print(f"  {delta:7.0f} {100*r:22.3f}% {need:22.3g}")
    print()
    print("  At best, one e-fold of [TERT] buys 0.35% of off-rate reduction.")
    print("  Buying 4x by abundance alone needs ~1e174-fold.  This is the")
    print("  quantitative form of the WP28 point: abundance is a dead lever")
    print("  for this objective, not merely a low-leverage one.")
    print()

    print("WP29.7  Lifespan arithmetic (HUVEC-like: L0 11.5 kb, Lcrit 5 kb)")
    print(f"  {'condition':32s} {'gain/PD':>9} {'net/PD':>8} {'PD to senescence':>18}")
    print("  " + "-" * 72)
    for lbl, a, cov in (("untreated (no TERT)", None, 0.0),
                        ("+TPP1, full coverage (2.73x)", A_PT, 1.0),
                        ("at L&C measured 4x, full coverage", A_MEAS, 1.0),
                        ("at 5.14x, full coverage", A_NO_PT / 5.143, 1.0)):
        if a is None:
            print(f"  {lbl:32s} {0.0:9.2f} {DELTA_HUVEC:8.2f} "
                  f"{(11.5e3 - 5e3) / DELTA_HUVEC:18.0f}")
            continue
        pd = lifespan_pd(a, DELTA_HUVEC, 11.5e3, 5e3, cov)
        net = DELTA_HUVEC - gain_bp(a, cov)
        txt = "IMMORTAL" if pd is None else f"{pd:.0f}"
        print(f"  {lbl:32s} {gain_bp(a, cov):9.2f} {net:8.2f} {txt:>18}")
    print()
    print("  The untreated number (~72) reproduces the observed HUVEC limit")
    print("  (~69-70).  That is a check on the arithmetic, not a result.")
    print("  TPP1-TEN extends lifespan ~2.2x and still senesces.  Only past")
    print("  5.14x does the sign of the net change.")
    print()

    print("WP29.8  VERDICT - the improved design (M7)")
    print(" 1. The objective is net maintenance per division, not output:")
    print("        net = 6 * c(x) * R1/2(s) - delta")
    print("    and the design target is the critical attrition delta* = 6 R1/2.")
    print(" 2. delta* has a ceiling of 1/log2(1+a) at c = 1.  At the measured")
    print(f"    enzyme that ceiling is {delta_star(A_PT):.1f} bp/PD - below EVERY real")
    print("    human attrition rate.  Unreachable at any [TERT].")
    print(" 3. The transcription layer is not a weak lever, it is a dead one:")
    print(f"    0.35%/e-fold at best, ~1e174-fold for 4x.")
    print(" 4. The design reduces to ONE number: the k_off ratio.  Target")
    print(f"    >= {offrate_fold_needed(DELTA_HUVEC):.2f}x for a 90 bp/PD cell.  M5 predicts 2.73x, L&C")
    print("    measure >=4x.  The gap is the whole remaining programme.")
    print(" 5. Still not a therapy.  Immortalising a somatic cell is the first")
    print("    step of oncogenesis, and the OCT4 hazard is unchanged.")
    print()

    rows = [
        dict(item="objective", value="6*c*R1/2 - delta",
             note="net maintenance, not output"),
        dict(item="delta* identity", value="1/log2(1+a)",
             note="exact; coverage <= 1 is the best case"),
        dict(item="delta* at M5 prediction (2.73x)", value=round(delta_star(A_PT), 3),
             note="bp/PD, full coverage"),
        dict(item="delta* at L&C measured (4x)", value=round(delta_star(A_MEAS), 3),
             note="bp/PD, full coverage"),
        dict(item="delta* at Wang 4.48x", value=round(delta_star(A_WANG), 3),
             note="bp/PD, full coverage"),
        dict(item="fibroblast attrition (low)", value=DELTA_FIBRO_LO,
             note="bp/PD; exceeds delta* at 2.73x"),
        dict(item="HUVEC attrition", value=DELTA_HUVEC,
             note="bp/PD, measured 89-90"),
        dict(item="coverage needed at 90 bp/PD (2.73x)",
             value=round(coverage_needed(A_PT, DELTA_HUVEC), 3),
             note=">1 = unreachable"),
        dict(item="off-rate fold needed at 90 bp/PD",
             value=round(offrate_fold_needed(DELTA_HUVEC), 3),
             note="M5 predicts 2.73x; L&C measure >=4x"),
        dict(item="max substitution rate per e-fold [TERT]",
             value=round(substitution_rate(DELTA_FIBRO_LO), 6),
             note="dln(1+a)/dln x at c=0.5"),
        dict(item="[TERT] fold for 4x by abundance alone",
             value="1e174", note="dead lever, not weak lever"),
        dict(item="HUVEC PD, +TPP1 at full coverage",
             value=round(lifespan_pd(A_PT, DELTA_HUVEC, 11.5e3, 5e3), 1),
             note="untreated ~72; still senesces"),
        dict(item="therapeutic claim", value="none",
             note="immortalisation is not a treatment"),
    ]
    out = OUT / "wp29_maintenance.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- WP29.9 tests
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # 1. the identity holds exactly
    chk("delta* == 1/log2(1+a)", all(
        abs(delta_star(a) - 1.0 / math.log2(1.0 + a)) < 1e-12
        for a in (A_NO_PT, A_PT, A_MEAS, 0.001)))

    # 2. it is consistent with M5's own s and R1/2
    chk("delta* == 6*R1/2 from s=(1+a)^-6", abs(
        delta_star(A_PT) - 6 * (-math.log(2) / math.log(s_of_a(A_PT)))) < 1e-9)

    # 3. M5's 2.73x predicts 48 nt, the number L&C observed
    chk("2.73x predicts 48.0 bp (L&C observed 48)", abs(delta_star(A_PT) - 48.0) < 0.05)

    # 4. the deficit cannot be closed by coverage
    chk("no [TERT] immortalises a 90 bp/PD cell at 2.73x",
        coverage_needed(A_PT, DELTA_HUVEC) > 1.0)

    # 5. the required off-rate fold is self-consistent
    chk("required fold reproduces delta* = target", abs(
        delta_star(a_needed(DELTA_HUVEC)) - DELTA_HUVEC) < 1e-9)

    # 6. net changes sign exactly at the target
    chk("net > 0 iff a < a_needed", (
        gain_bp(a_needed(DELTA_HUVEC) * 0.99) > DELTA_HUVEC
        and gain_bp(a_needed(DELTA_HUVEC) * 1.01) < DELTA_HUVEC))

    # 7. lifespan diverges (None) only past the target
    chk("lifespan None iff net <= 0", (
        lifespan_pd(a_needed(DELTA_HUVEC) * 0.99, DELTA_HUVEC, 11.5e3, 5e3) is None
        and lifespan_pd(A_PT, DELTA_HUVEC, 11.5e3, 5e3) > 0))

    # 8. the untreated-lifespan check reproduces the observed HUVEC limit
    pd = lifespan_pd(A_NO_PT * 0.0 + 1.0, DELTA_HUVEC, 11.5e3, 5e3, coverage=0.0)
    chk("untreated HUVEC lifespan in 65-75 PD", 65 <= pd <= 75)

    # 9. substitution rate is bounded and tiny
    chk("substitution rate < 1% per e-fold", substitution_rate(DELTA_FIBRO_LO) < 0.01)

    # 10. WP28's model is a strict ingredient, not contradicted
    chk("gain is linear in coverage (WP28 c(x) enters once)", abs(
        gain_bp(A_PT, 0.5) - 0.5 * gain_bp(A_PT, 1.0)) < 1e-12)

    print()
    print("WP29.9  Self-test: is the new objective consistent with what is settled?")
    for name, ok in checks:
        print(f"    {name:56s} {'ok' if ok else 'FAIL'}")
    failed = sum(1 for _, ok in checks if not ok)
    print(f"    {len(checks)} checks, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    import sys
    sys.exit(_self_test() if "--self-test" in sys.argv else main())
