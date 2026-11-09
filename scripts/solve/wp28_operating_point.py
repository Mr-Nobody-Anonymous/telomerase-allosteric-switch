#!/usr/bin/env python3
"""
WP28 - Improving M6: the transcription layer is saturating, not linear.

WHY THIS EXISTS
---------------
WP27 concluded "the enzyme knob wins by 4.3-11.5x per percent".  That statement
is true and it is also the wrong decision rule, because it compares a LOCAL
slope at one operating point.  Two things are wrong with stopping there:

  1. The enzyme elasticity is not a constant 4-12.  eps = -1/ln(s) DIVERGES as
     s -> 1.  Quoting "4.3-11.5x" is quoting two particular operating points and
     calling them a range.

  2. M6 assumed output = R1/2(s) * [TERT], i.e. a LINEAR dependence on TERT
     abundance.  Telomerase must first RECRUIT to a telomere.  The number of
     elongating complexes saturates:

         c(x) = x / (1 + x)        x = [TERT] / K_d      (recruitment occupancy)

         output = N_tel * c(x) * R1/2(s)                  <- M6'

     At x << 1 output is linear in x (WP27's regime).  At x >> 1 it is FLAT:
     more TERT buys nothing because every telomere already has one.

The consequence is the opposite of WP27's headline, and it is the useful result:

    PER PERCENT      the enzyme knob always wins  (eps_enzyme >= 4.27 > 1 >= eps_TERT)
    IN TOTAL RANGE   the winner depends on where you START

A TERT-negative somatic cell sits at x ~ 0, where the transcription knob spans
orders of magnitude and the measured enzyme knob spans 2.70x.  A cell already at
saturating coverage sits at x >> 1, where the transcription knob is dead and the
enzyme knob is the only one that moves.

SO: "which knob is better" is not a property of the knobs.  It is a property of
the operating point, and the design should name the operating point.

WHAT IS STILL NOT CLAIMED
-------------------------
Raising TERT to coverage is not a therapy.  It is the immortalisation step, and
the OCT4 hazard from WP27 does not go away.  This is a statement about which
MEASUREMENT to prioritise and which ASSAY can see what - not a treatment.

    python scripts/solve/wp28_operating_point.py
    -> results/solve/wp28_operating_point.csv
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# L&C anchors, identical to WP22-WP27 so the two work packages are comparable
S_NO_PT = 2 ** (-1 / 2.963)      # 0.791414
S_PT = 2 ** (-1 / 8.0)           # 0.917004


def R_half(s: float) -> float:
    """Median number of repeats added, from per-repeat survival s."""
    return -math.log(2.0) / math.log(s)


def eps_s(s: float) -> float:
    """d ln(R1/2) / d ln(s) = -1/ln(s).  NOT a constant; diverges as s -> 1."""
    return -1.0 / math.log(s)


def occupancy(x: float) -> float:
    """Fraction of telomeres carrying an elongating complex, x = [TERT]/K_d."""
    return x / (1.0 + x)


def eps_x(x: float) -> float:
    """d ln(output) / d ln(x) for M6'.  Equals 1 - c, bounded by 1."""
    return 1.0 / (1.0 + x)


def output(s: float, x: float, n_tel: float = 1.0) -> float:
    """M6': total repeats added per cell.  Reduces to M6 x N_tel when x << 1."""
    return n_tel * occupancy(x) * R_half(s)


def fold_to_match(target_ratio: float, x0: float) -> float | None:
    """Fold increase in [TERT] (from x0) needed to match an output ratio.

    Returns None when the required coverage exceeds 1 (i.e. unreachable)."""
    c1 = occupancy(x0) * target_ratio
    if c1 >= 1.0:
        return None
    x1 = c1 / (1.0 - c1)
    return x1 / x0


def main() -> int:
    L = "=" * 92
    print(L)
    print("WP28  Improving M6: the operating point decides which knob binds")
    print(L)
    print()

    # ------------------------------------------------------------ WP28.1
    print("WP28.1  The enzyme elasticity is not a range - it diverges")
    print(f"  {'s':>8} {'R1/2':>9} {'eps = -1/ln s':>16}   note")
    print("  " + "-" * 68)
    grid = [(S_NO_PT, "L&C no-TPP1 anchor"),
            (S_PT, "L&C +TPP1 anchor"),
            (0.95, ""), (0.98, ""),
            (0.99, "WP27 quoted this as the top of a '4-12x' range"),
            (0.999, "1 - s = 1e-3")]
    for s, note in grid:
        print(f"  {s:8.4f} {R_half(s):9.3f} {eps_s(s):16.2f}   {note}")
    print()
    print("  'A 1% gain in s buys 4.3-11.5%' is true only AT the two anchors.")
    print("  The same knob at s = 0.99 buys 100%.  Quoting the anchors as a")
    print("  range understates the enzyme layer and gives the false impression")
    print("  that its advantage is bounded.  WP27's '4-12x' is a local slope")
    print("  reported as if it were a range.")
    print()

    # ------------------------------------------------------------ WP28.2
    print("WP28.2  M6' - the transcription layer is saturating, not linear")
    print("      occupancy   c(x) = x/(1+x)        x = [TERT]/K_d")
    print("      output      = N_tel * c(x) * R1/2(s)")
    print()
    print(f"  {'x':>9} {'coverage c':>11} {'eps_x = 1/(1+x)':>16} {'output @ s0':>13}")
    print("  " + "-" * 54)
    for x in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
        print(f"  {x:9.2f} {occupancy(x):11.4f} {eps_x(x):16.4f} "
              f"{output(S_NO_PT, x):13.4f}")
    print()
    print("  Two regimes, and they behave nothing alike:")
    print("    x << 1   output is LINEAR in TERT; every extra molecule helps")
    print("    x >> 1   output is FLAT; every telomere already has an enzyme")
    print()
    print("  M6 (WP27) is the x << 1 limit of M6'.  It was not wrong, it was")
    print("  UNCONDITIONAL, and an unconditional model cannot say which regime")
    print("  a design is in.  That is the thing that needed fixing.")
    print()

    # ------------------------------------------------------------ WP28.3
    print("WP28.3  The comparison, made correctly")
    print("  Per percent, at ANY operating point:")
    print(f"    enzyme       eps_s  >= {eps_s(S_PT):.2f}   (minimum, at the +TPP1 anchor)")
    print(f"    transcription eps_x <= 1.00   (maximum, at x -> 0)")
    print("  => the enzyme knob wins per percent, everywhere, always.")
    print("     WP27 was right about that and this does not retract it.")
    print()
    rng_enzyme = R_half(S_PT) / R_half(S_NO_PT)
    print("  In total range, from a stated start:")
    print(f"    enzyme       measured range  {rng_enzyme:.2f}x  = "
          f"{math.log(rng_enzyme):.3f} ln units  (L&C no-TPP1 -> +TPP1)")
    for x0, lbl in ((0.01, "x0 = 0.01 (near-absent TERT)"),
                    (0.1, "x0 = 0.10"),
                    (1.0, "x0 = 1.00 (half coverage)")):
        span = math.log(1.0 / occupancy(x0))
        print(f"    transcription {lbl:32s} {span:6.3f} ln units to full coverage")
    print()
    print("  The enzyme range is a MEASURED number (it is the 2.73x M5 predicts")
    print("  and L&C measure, expressed as a median-length ratio).  The")
    print("  transcription spans above are GEOMETRIC - how far you can travel,")
    print("  not how far anything measured has travelled.  That distinction is")
    print("  the whole point and it is why the two columns are not the same kind")
    print("  of number.")
    print()

    # ------------------------------------------------------------ WP28.4
    print("WP28.4  Sizing: what a transcription-layer design has to deliver")
    print("  To match the enzyme layer's measured 2.70x, [TERT] must rise by:")
    print(f"  {'start x0':>9} {'coverage':>9} {'fold [TERT] needed':>20}")
    print("  " + "-" * 42)
    for x0 in (0.001, 0.01, 0.1, 0.5, 1.0, 2.0):
        f = fold_to_match(rng_enzyme, x0)
        txt = f"{f:20.2f}" if f is not None else f"{'UNREACHABLE':>20}"
        print(f"  {x0:9.3f} {occupancy(x0):9.4f} {txt}")
    print()
    print("  From near-zero TERT, 2.70x of output costs about 2.7x of [TERT]")
    print("  because the regime is linear.  At one-third coverage (x=0.5) the")
    print("  same 2.70x costs 18x of [TERT]; at half coverage (x=1) it cannot")
    print("  be bought at any [TERT], because 2.70x would need coverage 2.70.")
    print("  The")
    print("  transcription knob's price rises without limit; the enzyme knob has")
    print("  no such penalty, which is the real argument for prioritising it.")
    print()

    # ------------------------------------------------------------ WP28.5
    print("WP28.5  Why 26 work packages never touched this layer")
    print("  Every assay in this repository is in vitro: purified enzyme, fixed")
    print("  DNA, fixed [TERT].  In the model that is df/dx = 0 by construction -")
    print("  the assay cannot move x, so it has ZERO sensitivity to the")
    print("  transcription layer however long it runs.")
    print()
    print("  That is not an oversight in the experiments.  It is a structural")
    print("  property of the assay, and it is why the blind spot is legitimate:")
    print("  WP1-WP26 measure s correctly and say nothing about x because they")
    print("  cannot.  Any test of the transcription layer has to be a CELL")
    print("  experiment with TERT abundance as the manipulated variable.")
    print()

    # ------------------------------------------------------------ WP28.6
    print("WP28.6  Revised recommendation (this SUPERSEDES WP27's bottom line)")
    print("  1. Do not quote '4-12x'.  Quote eps = -1/ln(s) with the operating")
    print("     point attached, and note it diverges as s -> 1.")
    print("  2. Replace M6 with M6' (saturating).  M6 is its x << 1 limit.")
    print("  3. Name the operating point before comparing knobs.  From near-zero")
    print("     TERT the transcription layer is the binding constraint; from")
    print("     saturating coverage it is dead and the enzyme layer is all there")
    print("     is.")
    print("  4. The in vitro programme is still the right first investment - not")
    print("     because its per-percent leverage is 4-12x, but because it is the")
    print("     only knob this assay can see, and its price does not rise with")
    print("     the operating point the way transcription's does.")
    print("  5. Still no therapeutic claim.  Reaching coverage is immortalisation,")
    print("     and the OCT4 hazard is unchanged.")
    print()

    rows = [dict(item="enzyme elasticity at s_no_TPP1", value=round(eps_s(S_NO_PT), 4),
                 note="local slope, not a range"),
            dict(item="enzyme elasticity at s_+TPP1", value=round(eps_s(S_PT), 4),
                 note="local slope, not a range"),
            dict(item="enzyme elasticity at s=0.99", value=round(eps_s(0.99), 4),
                 note="diverges as s->1"),
            dict(item="max transcription elasticity (x->0)", value=1.0,
                 note="eps_x = 1/(1+x), bounded by 1"),
            dict(item="enzyme measured range", value=round(rng_enzyme, 4),
                 note="ln units %.4f" % math.log(rng_enzyme)),
            dict(item="output model", value="N_tel * c(x) * R1/2(s)",
                 note="M6' saturating; M6 is x<<1 limit"),
            dict(item="fold [TERT] to match 2.70x from x0=0.01",
                 value=round(fold_to_match(rng_enzyme, 0.01), 4), note="linear regime"),
            dict(item="fold [TERT] to match 2.70x from x0=1.0",
                 value="unreachable", note="needs c=2.7>1"),
            dict(item="in vitro assay sensitivity to x", value=0.0,
                 note="structural - assay holds [TERT] fixed"),
            dict(item="therapeutic claim", value="none", note="coverage = immortalisation"),
            ]
    out = OUT / "wp28_operating_point.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- WP28.7 tests
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # 1. eps_x is exactly 1 - c
    ok = all(abs(eps_x(x) - (1.0 - occupancy(x))) < 1e-15 for x in (0.01, 0.5, 3.0))
    chk("eps_x == 1 - coverage analytically", ok)

    # 2. enzyme elasticity strictly beats the transcription elasticity everywhere
    ok = all(eps_s(s) > eps_x(x) for s in (S_NO_PT, S_PT, 0.99) for x in (1e-6, 1.0, 1e3))
    chk("enzyme elasticity > transcription elasticity at all x", ok)

    # 3. M6' reduces to linear (M6) when x << 1: output/x -> R1/2
    x = 1e-8
    ok = abs(output(S_NO_PT, x) / x - R_half(S_NO_PT)) < 1e-6
    chk("M6' -> M6 (linear) as x -> 0", ok)

    # 4. output saturates: from x=1e4 tenfold more TERT adds < 0.1% of output
    o1, o2 = output(S_NO_PT, 1e4), output(S_NO_PT, 1e5)
    chk("output saturated at high x (<0.1% for 10x TERT)", (o2 / o1 - 1.0) < 1e-3)

    # 5. the fold_to_match sizing is self-consistent: applying it hits the target
    f = fold_to_match(2.7, 0.01)
    got = output(S_NO_PT, 0.01 * f) / output(S_NO_PT, 0.01)
    chk("fold_to_match hits its target ratio", abs(got - 2.7) < 1e-9)

    # 6. unreachable above coverage 1 is reported as None, not a huge number
    chk("unreachable target returns None", fold_to_match(2.7, 2.0) is None)

    # 7. monotonicity in both knobs
    ok = (output(S_PT, 1.0) > output(S_NO_PT, 1.0)) and (output(S_NO_PT, 2.0) > output(S_NO_PT, 1.0))
    chk("output monotone in s and in x", ok)

    # 8. the range identity: the enzyme range is exactly the R1/2 ratio,
    #    i.e. 8.000/2.963 = 2.69997, which is the "2.7x" quoted throughout.
    rng = R_half(S_PT) / R_half(S_NO_PT)
    chk("enzyme range == R1/2 ratio (2.69997)", abs(rng - 8.0 / 2.963) < 1e-12)

    print()
    print("WP28.7  Self-test: does the improved model break anything already decided?")
    for name, ok in checks:
        print(f"    {name:52s} {'ok' if ok else 'FAIL'}")
    failed = sum(1 for _, ok in checks if not ok)
    print(f"    {len(checks)} checks, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    import sys
    sys.exit(_self_test() if "--self-test" in sys.argv else main())
