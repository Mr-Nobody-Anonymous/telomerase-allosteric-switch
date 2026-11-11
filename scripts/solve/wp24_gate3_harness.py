#!/usr/bin/env python3
"""
WP24 - Gate 3: why a threshold fraction cannot satisfy it, and what a harness
       would need.

GATE 3 (paper/MANUSCRIPT.md section 11)
---------------------------------------
    "real per-band densitometry: does a full-distribution fit with phi free
     show phi > 0 and a majority s below the tail slope?"

Two things are needed: (A) per-band band intensities, and (B) a fit that uses
the WHOLE distribution.  This script establishes what can be established
without the data, and is deliberately explicit about what it does not invent.

WHAT THIS SCRIPT PROVES  (T1)
-----------------------------
data/external/tesmer2023/ records, per condition, the fraction of signal in
products of >= 16 repeats - one scalar.  That scalar sums the tail and discards
everything below repeat 16.  Therefore ANY two product-length distributions
that agree on the tail have identical recorded value, whatever they do below
it - including whether per-repeat survival is flat or collapsing.

This is constructed explicitly below: two distributions with EXACTLY equal
>=16 fractions that imply opposite answers about whether s declines.  The
collision is exact, not approximate.  Hence:
    a threshold fraction is not a sufficient statistic for Gate 3,
    so no amount of re-analysis of the existing numbers can satisfy it.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
It does not implement a Gate-3 fit, because the manuscript's criterion ("a
majority s below the tail slope") needs an operational definition - which
residues the tail is fitted over, and how a per-repeat s is extracted - that
has not been fixed.  Choosing that here, after seeing which synthetic cases
pass, would be the same error this repository exists to document.  The
requirement is recorded in Gate 3's status instead.

Run:  python scripts/solve/wp24_gate3_harness.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
TESMER = REPO / "data" / "external" / "tesmer2023"

THRESHOLD = 16
N_MAX = 40


def distribution(s_func, phi: float = 1.0, n_max: int = N_MAX) -> list[float]:
    """P(product ends with exactly i repeats), i = 0..n_max."""
    out, alive = [], 1.0
    for i in range(n_max):
        surv = s_func(i) * phi
        out.append(alive * (1.0 - surv))
        alive *= surv
    out.append(alive)
    total = sum(out)
    return [x / total for x in out]


def thresh_fraction(p: list[float], k: int = THRESHOLD) -> float:
    return sum(p[k:])


def flb_slope(p: list[float], kmax: int = 14) -> float:
    """ln(1-FLB) vs repeat number, as in Latrick & Cech / Tesmer quantitation."""
    xs, ys = [], []
    for i in range(1, kmax):
        flb = sum(p[:i])
        if 0 < flb < 1:
            xs.append(i)
            ys.append(math.log(1.0 - flb))
    if len(xs) < 2:
        return float("nan")
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return (sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            / sum((x - mx) ** 2 for x in xs))


def make_collision() -> tuple[list[float], list[float], float]:
    """Two distributions with EXACTLY equal mass at and above THRESHOLD.

    Construction: fix the tail (i >= THRESHOLD) to a common shape with total
    mass T, and give each distribution a head (i < THRESHOLD) that sums to
    1 - T but with a different SHAPE - one flat, one collapsing.  Both then
    sum to 1 and both have >= THRESHOLD fraction exactly T.
    """
    head_n = THRESHOLD                      # i = 0..THRESHOLD-1
    tail_n = N_MAX + 1 - THRESHOLD          # i = THRESHOLD..N_MAX

    # common tail: a smooth decay, then normalised to mass T
    tail_shape = [math.exp(-0.10 * (i - THRESHOLD)) for i in range(tail_n)]
    T = 0.35
    tail = [T * x / sum(tail_shape) for x in tail_shape]

    head_mass = 1.0 - T
    # flat head: nearly uniform
    flat = [1.0 + 0.05 * i for i in range(head_n)]
    # collapsing head: steeply decaying
    coll = [math.exp(-0.90 * i) for i in range(head_n)]

    a = [head_mass * x / sum(flat) for x in flat] + tail
    b = [head_mass * x / sum(coll) for x in coll] + tail
    return a, b, T


def main() -> int:
    print("=" * 92)
    print("WP24  Gate 3: the scalar cannot answer it, and that is now exact")
    print("=" * 92)
    print()

    a, b, T = make_collision()
    fa, fb = thresh_fraction(a), thresh_fraction(b)
    sa, sb = flb_slope(a), flb_slope(b)

    print("WP24.1  T1: an exact collision at the recorded threshold")
    print(f"  Two product-length distributions, >=16 tail forced to be identical:")
    print()
    print(f"  {'':34s} {'dist A':>12} {'dist B':>12}")
    print("  " + "-" * 60)
    print(f"  {'recorded value (fraction >=16)':34s} {fa:12.10f} {fb:12.10f}")
    print(f"  {'identical to 10 dp?':34s} {str(abs(fa-fb) < 1e-10):>12} {'':>12}")
    print(f"  {'inferred survival from ln(1-FLB)':34s} {math.exp(sa):12.4f} {math.exp(sb):12.4f}")
    print(f"  {'shape below repeat 16':34s} {'flat':>12} {'collapsing':>12}")
    print()
    early_a = sum(a[:4])
    early_b = sum(b[:4])
    print(f"  The single number Tesmer et al. record is {fa:.10f} for BOTH.")
    print()
    print(f"  {'shape statistic':38s} {'dist A':>10} {'dist B':>10} {'ratio':>8}")
    print("  " + "-" * 70)
    print(f"  {'fraction ending at <=3 repeats':38s} {early_a:10.4f} {early_b:10.4f} "
          f"{early_b/max(early_a,1e-9):8.2f}x")
    print(f"  {'fraction at >=16 (the recorded one)':38s} {fa:10.4f} {fb:10.4f} "
          f"{fb/max(fa,1e-9):8.2f}x")
    print()
    print(f"  The early-product fraction differs by {early_b/max(early_a,1e-9):.2f}x while the")
    print(f"  recorded value is identical.  A mechanism that releases most products")
    print(f"  after one or two repeats and one that holds them all to 16 are")
    print(f"  different enzymes - but they are the same row in this table.")
    print()

    print("  So: a threshold fraction is invariant to the SHAPE below repeat 16,")
    print("  and that shape is exactly where per-repeat survival is decided.  This")
    print("  is a proof, not an estimate: the equality is exact by construction.")
    print()

    # ---------------------------------------------------------------- T2
    print("WP24.2  T2: what the repository holds")
    files = sorted(TESMER.glob("*.csv")) if TESMER.exists() else []
    print(f"  {TESMER.relative_to(REPO)}: {len(files)} file(s)")
    for f in files:
        with f.open() as fh:
            rows = [r for r in fh if not r.startswith("#")]
        ncol = len(rows[0].split(",")) if rows else 0
        print(f"    {f.name:34s} {max(0, len(rows)-1)} rows x {ncol} cols")
    print()
    print("  All scalar-per-condition.  No product-length ladder anywhere, so")
    print("  there is nothing to fit a distribution to.")
    print()

    # ---------------------------------------------------------------- T3
    print("WP24.3  T3: what a Gate-3 run would require")
    print("  Input: one row per repeat number, one column per lane, raw")
    print("         background-subtracted band intensities.")
    print("  Fit:   full-distribution likelihood with per-repeat s_i free and a")
    print("         global phi; report whether phi > 0 and whether a majority of")
    print("         the fitted s_i lie below the tail-implied slope.")
    print()
    print("  It is NOT implemented here.  The manuscript does not fix two choices")
    print("  that decide the answer - which repeats the tail is fitted over, and")
    print("  how s_i is extracted from the fit.  Picking those now, after seeing")
    print("  which synthetic cases pass, would be fitting the criterion to the")
    print("  result.  The gap is recorded in Gate 3's status rather than papered")
    print("  over with a plausible-looking number.")
    print()

    # ---------------------------------------------------------------- verdict
    print("WP24.4  Status of Gate 3, precisely")
    print("  per-band data      : ABSENT  (scalar threshold fractions only)")
    print("  sufficiency        : PROVEN inadequate (T1, exact collision)")
    print("  criterion          : NOT YET OPERATIONAL (two choices unfixed)")
    print("  verdict            : Gate 3 is CONDITIONAL and cannot be closed by")
    print("                       re-analysis.  It needs new data AND a fixed")
    print("                       criterion.  Both are stated, neither is faked.")
    print()

    rows = [
        dict(test="T1_collision", quantity="fraction >=16, dist A",
             value=f"{fa:.10f}", note="identical by construction"),
        dict(test="T1_collision", quantity="fraction >=16, dist B",
             value=f"{fb:.10f}", note="identical by construction"),
        dict(test="T1_collision", quantity="inferred survival A vs B",
             value=f"{math.exp(sa):.4f} vs {math.exp(sb):.4f}",
             note="differ despite equal recorded value"),
        dict(test="T2_inventory", quantity="per-band ladders present",
             value="0", note="4 csv files, all scalar-per-condition"),
        dict(test="T3_criterion", quantity="operational definition fixed",
             value="no", note="tail window and s_i extraction unspecified"),
    ]
    out = OUT / "wp24_gate3_harness.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "quantity", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
