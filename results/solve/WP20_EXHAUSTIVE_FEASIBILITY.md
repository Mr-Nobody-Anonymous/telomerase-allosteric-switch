# WP20 — Exhaustive feasibility: every mechanism, swept

> **⚠ SUPERSEDED IN PART BY WP22.** Two things here are wrong:
> (i) `2.7` is a ratio of **R₁/₂ (median length)**, not of per-repeat survival `s`.
> L&C's own text gives R₁/₂ = 8 with POT1-TPP1, so the `s` ratio is **1.159**, not
> 2.700 — the "ceiling violation" below is a **units error**, not a contradiction.
> (ii) Wang's `s₀` (RRL-reconstituted enzyme) and L&C's 2.7 (HEK
> "super-telomerase") are **different enzymes** and must not be multiplied in one
> equation. WP22 scores them separately.
> What still stands: the **T2 sweep** and the **T3 Kd route**, both computed
> within a single system. See `WP22_RAP_MODEL_V2.md`.

WP19 gathered the literature numbers. This asks the question that follows:
**given those numbers, which mechanism combinations are actually possible?**
The answer is a computed feasible set, not an argument.

    python scripts/solve/wp20_exhaustive_feasibility.py
    -> results/solve/wp20_exhaustive_feasibility.csv

## The model

    phi = k_tr / (k_tr + k_off_tr)
    s   = (1 + a)^-6 · phi ,   a = k_off_ext / k_pol

`a0` is not free — the absolute anchor fixes it:

    Wang 2007: R_1/2 = 0.78 repeats  ->  s0 = 2^(-1/0.78) = 0.4112
    a0 = (s0/phi0)^(-1/6) - 1 = 0.05593

## T1 — A hard ceiling, and it is violated

Because `k_off ≥ 0`, we always have `a1 ≥ 0`, hence

    s1 = (1+a1)^-6 · phi1  ≤  phi1

so **`s1/s0 ≤ phi1/s0 = 0.99/0.4112 = 2.408`**, whatever the mechanism.

| Target | Value | vs ceiling 2.408 | |
|---|---|---|---|
| Wang 2007 absolute `s` ratio | 1.971 | below | feasible |
| L&C 2010 processivity fold | 2.700 ± 0.6 | **above** | **impossible** |

This is a bound on the *data*, not on a model. At `s0 = 0.4112` the identity
**cannot** produce a 2.7-fold change in per-repeat survival. So at least one of
the two published numbers is not per-repeat `s`.

**Sensitivity:** reaching 2.7× needs `s0 ≤ 0.3667`, i.e. `R_1/2 ≤ 0.691`
instead of 0.78 — an 11.4 % shift, well inside a densitometric read-off error.
This is a **cross-paper calibration tension, not a refutation**. It is the
clearest concrete reason the repo still needs real per-band densitometry
(Gate 3).

## T2 — Exhaustive sweep

`koff_drop` and `ktr_drop` each swept over 1→1000 (37 points, log grid) with
`phi1` fixed at 0.99:

| observable | feasible? | points |
|---|---|---|
| range covered by the sweep | — | 1.7368 → 2.4068 |
| Wang `s` ratio 1.971 ± 0.05 | yes | 74 |
| L&C 2.7 ± 0.6 | yes (within ±1 SD) | 1184 |

`ktr_drop` is **irrelevant** to every one of these: without an independent
`k_tr`, the two arms cannot be separated, so all `ktr_drop` values tie. That is
WP18's degeneracy, now confirmed by enumeration rather than by a single case.

## T3 — The Kd route, and what it is worth

Primer `Kd` falls **3.51×** (Wang 2007: 26.0 → 7.4 nM). If that acts by raising
`k_on`, then `a` falls with `Kd` while `k_off` is untouched:

    a1 = 0.05593/3.51 = 0.01594   s1 = 0.9004   ratio = 2.1897

Wang's target is 1.971, so this **overshoots by 11.1 %**. It is not a match —
but it is the closest single mechanism in the sweep, it uses one measured
number with no free parameters, and it never touches `k_off_ext`.

## T4 — Candidates against every measured item

| Candidate | `s` ratio | vs Wang 1.971 | resting `k_off` | φ |
|---|---|---|---|---|
| M1 strong (global `k_off` ↓) | 1.2607 | −0.710 | **FAIL** | **FAIL** |
| M1 weak (`k_on` ↑, Kd only) | 2.1897 | +0.219 | ok | ok |
| M2 (`k_tr` ↑ + `k_off_tr` ↓) | 1.7368 | −0.234 | ok | ok |
| M4 (`k_off_tr` ↓ only) | 1.7368 | −0.234 | ok | ok |
| M3 (`k_pol` ↑) | 1.2607 | −0.710 | ok | **FAIL** |

No single-rate candidate lands on 1.971. Wang's value sits **between** the
φ-only case (1.74) and the Kd-driven case (2.19) — the signature of both arms
contributing partially, which no one-rate model can express.

## T5 — What would break the degeneracy

With `a0` anchored, M2 and M4 both give exactly `s1 = phi1 = 0.99`. Their whole
difference lives in `k_tr`, which both `s` and φ integrate out:

    implied k_off_tr = k_tr (1 - phi)/phi
      no TPP1,  k_tr = 0.07  ->  0.05281 /min
      + TPP1,   k_tr = 0.07  ->  0.00071 /min
      + TPP1,   k_tr = 0.10  ->  0.00101 /min

A direct `k_tr` measurement (or stopped-flow `k_off_tr`) settles it. **More
processivity data cannot**, at any precision — that is now proven, not asserted.

## Verdict

1. The identity caps the achievable fold-change at **2.408×** given `s0 = 0.4112`.
   L&C's 2.7× exceeds it; Wang's 1.971× does not. One of them is not per-repeat `s`.
2. The measured 3.51× `Kd` change alone reaches 2.19× — closest single mechanism.
3. M1-strong dies on the resting `k_off`; M3 dies on φ and on the 35.8–42.1 Å reach.
4. **M2 and M4 remain exactly degenerate.** The decisive experiment is a direct
   `k_tr` measurement.
