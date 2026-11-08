# WP24 — Gate 3: what is blocked, exactly, and why re-analysis cannot fix it

Gate 3 (`paper/MANUSCRIPT.md` §11) asks for real per-band densitometry: a
full-distribution fit with φ free, showing φ > 0 and a majority of per-repeat
`s` below the tail slope.

    python scripts/solve/wp24_gate3_harness.py
    -> results/solve/wp24_gate3_harness.csv

## The proof: the recorded scalar cannot answer the question

`data/external/tesmer2023/` records, per condition, **the fraction of signal in
products ≥16 repeats** — one number. WP24 constructs two product-length
distributions that share that number **exactly**:

| | dist A | dist B | ratio |
|---|---|---|---|
| fraction at ≥16 (**the recorded one**) | 0.3500000000 | 0.3500000000 | **1.00×** |
| fraction ending at ≤3 repeats | 0.1270 | 0.6322 | **4.98×** |

Both sum to 1. A mechanism that releases most products after one or two repeats
and one that holds them all to 16 are **different enzymes** — and they are the
*identical row* in this table.

So the threshold fraction is invariant to the shape below repeat 16, which is
precisely where per-repeat survival is decided. **Gate 3 cannot be satisfied by
re-analysing the existing numbers at any level of care.** This is exact by
construction, not an estimate, and it is pinned by `tests/test_wp24_gate3.py`.

## What is missing, split into its two parts

| part | status |
|---|---|
| (A) per-band intensities | **ABSENT** — 4 CSVs, all scalar-per-condition |
| (B) an operational criterion | **NOT YET FIXED** |

(B) deserves to be stated plainly rather than hidden. The manuscript's phrase
"a majority `s` below the tail slope" needs two choices the manuscript does not
make: over which repeats the tail is fitted, and how a per-repeat `s` is
extracted from the fit. I deliberately did **not** fix those here. Choosing them
after seeing which synthetic cases pass would be fitting the criterion to the
result — the exact failure mode this repository exists to document. An earlier
draft of this script did precisely that and its own controls failed, which is
why it was rewritten rather than tuned.

## Status

Gate 3 remains **conditional**. It is no longer blocked by a missing tool — the
blocker is now stated precisely:

1. the repository has no product-length ladder, and
2. the criterion needs an operational definition fixed *before* data is seen.

Both are honest gaps, and neither is papered over with a plausible-looking
number.
