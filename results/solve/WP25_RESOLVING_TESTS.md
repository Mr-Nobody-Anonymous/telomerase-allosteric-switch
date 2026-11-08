# WP25 — What would settle gaps 2 and 3, with the numbers attached

Neither gap can be closed with data that exists. What *can* be done is to
convert "unresolved" into "one experiment, with a quantitative prediction for
each branch". Both are now enumerated.

    python scripts/solve/wp25_resolving_tests.py
    -> results/solve/wp25_resolving_tests.csv

## Gap 2 — M5 predicts 2.73×, L&C measured ≥4×

The ratio is **1.46×**. Only three things can produce a measured fold larger
than the distribution implies, and they are exhaustive:

| explanation | prediction | status |
|---|---|---|
| **E1** `k_tr` falls, so `t½` slows more | `k_tr` ratio **−0.003** | **EXCLUDED** |
| **E2** stalled complexes inflate `t½` | stall share 1.46× | open |
| **E3** ≥4× is a lower bound, not a value | true fold in [2.7, 4] | open |

**E1 is excluded numerically**: it requires `k_tr` to *fall* by ~300×, while L&C
fit it *rising* 0.07 → 0.10 /min. That is the only branch a calculation can kill,
and it is killed.

E2 and E3 cannot be separated from the reported numbers — and E3 is what L&C
themselves imply by writing "at least four-fold" and stating the data "could not
be fit to a single first-order constant".

**So the 1.46× gap is not evidence against M5.** It says the measured ≥4×
contains the true value 2.73× plus a stalling contribution.

### The resolving experiment

Add a **DNA trap** (or dilution) to the pulse-chase so released primer cannot
rebind and stalled complexes leave the bound pool:

| branch | predicted fold under trap |
|---|---|
| E2 — stalls inflate it | **falls to 2.73×** |
| E3 — ≥4× is a lower bound | **stays ≥4×** |

The assay already exists — L&C used an A5P chase; a trap is a reagent swap.

## Gap 3 — `k_tr` vs `k_off_tr`

`φ = k_tr/(k_tr + k_off_tr)` is two unknowns against one observation:

| `k_tr` (+PT) /min | `k_off_tr` implied /min | branch |
|---|---|---|
| 0.070 | 0.00071 | flat (M4) |
| 0.100 | 0.00101 | +43 % (M2) |
| 0.200 | 0.00202 | +186 % |

Every row reproduces φ = 0.99 **exactly**, so φ cannot separate M2 from M4 at any
precision — the WP18/WP20 result restated in M5's parameterisation.

**But this gap is not load-bearing.** φ left the product, so `s = (1+a)⁻⁶` never
sees `k_tr`. A direct `k_tr` measurement (or stopped-flow `k_off_tr`) under the
single-turnover condition L&C already used for their +2/+4 band assay would
settle it — different analysis of an existing assay, not new biology.

## What this changes

Neither gap is now an open-ended concern. Each has:

- a closed set of surviving explanations,
- the ones a calculation can eliminate already eliminated, and
- a stated experiment with a distinct numeric prediction per branch.

That is the most that inference can do. **What it cannot do is produce the data.**
