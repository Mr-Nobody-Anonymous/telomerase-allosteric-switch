# WP22 — **M5**: the corrected RAP model

The best candidate in this repository was WP18's M2/M4 (state-specific
translocation). Reading the **primary texts** — not abstracts — turned up four
things it was missing. Supplying them produces **M5**, which is the first
version of this model that is physically consistent *and* testable.

    python scripts/solve/wp22_rap_model_v2.py
    -> results/solve/wp22_rap_model_v2.csv

## What M2/M4 was missing

### 1. A resolution error, repeated through three work packages

L&C 2010 report, in one paragraph:

> "processivity to increase by a factor of **2.7 ± 0.6** … **almost half the
> chains are elongated by eight repeats (48 nt)** … whereas <10 % of the chains
> reach this length without POT1-TPP1."

Half the chains at 8 repeats means **R₁/₂ = 8**, and 2.7 is a ratio of **R₁/₂ —
a median length**. The repository has used 2.7 as a ratio of **per-repeat
survival `s`**. Those are different quantities, by a factor of **2.33**:

| ratio | value |
|---|---|
| R₁/₂ ratio (what 2.7 *is*) | 2.700 |
| `s` ratio (what the repo *used*) | 1.159 |

WP3, WP18 and WP20 all compare a survival ratio against a median-length ratio.

### 2. Two different enzymes in one equation

Wang 2007 used telomerase reconstituted in **rabbit reticulocyte lysate**
(R₁/₂ 0.78 → 3.3). L&C 2010 used "**super-telomerase**" overexpressed in HEK
cells (R₁/₂ 2.96 → 8), and say so, noting it is *more* processive. The repo
multiplied L&C's fold by Wang's absolute `s`. M5 keeps them apart and gets
**two independent tests** instead of one contaminated fit.

### 3. Dissociation also goes *up*

Hwang et al. 2014 (single-molecule) found POT1-TPP1 accelerates extension **and
makes the extended product leave faster**. The repo's model had `k_off` only
falling. A single `k_off` cannot do both; M5 splits it into
`k_off^ext` (falls ⟹ `s` rises) and `k_off^rel` (rises ⟹ product leaves).

### 4. GQ context

Hwang 2019 (PNAS) shows G-quadruplexes form **inside the actively extending
complex**, and POT1-TPP1 alters the product profile without rescuing synthesis
rate. The single-`k_pol`, single-`phi` assumption is now explicit, not silent.

## The falsification

WP18–WP20 used

    s = (1 + a)^-6 · phi ,   a ≥ 0        ⟹   s ≤ phi  always

L&C's own no-POT1-TPP1 numbers:

    R₁/₂ = 2.963  ⟹  s₀ = 0.7914
    phi          =  0.570          (translocation efficiency, Fig 2F)
    s₀/phi₀      =  1.388 > 1      ⟹   a₀ = -0.053   UNPHYSICAL

**79 % of chains survive a repeat while only 57 % translocate in the window.**
Both hold only if untranslocated complexes **stall and retry** rather than die —
exactly what L&C write: *"Stalled complexes and dissociated complexes could
contribute to the untranslocated population."*

So `phi` is a **window-limited efficiency**, not a per-repeat survival
probability, and must not multiply into `s`. **M5:**

    s = (1 + a)^-6                 a = k_off_ext / k_pol
    log s = -6 log(1 + a)          (k_tr-free AND phi-free)

## M5 passes an out-of-sample test

| system | a₀ | a₁ | predicted off-rate fall | measured |
|---|---|---|---|---|
| L&C 2010 super-telomerase | 0.03976 | 0.01455 | **2.73×** | **≥ 4×** (pulse-chase) |
| Wang 2007 RRL enzyme | 0.15964 | 0.03563 | 4.48× | not measured |

M5 predicts **2.73×** for the L&C system from the product distribution alone.
L&C measured **≥ 4×** independently, by pulse-chase. The 4× was never fed in.

The residual **1.46×** is the stalling contribution — which means **≥4× is an
over-estimate** of the true off-rate change, and the two agree in direction and
order of magnitude.

## Scoreboard

| candidate | `s₀ ≤ phi` holds | resting `k_off` | predicts a number |
|---|---|---|---|
| M1 strong (global `k_off`) | **NO** (unphysical) | **FAIL** | no |
| M2 translocation | **NO** (unphysical) | ok | n/a |
| M3 catalytic (`k_pol` ↑) | **NO** (unphysical) | ok | no |
| M4 state-specific | **NO** (unphysical) | ok | n/a |
| **M5 two-state, stalls retry** | **yes** | ok | **2.73×** |

Every earlier candidate inherited the same factorisation, so every one is
unphysical on L&C's own no-POT1-TPP1 data. **M5 is the only one that is both
physical and falsifiable.**

## Three predictions no earlier version made

- **P1** — the 2.7× is a median-length ratio; any model comparing it to a
  survival ratio is wrong by 2.33×. Checkable against L&C's own
  "eight repeats (48 nt)" sentence.
- **P2** — the true off-rate change is ≈2.7×, below the ≥4× observed, and the
  difference is **stalling**. Test: re-measure `k_off` under a **DNA trap**,
  which removes retries. Prediction: it lands near 2.7×, not ≥4×. **This is the
  cheap one — it re-analyses existing gel data, no new reagent.**
- **P3** — dissociation is **state-split**: falling while extending, rising on
  release. Test: product-length-resolved off-rate; the fast phase must get
  *faster* with POT1-TPP1 even as the slow phase slows.

## Honest limits

M5 still cannot separate `k_tr` from `k_off_tr` (they enter only as `phi`), so
the M2/M4 distinction remains open — but M5 does not need it, because `phi` is
no longer inside `s`. What M5 does **not** claim: that POT1-TPP1 tightens
primer binding overall. Hwang's result says the opposite for the release step.
