# WP44 — All six blocked items, attempted

    python scripts/solve/wp44_attempt_all.py
    python scripts/solve/wp44_attempt_all.py --self-test    -> 31 checks, 0 failed

The instruction was to try to solve every item rather than accept the blocked
list. Doing so changed the status of five of the six — and **not in one
direction**. One was never load-bearing, one turned out to be vacuous on its own
terms, one gained a real (pre-existing) test, and three survive with narrower
blockers.

| # | item | before | after the attempt |
|---|---|---|---|
| 1 | measure `n` | blocked | **never load-bearing** — the gate is `n`-free |
| 2 | second-enzyme fold | blocked | **blocked** — panel has no rates; adds a mechanism constraint |
| 3 | P2 as two-chase | design exists | **the design is vacuous** — it is one relation written twice |
| 4 | TEL-patch affinity | blocked | **blocked, demonstrably** — 0 calibration ligands |
| 5 | Gate 3 densitometry | unreachable | **unreachable, now shown** — 2.3× spread from one number |
| 6 | Y176A panel | designed, not run | **blocked** — the assay is blind to attempts |
| — | *bonus* | — | **new pass**: published, model-free evidence for processive action |

---

## 1. `n` was never load-bearing

`n` enters M5 **only** through `a = s^(−1/n) − 1`, which *is* the identity
`δ* = n·R½` that C1 retracted. So I recomputed `R½` from the measured `s₀` all
the way through the exponent:

| n | 1 | 3 | 6 | 12 | 24 | 64 |
|---|---|---|---|---|---|---|
| `R½` | 7.530 | 7.903 | **8.000** | 8.050 | 8.074 | 8.090 |

**7.4 % across a 64-fold range of `n`.** The exponent cancels to first order.
`R½`, `δ* = 6·R½`, and the gate are all `n`-free, which is why WP37's
exponent-free definition was the right move and why WP29's `δ* = n·R½` was the
wrong question.

**Status: not a gap.** "No measurement of `n`" was a hole in an alternative
parameterisation, not in the argument. That is the single biggest thing this
work package changes.

---

## 2. The second-enzyme fold — blocked, now for a measured reason

The Tesmer panel's three readouts are **0.052, 0.170, 0.731** — all unitless
activity fractions. No rate, no affinity, no per-repeat survival. There is
nothing to fit a `k_off` to and no ratio that yields one.

What it *does* carry is a **mechanism constraint**:

| | processivity | activity | activity rank |
|---|---|---|---|
| **F776E** | 0.013 | 0.418 | **5 of 9** |
| mutant mean | 0.228 | 0.418 | — |

F776E sits mid-panel on activity while processivity is at the floor. A mutation
that leaves catalysis average while destroying processivity **cannot act on a
shared catalytic step**. That constrains M5's `a = k_off/k_pol` — the two are
separable, at least for this residue — and it is a published, independent fact.
It is still not a fold.

> **A naming error caught here.** My first draft named **H762E** as the
> discriminating mutant. It is not: H762E is high on *both* (0.416 processivity,
> 0.680 activity). The correct one is F776E.

---

## 3. The two-chase inversion is an identity — this is the uncomfortable one

`m = κ(1+τ)` and `c = κ(1−τ)` are **one relation written twice**. So `m/c`
returns `τ` by construction:

| τ in | m | c | τ out | error |
|---|---|---|---|---|
| 0.1833 | 3.1954 | 2.2054 | 0.1833 | 8.3e-17 |
| 0.3500 | 3.5875 | 1.7273 | 0.3500 | 5.6e-17 |

The round-trip to 1e-17 is **algebra, not biology**. It would succeed whether or
not any stalling exists. This is precisely the **C1 pattern**: a quantity
recovered from a construction that already contains it.

A real two-chase needs an **independent count of the stalled population** — the
trap-resistant fraction, not a second chase of the same survival. Then
`κ = c/(1−S)` must equal `m/(1+S)`, and disagreement refutes the model. WP37
called this "the design" and it was not one.

---

## 4. Affinity — blocked, and now demonstrably so

WP39's `−0.046` per heavy atom is a soft potential around one conformer. It has
**no solvent, no entropy, no flexibility, and zero ligands of known affinity** to
calibrate against. Converting a score to a `Kd` needs a calibration curve, and
there is nothing to build one from. Irreducibly wet lab or FEP.

---

## 5. Gate 3 — unreachable, now shown

Both plausible conventions hit the same reported fraction (0.485) at very
different survivals:

| convention | s | R½ |
|---|---|---|
| chain-count | 0.9529 | 2.395 |
| signal-mass | 0.8945 | 1.036 |

**2.3× spread in `R½` from a single number.** A threshold statistic is not
invertible to a per-repeat survival. This was asserted in WP24; it is now
computed.

---

## 6. The Y176A panel needs an assay that is blind to its own question

Zhao 2011's *in vivo* method is **BrdU/CsCl density**. It measures *how many
repeats were added*. It cannot measure *how many attempts were made*. A mutant
panel read through it inherits that blindness, so it cannot answer what it was
designed to answer.

### But the same paper contains a model-free mode discriminator

Partial inhibition (GRN163L, 0.5 µM ≈ 50 %, 1.0 µM ≈ 70 %):

| inhibition | extended fraction | product length |
|---|---|---|
| 0 % | 1.00 | ~60 nt |
| 50 % | 0.50 | ~60 nt |
| 70 % | 0.30 | ~60 nt |

A **distributive** enzyme spreads synthesis over successive rounds, so fewer
active molecules means **shorter** products. A **processive** one loses whole
ends and keeps length constant.

**Length is constant; the fraction falls.** That is published, pre-existing,
model-free positive evidence for processive action — and it removes an
alternative rather than merely failing to refute one.

It does **not** test M5's exponent. It is a constraint on the mechanism space.

---

## Bonus: the stall renewal does not break the framework

`s_obs = s_true / (1 − τ(1 − s_true))` is itself a geometric survival with a
rescaled per-repeat probability, so `R½` stays definable:

| s_true | s_obs (τ = 0.18) | apparent R½ |
|---|---|---|
| 0.700 | 0.7407 | 0.385 |
| 0.800 | 0.8304 | 0.622 |
| 0.917 | 0.9312 | 1.620 |
| 0.950 | 0.9588 | 2.745 |

**Caveat, stated:** this holds only if the retry is **memoryless**. If a stall
raises the probability of stalling again, no single `R½` exists and every
`R½`-based statement here becomes an effective value. Nobody had checked this.

---

## What this changes in the ledger

| row | before | after |
|---|---|---|
| `n` | unmeasured, blocking | **not load-bearing** — remove from the ledger |
| P2 design | exact, pre-registered | **vacuous** — needs an independent stall count |
| evidence | one fit datum | one fit datum **+ model-free mode evidence** |
| second enzyme | no fold | no fold; **+ F776E separability constraint** |

The programme is *smaller* and *better specified* than it was. Two things I
previously presented as strengths — the two-chase inversion and the blocked-`n`
caveat — were respectively an identity and a non-issue.

## Still genuinely blocked

1. A measured `k_off` fold in a second system.
2. An independent count of the stalled population.
3. Any `Kd` or `ΔG` for any compound.
4. Gate 3 per-band data.
5. An assay that counts attempts, for the Y176A panel.

All five need a bench. None is closed by more analysis here, and this work
package is the evidence for that claim rather than an assertion of it.
