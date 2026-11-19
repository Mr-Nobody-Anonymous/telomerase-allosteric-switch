# Preregistration — A1: is `s_true = φ`?

**Registered before execution.** This document is written so that it cannot be
adjusted after the data arrive. If the result is ambiguous, the outcome is
"ambiguous", not "adjusted".

- **Test ID:** A1
- **Gap class:** G2 (measurement) — issue [#16](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/16)
- **Informational value:** EIG **0.4547 bits**, rank 1 of 5 (WP46)
- **Decides:** whether **R12** (M5 refuted) and **R13** (gate anchor defect) stand
- **Registered:** WP46, commit on branch `arena/3d3bb488-telomerase-allosteric-switch`

---

## 1. The claim under test

WP45 derived `τ₀ = 0.6506` from two of Latrick & Cech's own published numbers:

```
s_obs(−PT) = 0.7914        R½ = 2.963, from the product distribution
φ(−PT)     = 0.570          translocation efficiency, Fig 2F

τ₀ = (s_obs − φ) / (s_obs (1 − φ)) = 0.6506
```

This is only valid if **φ is the per-repeat success probability that the stall
later inflates** — i.e. that translocation efficiency and extension efficiency are
the *same* per-repeat event, with the stall acting only on extension.

**That is the assumption. It has never been tested. Everything in R12 and R13 rests
on it.**

## 2. Prediction, stated numerically in advance

If A1 holds, then a preparation with reduced translocation efficiency (via
competitor DNA / diffusive-access perturbation) should show:

| quantity | predicted under A1 | prediction if A1 is FALSE |
|---|---|---|
| `φ` | falls | falls |
| `s_obs` | falls **sub-proportionally** — the stall partly compensates | falls **proportionally** |
| `τ` recovered from the pair | **constant** at 0.6506 ± 0.08 | drifts with `φ` |
| `R½` | falls less than `φ` does | tracks `φ` |

**The discriminating statistic is the constancy of τ across a φ series.** Under A1
it is an invariant; if φ and `s_true` are different quantities, τ computed this way
is meaningless and will drift.

### The numeric prediction

For each point on the φ ladder, the two readings give these `s_obs` values. Under
A1 the observed `s_obs` must sit on the middle column; if A1 is false it tracks the
last column instead:

| φ (set) | `s_obs` predicted under A1 | τ recovered | `s_obs` if A1 false (τ→0) | separation |
|---|---|---|---|---|
| 0.300 | 0.5509 | 0.6506 | 0.3000 | **0.251** |
| 0.450 | 0.7007 | 0.6506 | 0.4500 | 0.251 |
| 0.570 | 0.7914 | 0.6506 | 0.5700 | 0.221 |
| 0.700 | 0.8698 | 0.6506 | 0.7000 | 0.170 |
| 0.850 | 0.9419 | 0.6506 | 0.8500 | 0.092 |

*Computed from the renewal at τ = 0.6506, not written from memory — the first
draft of this table contained plausible values that were wrong by up to 0.1.*

The separation is **largest at low φ** (0.251 at φ = 0.30) and shrinks to 0.092 as
φ rises, so **the discriminating power is in the low-φ points.** A ladder spanning
only 0.70–0.90 tests almost nothing — stated here so the range cannot be chosen
after the fact to suit the result.

## 3. Design

- **System:** *Tetrahymena* telomerase + TPP1, as in Latrick & Cech 2010 (matched
  to the source of `s_obs` and `φ`, not a substitute).
- **Perturbation ladder:** ≥5 points spanning `φ` ≈ 0.30–0.90, titrated with
  competitor primer or an equivalent diffusive-access modifier.
- **Per point, measure both** `φ` (translocation, Fig 2F method) and the product
  distribution (for `s_obs`). Measuring only one is not a test.
- **Replicates:** n ≥ 3 independent preparations per point.

## 4. Stopping rule

- **n is fixed at 3 per point.** No optional continuation, no "collect until
  significant".
- If a point fails QC (signal < 3× background, or φ undefined), **drop the point
  and report the drop**; do not replace it with a better-looking one.

## 5. Analysis, fixed in advance

1. Compute τ at each point from the pair, using the WP45 formula
   (`wp45_creative.tau_from_pair`).
2. **Primary statistic:** the coefficient of variation of τ across the ladder.
3. **PASS (A1 holds):** CoV(τ) < 0.15 **and** mean τ within 0.08 of 0.6506.
4. **FAIL (A1 false):** CoV(τ) > 0.30, or τ correlates with φ at |r| > 0.7.
5. **AMBIGUOUS (report as such):** anything else. **Ambiguous must not be reported
   as a pass.**

## 6. Outcomes and what each does to the programme

| outcome | consequence |
|---|---|
| **PASS** | R12 and R13 stand. M5 stays refuted. δ\* is 7.9 bp/PD, not 48. The plan moves to B1. |
| **FAIL** | **R12 and R13 are withdrawn.** M5 returns as a viable mechanism and `P(any)` must be recomputed. This is a large, publishable reversal and must be reported as one — not quietly. |
| **AMBIGUOUS** | Nothing changes and nothing is claimed. A1 stays open and the programme does not pretend to have closed it. |

## 7. Independent replication slot

One full repeat of the ladder at a **second site**, with the analysis script run by
someone who did not write this document. **Booked before the first run.** If no
second site is available, that is a limitation to state up front, not to discover
later.

## 8. What this test cannot settle

- It does not test **A2** (memoryless retry). That is a separate preregistration.
- It does not produce any binding datum, so it does not touch **G1**.
- A PASS confirms the *consistency* of the WP45 derivation; it does not prove M5 is
  the only mechanism M1–M6 could act through.

---

*Register entry: this file is scanned by `tests/test_open_findings.py`. Removing
the prediction table or the stopping rule breaks the build.*
