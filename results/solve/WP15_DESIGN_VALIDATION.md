# WP15 — Experimental design validation package

Script: `wp15_design_validation.py` → `wp15_bounds.csv`, `wp15_active_fraction.csv`,
`wp15_trap_readouts.csv`, `wp15_hierarchical_power.csv`.

Three priorities as requested: active-fraction calibration, trap validation, hierarchical power.
One thing was added because it was available and removes a caveat: a proof replacing "tested
profiles".

---

## 15.0 The heterogeneity caveat is discharged by proof, not by more sampling

The audit's caution was correct: the tested profiles did not establish universality. But the claim
can be proven, so it no longer needs to rest on examples.

With `uᵢ = log(1 + aᵢ)` and `Σuᵢ = L = −log A₀`:

```
log R(f) = n·log(f) − Σᵢ h(uᵢ),     h(u) = log(1 + (f−1)e^(−u))
```

`h` is convex in `u` (it is a log-sum-exp). So `−Σh` is maximised at *uniform* `u` and minimised at
the most *concentrated* feasible vector. Both extremes have closed forms:

```
R_max = f^n · exp(−n·h(L/n))        [all steps equal]
R_min = f   · exp(−h(L))            [all survival in ONE step]
```

| A₀ | R_min | R_max | spread over ALL profiles |
|---|---|---|---|
| 0.99 | 1.0076 | 1.0076 | 0.00% |
| 0.95 | 1.0390 | 1.0392 | 0.02% |
| 0.90 | 1.0811 | 1.0820 | 0.09% |
| 0.80 | 1.1765 | 1.1813 | 0.41% |
| **0.70** | **1.2903** | **1.3041** | **1.07%** |
| 0.50 | 1.6000 | 1.6690 | 4.31% |

240,000 random profiles (5 Dirichlet concentrations + heavy log-normal) all fell inside the
predicted interval.

**So at A₀ = 0.70 a 4-fold fall in extension-phase k_off gives R ∈ [1.2903, 1.3041] for *every*
within-population heterogeneity, not just the tested ones.** The caveat is discharged for
per-position heterogeneity. It is **not** discharged for population mixtures, which are a different
model — that is 15.1, and it is where the real problem lives.

## 15.1 Active fraction w — and the failure mode nobody flagged

`A_obs = Σⱼ wⱼAⱼ` is one equation in `{wⱼ, Aⱼ}`, so **w is not identifiable from A**. Confirmed: the
audit's "circular reasoning" objection holds.

But the danger is the opposite of what the earlier analysis assumed. It is not a low `A_obs`. It is
a **low `A_obs` that passes the gate while containing an unmeasurable subpopulation**:

| mixture | w | A_obs | gate | R_obs | n_obs | n needed for worst subpopulation |
|---|---|---|---|---|---|---|
| 50/50 of A = 0.30 and **0.95** | 0.50 | 0.625 | **PASS** | 2.52 | 3 | **26** |
| 50/50 of A = 0.30 and **0.99** | 0.50 | 0.645 | **PASS** | 2.41 | 3 | **688** |
| 60/40 of A = 0.20 and **0.95** | 0.60 | 0.500 | **PASS** | 3.83 | 3 | **26** |
| 50/50 of A = 0.30 and 0.70 (both powered) | 0.50 | 0.500 | PASS | 3.40 | 3 | 3 |
| 50/50 of A = 0.70 and dead | 0.50 | 0.850 | fail | 1.35 | 3 | 3 |

Rows 1–3 pass the gate and look well powered at n = 3, while half the enzyme sits in a regime that
would need 26–688 replicates to measure. **The aggregate is powered; the population is not.**

Methods for w, ranked:
1. **active-site titration** against a purified standard — direct, gives a CI, calibration-dependent;
2. **TPP1 binding stoichiometry** — gives an upper bound on w;
3. product-length curvature across a dilution series — a consistency check only, because it shares
   information with A and is therefore partly circular.

**Required:** w with a 95% CI and the method. If the CI includes values below 0.8, the feasibility
gate is not applied.

## 15.2 Trap validation — two methods, five readouts, and method-dependence

Acceptance: every readout within ±10% of the no-trap control, **in both methods**.

| readout | dilution (100×) | excess DNA | heparin |
|---|---|---|---|
| basal nucleotide addition | 98% | 95% | 75% |
| TPP1 binding (K_d) | 97% | 94% | 85% |
| active concentration | 95% | 97% | 90% |
| primer dissociation | 10% | 5% | 60% |
| product-length distribution (s) | 96% | 93% | 80% |

*(Values are illustrative — the numbers show the structure of the requirement, not measurements.)*

The distinction that matters: dilution and excess-DNA trap remove rebinding without binding the
enzyme; **heparin binds the enzyme** and is the most likely to change k_pol or K_d. Preference order
is dilution > excess DNA ≫ heparin.

Because the trap is applied to both arms, its systematic bias does **not** cancel in the ratio. And
if the two methods disagree beyond their CIs, the inferred dissociation is method-dependent and
neither is reportable as `k_tr/k_off_tr` — the disagreement itself is then the result.

## 15.3 Hierarchical power — purification is the unit

Model: paired within purification, with a purification-level random slope on the effect:

```
log A_trt[p] − log A_ctrl[p] = mu + b_p + e_p,   b_p ~ N(0, tau²),  e_p ~ N(0, 2σ²)
⇒ Var = (tau² + 2σ²)/P over P purifications
```

| tau (prep SD on the effect) | P at A₀ = 0.70 | P at A₀ = 0.90 |
|---|---|---|
| 0.00 | 2 | 7 |
| 0.05 | 2 | 10 (1.4×) |
| 0.10 | 2 | 20 (2.9×) |
| 0.20 | 6 (3.0×) | 59 (8.4×) |

At A₀ = 0.70 the effect is large enough that prep noise barely matters. **At A₀ = 0.90 it is
decisive: tau = 0.10 costs ~3× and tau = 0.20 costs ~8×.** Extra reactions per preparation buy
nothing once tau > 0 — only more purifications reduce error.

WP13 and WP14 both sized on σ alone. That understates the requirement by up to ~8× at exactly the
operating point where the study is otherwise believed feasible.

## 15.4 Stopping rules

| condition | action |
|---|---|
| w CI includes values < 0.8 | **STOP** — preparation characterisation incomplete |
| tau cannot be estimated from a pilot | **STOP** — study cannot be sized |
| trap readouts disagree between methods | **STOP for channel B**; report A only |
| any trap readout outside ±10% | **STOP for channel B**; report A only |
| A_ctrl ≥ 0.90 at the pilot | **STOP** — no power; do not run the comparison |
| all checks pass | PROCEED; apply the WP14 decision tree |

Every branch terminates in a reportable outcome, including STOP. A study that stops for a stated
reason is a result; one that proceeds past a failed check is a hazard.

---

## Status

| question | conclusion |
|---|---|
| Does A ≤ 0.70 survive position-dependent k_pol? | **yes — proven**, not merely tested: R ∈ [1.2903, 1.3041] over all profiles |
| Does it survive population mixtures? | no |
| Were the sample sizes reliable? | no — corrected twice (exact-t, then hierarchical) |
| Can a trap perturb the discriminator? | yes; bias does not cancel across arms |
| Can w be inferred from A? | **no** — one equation, two unknowns |
| Is the mechanism identified? | **no** |

**WP15 is the last work package the mathematics can justify.** Everything further requires
measurement. The three pilot measurements in priority order are: **A on untreated enzyme** (one
reaction — sets the operating point), **w with a CI** (sets whether the gate is interpretable), and
**tau** (sets whether the study can be sized at all).
