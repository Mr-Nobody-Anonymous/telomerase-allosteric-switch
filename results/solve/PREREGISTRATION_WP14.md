# Preregistration — TPP1 retention/allostery assay, WP14 revision

**Status:** supersedes the assay section of WP13. Written before any data are collected.
**Revised by WP15** — see `WP15_DESIGN_VALIDATION.md`. Three further changes are incorporated
below: the heterogeneity caveat is now a proof, w must be measured independently, and power is
computed hierarchically over purifications.
**Scope:** the two-channel (A, B) experiment only. Says nothing about the structural work in
WP1–WP9, which is already complete and reported.

**Why a revision.** WP14 found three defects in the WP13 design that would have produced a
misleading study if carried forward. They are recorded in §0 rather than quietly fixed.

---

## §0. What changed from WP13, and why

| WP13 said | WP14 found | status |
|---|---|---|
| n = 686 / 26 / 6 at CV 5% (A = 0.99 / 0.95 / 0.90) | closed form is a z approximation; exact two-sample t needs **686 / 28 / 8**. At A = 0.90, n = 6 gives **69% power, not 80%** | corrected |
| n = 1 is sufficient when A ≤ 0.70 | n = 1 is not a runnable design; **n = 3** is the floor to estimate variance | corrected |
| feasibility gate: A ≤ 0.70 | valid only if the active fraction **w ≥ 0.8**; a 50% inert preparation reads A = 0.815 and is wrongly rejected while containing 50% fully responsive enzyme | corrected |
| A near 1 is "uninformative" | confirmed, and now quantified: at A = 0.99 a null is consistent with a 99-fold fall in k_off_ext | retained |
| trap is advisable | trap moves **both** baseline A and effect size, in the same direction, and does not cancel across arms | strengthened |
| A ≤ 0.70 survives "tested profiles" | **proven** over all within-population profiles: R ∈ [1.2903, 1.3041] at A₀ = 0.70 | upgraded |
| w ≥ 0.8 gate | the hazard is a *low* A_obs that passes the gate while containing a subpopulation at A ≥ 0.95 needing n = 26–688 | corrected |
| n from a z-formula on σ | must be computed hierarchically: tau = 0.10 raises P from 7 to 20 at A₀ = 0.90 | corrected |

---

## §1. Replicate unit — this governs everything

**The replicate is an independent biological preparation (a separate purification/transfection).**
Technical replicates — repeat loadings of the same reaction, repeat gels, repeat densitometry — are
recorded for QC but **never counted toward n**.

Rationale: the variance that limits this design is preparation-to-preparation (activity titre,
assembly efficiency, TPP1 stoichiometry), not pipetting. Counting technical replicates as n would
inflate the apparent power by a factor we cannot estimate in advance.

**Consequence:** all sample sizes below are in preparations. A study with 8 preparations that are
split into 3 reactions each has n = 8, not n = 24.

---

## §2. Analysis model (fixed in advance)

Two arms: control (no TPP1) and treated (saturating TPP1–POT1). Both arms use the *same*
preparation, paired.

- Primary estimate: **A = s/φ**, with s from a full product-length distribution fitted by
  likelihood (never a length cut-off) and φ from a single-round translocation assay on the same
  material.
- Secondary estimate: **B = φ/(1−φ)**, reported as a **lower bound** from the exact
  (Clopper–Pearson) lower limit of φ whenever φ̂ > 0.95.
- Test: **two-sided, α = 0.05**, on log(A) with a paired design (same preparation across arms).
  Paired testing is the default because it removes between-preparation variance; the unpaired
  sample sizes in §3 are the conservative case.
- Effect size threshold for "detected": |A ratio − 1| ≥ 0.05 **and** 95% CI excluding 1.

---

## §2b. Active fraction (WP15)

**w must be measured independently of A.** `A_obs = ΣⱼwⱼAⱼ` is one equation in `{wⱼ, Aⱼ}`, so w is
not identifiable from A. Primary method: active-site titration with a 95% CI. Secondary, as an
upper bound: TPP1 binding stoichiometry. Product-length curvature is a consistency check only
(partly circular).

Report w with its CI. **If the CI includes values below 0.8, the feasibility gate is not applied.**

A low A_obs is not by itself reassuring: a 50/50 mixture of A = 0.30 and A = 0.99 reads A_obs = 0.645,
passes the gate, and looks powered at n = 3 while half the enzyme needs n = 688 to measure.

## §3. Sample size (exact two-sample t, 80% power, CV 5% on A)

| baseline A | R if k_off_ext falls 4× | n per arm (exact) | power |
|---|---|---|---|
| 0.99 | 1.008 | 686 | 0.80 |
| 0.95 | 1.039 | 28 | 0.80 |
| 0.90 | 1.082 | 8 | 0.83 |
| ≤ 0.70 | ≥ 1.30 | 3 (floor) | ≥ 0.99 |

CV 10% multiplies these by 4. **If CV is not known from a pilot, the study cannot be sized and must
not be preregistered with a number** — the pilot comes first.

### §3b. Hierarchical correction — purification is the unit

The table above assumes a single variance component. The correct model is paired within
purification with a purification-level random slope on the effect:

`Var(effect) = (tau² + 2σ²)/P`, where tau is the preparation-to-preparation SD of the *effect*.

At A₀ = 0.90: tau = 0.05 → P = 10 (1.4×), tau = 0.10 → P = 20 (2.9×), tau = 0.20 → P = 59 (8.4×).
Additional reactions per preparation do **not** reduce error once tau > 0.

**Both σ and tau must be reported from the pilot.** Sizing on σ alone understates the requirement
by up to ~8× at the operating point where the study is otherwise judged feasible.

---

## §4. Mandatory controls, each with a stated pass criterion

| control | why | pass criterion fixed in advance |
|---|---|---|
| **Active-fraction titration** | the A ≤ 0.70 gate is invalid if w < 0.8 | w reported with 95% CI; proceeds only if the CI excludes w < 0.8 |
| **Absorbing trap validation, 2 methods** | B is uninterpretable without irreversible loss | five readouts (basal addition, TPP1 K_d, active concentration, primer dissociation, s) each within ±10%, in **both** dilution and excess-DNA trap; heparin not accepted as a substitute |
| **TRAP-wedge/interface mutant** (H752A, H754A) | separates TPP1-driven from cycle-driven effects | mutant reports its own A and φ; no pooling with WT |
| **Trap-free paired arm** | detects whether the trap itself creates the effect | A ratio must agree in sign in both arms |

If the trap validation fails, **the B-channel conclusion is withheld** and only A is reported.
If w < 0.8, **the feasibility gate is not applied** and the study reports "preparation
characterisation incomplete" rather than a mechanism.

---

## §5. Decision rule — including every route to INCONCLUSIVE

Applied in order. The first matching row governs.

| # | condition | conclusion |
|---|---|---|
| 1 | w CI includes values < 0.8 | **INCONCLUSIVE — preparation**; A gate not interpretable |
| 1b | tau not estimable from a pilot | **INCONCLUSIVE — cannot size**; do not run |
| 2 | trap validation failed | report A only; **B withheld** |
| 3 | A_ctrl ≥ 0.90 | **INCONCLUSIVE — no power**; explicitly *not* "no effect" |
| 4 | A ratio CI excludes 1, \|ratio−1\| ≥ 0.05 | **extension-phase retention detected** (quantify fold-change) |
| 5 | A ratio CI includes 1, A_ctrl ≤ 0.70, w ≥ 0.8, CI narrow | **extension-phase retention excluded**; report the CI-derived upper bound |
| 6 | B lower bound > control B upper bound **and** trap validated | **translocation-phase effect detected** |
| 7 | φ difference within resolution, both < 0.95 | **translocation-phase effect excluded** |
| 8 | φ_trt > 0.95 with CI touching 1.0 | **INCONCLUSIVE on B**; report bound only |
| 9 | anything else | **INCONCLUSIVE** — do not classify |

**Rows 1, 2, 3, 8, 9 are required to return inconclusive.** A study that reports a mechanism while
any of them applies is out of compliance with this document.

---

## §6. What would falsify the design itself

Recorded so the design can fail honestly:

- **If the A-channel effect is undetectable at A_ctrl ≤ 0.70 with the hierarchical n**, the s/φ
  channel is not a usable discriminator for human telomerase and the programme should abandon it
  rather than escalate n.
- **If the trap cannot be validated within 10%** on all five readouts in both methods, the B channel
  is not measurable by this route and should be dropped, not reported as a bound.
- **If the two trap methods disagree beyond their CIs**, the inferred dissociation is
  method-dependent; report the disagreement as the result and claim neither value.
- **If paired and unpaired analyses disagree in sign**, the between-preparation variance is driving
  the result and neither is reportable.

---

## §7. Out of scope

This preregistration does not test, and makes no claim about: the catalytic-step (k_pol) allosteric
route; long-timescale conformational dynamics; cellular telomere maintenance; or the RPA
processivity pathway. Those require different designs.

---

## §8. The honest summary of what the mathematics bought

Fifteen work packages established that the discriminator is *sound* — the identity is exact, and
the invariance to within-population heterogeneity is now **proved**, not sampled — and that it is
*conditional*: it has power only in an operating regime nobody has measured, on preparations whose
active fraction nobody has reported, at a sample size that depends on a variance component nobody
has estimated.

Three pilot measurements therefore precede any mechanism test, in order:
1. **A on untreated enzyme** — one reaction; sets the operating point.
2. **w with a 95% CI** — sets whether the gate is interpretable at all.
3. **tau across ≥3 purifications** — sets whether the study can be sized.

**That is a result about experimental design, not about TPP1.** The mechanism remains unresolved,
and this document does not change that.
