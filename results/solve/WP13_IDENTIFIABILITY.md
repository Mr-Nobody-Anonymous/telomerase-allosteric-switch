# WP13 — Experimental identifiability and assay operating range

Status: **complete**. Mechanism: **unresolved**, and now with a measured reason.

Script: `wp13_identifiability.py` → `wp13_sensitivity_map.csv`, `wp13_sample_sizes.csv`,
`wp13_intervals.csv`.

The question WP13 answers is not "which mechanism" but **"can the proposed experiment distinguish
them at realistic human telomerase rates?"** The answer is: only inside an operating regime that
nobody has measured, and the first experiment should measure that regime rather than the mechanism.

---

## 13.1 Sensitivity map

Channel A responds only to `a = k_off_ext / k_pol`. Nothing else enters it. Selected rows
(reported value = A(treated)/A(control) after the stated fold-drop in `k_off_ext`):

| k_pol | k_off_ext | a | A | R (f=2) | R (f=4) | R (f=10) |
|---|---|---|---|---|---|---|
| 0.02 | 0.012 | 0.600 | 0.060 | 3.48 | 7.25 | 11.8 |
| 0.05 | 0.003 | 0.060 | 0.705 | 1.19 | 1.30 | 1.37 |
| 0.10 | 0.003 | 0.030 | 0.838 | 1.09 | 1.14 | 1.17 |
| 0.10 | 0.012 | 0.120 | 0.507 | 1.39 | 1.65 | 1.84 |
| 0.20 | 0.003 | 0.015 | 0.915 | 1.05 | 1.07 | 1.08 |
| 0.50 | 0.003 | 0.006 | 0.965 | 1.018 | 1.027 | 1.033 |
| 2.00 | 0.003 | 0.0015 | 0.991 | 1.005 | 1.007 | 1.008 |
| 2.00 | 0.012 | 0.006 | 0.965 | 1.018 | 1.027 | 1.033 |

**Because A depends only on `a`, one measured value of A places the assay on this map without any
separate knowledge of k_pol.** That is what makes the regime measurement cheap.

## 13.2 Sample sizes

| A | R (f=4) | n (WP13 z-form) | n (exact t, WP14) | power at WP13 n | verdict |
|---|---|---|---|---|---|
| 0.990 | 1.0076 | 686 | 686 | 0.80 | infeasible |
| 0.950 | 1.0392 | 26 | **28** | 0.78 | hard |
| 0.900 | 1.0820 | 6 | **8** | **0.69** | feasible |
| 0.700 | 1.3041 | 1 | **3** (floor) | n/a | feasible |
| 0.500 | 1.6690 | 1 | **2** | n/a | feasible |

> **⚠️ Corrected by WP14.1.** The z-formula understates n at moderate effect sizes: at A = 0.90 it
> gives n = 6, which reaches only 69% power, not 80%. It also returns n = 1 at A ≤ 0.70, which is
> not a runnable design — 3 is the floor for estimating variance. Use the exact-t column.
>
> At A = 0.99 the study is not merely insensitive — it is infeasible ("infeasible" here means
> ~700 independent preparations; ~2700 at CV 10%).

## 13.3 Intervals — and the correction you asked for

Exact (Clopper–Pearson) intervals, and why the Wald interval must not be used near φ = 1:

| k/n | φ̂ | Wald 95% | Clopper–Pearson | B point | B lower bound |
|---|---|---|---|---|---|
| 27/30 | 0.900 | [0.793, 1.007] | [0.735, 0.979] | 9.00 | 2.77 |
| 29/30 | 0.967 | [0.902, 1.031] | [0.828, 0.999] | 29.0 | 4.81 |
| 30/30 | 1.000 | [1.000, 1.000] | [0.884, 1.000] | ∞ | **7.64** |
| 99/100 | 0.990 | [0.970, 1.010] | [0.946, 0.9997] | 99.0 | 17.4 |
| 992/1000 | 0.992 | [0.986, 0.998] | [0.984, 0.997] | 124 | 62.7 |

At 30/30 the Wald interval asserts **infinite precision** and B = ∞. Clopper–Pearson gives
[0.884, 1.000], so the honest statement is **B ≥ 7.6**. Reporting a point estimate of φ near 1 is
reporting an artefact of the interval, not the data.

**Rule adopted:** report B as a lower bound from the Clopper–Pearson lower limit of φ whenever
φ̂ > 0.95 and n < 200.

## 13.4 Decision rule, with an explicit INCONCLUSIVE branch

| condition | verdict | note |
|---|---|---|
| A ratio CI excludes 1 **and** \|ratio−1\| ≥ 0.05 | extension-phase retention **detected** | quantify fold-change in k_off_ext; do not call it k_pol |
| A ratio CI includes 1 **and** A_ctrl ≤ 0.70 | extension-phase retention **excluded** | report the upper bound implied by the CI edge |
| A ratio CI includes 1 **and** A_ctrl > 0.90 | **inconclusive** on A | the assay has no power here; report the experiment as uninformative |
| B lower bound exceeds control's B upper bound | translocation-phase effect **detected** | direction ambiguous: k_tr↑ or k_off_tr↓ |
| φ difference CI within ±resolution, both < 0.95 | translocation-phase effect **excluded** | valid only below saturation |
| φ_trt > 0.95 with CI touching 1.0 | **inconclusive** on B | report B as a bound only |

## 13.5 The first experiment measures the operating regime

Measure **A in the untreated condition, with no TPP1 at all** — one reaction, `s` from a full
product-length distribution and `φ` from single-round translocation in the same material.

- `A_ctrl ≤ 0.70` **and active fraction w ≥ 0.8** → channel A has power; proceed to the ±TPP1
  comparison. (The w ≥ 0.8 condition was added by WP14: a 50% inert preparation reads A = 0.815 and
  would be wrongly rejected while containing 50% fully responsive enzyme.)
- `A_ctrl ≥ 0.90` → channel A has no power; move to channel B or to a different readout.

This is a **pre-registration input**, not a result: it decides which study is worth doing, and it
costs one reaction.

### Correction I owe you on resting vs extension-phase k_off

You were right to separate them, and it matters more than "different states":

| k_pol | A using resting 0.012 | A using ext-phase 0.003 | R(f=4) resting | R(f=4) ext |
|---|---|---|---|---|
| 0.05 | 0.275 | 0.705 | 2.56 | 1.30 |
| 0.10 | 0.507 | 0.838 | 1.65 | 1.14 |
| 0.50 | 0.867 | 0.965 | 1.11 | 1.03 |
| 2.00 | 0.965 | 0.991 | 1.03 | 1.01 |

At k_pol = 0.1 the resting value predicts a 1.65× effect where the extension-phase value gives
1.14×. A design sized on the resting number is underpowered by roughly an order of magnitude in
replicates. **Resting k_off is a baseline and a control; it must not be substituted for
k_off_ext in a power calculation.**

## 13.6 Can the anatomy be validated? Partly — and one assumption is load-bearing

Extending the scheme with (i) dissociation from the closed complex before chemistry and (ii) a
dead-end realignment state (the P2 pause signal), 24 rate parameters, symbolic solve:

```
∂A/∂k_on  = 0                 ∂A/∂k_tr   = 0        ∂A/∂k_off_tr = 0
∂A/∂k✳_pol = 0                 ∂A/∂k_cl   = 0
∂B/∂k_pol = 0                  ∂B/∂k✳_pol = 0        ∂B/∂k_on = 0
```

The two-channel anatomy **is exactly invariant** to rate constants on the other side of each
channel — absent cooperative branch points. A channel-commitment test (extend vs translocate from
the same post-repeat state) shows commitment fractions sum to 1 exactly, so the anatomy is *a
decomposition, not a kinetic parameterisation*: it does not assume independence, it survives its
violation on the main branch. The invariance breaks only when a rate is shared across both phases
(the documented dG1 pause signal is exactly such a term).

**But one assumption is not free and the assay must enforce it: dissociation must be absorbing.**

| rebinding (k_on·[P]) | s (1 repeat) | s (5 repeats) | A |
|---|---|---|---|
| 0.000 | 0.8031 | 0.3340 | 0.8375 |
| 0.010 | 0.8293 | 0.3922 | 0.8375 |
| 0.100 | 0.8364 | 0.4094 | 0.8375 |
| 1.000 | 0.8374 | 0.4117 | 0.8375 |

Rebinding inflates the apparent translocation term without touching A. So the channels stay
separate, but **B becomes an upper bound** on k_tr/k_off_tr. The fix is procedural, not
statistical: pulse-chase, dilution, or a DNA trap, so loss is irreversible. Latrick & Cech used
exactly this design.

Three further assumptions, none free:
1. phases are strictly sequential (no extension concurrent with translocation);
2. TPP1 does not act on rebinding — otherwise k_off is not the whole story, which is precisely the
   k_off vs k_on distinction the manuscript insists on;
3. every repeat has the same (k_pol, k_off, k_tr) — otherwise A and B are population averages.

These belong in the preregistration, not in a results section.

---

## Milestone

| item | status |
|---|---|
| WP11–WP12 audit | complete |
| WP13 identifiability and operating range | complete, **sample sizes corrected by WP14** |
| mechanism | **unresolved** |
| next experimental step | measure A in the untreated condition (one reaction, no TPP1) |
