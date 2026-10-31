# Kinetic model comparison — what can be asked, and what the algebra already decides

**No model was fitted to the published measurements.** This document compares five pre-registered
product-length models (`PREREGISTRATION.md` §6) in two ways that do not consume the holdout and do
not tune anything to the data: an algebraic identity (G1 ≡ G5), and an information-theoretic
calculation of how much data would be needed to separate the families if per-band data existed.

Reproduce: `python scripts/affinity_processivity_analysis.py` (≈2 s).
Outputs: `identifiability_results.csv`, `scalar_observable_limits.csv`,
`replicates_needed.csv`, `identifiability_collapse.json`.

---

## 1. The models

Let *q*(*n*) be the probability that a product has length *n* repeats, *S*(*k*) = Σ_{n≥k} *q*(*n*) the
survival function (so *S*(1) = 1 in every model), and *s* the per-repeat survival. The observable in
the published data is a single number, *P*(≥16) = *S*(16) = *s*¹⁵ for the constant-hazard model.

| ID | Model | Survival *S(n)* | Free parameters | Curvature it imposes on log *q*(*n*) |
|---|---|---|---|---|
| **G1** | constant-hazard retention (**affinity**) | *S*(*k*) = *s*^(*k*−1) | 1 | straight line (geometric) |
| **G2** | heterogeneous population (two species) | *w·s₁^(*k*−1) + (1−w)·s₂^(*k*−1) | 3 | concave (mixture of slopes) |
| **G3** | repeat-dependent extension | Π_{j=2..k} *s_j*, *s_j* = *s₀·ρ^(*j*−2) | 2 | convex (hazard grows with length) |
| **G4** | retention with a reversible pause | two-state (E ⇄ P, dissociation only from P) | 3 | excess of long products |
| **G5** | conformational (**k_ext** changes, **k_off** does not) | same form, with *k_ext*′ > *k_ext* | 2 | identical to G1 by construction |

In all five, *s* is a ratio of an extension rate to a dissociation rate. The models differ only in
*which* rate varies, and with what dependence on *n*.

## 2. The identity that settles half the question: G1 ≡ G5

G5 says TPP1 raises *k_ext*; G1 says TPP1 lowers *k_off*. Both enter the per-repeat survival
exclusively through

> *s* = *k_ext* / (*k_ext* + *k_off*).

Consequently, for **every** pair (*k_ext*, *k_off*) there is a pair (*k_ext*′, *k_off*) that yields
the same *s*, and therefore **the same product-length distribution and the same processivity,
exactly**. Verified numerically for four widely separated pairs (`identifiability_collapse.json`,
`"identical": true`, `n_distinct_P16 = 1`).

**Therefore: no measurement of the product-length distribution, at any resolution or depth, can
distinguish "affinity" from "conformational" as defined here.** The distinction is not
empirically hard; it is unobservable, because the two models use the same observable through the
same parameter.

Two things follow, and they are the most important results in this document:

1. **The product-length distribution cannot answer the mechanism question** as the models were
   originally framed. It can separate *constant hazard* from *heterogeneous* or
   *repeat-dependent* or *pausing* — real and useful distinctions about the enzyme's kinetic
   architecture — but it cannot separate the two causal stories labelled "affinity" and
   "conformational". The request that motivated this investigation (compare distributions to
   decide between them) is therefore answered in the negative by algebra, before any data are
   collected. **This is a genuine result and it is cheap.**
2. **A second, independent observable is required.** The two models make different commitments
   about *what else* moves when TPP1 binds: G1/G5-C lowers the apparent primer *K*m while leaving
   the maximum plateau unchanged; G5-B raises the plateau (*k_ext*) while leaving *K*m unchanged.
   That is the primer-titration signature (Panel 3 of the figure, and
   `decisive_experiment.md` §1). It is a different experiment, not a finer version of the same
   one.

## 3. How much per-band data would separate the *other* pairs

Using the Kullback–Leibler divergence between model-implied distributions, each constrained to
reproduce the observed *P*(≥16) = 0.49, evaluated at a densitometry error of 20% per band, a
quantification floor of 10⁻³, and restricted to the 50 bands a gel can resolve
(`replicates_needed.csv`):

| pair | informative bands | replicates needed (median) |
|---|---|---|
| G1/G5 vs G2 heterogeneous | 42 | <1 |
| G1/G5 vs G3 repeat-dependent | 31 | <1 |
| G1/G5 vs G4 reversible pause | 0.5 | **≈1.7** |
| G2 vs G3 | 43 | <1 |
| G2 vs G4 | 40 | <1 |
| G3 vs G4 | 26 | <1 |

**The "fewer than one replicate" entries are not a usable measurement.** They say that the two
models differ by far more than 20% in dozens of bands — which is true of the *specific member* of
each family chosen, and false in general, because the families are **nested**: G2 with *w* → 1
*is* G1, and G3 with *ρ* → 1 *is* G1. A comparison between the extreme members of two families is
easy; a comparison between a family and its own limit is impossible. Reporting the first number as
if it answered the question would be a category error of exactly the kind this investigation was
set up to avoid.

The one entry that is informative without qualification is **G1/G5 vs G4 (≈1.7 replicates)**:
pausing produces a shape that a constant hazard cannot mimic at *any* parameter setting, and it
needs about two good densitometry profiles to detect. That is a realistic, useful experiment.

The honest summary of the counting exercise is therefore: **distribution shape can separate
constant-hazard from the other architectures with one good gel; it cannot separate affinity from
conformational regulation at any depth** (§2).

## 4. The scalar-observable limit — why the published data cannot decide anything

`scalar_observable_limits.csv`: with one number per replicate, only a one-parameter model is
identified. G2 (3 parameters), G3 (2) and G4 (3) are all **exactly** consistent with the observed
scalar for some parameter setting — the KL divergence from G1 is 0.0 for all of them at the
observed value. **This is not a statement about how good the data are; it is a statement about how
many numbers they contain.** No amount of statistical sophistication recovers four parameters from
one observation per replicate.

## 5. What this changes about the original plan

The requested comparison — geometric vs mixture vs repeat-dependent vs heterogeneous, by held-out
predictive performance — is well posed and would be informative about **kinetic architecture**. It
is **not** a test of affinity vs allostery. Running it and concluding "affinity wins" would be the
exact error `PREREGISTRATION.md` §4.6 was written to prevent.

Consequently this work package delivers:

* the formal model statements (§1),
* the identity that rules out the originally proposed discriminator (§2),
* the data requirement for the discriminations that *are* possible (§3),
* and the identification of the observable that *is* discriminating (§2, second consequence).

It does **not** deliver a winning model, and it must not be cited as one.

## 6. Verdict

| question | answer |
|---|---|
| Can the product-length distribution separate affinity from conformational regulation? | **No — algebraically impossible** as the models are defined (G1 ≡ G5) |
| Can it separate constant-hazard from heterogeneous / repeat-dependent / pausing? | **Yes**, with one good densitometry profile, subject to the gel's dynamic range |
| Can the published data support *any* model comparison? | **No** — one scalar per replicate, 2–3 free parameters |
| Is a model winning? | **No model wins.** Nothing here ranks affinity above allostery or vice versa |
