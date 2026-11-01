# Mechanism classification — what processivity can and cannot identify

**Status:** analysis document. No mechanism is proven. This file records which mechanistic
questions telomerase processivity can answer, which it provably cannot, and which measurement
would settle each open case.

> ## Thesis
>
> **An identifiable effective kinetic parameter is not necessarily an identifiable molecular
> mechanism.** Mechanistic interpretation requires independent measurements that constrain the
> underlying physical rates.

Reproduce every number here: `python scripts/mechanism_classification_checks.py`
Guards: `tests/test_mechanism_classification.py` (25 tests) — the document is frozen by
`MECHANISM_CLASSIFICATION.sha256`, so editing it requires regenerating the hash deliberately.

**Committed on:** `arena/ab9d1ce3-telomerase-allosteric-switch`, based on `main` = `9e4cd39`.

---

## 0. What is established, and how strongly

| level | finding | status |
|---|---|---|
| **exact algebra** | Extension-rate enhancement (B) and affinity-mediated retention (C1) are **identical at any single primer concentration once *s* is matched** — the distribution is a function of *s* alone, so no processivity observable separates them there. Across a titration they differ only in the plateau, i.e. only through *k*_loss, the worst-determined parameter | proven (symbolic + numerical) |
| **exact algebra** | Every rate set is invariant under scaling θ → λθ: **no absolute kinetic rate is recoverable from processivity**, at any precision | proven (symbolic, simplifies to 0) |
| **exact algebra** | The C1a (*k*_off ↓) and C1b (*k*_on ↑) variants are **not** identical at matched *K*_d: the loss term is *R* = *k*_off²/(*k*_off + [P]·*k*_on), so a matched-*K*_d pair differs by exactly a factor λ — which is algebraically confounded with *k*_loss. Direct binding, *k*_off and *k*_on reported separately, is what separates them | proven (symbolic + numerical) |
| **exact algebra** | B vs C2 has **exactly one** crossing in primer concentration, linear in [P] | proven (degree 1) |
| **structural identifiability** | Processivity vs [P] has **rank 3 of 4**; direct binding is the only measurement that adds genuine independent rank | computed (complex-step Jacobian) |
| **simulation-derived** | The **tail-slope estimator** can manufacture an apparent concentration-independent loss term from zero true loss | simulated, **not yet validated on real bands** |
| **unresolved biology** | What microscopic action TPP1 actually changes | **no experiment yet distinguishes** |

The **tail-slope finding** (row 6) is exact as a statement about a mixture, but its *practical*
relevance is simulation-derived: it is promoted or retired only by fitting real per-band
densitometry (gate 3 below). The confidence tiers above are ordered —
exact algebra is not revisable by measurement, structural identifiability is not revisable by
measurement either but is conditional on the model class, and everything simulation-derived is
provisional until it meets real bands.

### 0.1 Corrections made to this document before publication

The repository's standard is that a retraction is published, not hidden. Three errors were found
in the first draft of this document by writing its own checks, and are recorded here:

| error | how it was caught | correction |
|---|---|---|
| claimed C1a and C1b are **identical** at matched *K*_d | the check that was supposed to verify it returned *false*; deriving *R* in closed form showed the two differ by an exact factor λ | §2 rewritten; the confound is with *k*_loss, not an identity |
| displayed *P*\* with the **opposite sign** in the numerator | comparing the displayed formula against `sympy.solve` and against a numeric scan: the displayed version gave *P*\* > 0 where no crossing exists | §3.1 sign corrected and the numeric confirmation table added |
| claimed the tail slope shifts *s* by ≈ φ(1−*s*), and by −0.040 at saturation | computing r(*n*) directly: the bias is 50–90× larger, is set by (*s*₂ − *s*), and at saturation is −0.0014, not −0.040 | §5 rewritten around the window-dependence result |

Each of the three now has a regression test.

---

## 1. Microscopic versus effective kinetics

The same word is used for four different things in the telomerase literature, including in this
repository's earlier write-ups. They must be kept apart.

| level | question | example |
|---|---|---|
| **molecular** | what physically happens to the primer? | dissociates and rebinds; enters a committed non-binding state; stays bound and is extended |
| **microscopic rate** | which physical rate constant changes? | *k*_off, *k*_on, *k*_commit, *k*_ext |
| **effective rate** | what does the assay measure? | *k*_loss (loss that is not recoverable by raising [P]), *R* = *k*_off·*K*_d/(*K*_d+[P]) |
| **observable** | what is actually recorded? | processivity distribution; primer-binding traces; single-turnover incorporation |

**Rule.** No claim in this repository may name a molecular event using an effective parameter.
*k*_loss in particular is an *effective* quantity: it means "escape that more primer cannot fix."
Whether it corresponds to a committed molecular state is a separate, experimental question (§4).

### 1.1 The effective-survival model

Per repeat, the enzyme either extends or escapes. With primer binding in equilibrium:

> **s = k_ext / (k_ext + k_loss + R),  R = k_off·K_d/(K_d + [P])**

The processivity observable is a function of *s* alone: P(≥ *n*) = *s*^(*n*−1) for a single
population.

### 1.2 Scale invariance (proven)

Substituting θ → λθ for θ = (*k*_ext, *k*_loss, *k*_off, *k*_on) leaves *s* unchanged. SymPy
returns **exactly 0** for `s(λθ) − s(θ)`, for all λ > 0.

**Consequence.** Processivity identifies only *ratios*. The identifiable combinations are

* *k*_loss/*k*_ext  (from the [P] → ∞ plateau)
* *k*_off/*k*_ext  and  *k*_on/*k*_ext  (from the shape of the [P] dependence)

and nothing else. **No absolute rate — not *k*_ext, not *k*_off, not *k*_on, not *k*_loss — is
recoverable from processivity at any precision.** This is why an analysis can look quantitatively
convincing and still establish no molecular rate.

---

## 2. C1a / C1b: affinity is not one hypothesis

*K*_d = *k*_off/*k*_on, so "TPP1 lowers the primer *K*_d" splits into three distinct molecular
claims. They are **not** interchangeable at the processivity level, and the earlier draft of this
document said they were. The correction follows from rewriting the loss term.

With *R* = *k*_off·*K*_d/(*K*_d + [P]):

> **R = k_off² / (k_off + [P]·k_on)**  (verified identical, symbolically)

Three consequences, each verified symbolically and numerically:

1. **[P] → 0:** *R* → *k*_off, independent of *K*_d and of *k*_on. At low primer the effective loss
   *is* the off-rate. This is the only place the *k*_off/*k*_on split is directly visible.
2. **Matched *K*_d, both rates scaled by λ:** *R* → λ·*R*. The curve keeps its shape in [P] and
   changes amplitude. So a *k*_off-only change and a *k*_on-only change at matched *K*_d differ by
   an **exact multiplicative factor**, not by nothing.
3. Because *s* = *k*_ext/(*k*_ext + *k*_loss + λ*R*), that factor λ is **algebraically indistinguish-
   able from a change in *k*_loss**. The confound is with *k*_loss specifically — the parameter the
   C2 question rests on, and the worst-determined one in §8.2.

| variant | microscopic change | *K*_d | effect on the loss term | processivity prediction |
|---|---|---|---|---|
| **C1a** | *k*_off ↓ | ↓ | *R* ∝ *k*_off at low [P]; whole curve scales | lower effective loss, all [P] |
| **C1b** | *k*_on ↑ | ↓ | only the *K*_d-saturating part moves | lower loss mainly at low [P] |
| **C1c** | both | ↓ | both effects | as above |

Worked example (wild type *k*_off = 2.0, *k*_on = 1.0, *k*_ext = 2.0, *k*_loss = 0.5). C1a
(*k*_off = 1.0) and C1b (*k*_on = 2.0) both give *K*_d = 1.0, yet *s* differs at every primer
concentration, converging only as [P] → ∞:

| | *K*_d | *s*([P]=0.1) | *s*([P]=1) | *s*([P]=100) |
|---|---|---|---|---|
| wild type | 2.0 | 0.454 | 0.522 | 0.788 |
| C1a (*k*_off ↓) | 1.0 | **0.587** | **0.667** | 0.797 |
| C1b (*k*_on ↑) | 1.0 | **0.463** | **0.571** | 0.794 |

**Therefore:** "TPP1 lowers the primer *K*_d" is under-specified. Processivity *can* in principle
separate C1a from C1b, but only through an amplitude change that trades off against *k*_loss*
(§8.2: ~49% relative SE even with direct binding). Distinguishing them therefore needs **direct
primer-binding kinetics that report *k*_off and *k*_on separately**, not *K*_d alone — and the
low-primer asymptote is where the separation is cleanest, and where processivity is measured worst
(short products, few extension events). That tension is a design constraint, not a footnote.

---

## 3. The B–C2 crossing, and its parameter domain

### 3.1 Derivation

With α = *k*_ext multiplier and β = *k*_loss divisor:

> *s*_B  = α*k*_ext / (α*k*_ext + *k*_loss + *R*)
> *s*_C2 = *k*_ext / (*k*_ext + *k*_loss/β + *R*)

Cross-multiplying, the α*k*_ext·*k*_ext terms cancel, leaving an equation **linear in *R***:

> *R*·(α − 1) = *k*_loss·(1 − α/β)

Substituting *R* = *k*_off·*K*_d/(*K*_d + [P]) gives a result **linear in [P]**:

> **[P]\* = *k*_off·(−αβ·*k*_off − α·*k*_L + β·*k*_L + β·*k*_off) / [ *k*_L·*k*_on·(α − β) ]**

### 3.2 Three cases, all proven symbolically

| case | condition | meaning |
|---|---|---|
| **positive crossing** | *P*\* > 0 (equivalently 0 < *R*\* < *k*_off) | B and C2 swap dominance at an accessible [P] |
| **negative crossing** | *P*\* < 0, i.e. *R*\* > *k*_off | **no** crossing at any physical [P]: one mechanism dominates at every primer concentration |
| **degenerate** | α = β | numerator reduces to α·*k*_off²(α−1), independent of [P] — crossing only at [P] → ∞ |

**Linearity ⇒ uniqueness.** At most one crossing exists for all α, β, *k*_L > 0, α ≠ β. This is a
proof, not a scan: a linear equation has one root.

**Numerical confirmation of the closed form** (*k*_ext = 2, *k*_L = 0.5, *k*_off = 2, *k*_on = 1),
scanning [P] over 10⁻³–10⁴ with 20 001 points:

| α | β | *P*\* | *s*_B − *s*_C2 at *P*\* | sign changes in [P] |
|---|---|---|---|---|
| 0.5 | 1.0 | −10.00 | n/a | **0** |
| 2.0 | 1.0 | −10.00 | n/a | **0** |
| 30.0 | 3.42 | −31.85 | n/a | **0** |
| 2.0 | 3.42 | **17.27** | 1.1×10⁻¹⁶ | **1** |
| 2.0 | 10.0 | **8.00** | 2.2×10⁻¹⁶ | **1** |
| 5.0 | 10.0 | **62.00** | 1.1×10⁻¹⁶ | **1** |

The closed form satisfies *s*_B = *s*_C2 to machine precision exactly when *P*\* > 0, and the number
of sign changes equals 1 if and only if *P*\* > 0. This table is a regression test, not decoration:
an earlier draft of this document displayed the numerator with the opposite sign, which turns every
positive case into a negative one. The script fails if that happens again.

### 3.3 What the domain caveat means

A negative [P]\* is not a mechanism switch. It says one model dominates the whole accessible
range — a *different* result from "they cross," and it must be reported as such. The sign of
*P*\* is therefore a **precondition** of any crossing statement, and must be computed before one is
made.

The earlier claim in this project that B and C2 are "indistinguishable" was **wrong** — but only
*when P*\* > 0*. In that case they coincide at exactly one concentration and differ everywhere else,
so they are **distinguishable in principle, though the size of the difference is small** (§8.3).
When *P*\* < 0 there is no coincidence point in the physical domain at all, and the correct report
is "one mechanism dominates throughout," which is a weaker statement than either "they cross" or
"they are indistinguishable." Conflating the three cases is the error this section exists to
prevent.

---

## 4. The concentration-independent effective loss term

### 4.1 What the [P] → ∞ limit establishes

> *s*_∞ = *k*_ext / (*k*_ext + *k*_loss)  ⇒  *k*_loss/*k*_ext = (1 − *s*_∞)/*s*_∞

A plateau below 1 therefore establishes **a loss channel that more primer cannot rescue**. Since
*s*_∞ = 1 would mean infinite processivity and every measured telomerase processivity is finite,
**such a term necessarily exists in this model** — that is a consequence of the equations plus one
uncontroversial observation, not a finding. Its *existence* is therefore not the empirical
question.

The empirical questions are: is the term **reproducible**, how **large** is it, and does it have
an **independently measurable molecular identity**.

### 4.2 Renaming

This is the **concentration-independent effective loss test**. It is *not* a committed-loss test.
A plateau below 1 is equally consistent with:

| alternative | where it sits in the model | distribution signature |
|---|---|---|
| concentration-independent dissociation | *k*_loss | raises plateau; no curvature |
| unmodelled inactive enzyme fraction φ | separate multiplicative parameter | scales intercept; **slope unchanged** |
| competing nonproductive conformational state | φ or *k*_loss, state-dependent | intercept and/or slope |
| heterogeneous primer-bound complexes | mixture of *s* | **curvature** in log *q*(*n*) |

Claiming a committed molecular intermediate requires independent evidence of that state or its
transition kinetics. Saturation data cannot supply it.

### 4.3 Whether C2 is a third mechanism at all

*k*_loss is effective. If irreversible loss simply *is* dissociation that fails to rebind, then
**C2 is a coarse-grained restatement of C1 plus local geometry and [P]**, not an independent
mechanism — and the three-model problem is a two-model problem with two parameterisations.

C2 survives as a distinct hypothesis only if irreversible loss is a real committed step. That is
an open question, and it is the first thing to settle, because it determines how many models the
programme actually has.

---

## 5. Operational saturation, and why the tail slope is forbidden

### 5.1 Saturation, defined on the response rather than on occupancy

Not "10 × *K*_d". Saturation is established when:

1. **three consecutive doublings of [P]** each change the log fold-stimulation by less than a
   prespecified margin δ;
2. **TOST equivalence** holds — the 90% CI of each change lies entirely inside ±δ (CV alone is not
   an equivalence margin, and "flat within the CV" is also satisfied by a dead assay);
3. a **sensitivity control** run at the plateau concentration still resolves δ (spike-in of a
   known processivity difference);
4. the same plateau is reached **from both directions** of approach.

**Declared failure mode:** if no plateau is bracketed with the CI inside ±δ within the accessible
[P] range, the result is *"saturation not established"* and the C2 test is **inapplicable** — not
extrapolated.

### 5.2 The tail slope is retired as a basis for any mechanistic claim

For a two-species mixture, *q*(*n*) ∝ (1−φ)*s*^(*n*−1) + φ(1−*s*₂)*s*₂^(*n*−1), the adjacent-band
ratio is closed-form:

> **r(*n*) = *s* · [1 + κρ^n] / [1 + κρ^(*n*−1)],  κ = φ(1−*s*₂)/[(1−φ)(1−*s*)],  ρ = *s*₂/*s***

so the tail slope converges to **max(*s*, *s*₂)** — the *longest-lived* species, not the majority.
Measured over the long window *n* = 40–80, where convergence is essentially complete:

| true majority *s* | φ | *s*₂ | tail slope | apparent Δ*s* | φ(1−*s*) would predict |
|---|---|---|---|---|---|
| 0.60 | 0.01 | 0.95 | **0.9500** | **+0.350** | 0.004 |
| 0.70 | 0.01 | 0.90 | **0.8997** | **+0.200** | 0.003 |
| 0.85 | 0.05 | 0.95 | **0.9382** | **+0.088** | 0.0075 |
| 0.999 | 0.01 | 0.95 | 0.9976 | −0.0014 | 0.00001 |

Four consequences:

* the bias is **50–90× larger** than the naive φ(1−*s*) estimate, and is set by (*s*₂ − *s*), not by φ;
* at saturation (*s* → 1) the bias essentially **vanishes** (−0.0014 at *s* = 0.999). The earlier draft
  of this table asserted −0.040 here; that number is not reproducible at any window. The correct
  saturation statement is that the *finite-window* reading is biased slightly **downward**, because a
  slower contaminant has not yet decayed out — see the window sweep below;
* the bias **depends on which length window is read**, because r(*n*) varies with *n* by construction;
* therefore the tail slope is **not a function of the distribution at all** — it is a function of the
  distribution *and the length window*, and no experiment can report "the" tail slope without
  reporting both.

**The window effect is the decisive objection.** Same underlying truth, five read windows:

| true *s* | φ | *s*₂ | *n*=3–10 | 10–16 | 16–25 | 25–40 | 40–80 | spread |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 0.01 | 0.95 | 0.610 | 0.707 | 0.909 | 0.950 | 0.950 | **0.340** |
| 0.70 | 0.01 | 0.90 | 0.703 | 0.716 | 0.771 | 0.875 | 0.900 | **0.196** |
| 0.85 | 0.05 | 0.95 | 0.853 | 0.857 | 0.864 | 0.889 | 0.938 | **0.085** |
| 0.999 | 0.01 | 0.95 | 0.986 | 0.989 | 0.991 | 0.994 | 0.998 | **0.012** |

For a majority processivity of *s* = 0.60, reading the near window (*n* = 3–10) returns 0.610 and
reading the far window returns 0.950. That is an implied *k*_loss/*k*_ext of 0.639 versus 0.053 —
a **12-fold** error in the quantity a C2 claim is made of, from the same data. A C2 claim is a
claim about the residual loss after saturation, and that is exactly the quantity this estimator
cannot pin down.

**Fabricated C2 from zero true loss:** with *k*_loss = 0 (true plateau 1.0), the tail-slope
estimator reads 0.998 and infers *k*_loss/*k*_ext ≈ **0.017–0.036** from 0.5–2% contamination,
while an earlier near-window reading of the same distributions implied **0.05–0.10**.

A full-distribution fit does not make this error, because φ is fitted as a separate parameter
(recovered to ±0.02 in simulation). **Therefore: every mechanistic conclusion in this project must
come from a full-distribution fit with φ free, never from a tail slope.**

### 5.2a Two further exact results about the mixture ratio

Writing *u* = κρ^(n−1), the ratio r(*n*) = *s*(1+*u*ρ)/(1+*u), so

> **dr/du = s(ρ−1)/(1+u)²**   (verified symbolically)

Because *u* is monotone in *n*, r is monotone in *n* and therefore has no interior extremum. Two
consequences follow, and both bound how strongly the bias result may be stated:

**Interpolation bound.** For every *n*, φ and *s*₂,

> **min(s, s₂) ≤ r(n) ≤ max(s, s₂)**

so the adjacent-band ratio can never report a processivity beyond the most extreme subpopulation
actually present. (Worst relative violation 1.0×10⁻¹³ over 5 000 random parameter sets — numerical
agreement with exactness.) **r(n) → max(s, s₂)** as *n* → ∞, which is the limit stated at the top of
this subsection.

**Onset length.** The two components of *q*(*n*) are equal, (1−φ)(1−*s*)*s*^(n−1) =
φ(1−*s*₂)*s*₂^(n−1), at

> **n\* = 1 + ln(1/κ)/ln(ρ)**

Below *n*\* the majority controls the per-band signal; above it, the minority does. Checked against
the numerically located midpoint of r(*n*) in eight cases: the closed form lands within one band
every time (e.g. *s* = 0.60, φ = 0.01, *s*₂ = 0.95 → *n*\* = 15.52, midpoint at *n* = 15).

**Why these two matter operationally.** The bound says the damage is finite — a real assay cannot
report a processivity no molecule has, so the window artefact cannot exceed |*s*₂ − *s*|. The onset
length says *where* the artefact begins, and it is computable from quantities an experimenter can
estimate. If *n*\* lies inside the readable band range, the tail slope is measuring the minority; if
it lies beyond the last readable band, the tail slope is measuring the majority and the artefact is
harmless. That is a testable, dataset-specific question, and it is the one gate 3 asks.

### 5.3 Equivalence margin δ (open decision)

Paired TOST design, 80% power, two-sided, Bonferroni over three doublings:

| δ | CV 0.345 (published) | CV 0.15 (optimised assay) |
|---|---|---|
| 0.10 log | 74 (1 test) / **105** (3 tests) | 14 / **20** |
| 0.20 log | 19 / 27 | 4 / 5 |

δ must be set from the **smallest biologically meaningful processivity difference**, not from
convenience. Calibration from the existing mutant panel: adjacent published processivity ratios
differ by 0.05–0.38 log, so δ = 0.20 would fail to separate neighbouring mutants that the field
already treats as distinct. **δ remains undecided and is a decision for the experiment owner.**

---

## 6. Chemistry versus repeat-cycle kinetics

Single-turnover incorporation does **not** automatically report *k*_ext as defined in §1.1.

| quantity | measurement | rests on |
|---|---|---|
| *k*_chem | single-turnover, pre-formed complex, saturating dNTP, quench ≪ one cycle | chemistry is irreversible |
| *k*_rep | product-length advance per unit time at steady state | rate-limiting step constant across repeats |
| *k*_transloc | difference between the two — **not a direct observable** | only two steps separate them |

**Rank consequence.** Adding single-turnover chemistry to the design **adds zero identifiability
unless one assumes the measured step is the model's extension step.** With a free coupling
parameter η between them, the Jacobian gains one observable and one parameter and the rank gain is
zero. The "third independent scale" is an assumption, and it needs its own test: a dNTP-chase
experiment confirming that the single-turnover step is the processivity-rate-limiting step.

Any residual discrepancy after that is **unassigned**, not "translocation".

---

## 7. Alternative explanations for each observable

| observable | what it is usually taken to mean | alternatives that produce the same signal |
|---|---|---|
| P(≥16) fraction alone | processivity of the population | scaled by inactive fraction φ; a 1% long-lived subpopulation raises it more than a uniform improvement |
| tail slope | the population's survival probability | reports max(*s*, *s*₂); biased by length-dependent detection exp(−*n*/L); biased by any long-lived contaminant |
| plateau of *s* vs [P] | concentration-independent loss (C2) | inactive fraction; slow contaminant; competing nonproductive state |
| TPP1 stimulation persisting at high [P] | a signal, not affinity | a genuine committed state (C2) — same observation, different mechanism |
| dNTP × TPP1 interaction | TPP1 acts on extension | nucleotide-dependent pausing with TPP1 acting elsewhere (simulated interaction 1.62 with no extension change) |
| single-turnover rate change | TPP1 changes chemistry | a changed translocation step or nucleotide-dependent pause, if the coupling assumption fails |

### 7.1 The dNTP × TPP1 experiment is rejected

Recorded because it was proposed and then falsified as a design, in this repository:

* **not mechanism-specific:** a nucleotide-dependent pausing model with TPP1 acting on a different
  step reproduces an interaction (1.62) without any change in extension efficiency;
* **infeasible in its informative regime:** the interaction is large only below ~0.75 × *K*_m, where
  P(≥16) = 2.4×10⁻⁹ at 0.1 × *K*_m — four orders of magnitude below a gel's quantification floor;
* inside the experimentally usable window (0.75–100 × *K*_m) the interaction is **0.84–0.85**, a
  10–15% change, below the 15–35% CV of the published assay.

---

## 8. Identifiability, scale degeneracies, and release gates

### 8.1 Structural identifiability (Jacobian rank, complex-step derivatives)

| design | parameters | rank | smallest relative σ |
|---|---|---|---|
| processivity vs [P] | 4 | **3 / 4** | 5.8×10⁻¹⁷ |
| + single-turnover, **assuming** *k*_ext ≡ *k*_chem | 4 | 4 / 4 | 0.087 |
| + single-turnover, **free** coupling η | 5 | **4 / 5** | 5.8×10⁻¹⁷ |
| + direct *k*_on, *k*_off | 4 | **4 / 4** | 0.20 |

Rank is grid-independent (5.8–1.1×10⁻¹⁶ across three primer grids): the deficiency is
**structural**, not a sampling artefact.

**Conclusion: direct primer-binding kinetics are the strongest independent measurement in the
programme.** They are the only one that adds genuine rank.

### 8.2 Practical identifiability (full rank is not sufficient)

All condition numbers below are from a **log-space** Fisher matrix with a 15% assay CV on *s* at 12
primer concentrations; condition numbers from different parametrisations are not comparable, so
the parametrisation is stated in every row.

**Ratio form** — 3 parameters, *k*_ext absorbed:

| ratio | estimate | relative SE |
|---|---|---|
| *k*_loss/*k*_ext | 0.250 | 44% |
| *k*_off/*k*_ext | 2.000 | 12% |
| *k*_on/*k*_ext | 0.500 | 60% |

Condition number **60.3**.

**Absolute four-rate form** — *k*_ext, *k*_loss, *k*_off, *k*_on:

| design | condition number | *k*_ext | *k*_loss | *k*_off | *k*_on |
|---|---|---|---|---|---|
| processivity vs [P] | 2.4×10⁸ | 31% | 34% | 29% | 67% |
| + direct binding at 10% CV | **7.5** | 22% | **49%** | 9.9% | 9.9% |

Adding direct binding fixes *k*_off and *k*_on to the precision it is measured at, and collapses the
condition number by seven orders of magnitude. It does **not** fix *k*_loss: its relative SE goes
**34.3% → 49.0%**, i.e. it gets *worse*, because binding removes the correlations that were
previously absorbing part of that variance. *k*_loss is the worst-determined quantity in every
design tested, and it is precisely the parameter the C2 question rests on. Direct binding is
therefore necessary but **not sufficient** for the C2 claim; nothing in this document may be read
as saying the binding experiment alone settles it.

### 8.3 Non-identifiability conditions, declared in advance

| condition | consequence |
|---|---|
| matched *s*, B vs C1 | **exactly non-identifiable by any processivity observable** |
| matched *K*_d, C1a vs C1b | differ by a factor λ on the loss term, which is confounded with *k*_loss (49% SE); needs *k*_off and *k*_on reported separately |
| B vs C2 within ±(2 × CV) of *s* | **non-identifiable in practice** (differences of ~0.10 in *s* need CV ≤ 0.05) |
| free coupling η, single-turnover only | rank-deficient; the scale is an assumption, not a measurement |
| φ poorly determined in the full-distribution fit | the assay lacks low-*n* dynamic range; C2 conclusions are unsupported |

If a measurement cannot separate the remaining models, the deliverable is a **non-identifiability
report naming the measurement that would**, not a selected winner.

### 8.4 Release gates

| gate | status | remaining action |
|---|---|---|
| 1. Symbolic verification (invariance, crossing, domains) | **satisfied** | independent third-party reproduction |
| 2. Structural identifiability (Jacobian rank) | **satisfied** | distinguish structural from practical (done, §8.2) |
| 3. Observation-model validation | **conditional** | real per-band densitometry: does a full-distribution fit on existing bands show φ > 0 and a majority *s* below the tail slope? |
| 4. Kinetic mapping (*k*_chem ≡ *k*_ext?) | **experimental** | dNTP-chase test; no analysis can substitute |
| 5. Equivalence design (δ, power) | **pending decision** | biological judgement on the smallest meaningful effect |

**Gate 3 is the next real advance.** It requires no new experiment — only densitometry of an
existing gel — and it can fail in both directions: if φ ≈ 0, Finding 3 has no practical importance
for that dataset and should be reported as such.

---

## 9. What this document does not claim

* It does not claim TPP1 changes primer retention, extension chemistry, or a committed-loss state.
* It does not claim a mechanical or allosteric component exists — or is absent.
* It does not claim C2 is a distinct molecular mechanism; that is unresolved (§4.3).
* It does not promote affinity/anchoring above conformational regulation; **within the stated
  model they are the same function of the observable** (§0), which is a statement about
  identifiability, not about biology.

The one substantive claim is negative and general: **processivity distributions alone cannot
identify TPP1's mechanism, and the framework above specifies which independent measurements can.**
