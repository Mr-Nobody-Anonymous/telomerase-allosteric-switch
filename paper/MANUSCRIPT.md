# What Processivity Distributions Cannot Tell Us About TPP1

**An identifiability analysis of telomerase repeat-addition processivity, with a corrected treatment of heterogeneous populations**

**Beamlak Sintayehu**

2026-11-01

---

## Abstract

TPP1 stimulates the repeat-addition processivity (RAP) of human telomerase, and the effect requires the TERT TEN domain, the IFD–TRAP motif, and the TPP1 TEL patch. What TPP1 changes *mechanically* is unresolved. Here we ask a narrower and answerable question: what can a processivity measurement — a product-length distribution — identify about that mechanism?

We show that it identifies less than is commonly assumed, and we state exactly where the boundary lies. Three results are exact within a stated kinetic model class. First, an extension-rate mechanism (B) and an affinity-mediated-retention mechanism (C1) are **the same function of the observable**: they act only through a single effective survival probability *s*, so under matched *s* they produce identical distributions, and a distribution measured at one primer concentration carries no mechanistic information beyond *s*. Second, no absolute kinetic rate is recoverable from processivity at any precision: every observable is invariant under scaling all rates by a common factor. Third, a processivity design yields a structural Jacobian of rank 3 for 4 parameters; direct primer-binding kinetics is the only independent measurement that adds rank.

We also correct a claim that we and others have made: the two ways of lowering the primer dissociation constant — reducing *k*_off (C1a) or raising *k*_on (C1b) — are **not** equivalent at matched *K*_d. The concentration-dependent loss term rewrites exactly as *R* = *k*_off²/(*k*_off + [P]*k*_on), so a matched-*K*_d pair differs by an exact factor λ, and that factor is algebraically confounded with the concentration-independent loss *k*_loss. Reporting *K*_d alone is therefore under-specified.

Finally, we treat heterogeneous populations, the main practical obstacle. For a two-species mixture the adjacent-band ratio is *r*(*n*) = *s*(1+κρⁿ)/(1+κρⁿ⁻¹) with κ = φ(1−*s*₂)/[(1−φ)(1−*s*)]. Three consequences follow exactly: *r* is monotone in *n* and always bounded by min(*s*,*s*₂) and max(*s*,*s*₂); the tail slope converges to the longest-lived species, not the majority, so the bias is positive for a more processive contaminant and negative for a less processive one; and the minority takes over the signal at a computable length *n*\* = 1 + ln(1/κ)/ln(ρ). Because *n*\* can fall inside a readable gel range, the *same* underlying distribution yields an apparent *s* anywhere from 0.61 to 0.95 depending only on which band window is read — a 12-fold change in the implied *k*_loss/*k*_ext. We therefore forbid the tail slope as a basis for any mechanistic claim and require a full-distribution fit with the mixture weight free. An earlier approximation, Δ*s* ≈ φ(1−*s*), is retired: it underestimates the bias by 50–90× and gets the sign wrong.

These results are mathematical, not biological. Within the stated model class the derivations are exact; whether they matter for a given experiment is an empirical question that requires per-band densitometry of an existing gel, not another simulation. We specify that test, together with the equivalence design and the release gates a mechanistic claim would have to pass. **No mechanism is demonstrated here.** The contribution is a boundary: it tells you what evidence would not settle the question, before that evidence is collected.

---

## 1. The question, and why it is narrower than it looks

Human telomerase extends a chromosome end by one telomeric repeat, then must realign the DNA product on its RNA template and extend again. The fraction of enzyme that completes a further repeat rather than releasing the DNA is its repeat-addition processivity. POT1–TPP1 raises it *in vitro* [1], and the stimulation requires contacts among the TERT TEN domain, the IFD–TRAP motif and the TPP1 TEL patch [2,3,4]. Structures of the telomerase–TPP1 complex place TPP1 at the TEN–IFD–TRAP surface and do not show TPP1 contacting DNA [3,4]. Recently, RPA has been identified as a second, independent processivity factor acting at a distinct surface [5].

The field-level question — does TPP1 act through a local structural conduit, broad conformational allostery, or affinity/anchoring — has been difficult to settle because the available perturbations do not separate these possibilities. Kinetic work has argued for reduced primer dissociation and aided translocation [6], and the structural work is consistent with stabilisation of the DNA anchor site [3,4]. Neither establishes that a *mechanical* signal exists, as distinct from a change in how well the enzyme holds its primer.

We therefore replace the mechanism question with one a measurement can address:

> **Does TPP1 change how tightly telomerase retains the primer, how quickly it extends the primer, or both — and can a processivity distribution distinguish these?**

The answer we reach is: the third part is no, under conditions we make explicit. That is a negative result about inference, and it is the point of the paper. A negative result about identifiability is useful only if it is precise, so we state the model class, prove the results within it, and then state exactly what would be needed to go beyond it (§8).

**A distinction we will keep throughout.** An exact result within a model is not an experimentally established biological fact. Every theorem below is conditional on a kinetic scheme; every number is reported with the assumption that produced it. Where the two could be confused, we say which one we are asserting.

---

## 2. Model class and notation

We use the smallest model that contains both mechanisms under discussion. Per repeat, the enzyme either extends the primer or loses it, and primer binding is treated as an equilibrium:

> **s = k_ext / (k_ext + k_loss + R),  R = k_off·K_d/(K_d + [P])**   (1)

with *K*_d = *k*_off/*k*_on, and the observable *P*(≥ *n*) = *s*^(n−1) for a single population. *R* is the primer-dependent loss: at low [P] enzyme is mostly unbound, at high [P] binding is saturated. *k*_loss is the residual loss that more primer cannot rescue — an **effective** quantity, not a molecular state.

Four comments on the model, because they bound every claim that follows.

1. It is a *phenomenological* per-repeat model. It does not assert that escape happens at a particular step, only that a net per-repeat survival exists.
2. *k*_ext is the effective extension rate for one repeat. Its relation to a single-turnover chemical rate *k*_chem is not assumed (§7).
3. *k*_loss is defined operationally, as what is left over at saturating primer. Calling it a "committed" state is an extra claim, and we do not make it.
4. The distribution is treated as exponential in *n*. Real gels depart from this for reasons unrelated to mechanism, principally heterogeneity (§5) and length-dependent detection (§6).

Terminology, used consistently: **molecular event** (what happens to the primer), **microscopic rate** (which rate constant changes), **effective rate** (what the assay measures, e.g. *k*_loss), **observable** (what is recorded). A claim may not name a molecular event using an effective parameter. Most of the confusion in this literature, including in our own earlier writing, comes from that substitution.

---

## 3. Model B and model C1 are the same function of the observable

Consider two mechanisms that are usually treated as competitors.

- **Model B — extension-rate enhancement.** TPP1 raises the per-repeat extension rate: *k*_ext → α*k*_ext.
- **Model C1 — affinity-mediated retention.** TPP1 strengthens primer binding: *k*_off → *k*_off/α, at unchanged *k*_ext and *k*_loss.

Both act on the survival probability through the same denominator of Eq. 1. Consequently:

**Result 1.** The distribution depends on the mechanism *only through the scalar s*: *P*(≥ *n*) =
*s*^(n−1). Consequently, at any given primer concentration, a B-parameter set and a C1-parameter set
that produce the same *s* produce **identical** data, and no measurement of the distribution can
separate them — at any precision and any sample size. A single-concentration processivity
measurement is therefore a measurement of *s*, and of nothing else.

The scope of this statement deserves care, because it is easy to overstate in two different
directions.

*It is not a claim that B and C1 are the same model.* They are not. With one fixed parameter set
each they agree at *no* primer concentration: for α = 2 at *k*_ext = 2, *k*_loss = 0.5, the B curve
runs from 0.616 at [P] = 0.01 to a plateau of 0.889, while the C1 curve with *k*_off halved runs
from 0.573 to a plateau of 0.800. What is identical is the *conditional* statement: fix *s*, and the
distribution is fixed.

*Nor is it a claim that a primer titration is useless.* A titration measures the function *s*(P), and
B and C1 differ in its **plateau** — B raises it, C1 leaves it unchanged (Fig. 1a). So a titration can
in principle separate them, but only through the plateau, which means only through *k*_loss — the
worst-determined quantity in the design (§7: 34–49% relative SE, and direct binding does not improve
it). The separation is therefore real in principle and fragile in practice, which is a more accurate
statement than either "indistinguishable" or "distinguishable" alone.

To be explicit about the level at which this is being asserted: it is a statement about the model
equations, not about biology. B and C1 differ in which microscopic rate is perturbed, and that
difference is real; what the equations say is that a distribution read at one concentration reports
only the composite *s*, and therefore cannot report the difference. Figure 1 shows the two halves of
the argument — panel b the conditional identity, panel a the plateau that a titration would have to
resolve.

**Result 2 (scale invariance).** *s* is invariant under θ → λθ applied to *k*_ext, *k*_loss, *k*_off and *k*_on simultaneously. Symbolically, *s*(λθ) − *s*(θ) = 0 exactly. Therefore **no absolute kinetic rate is recoverable from processivity**, at any precision — only ratios. Any reported "rate" from a processivity fit is a rate in units of one arbitrary scale.

Result 2 is why the practical design must be phrased in ratios (§7), and why a fit that returns four absolute rates has not measured four things.

**What breaks the equivalence.** This identity holds because both mechanisms enter through the same
denominator. A mechanism that changes the *plateau* — the concentration-independent loss — is not absorbed. Writing the plateau as *s*_∞ = *k*_ext/(*k*_ext + *k*_loss), a mechanism that lowers *k*_loss raises *s*_∞, and no choice of *k*_off reproduces that. We call this **C2**. It is the only class of mechanism this observable can detect without an additional measurement, and it is exactly the class about which we can say least (§5, §8).

---

## 4. Affinity is not one hypothesis

*K*_d = *k*_off/*k*_on can fall because *k*_off falls (C1a), *k*_on rises (C1b), or both. It is tempting to conclude that these are equivalent at matched *K*_d, and we initially assumed they were. They are not, and the reason matters.

Substituting *K*_d = *k*_off/*k*_on into *R*:

> **R = k_off·K_d/(K_d + [P]) = k_off²/(k_off + [P]·k_on)**   (2)

verified symbolically. Three consequences:

1. **As [P] → 0, R → k*_off*, independent of *K*_d and of *k*_on.** At low primer the effective loss *is* the off-rate.
2. **At matched K*_d*, scaling both rates by λ scales R by λ exactly:** *R*(λθ) = λ*R*(θ). The curve keeps its shape in [P] and changes amplitude.
3. Because *s* = *k*_ext/(*k*_ext + *k*_loss + λ*R*), that λ is **algebraically indistinguishable from a change in *k*_loss**.

So a matched-*K*_d pair is not identical; it differs by an exact multiplicative factor, and the factor is confounded with precisely the parameter the C2 question rests on. Figure 2 shows the decomposition: *K*_d fixes the *shape* of the loss term (panel b, where C1a and C1b coincide), while *k*_off fixes its *size* (panel c, where they differ by λ).

A worked example at *k*_ext = 2.0, *k*_loss = 0.5 makes the size of the effect concrete:

| | *K*_d | *s*([P]=0.1) | *s*(1) | *s*(100) |
|---|---|---|---|---|
| wild type (*k*_off 2.0, *k*_on 1.0) | 2.0 | 0.454 | 0.522 | 0.788 |
| C1a (*k*_off 1.0) | 1.0 | **0.587** | **0.667** | 0.797 |
| C1b (*k*_on 2.0) | 1.0 | **0.463** | **0.571** | 0.794 |

Both variants halve *K*_d, and at low primer they differ by 0.13 in *s* — comparable to the whole B-versus-C2 signal. The separation is cleanest at low primer, which is where processivity is measured worst (short products, few extension events), and they converge only as [P] → ∞.

**Consequence for reporting.** "TPP1 lowers the primer *K*_d" is under-specified. The resolved statement requires *k*_off and *k*_on separately, from direct binding kinetics. That is a stronger requirement than the field's usual practice of reporting an affinity, and this section is the reason for it.

---

## 5. Heterogeneous populations: a closed form, and what it does to the tail slope

Real telomerase preparations are heterogeneous. Some fraction φ of product molecules may come from a subpopulation with a different survival probability *s*₂ — dead or partially active enzyme, a damaged fraction, a second assembly state, or a contaminant activity. It is standard to summarise the resulting distribution by its **tail slope**: the ratio of adjacent band intensities at long product lengths.

For a two-species mixture with *q*(*n*) ∝ (1−φ)(1−*s*)*s*^(n−1) + φ(1−*s*₂)*s*₂^(n−1), the adjacent-band ratio is exactly

> **r(n) = s · [1 + κρⁿ] / [1 + κρ^(n−1)],  κ = φ(1−s₂)/[(1−φ)(1−s)],  ρ = s₂/s**   (3)

with a maximum absolute deviation from direct evaluation of 4.4×10⁻¹⁶ over *n* = 1…400 in our numerical check.

Equation 3 yields three exact results, none of which requires simulation.

**Result 3 (interpolation bound).** With *u* = κρ^(n−1), Eq. 3 is *r* = *s*(1 + *u*ρ)/(1 + *u*), so *dr*/*du* = *s*(ρ−1)/(1+*u*)². Hence *r* is monotone in *u*, *u* is monotone in *n*, and

> **min(s, s₂) ≤ r(n) ≤ max(s, s₂)  for every n and φ.**

The tail-slope estimator can never report a processivity beyond the most extreme subpopulation actually present. This matters for how strongly the bias result may be stated: it bounds the damage. A worst-case relative violation of the bound of 1.0×10⁻¹³ across 5000 random parameter sets was numerically consistent with exactness.

**Result 4 (the tail slope converges to the longest-lived species).** Taking *n* → ∞ in Eq. 3, *r*(*n*) → ρ*s* = *s*₂ if *s*₂ > *s*, and → *s* otherwise. The estimator reports **max(*s*, *s*₂)** — the *longest-lived* species, not the majority one. Consequently the bias in apparent processivity is

- **positive** when the contaminating population is more processive than the majority;
- **negative** when it is less processive;
- **dependent on the product-length window** read, because *r*(*n*) varies with *n* by construction.

**Result 5 (onset length).** The two components of *q*(*n*) are equal when (1−φ)(1−*s*)*s*^(n−1) = φ(1−*s*₂)*s*₂^(n−1), i.e. at

> **n\* = 1 + ln(1/κ)/ln(ρ).**   (4)

Below *n*\* the majority controls the per-band signal; above it the minority does. Equation 4 matches the numerically located midpoint of *r*(*n*) to within one band in every case tested (e.g. *s* = 0.60, φ = 0.01, *s*₂ = 0.95 → *n*\* = 15.5, midpoint at *n* = 15).

**Retirement of a previous approximation.** We previously used Δ*s* ≈ φ(1−*s*) for the bias. That approximation is wrong and is withdrawn. It omits the (1−*s*₂) weighting, it is linear in φ where the true bias is not, and it is set by (*s*₂ − *s*) rather than by φ. In the cases below it is **50–90× too small**, and it cannot change sign, whereas the true bias does:

| majority *s* | φ | *s*₂ | tail slope (n = 40–80) | apparent Δ*s* | φ(1−*s*) predicts |
|---|---|---|---|---|---|
| 0.60 | 0.01 | 0.95 | 0.9500 | **+0.350** | 0.004 |
| 0.70 | 0.01 | 0.90 | 0.8997 | **+0.200** | 0.003 |
| 0.85 | 0.05 | 0.95 | 0.9382 | **+0.088** | 0.0075 |
| 0.999 | 0.01 | 0.95 | 0.9976 | −0.0014 | 0.00001 |

---

## 6. The result that changes practice: the answer depends on the window

Equation 3 makes *r*(*n*) a function of *n*. It is therefore not a property of a distribution at all: it is a property of a distribution **and a length window**. Two labs reading the same gel over different band ranges will report different processivities, and both can be internally consistent.

With a true majority *s* = 0.60, φ = 0.01 and *s*₂ = 0.95:

| true *s* | φ | *s*₂ | *n*=3–10 | 10–16 | 16–25 | 25–40 | 40–80 | spread |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 0.01 | 0.95 | 0.610 | 0.707 | 0.909 | 0.950 | 0.950 | **0.340** |
| 0.70 | 0.01 | 0.90 | 0.703 | 0.716 | 0.771 | 0.875 | 0.900 | **0.196** |
| 0.85 | 0.05 | 0.95 | 0.853 | 0.857 | 0.864 | 0.889 | 0.938 | **0.085** |
| 0.999 | 0.01 | 0.95 | 0.986 | 0.989 | 0.991 | 0.994 | 0.998 | **0.012** |

The first row is the practical problem. Reading the near window returns *s* = 0.610 and the far window *s* = 0.950. Inferred *k*_loss/*k*_ext = (1−*s*)/*s* is 0.639 versus 0.053 — a **12-fold** difference in the quantity a mechanism claim would be made of, from identical data (Fig. 5b).

The same mechanism fabricates a concentration-independent loss term where none exists. With *k*_loss = 0 (true plateau 1.0) and 0.5–2% contamination, the tail-slope estimator reads 0.966–0.983 and infers *k*_loss/*k*_ext ≈ **0.017–0.036**; read over a near window the same distributions imply 0.05–0.10. A C2 claim is precisely a claim about residual loss after saturation, so this estimator is disqualified from making one.

**Therefore: the tail slope is forbidden as a basis for any mechanistic claim, not merely deprioritised.** Every C2 conclusion must come from a full-distribution fit in which φ is a free parameter. Such a fit does not make this error, because the mixture weight and the majority *s* are estimated jointly rather than read off one window (§8, gate 3). Whether a real experiment is in the regime where this matters is an empirical question: it depends on whether φ > 0 at a level where *n*\* falls inside the readable band range. That is testable on existing data (§10).

**A saturation subtlety.** At *s* → 1 the bias from a slower contaminant is small and negative (−0.0014 at *s* = 0.999) but the *finite-window* reading is biased slightly downward, because the slower species has not yet decayed out. An earlier draft of this work asserted a −0.040 bias at saturation; that number is not reproducible at any window and is withdrawn here.

---

## 7. Identifiability: structural rank, practical precision, and equivalence

**Structural rank.** Differentiating *s* with respect to (*k*_ext, *k*_loss, *k*_off, *k*_on) at 12 primer concentrations using complex-step derivatives (finite differences give unreliable rank here), the Jacobian has rank **3 of 4**, with a smallest relative singular value of 3.9×10⁻¹⁷. The deficiency is structural, not a sampling artefact. Adding single-turnover chemistry increases rank to 4 **only if one assumes the measured step is the model's extension step** (*k*_chem ≡ *k*_ext); with a free coupling parameter η between them, the design gains one observable and one parameter and returns rank **4 of 5**. Adding direct *k*_on and *k*_off measurements gives rank 4 of 4 with a smallest relative singular value of 0.20. **Direct primer-binding kinetics is the only independent measurement in this design space that adds genuine rank.**

**Practical precision.** Full rank is necessary, not sufficient. From a log-space Fisher matrix with 15% assay CV on *s* at 12 primer concentrations, the identifiable ratio form gives condition number 60.3 and relative standard errors of 44% (*k*_loss/*k*_ext), 12% (*k*_off/*k*_ext) and 60% (*k*_on/*k*_ext). In the absolute four-rate form:

| design | condition number | *k*_ext | *k*_loss | *k*_off | *k*_on |
|---|---|---|---|---|---|
| processivity vs [P] | 2.4×10⁸ | 31% | 34% | 29% | 67% |
| + direct binding at 10% CV | **7.5** | 22% | **49%** | 9.9% | 9.9% |

Direct binding collapses the condition number by seven orders of magnitude and fixes *k*_off and *k*_on at the precision they are measured with. It does **not** fix *k*_loss: the relative standard error moves 34.3% → **49.0%**, i.e. it gets worse, because binding removes the correlations that were previously absorbing part of that variance. *k*_loss is the worst-determined quantity in every design we tested, and it is exactly the parameter the C2 claim depends on. **Direct binding is necessary but not sufficient for a C2 claim**, and nothing here should be read as saying the binding experiment settles it.

**Equivalence design.** Saturation cannot be established by a null result. It requires replicated measurements, a prespecified maximum acceptable change δ, and confidence intervals lying inside ±δ (TOST) — "flat within the CV" is also satisfied by a dead assay. With paired TOST at 80% power, two-sided, Bonferroni over three doublings:

| δ (log) | CV 0.345 (published) | CV 0.15 (optimised assay) |
|---|---|---|
| 0.10 | 74 (1 test) / **105** (3 tests) | 14 / **20** |
| 0.20 | 19 / 27 | 4 / 5 |

δ must come from the smallest biologically meaningful processivity difference, not from convenience. Adjacent published processivity ratios differ by 0.05–0.38 log, so δ = 0.20 would fail to separate neighbouring mutants the field already treats as distinct. **We do not fix δ**; it is a decision for the experimenter, and the power calculation above is what it costs.

---

## 8. What this does and does not establish

**Exact, within the stated model class.** B ≡ C1 under matched *s* (Result 1); scale invariance (Result 2); the rewrite of the loss term and the λ-confound (Eq. 2); the closed form for *r*(*n*) and the interpolation bound, longest-lived limit and onset length (Results 3–5); the degree-1 crossing between B and C2, which means at most one crossing for all α, β > 0, α ≠ β, with *P*\* = *k*_off(−αβ*k*_off − α*k*_L + β*k*_L + β*k*_off)/[*k*_L*k*_on(α−β)]. The sign of *P*\* is a precondition of any crossing statement: when *P*\* < 0 there is no coincidence point in the physical domain and the correct report is "one mechanism dominates throughout," which is a different claim from either "they cross" or "they are indistinguishable."

**Computed, and conditional on the model.** The rank and precision tables. They are arithmetic consequences of Eq. 1, not measurements.

**Simulation-derived.** That the mixture bias is *practically* relevant — that real preparations have φ and *s*₂ in the range where *n*\* falls inside a readable window. Eq. 3 is exact; the claim that it describes a real gel is not yet tested.

**Not established, and not claimed.** We do not claim TPP1 changes primer retention, extension chemistry, or any committed state. We do not claim a mechanical or allosteric component exists, or that it does not. We do not claim C2 is a molecular species rather than an effective residual loss. We do not promote affinity/anchoring above conformational regulation: within this model they are the same function of the observable, which is a statement about identifiability, not about biology.

**The boundary, stated once.** A mathematically exact result within a model is not an experimentally established biological fact. Results 1–5 are of the first kind. Turning any of them into the second kind requires an independent measurement that constrains the underlying physical rates, and the analysis above specifies which measurements can do that and which cannot.

---

## 9. Relationship to prior work

The kinetic case that POT1–TPP1 slows primer dissociation and aids translocation rests on primer-extension and dissociation measurements [6], and the structural case places TPP1 at the TEN–IFD–TRAP surface, stabilising the DNA anchor site [3,4]. Our Result 1 is consistent with both: it explains why those two descriptions can both be right while neither identifies a mechanism, since a change in retention and a change in extension rate are absorbed into the same parameter by a processivity measurement.

Two later developments frame the result. RPA stimulates processivity through a surface distinct from TPP1's [5], which is a second independent modulation of the same observable and therefore subject to the same identifiability limits — a design that cannot separate B from C1 for TPP1 cannot separate it for RPA either. And recent cryo-EM work has begun to resolve nucleotide-bound states across the repeat-addition cycle [7], which is the kind of structural information that could in principle distinguish retention from catalysis, provided it is coupled to kinetics rather than inferred from occupancy.

To our knowledge the mixture closed form (Eq. 3) and the window dependence of the tail-slope estimator are not standard in this literature, although the underlying statistics of finite mixtures are classical; our contribution is the specific consequence for RAP assays and the demonstration that it disqualifies a widely used summary statistic. We would welcome being shown otherwise — a prior derivation of Eq. 3 in this context would improve the paper and would not weaken it.

---

## 10. The next experiment is not a simulation

The analysis above is worthless if it does not change what someone measures. It changes three things.

**First, an existing gel can test the practical claim.** Result 4 predicts that a full-distribution fit with φ free will place the majority *s* below the tail slope whenever a long-lived contaminant is present, and Eq. 4 predicts where the crossover occurs. One gel with recoverable per-band intensities — densitometry, not a new experiment — can test this. It can also **falsify** it: if φ ≈ 0, the mixture bias has no practical importance for that dataset and should be reported as such rather than kept as a caution. Gate 3 in §11 is that test. It is the next real advance, and it requires no simulation.

**Second, saturation must be shown, not assumed.** With both primer and TPP1 varied over a range established by TOST against a prespecified δ and a positive control proving the assay can still resolve a change, the models separate as follows: if TPP1 stimulation persists at saturating primer, C1 is rejected and the effect is on the concentration-independent term; if it vanishes under verified saturation, B and C2 are rejected.

**Third, direct binding kinetics is the measurement that adds rank**, and it must report *k*_off and *k*_on separately rather than as *K*_d (§4), with the caveat that it does not by itself rescue *k*_loss (§7).

A dNTP × TPP1 titration was considered and **rejected** as a route to this question. In a model with nucleotide-dependent pausing, an interaction of 1.62 can appear with no change in extension rate, so the readout does not discriminate; the interaction in the feasible concentration window was 1.00/0.85/0.84, and extending the range further drove nucleotide occupancy below the gel floor.

---

## 11. Release gates for any mechanistic claim

We state these so that a future claim — by us or by others — can be checked against them rather than argued about.

| gate | status | what it requires |
|---|---|---|
| 1. Symbolic verification | **satisfied** | invariance, crossing, and parameter domains proved, not sampled |
| 2. Structural identifiability | **satisfied** | rank deficiency established with the design stated; structural separated from practical |
| 3. Observation-model validation | **conditional** | real per-band densitometry: does a full-distribution fit with φ free show φ > 0 and a majority *s* below the tail slope? |
| 4. Kinetic mapping *k*_chem ≡ *k*_ext | **experimental** | a dNTP-chase test confirming the single-turnover step is rate-limiting for processivity |
| 5. Equivalence margin δ | **pending** | a biological judgement on the smallest meaningful effect, not a statistical convenience |

Gates 1–2 are properties of the analysis and are met. Gates 3–5 are properties of the evidence and are not. **No mechanistic claim should be published on the strength of gates 1–2 alone**, which is the failure mode this paper is written against.

---

## 12. Reproducibility

Every number in this paper is regenerated by one command:

```bash
python scripts/mechanism_classification_checks.py   # ~2 s; writes classification_checks.json + .csv
python paper/make_figures.py                        # regenerates all six figures
pytest tests/test_mechanism_classification.py -q    # 25 guards on the numbers and the claims
```

The checks script is self-contained (numpy + sympy) and writes machine-readable output; the test module fails if any quoted number, or any retracted claim, reappears. Source data used for the published processivity values are deposited with SHA-256 hashes and provenance notes in `data/external/tesmer2023/`. The analysis document underlying this paper is `results/causality/MECHANISM_CLASSIFICATION.md`, frozen by hash, with its own retraction record (§0.1 there).

Repository: `https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch`

We note plainly that this repository also contains earlier claims that were subsequently retracted — a "one continuous IFD–TRAP helix," a K757A separation-of-function "breakthrough," an allosteric-coupling result that proved inconclusive. They are retained rather than deleted, and the corrections are recorded alongside them. The same standard is applied here: the three errors in this paper's own first draft are listed in §13 rather than quietly fixed.

---

## 13. Corrections made to this work before submission

| error | how it was caught | correction |
|---|---|---|
| claimed C1a and C1b are identical at matched *K*_d | the check written to verify it returned *false*; deriving *R* in closed form showed an exact factor λ | §4; the confound is with *k*_loss, not an identity |
| displayed *P*\* with the opposite sign | comparison against `sympy.solve` and a 20 001-point numeric scan; the displayed version put *P*\* > 0 where no crossing exists | §8; numeric confirmation table added (residual 1.1×10⁻¹⁶) |
| claimed a −0.040 saturation bias from the mixture | direct evaluation of *r*(*n*) gave −0.0014; the bias is window-dependent, not a single number | §6; replaced by the window sweep and the 12-fold result |

The Δ*s* ≈ φ(1−*s*) approximation (§5) is likewise withdrawn. Each correction now has a regression test.

---

## References

1. Wang F, Podell ER, Zaug AJ, Yang Y, Baciu P, Cech TR, Lei M. The POT1–TPP1 telomere complex is a telomerase processivity factor. *Nature* **445**, 506–510 (2007). doi:10.1038/nature05454
2. Padmanaban S, Tesmer VM, Nandakumar J. Interaction hub critical for telomerase recruitment and primer-template handling for catalysis. *Life Science Alliance* **6**, e202201727 (2023). doi:10.26508/lsa.202201727
3. Sekne Z, Ghanim GE, van Roon A-MM, Nguyen THD. Structural basis of human telomerase recruitment by TPP1–POT1. *Science* **375**, 1173–1176 (2022). doi:10.1126/science.abn6840
4. Liu B, He Y, Wang Y, Song H, Zhou ZH, Feigon J. Structure of active human telomerase with telomere shelterin protein TPP1. *Nature* **604**, 578–583 (2022). doi:10.1038/s41586-022-04582-8
5. Agrawal S, Lin X, Susvirkar V, O'Connor MS, Chavez BL, Tholkes VR, Tauber GP, He Q, Abe K, Huang X, Lim CJ. Human RPA is an essential telomerase processivity factor for maintaining telomeres. *Science* **390**, 495–502 (2025). doi:10.1126/science.ads5297
6. Latrick CM, Cech TR. POT1–TPP1 enhances telomerase processivity by slowing primer dissociation and aiding translocation. *EMBO Journal* **29**, 924–933 (2010). doi:10.1038/emboj.2009.409
7. Structures of nucleotide-bound human telomerase at several steps of its telomeric DNA repeat addition cycle. *Nature Communications* (2026). doi:10.1038/s41467-026-68560-8

*Reference 7 is cited for the general claim that nucleotide-bound states across the repeat-addition cycle have been resolved by cryo-EM; readers should consult the primary article for its specific findings. Reference 2 is the source of the processivity and activity values used in the repository's deposited data and is referred to there as "Tesmer 2023"; the first author is Padmanaban.*

---

## Figure legends

**Figure 1. What a processivity distribution can and cannot identify.**
(a) Effective survival *s* versus primer concentration for wild type and for three mechanisms: extension-rate enhancement (B), affinity-mediated retention (C1), and a reduction in concentration-independent loss (C2). B and C2 raise the plateau *s*_∞; C1 leaves it at exactly 0.800, because it acts only on the primer-dependent term. (b) Two unrelated mechanisms tuned to the same *s* = 0.75 produce indistinguishable distributions: *P*(≥ *n*) = *s*^(n−1) superimposes exactly, with max|Δ*P*(≥*n*)| = 0.

**Figure 2. Affinity is not one hypothesis.** (a) Wild type and the two affinity variants C1a (*k*_off ↓) and C1b (*k*_on ↑) tuned to the same *K*_d = 1. They are not equivalent. (b) Normalising the loss term by *k*_off, *R*/*k*_off = *K*_d/(*K*_d+[P]) depends on *K*_d alone, so C1a and C1b coincide: *K*_d fixes the shape. (c) In absolute terms *R* = *k*_off*K*_d/(*K*_d+[P]) scales with *k*_off, so the two differ by an exact factor λ = 2 — which is algebraically confounded with *k*_loss.

**Figure 3. The B-versus-C2 crossing, and its sign.** Difference in survival *s*_B − *s*_C2 versus primer concentration. (a) For α = 2, β = 3.42 the difference changes sign at *P*\* = 17.3, and the closed form satisfies *s*_B = *s*_C2 to 1.1×10⁻¹⁶ there. (b) For α = 30, β = 3.42 the solution is *P*\* = −31.9: there is no root in the physical domain and one mechanism dominates throughout. Reporting (b) as a crossing, or as indistinguishability, are both errors; the sign of *P*\* decides which regime applies.

**Figure 4. The mixture ratio and its three regimes.** (a) Adjacent-band ratio *r*(*n*) from Eq. 3 for four mixtures; each curve is monotone and lies between min(*s*,*s*₂) and max(*s*,*s*₂) (dotted lines). (b) The same distribution (*s* = 0.60, φ = 0.01, *s*₂ = 0.95) read over three windows returns 0.61, 0.91 and 0.95. The minority component overtakes the majority at *n*\* = 15.5 (Eq. 4).

**Figure 5. Window dependence.** (a) Reported tail slope versus read window for four mixtures. For a true *s* = 0.60 the reported value runs from 0.61 to 0.95. (b) The implied *k*_loss/*k*_ext for the same distribution: 0.630 from the near window and 0.053 from the far one, a 12-fold spread from identical data.

**Figure 6. Identifiability.** (a) Relative singular values of the processivity Jacobian: rank 3 of 4 for processivity alone (drop below the 10⁻⁸ cut-off at index 4); rank 4 when single-turnover is added under the assumption *k*_chem ≡ *k*_ext, and when direct binding is added. (b) Relative standard errors per parameter. Direct binding fixes *k*_off and *k*_on at 9.9% and reduces the condition number by seven orders of magnitude, but *k*_loss worsens from 34.3% to 49.0% — the parameter a C2 claim depends on is the one binding does not rescue. (c) Replicates per arm for a paired TOST equivalence design at 80% power, Bonferroni-corrected over three primer doublings, as a function of assay CV for δ = 0.10 and 0.20 log. At the published CV of 0.345, δ = 0.10 requires 105 replicates per arm.
