# PREREGISTRATION — Affinity–Processivity Causality Investigation

**Status:** frozen before any model fitting. Commit this file; do not edit it after the first fit
run. Corrections go in `PREREGISTRATION_AMENDMENTS.md` with a date and a reason, never as silent
edits here.

**Committed on branch:** `arena/f2375cbf-telomerase-allosteric-switch`, based on `main` = `94a857d`.

---

## 1. The single question

> Does a change in primer retention / DNA-binding affinity quantitatively account for the
> TPP1-dependent increase in telomerase processivity, without requiring a mechanical allosteric
> signal?

This is a narrower question than "what is the mechanism", and it is the one the evidence can
address. Three outcomes are pre-declared as acceptable: **affinity sufficient**, **affinity
insufficient**, or **not identifiable with available data**. The third is a result, not a failure.

## 2. Frozen directional predictions (holdout = Sekne et al. 2022)

These predictions are written down **before** any Sekne figure is examined for this investigation,
and they are directional on purpose: Sekne provides qualitative rescue/no-rescue observations, and
direction is all that can be tested without inventing precision the source does not have.

Models: **C** = affinity/anchoring (retention), **B** = conformational/translocation-efficiency,
**A** = local conduit (760–766; already contradicted, included for completeness).

| # | Held-out observation (Sekne 2022, Fig 4E–G) | C predicts | B predicts | A predicts |
|---|---|---|---|---|
| P1 | Anchor-site mutants (Y176A/Q177A, K757A/F759A, Q794A) have **reduced** processivity with TPT vs WT | reduced — the anchor is the retention element | reduced, no stronger claim | reduced |
| P2 | Adding TPT **rescues** processivity of each anchor mutant (mutant+TPT > mutant alone) | **rescued** — mass action compensates weaker binding | rescue is not predicted by a pathway model (it is not forbidden, but it is not expected) | **not rescued** — a broken conduit cannot be bypassed by more ligand |
| P3 | The magnitude of rescue tracks the magnitude of the DNA-contact defect | **monotone concordance** (stronger anchor defect → larger rescue headroom) | no ordering predicted | no ordering predicted |
| P4 | Mutants that keep TPP1 binding but damage DNA contact are still stimulable | **yes** (stimulation needs binding, not contact) | yes | no |
| P5 | Mutants whose defect is a folding lesion (F776E, L786E in the Tesmer panel) are not informative about either model | silent | silent | silent |
| P6 | Any observation in Sekne that separates retention from translocation rate | none expected in this dataset | none expected | none expected |

**P3 is the only prediction that could actively favour C over B.** P1, P2 and P4 are all
consistent with C and merely unhelpful for B; that asymmetry is recorded now so that a later
report cannot present "C predicted rescue and B did not" as if it were a symmetric test.

**Predictions are not revised after seeing Sekne.** If a prediction fails, it is reported as failed.

## 3. Data provenance and the verification gate

Provenance is frozen in `data/external/tesmer2023/SOURCE_MANIFEST.csv` and explained in that
directory's `README.md`. Summary:

| item | value |
|---|---|
| Development dataset | Tesmer et al. 2023, *Life Sci Alliance* 6:e202201727, **Fig 3B/3C/3E** source data |
| Upstream file | `LSA-2022-01727_SdataF3.1.pdf`, PMC10055720 |
| Upstream SHA-256 | **not obtainable** — direct download is blocked in this environment; recorded as `NO`, not left blank |
| Local SHA-256 | recorded per file in `SOURCE_MANIFEST.csv` (guards against later silent edits) |
| Measurement type | per-experiment derived scalars: processivity fraction (≥16 repeats); activity ratio normalised to WT **within experiment** |
| Independent experiments | **n = 3** per mutant, both readouts |
| Technical replicates in file | **none** |
| Held-out dataset | Sekne et al. 2022, *Science* 375:1173, Fig 4E–G — **qualitative direction only** |
| Held-out construction difference | ΔTEN *trans*-complementation with 500 nM POT1–TPP1 vs reconstituted telomerase ± TPT complex |

**Transcription verified by reproduction, not by hash.** Column means of the transcribed CSVs
recover H762E = 85.6% processivity (text: 86%), D788R = 64.3% (64%), H762E = 68.0% activity (68%),
D788R = 23.5% (24%). Four independently stated published numbers are reproduced from the file.
This is automated in `tests/test_causality_preregistration.py`, which allows a tolerance of one
percentage point. **One discrepancy is recorded rather than smoothed over**: the mean of the D788R
activity rows is 23.5%, while the paper's text says 24%. H762E reproduces exactly (68.0% vs "68%").
The 0.5-point gap is consistent with the authors summarising the source values differently (or
rounding from a slightly different mean), and it is why the preregistration treats published
percentages as ordinal, not exact.

**Exclusions, frozen now:**

| exclusion | reason |
|---|---|
| `F776E` experiment 3 processivity = −0.33789 | negative fraction is physically impossible (≥16-repeat signal at/below background); mutant flagged `processivity_undefined_at_low_signal`, excluded from quantitative comparison, retained as "essentially dead" qualitatively |
| `L786E` | same folding-lesion class; reported, not used to discriminate models |
| activity *absolute* magnitudes across experiments | Fig 3C is normalised to WT within each experiment; absolute activity is therefore not recoverable from the source and must not be claimed |
| any Sekne number not stated in the text | figures are not digitised for quantitative use |

## 4. Inferential rules (fixed in advance)

**4.1 Unit of replication.** The independent experiment. Not gel lanes, not bands, not MD frames.
n = 3 everywhere. No analysis may report a per-band or per-pixel n.

**4.2 Two response variables, both required.** For every condition:

* **R1 — processivity ratio:** paired within-experiment ratio mutant/WT (the quantity the model
  explains), and
* **R2 — absolute retained product mass:** the total ≥16-repeat signal relative to WT in the same
  experiment, unnormalised.

ΔRAP is a **derived** quantity: ΔRAP = (processivity with POT1–TPP1) ÷ (processivity without
POT1–TPP1), computed per experiment, then summarised. The source data contain **only the +PT
condition**, so ΔRAP **cannot be computed from them at all**. This is stated here so that no later
report implies otherwise: the −PT lanes are in the gels (Fig 3A) but not in the numeric source file.
A model that fits R1 while missing R2 is reported as failing.

**4.3 Precision limits, computed from the development data's structure (not from model results).**

Median CV across mutants: **processivity 34.5%**, **activity 90.5%**, excluding WT (fixed at 1 by
normalisation) and F776E (flagged, §3). Minimum detectable ratio at n = 3, two-sided α = 0.05, 80% power:

| readout | CV | n per group | smallest detectable ratio at n = 3 |
|---|---|---|---|
| processivity | 34.5% | 3 | **≈2.2×** |
| activity | 90.5% | 3 | **≈7.9×** |
| processivity, if CV were 15% | 15% | 3 | ≈1.4× |

So the existing data can resolve a 2.2-fold difference in processivity and essentially nothing
finer, and an 8-fold difference in activity. **Published claims of "64% vs 100%" or "68% vs 100%"
are differences well inside this blur** and must be reported as ordinal, not quantitative.

**4.4 Non-identifiability rule (binding).** If the best-fitting retention model and the
best-fitting conformational model (i) achieve indistinguishable held-out predictive performance,
or (ii) differ by less than the measurement error of the available data, the conclusion is
**NON-IDENTIFIABILITY**. No model is declared a winner. The report must then state which specific
measurement (replicate count, dynamic range, or per-band resolution) would break the tie, and by
how much.

**4.5 No post-hoc winner.** A model may not be preferred because it has lower in-sample fit.
Selection requires either (a) better held-out predictive performance, or (b) a pre-declared
directional prediction (P1–P6) that it alone makes and that is confirmed.

**4.6 A geometric distribution is not evidence for affinity.** A constant per-repeat hazard
produces a geometric product-length distribution. So do heterogeneous enzyme populations, an
unmodelled population of inactive molecules, and repeat-dependent extension probabilities that
average out. Any report that treats geometric shape as a fingerprint of retention has failed this
preregistration.

**4.7 Language rules.** "Consistent with" is not "demonstrates". Model C may be described as
best-supported and must not be described as a proven mechanism. Nothing in this investigation can
promote C beyond "affinity/anchoring model".

## 5. Feasibility verdict (gate G0) — recorded before modelling

The verification in §3 surfaced two findings that constrain the whole plan. Both come from the
structure of the data, not from any model result, so they do not consume the holdout.

**F1. There is no product-length distribution in the numeric data.** Each replicate contributes one
scalar per mutant (fraction ≥16 repeats). The five distribution models in §6 require per-band
intensities versus product length. Those exist only inside the gel images
(`SdataF3.2.tif`, `SdataF3.3.tif`; 8.1 and 7.5 MB), which cannot be retrieved here. **Distribution
modelling is blocked on data.**

**F2. n = 3 with CVs of 35–91%** cannot constrain a model with two or more free parameters per
mutant group. Any fit would be under-determined, and reporting a best-fit would be reporting
noise.

**Consequence — two pre-declared arms, and which one is active is decided by data availability,
not by results:**

* **Arm A (distribution data obtained).** Fit G1–G5 to per-band intensities; validate
  directionally on Sekne P1–P6; report held-out performance. Only in this arm can "affinity
  sufficient / insufficient" be concluded.
* **Arm B (no distribution data — the current state).** No fit is performed and no model is
  declared a winner. Deliverables become: (i) formal model statements and analytic expectations
  for each distribution family; (ii) an **identifiability analysis** — simulated *hypothetical*
  datasets under each model, used only to ask how much data would be needed to tell them apart;
  (iii) the decisive-experiment specification. Arm B may conclude "not identifiable with existing
  data" and may recommend a measurement; it may **never** conclude that affinity is or is not the
  mechanism.
* **Transition rule.** If per-band data arrive (from the authors' source data, from gel images, or
  from a new experiment), Arm A activates and this file is not edited — the arm change is recorded
  as an amendment with the data's provenance.

## 6. The competing models (pre-specified; two response variables each)

Let *p(n)* be the probability that a product reaches ≥ *n* repeats; *s* = per-repeat survival; the
observed processivity metric is the fraction of signal ≥ 16 repeats.

| ID | Model | Form | Free parameters | What it implies for the distribution |
|---|---|---|---|---|
| **G1** | Constant-hazard retention (affinity) | p(n) = s^n; s = k_ext/(k_ext + k_off) | 1 per condition (s) | strict geometric, straight line on log-scale |
| **G2** | Heterogeneous population | p(n) = Σ w_i s_i^n, i = 1,2 | 3 | concave on log-scale (mixture of slopes) |
| **G3** | Repeat-dependent extension | s_n varies with n (e.g. translocation efficiency decaying with product length) | 2 | convex on log-scale |
| **G4** | Retention with reversible pause | two-state: primer can re-engage after a pause, so dissociation is not absorbing | 2 | excess at long products vs G1 |
| **G5** | Conformational (TPP1 changes k_ext, not k_off) | s = k_ext(TPP1)/(k_ext(TPP1) + k_off), with TPP1 entering the numerator | 2 | *shifts* the distribution's scale; in the simplest case indistinguishable from G1 by shape alone |

**G5 is the model that makes the question answerable, and it is analytically awkward: changing
k_ext and changing k_off both raise s.** They differ in what else moves — the total product mass,
the fraction of enzyme that ever elongates, and the dependence on primer concentration — not in the
geometric shape. This is recorded now because it is the most likely source of a false "affinity
wins" conclusion later.

## 7. Deliverables and what each may claim

| file | may claim | may not claim |
|---|---|---|
| `results/causality/rap_reanalysis.md` | what the source data do and do not contain; precision limits; ordinal re-reading of the mutant panel | any ΔRAP value (not in data) |
| `results/causality/kinetic_model_comparison.md` | formal model definitions; analytic shape differences; identifiability | a winning model without held-out evidence |
| `results/causality/held_out_predictions.md` | P1–P6 tested as TRUE / FALSE / UNTESTABLE | quantitative Sekne numbers |
| `results/causality/decisive_experiment.md` | the smallest measurement that separates retention from conformational regulation, with expected outcomes under each model | that such an experiment has been done |
| `results/causality/figures/figure_causality_summary.png` | three panels as specified by the requester | extrapolation beyond the panels' data |

## 8. What would change the plan

* **New data** (per-band intensities, −PT lanes, or a new replicate set): Arm B → Arm A, by
  amendment.
* **A failed prediction P1–P6**: reported as a failure; if P3 fails, the quantitative form of C is
  weakened while its qualitative form survives, and the report must say exactly that.
* **A finding that G1 and G5 are analytically identical under the observable**: this becomes the
  headline result (non-identifiability), not a reason to pick one.

## 9. What this preregistration refuses to do

* Fit a model and call the fit validation.
* Treat n = 3 technical or digitised values as independent biological replicates.
* Declare affinity causal because a retention model fits.
* Declare allostery absent because no conformational model can be fitted.
* Digitise bar heights and present them as measurements.
* Produce a winner from data that cannot separate the competitors.
