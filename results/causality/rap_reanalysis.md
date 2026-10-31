# RAP re-analysis — what the existing data support, and what they do not contain

**Status:** descriptive only. **No model was fitted.** This document re-analyses the numbers that
exist, states their precision, and records the two quantities the source omits.

**Source:** `data/external/tesmer2023/` (Tesmer et al. 2023, Fig 3B/3C/3E source data).
Provenance, SHA-256 and replicate semantics: `data/external/tesmer2023/README.md` and
`SOURCE_MANIFEST.csv`. Rules governing this analysis: `PREREGISTRATION.md` §§3–4.

---

## 1. The measurements, re-derived from the source file

n = 3 independent experiments per genotype; 95% CI = t(0.975, 2) × SEM. R1 = paired
within-experiment processivity ratio vs WT. **R2 (absolute retained product mass) is not in the
source data** — see §4.

| genotype | P(≥16) mean | 95% CI | R1 processivity ratio | 95% CI | activity ratio | 95% CI |
|---|---|---|---|---|---|---|
| WT | 0.485 | 0.485–0.485 | 1.00 | — | 1.00 | — |
| K749E | 0.279 | −0.020–0.578 | 0.56 | — | 0.567 | −0.696–1.830 |
| R756E | 0.287 | 0.086–0.489 | 0.58 | — | 0.466 | −0.706–1.639 |
| K757E | 0.197 | 0.159–0.234 | 0.42 | — | 0.427 | −0.297–1.152 |
| F759E | 0.117 | −0.067–0.301 | 0.24 | — | 0.386 | −0.324–1.096 |
| K760E | 0.264 | 0.163–0.364 | 0.55 | — | 0.333 | −0.128–0.793 |
| **H762E** | **0.416** | 0.137–0.694 | **0.85** | — | **0.680** | −0.864–2.224 |
| F776E | 0.013 | −0.748–0.773 | **undefined** | — | 0.418 | −0.784–1.621 |
| L786E | 0.171 | 0.024–0.318 | 0.35 | — | 0.248 | −0.459–0.955 |
| **D788R** | **0.312** | −0.043–0.667 | **0.62** | — | **0.235** | −0.573–1.039 |

CIs are wide because n = 3 and the between-experiment spread is large. Every mutant's activity CI
includes 1.0 (no significant activity difference from WT) **and** several include values above
1.0; every mutant's processivity CI overlaps at least one other mutant's. Numeric values:
`real_data_summary.csv`.

## 2. Precision: what this dataset can and cannot resolve

From the data's structure, not from any model (frozen in `PREREGISTRATION.md` §4.3):

| readout | median CV | minimum detectable ratio at n = 3 (α=0.05, 80% power) |
|---|---|---|
| processivity | 34.5% | **≈2.2×** |
| activity | 90.5% | **≈7.9×** |

Consequences that must be respected downstream:

* The widely quoted contrasts **"86% of WT processivity" (H762E)** and **"64%" (D788R)** are
  differences of 14 and 36 percentage points inside a blur of roughly ±30 points. They are
  **ordinal** statements ("mild", "intermediate"), not measurements with 1–2% precision. Any model
  calibrated to them at face value would be fitting noise.
* The activity readout cannot separate any mutant from WT, or from any other mutant. **Activity is
  therefore not usable as a quantitative constraint on a kinetic model** in this dataset; it can
  only order mutants (H762E > K749E ≈ R756E > K757E > F776E ≈ F759E ≈ K760E > L786E > D788R).
* Only two contrasts survive at the detectable effect size: **H762E vs the anchor mutants**
  (R1 0.85 vs 0.24–0.58) and **F776E vs everything** (dead).

## 3. The one replicated ordinal pattern

Ordering by R1, with mutants grouped where their CIs overlap:

1. **H762E ≈ WT** (0.85, CI spans 0.14–0.69 on the absolute scale)
2. **K749E, R756E, D788R, K760E** (0.55–0.62) — template contacts (K749/R756), TEN–IFD interface
   (D788), and the unassigned K760
3. **K757E, L786E** (0.35–0.42) — DNA backbone contact, and a core residue
4. **F759E** (0.24) — DNA base stacking

This ordering is the only structure in the dataset that a model could be asked to reproduce, and
it is a weak target: one mutant (H762E) sits at WT and the rest span 0.24–0.62 with overlapping
CIs.

## 4. What is missing — the two quantities that decide the question

**(a) There is no product-length distribution.** Each replicate yields a single scalar (fraction
≥16 repeats). The models in `kinetic_model_comparison.md` need per-band intensities as a function
of length. Those exist only in the gel images (`SdataF3.2.tif`, `SdataF3.3.tif`), which cannot be
retrieved here. **The distribution analysis is blocked on data, not on method.**

**(b) There is no ΔRAP.** ΔRAP requires the −POT1–TPP1 lanes. The numeric source file contains
**only the +PT lanes**, so ΔRAP is **not computable from this dataset at all**. The −PT lanes are
in the gel image (Fig 3A) but not in the numbers. Any statement about "fold-stimulation" or
"stimulation index" from this source would be fabricated.

**(c) There is no absolute retained mass (R2).** Fig 3C is normalised to WT within each experiment,
so inter-experiment intensity differences — the only route to absolute mass — are removed by
construction.

**Why (c) matters, concretely.** D788R has R1 = 0.62 but an activity ratio of 0.235. A model that
explains the R1 alone can say "D788R is stimulated normally". A model that must also explain why
the absolute stimulated product mass is a quarter of WT has to account for the *catalytic* defect
as well. Without R2, that constraint is unavailable, and any model comparison is correspondingly
weaker. This is the strongest argument for re-quantifying the gels rather than re-analysing the
tables.

## 5. Verdict

* The source data are **sufficient to state the problem precisely** and **insufficient to solve
  it**: one scalar per replicate, three replicates, CVs of 35–90%.
* Nothing here supports a claim that affinity does or does not explain the phenotype. The most
  that can be said is that the **H762E contrast is real at this dataset's resolution** and that
  every other comparison is at or below it.
* Two specific, cheap data recoveries would change this: **per-band densitometry** of the existing
  gels (unlocks the distribution models) and **the −PT lanes** (unlocks ΔRAP). Neither requires a
  new experiment; both require the images or the original data.
