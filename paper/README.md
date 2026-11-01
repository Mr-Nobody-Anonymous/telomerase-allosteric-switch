# Manuscript package

**What Processivity Distributions Cannot Tell Us About TPP1**
*An identifiability analysis of telomerase repeat-addition processivity, with a corrected treatment of heterogeneous populations*

**Author:** Beamlak Sintayehu · **Date:** 2026-11-01

---

## What is here

| file | what it is |
|---|---|
| [`MANUSCRIPT.md`](MANUSCRIPT.md) | the full paper, in Markdown, with all six figure legends and the reference list |
| `make_figures.py` | regenerates every figure from the same equations the checks script verifies |
| `figures/` | `fig1`–`fig6` as 300 dpi PNG and vector PDF, plus `MANIFEST.json` |
| [`COVER_LETTER.md`](COVER_LETTER.md) | covering note for correspondence, with an honest statement of scope |

## Regenerate everything

```bash
python scripts/mechanism_classification_checks.py   # every quoted number  (~2 s)
python paper/make_figures.py                        # all six figures
pytest tests/test_paper_manuscript.py -q            # 29 guards on the text
```

`tests/test_paper_manuscript.py` fails if any number in the manuscript stops matching the
computation, if a figure is missing, if a reference loses its DOI, or if any of the three retracted
claims reappears in the text.

## The paper in one paragraph

TPP1 stimulates telomerase repeat-addition processivity, and what it changes mechanically is
unresolved. This paper asks what a processivity distribution can *identify* about that mechanism,
and answers: less than is commonly assumed. Three results are exact within a stated model class —
an extension-rate mechanism and a retention mechanism are the same function of the observable, no
absolute rate is recoverable at any precision, and a processivity design has structural rank 3 of 4.
A fourth result corrects an error in the previous draft of this work: reducing *k*_off and raising
*k*_on at matched *K*_d are **not** equivalent; the loss term rewrites as *R* = *k*_off²/(*k*_off +
[P]*k*_on), so they differ by an exact factor λ that is confounded with *k*_loss. The paper then
gives a closed form for the adjacent-band ratio of a two-species mixture, proves it is bounded by
the most extreme subpopulation, shows the tail slope converges to the **longest-lived** species, and
derives the length at which a minority contaminant takes over the signal. That last result has a
practical consequence: the same gel yields an apparent processivity anywhere from 0.61 to 0.95
depending only on which band window is read, a 12-fold change in the implied loss ratio. The tail
slope is therefore forbidden as the basis for any mechanistic claim.

## What the paper does not claim

No mechanism is demonstrated. Every result is conditional on a kinetic scheme, and the paper states
that boundary explicitly (§8): *a mathematically exact result within a model is not an
experimentally established biological fact.* The steps that would convert the mathematics into
evidence are listed as release gates (§11), and the next one — real per-band densitometry of an
existing gel — needs no new experiment and can falsify the practical claim.

## Status and provenance

The analysis is committed; the empirical test is not done. Source values for the published
processivity data used in the repository are in [`data/external/tesmer2023/`](../data/external/tesmer2023/)
with SHA-256 provenance. The underlying analysis document is
[`results/causality/MECHANISM_CLASSIFICATION.md`](../results/causality/MECHANISM_CLASSIFICATION.md),
frozen by hash and carrying its own retraction record.

*This package is not a journal submission and has not been peer reviewed. It is deposited in the
repository so that the arithmetic can be checked independently before anyone acts on it.*
