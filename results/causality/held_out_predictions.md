# Held-out directional validation against Sekne et al. 2022

**Test set used once, qualitatively, exactly as frozen.** Predictions P1–P6 were committed in
`PREREGISTRATION.md` §2 (hash `050f8d1f…`, guarded by `tests/test_causality_preregistration.py`)
before any Sekne data were examined for this investigation. Sekne was never a fitting target; the
development data are Tesmer only.

**Source consulted:** Sekne et al., *Science* 375:1173 (2022), text of Fig 4E–G and Results.
No figure was digitised; no number is used that is not stated in the text or in the paper's
sentence-level summary ("addition of the TPT complex rescued the RAP defects").

---

## 1. Prediction-by-prediction

| # | Prediction (frozen) | Held-out observation | Verdict |
|---|---|---|---|
| **P1** | anchor-site mutants have reduced processivity with TPT vs WT | Y176A/Q177A, K757A/F759A, Q794A all "significantly reduced" activity and RAP | **TRUE** — confirms C, does not discriminate (all three models predict a defect) |
| **P2** | C: TPT **rescues**; B: rescue not predicted; A: not rescued | "Addition of the TPT complex **rescued** the RAP defects in these mutants" | **TRUE for C.** B: unhelpful rather than false (it does not forbid rescue). **A: FALSE** — a broken conduit is not bypassed by more ligand. |
| **P3** | magnitude of rescue tracks magnitude of the DNA-contact defect | Sekne reports rescue for all three mutants but does **not** report per-mutant rescue magnitudes; the figure is a bar chart with SEM in triplicate | **UNTESTABLE** — the paper states direction, not ordering. Testing this would require digitising bars against a stated normalisation, which §3 of the preregistration forbids. |
| **P4** | mutants that keep TPP1 binding but damage DNA contact are still stimulable | The anchor mutants are DNA-contact lesions, not interface lesions, and TPT still raises their RAP | **TRUE for C** (stimulation needs TPP1 binding, not an intact DNA contact) |
| **P5** | folding-lesion mutants are uninformative about either model | Not in Sekne's panel (they are Tesmer's F776E/L786E, development data) | **NOT APPLICABLE** — correctly declared uninformative in advance |
| **P6** | no observation in this dataset separates retention from translocation rate | Confirmed: Sekne reports no primer titration, no single-molecule measurement, no active-site conformational readout | **TRUE** — and this is a statement about the *absence* of the discriminating data, not support for either model |

**Score: 3 TRUE (P1, P2-C, P4), 1 FALSE (P2-A), 1 UNTESTABLE (P3), 1 NOT APPLICABLE (P5),
1 TRUE-as-an-absence (P6).**

## 2. What this does and does not establish

**Established:**

* The rescue observation is real, is replicated across three independent mutant classes in an
  independent lab, and **contradicts Model A** — the conduit cannot be rescued by ligand, because a
  conduit is a physical connection and adding more TPP1 does not restore a broken one.
* The observation is **exactly what an affinity/compensation account predicts** and is not
  predicted by a pathway account. That is the only asymmetric prediction in the set (P2).
* **No existing observation separates affinity from conformational regulation** (P6). Sekne
  supplies none of the observables that would.

**Not established:**

* **The asymmetric result is one observation in one direction.** A pathway model can survive P2 by
  adding a compensatory step ("more TPP1 stabilises TEN enough to offset the lost contact") — that
  is no longer a *conduit*, but it is still not pure affinity. P2 kills the 760–766 conduit, not
  every possible allosteric model.
* **P3 — the only prediction that could have actively favoured C over B — is untestable in this
  dataset.** So C does not win on predictive asymmetry beyond P2; it survives, and it is the only
  model that predicted the direction of an observation before it was seen.
* Sekne's ΔTEN *trans*-complementation, 500 nM POT1–TPP1 and reconstituted-vs-TPT difference from
  Tesmer's system means **no magnitude is comparable across the two datasets**. Only signs and
  orderings were compared, per the non-comparability clause.

## 3. The prediction that failed, reported as a failure

**P2 for Model A failed.** Model A is already contradicted on structural grounds
(`../mechanism/model_comparison.md`, `../mechanism/structural_validation.md`). The held-out test
fails it a second time, on independent experimental grounds, and in the direction the
preregistration committed to in advance. This is the cleanest result in the causality
investigation: **it was predicted, it was checked against an untouched dataset, and the model
lost.**

## 4. Verdict

* **Model A (local conduit): rejected**, twice, once structurally and once experimentally, the
  second time by a pre-registered directional prediction.
* **Model C (affinity/anchoring): supported, not confirmed.** It predicted the rescue direction
  correctly. It was not tested against the one prediction that could have distinguished it from
  Model B (P3), because that observation does not exist in the source.
* **Model B (conformational): untouched.** Nothing in this test set bears on it. The correct
  description remains *unresolved*, not *excluded*.
