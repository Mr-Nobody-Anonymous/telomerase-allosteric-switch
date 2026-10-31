# Work package 4 — decision matrix for the next experiments

The purpose of this matrix is to rank candidate experiments by their ability to **change our
mind**, not by how interesting they are or how well they show off existing machinery. An
experiment is only worth doing if the three models in `model_comparison.md` predict different
outcomes for it and if the difference is larger than the noise of the assay.

Columns: what is perturbed · what is measured · what each model predicts · whether the
prediction is sharp · whether it can be done with existing reagents · decision.

## 1. Ranking

| # | Experiment | Perturbation | Readout | A predicts | B predicts | C predicts | Discriminates? | Cost | Decision |
|---|---|---|---|---|---|---|---|---|---|
| **E1** | **ΔRAP re-analysis of published gels** | none — re-quantify existing processivity gels | fold-stimulation with/without POT1–TPP1, per mutant | the stimulation requires the conduit, so mutants in/near 766 lose ΔRAP even if basal processivity survives | a signal mutant should keep ΔRAP while losing catalytic rate | ΔRAP tracks DNA affinity: mutants with weakened anchor lose absolute processivity but keep or gain fold-stimulation (rescue) | **yes** — this is the only readout that separates "being stimulated" from "being good" | zero (existing data) | **do first** |
| **E2** | **Primer-affinity phenocopy** | raise primer concentration / use a higher-affinity primer, no TPP1 | RAP of TEN `PLYQ` and β17 mutants | no effect: the defect is in the signal, not the affinity | no effect | **phenocopies TPP1**; the mutant should recover RAP without TPP1 | **yes** — clean and binary | low (existing constructs) | **do second** |
| **E3** | **H762E as the required negative control in any conduit panel** | H762E against a fresh cohort of conduit candidates | processivity, ΔRAP, TPP1 binding | H762E must be defective if 762 is in the conduit | silent | silent (surface residue, no partner) | **no** — it is a control, not a test | zero | **mandatory control, never counted as evidence for A** |
| **E4** | TPP1 separation-of-function panel | TPP1 OB-surface variants (TEL patch / E169-K170-L95 region) that bind TERT normally | TPP1 binding (BLI/ITC) + ΔRAP | a TPP1 mutant that binds but does not stimulate revives A | same — a binding-competent, stimulation-dead TPP1 revives B | no such mutant should exist | **yes, if one is found** | medium (protein engineering) | do after E1/E2 |
| **E5** | Conformational readout at the active site | ± TPP1, ± DNA | HDX-MS or smFRET between IFD-TRAP and the catalytic site | local change confined to the anchor region | change at the active site, TPP1-dependent | no TPP1-dependent change beyond TEN stabilisation | **yes** for B | high (new instrumentation) | only if E1–E4 say "signal" |
| **E6** | Conduit mutagenesis sweep (760–766) | alanine/charge-reversal panel | activity, ΔRAP | graded loss from 766 towards 760 | non-specific | no pattern; defects track nucleic-acid proximity | **weak** — already confounded by K760E vs H762E | medium | **do not do** |

## 2. Why the panel changed

* **H762E is now a required negative control, not a candidate.** It is inside the originally
  proposed conduit and is near-WT. Any future conduit panel that reports defects at 761, 763,
  764, or 765 without an H762E comparator is uninterpretable, because the comparator proves the
  stretch tolerates substitution.
* **762 is removed from the conduit hypothesis entirely.** WP1 shows 762 is 13 Å from TPP1 and
  760–765 make no contact at all; WP2 shows the mutation is the mildest in its panel. Keeping it
  in a candidate list would be motivated reasoning.
* **L766 remains the only defensible candidate residue in the window**, because it is the only
  one at the interface — and it sits exactly at the strand→turn transition, so its phenotype will
  always be confounded between "interface" and "local structure". It should be treated as a
  boundary marker, not as a pathway node.
* **The original 40-second workflow (build a network, read off a coupling change) is dropped.**
  WP3 §4 explains why: the effect size equals a short-range positive control, and the
  bound/unbound axis does not exist in the coordinates.

## 3. Statistical design (applies to E1–E5)

* **Unit of replication: independent transfections (or independent protein preparations), not
  gel lanes and not MD frames.** Both published datasets use n = 3 independent experiments with
  SEM; that is the standard to match, not to exceed.
* **Effect size to detect.** The published phenotypes are large (24% vs 100% activity; 86% vs
  ~100% processivity), so a difference of ≥1.5× in ΔRAP is the practical threshold. With
  triplicate SEMs of ~10% of the mean (as reported), n = 3 resolves a 1.5× difference only if
  the SEMs are at the optimistic end; **n = 4–5 with the same assay is the minimum defensible
  design**, and the analysis should be a ratio of paired measurements within experiment, not a
  comparison of two independently normalised means.
* **Normalisation.** Report ΔRAP as (processivity +POT1–TPP1)/(processivity −POT1–TPP1) per
  experiment, so that the mutant's own baseline cancels. This is what makes the readout
  separable from basal activity — the quantity that both published papers leave implicit.
* **What counts as a failure to replicate.** If H762E fails to reproduce as near-WT in a new
  cohort at n = 4, the whole comparison collapses and the correct response is to stop and
  re-examine the assay, not to reinterpret the model.
* **Blinding.** Gels should be quantified blind to genotype; the processivity metric is a
  ratio of band intensities and is sensitive to band-assignment choices.

## 4. Decision rule

Run E1 and E2 first. Then:

* If **E2** shows that raising primer affinity phenocopies TPP1 for anchor-site mutants, Model C
  is confirmed at the level of mechanism-as-affinity and there is no signal to find. **Stop.**
  The remaining question becomes a thermodynamic one (what state does TPP1 stabilise), which is
  a structural-biology question, not a signalling one.
* If **E2** shows no phenocopy and **E1** shows a mutant that keeps basal processivity but loses
  ΔRAP, a signal exists and Model A/B is alive — go to E4 and E5, in that order.
* If **E1** and **E2** are both ambiguous (differences below the assay's resolution), the honest
  conclusion is that the current assays cannot answer the question and E4/E5 are the only route.

## 5. What this matrix deliberately excludes

* **More network models / ANM / coarse-grained dynamics on the same two deposits.** They cannot
  discriminate (WP3 §4) and the deposits do not contain the bound/unbound pair they assume.
* **More TERT point mutants in the IFD-TRAP.** The phenotype space is already saturated with
  confounded mutants; another one adds no information.
* **A 40-second trajectory-based animation.** It would be a presentation artefact, not a
  measurement, and the repository's existing 3 ps trajectories cannot support any claim about a
  process that takes seconds.
* **Any claim ranked by how convincing a figure looks.** The ranking above is by discriminating
  power and cost, and the top two entries require no new simulation and almost no new
  experimental work.
