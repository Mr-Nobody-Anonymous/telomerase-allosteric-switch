# Work package 3 — the three models, scored against the existing evidence

This document does what the earlier analysis in this repository did not: it states each model
before looking at the data, lists what each model predicts, and then checks the predictions
against facts that were measured independently of the model. **No new simulation was run for
this document.** The point is precisely that the models can be ranked without one, and that the
existing mechanical calculation cannot do the ranking.

---

## 1. The models, stated so that they can fail

**Model A — local structural conduit.**
TPP1 binding at the IFD-TRAP interface converts the anchor-site region into a state that holds
the primer longer, and the information travels along the polypeptide from the TPP1 footprint
towards the DNA anchor (the earlier version named 760–766 as the proximal segment).
*Required:* residues between the TPP1 footprint and the anchor must be at the interface or
mechanically contiguous with it, and mutating them must break the TPP1-dependent stimulation.

**Model B — broad conformational allostery.**
TPP1 binding propagates through the TERT body to the active site and changes the catalytic
cycle itself (e.g. translocation efficiency), without the change being local to the anchor.
*Required:* a measurable conformational difference at the active site that depends on TPP1 and
correlates with processivity, in a region not in contact with TPP1.

**Model C — affinity/anchoring, no propagation required.**
TPP1 does not send a signal. It binds TERT and stabilises the TEN domain in a state that holds
the DNA, raising the effective affinity of the enzyme for the primer; processivity rises because
dissociation becomes rarer. Any "allostery" is a thermodynamic shift in the TEN ensemble, not a
propagating pathway.
*Required:* (i) TPP1 must not contact DNA; (ii) the effect must be substitutable by raising DNA
affinity; (iii) mutations that weaken DNA binding should be *rescued* by more TPP1, not
unaffected by it.

---

## 2. The observations to be predicted

Facts from WP1 (recomputed from coordinates) and WP2 (published experiments):

* **O1** TPP1 makes no contact with DNA (28.7 Å) or hTR (14.2 Å); it contacts the TEN helices,
  IFD β17, and IFD helix α25.
* **O2** TPP1's closest TERT contact is R774 (2.3 Å) on helix α25; the DNA anchor residues
  K757/F759 lie on a **β-strand**, and 768–771 is a loop with P771. There is no continuous helix
  between the anchor and the footprint.
* **O3** Of residues 760–766, only L766 contacts TPP1; 760–765 are 12–18 Å away.
* **O4** H762E (in the middle of the named conduit) keeps 86% of WT processivity with
  POT1–TPP1 and 68% activity. K760E, three residues away, is severely defective.
* **O5** Adding the TPT complex **rescues** the RAP of anchor-site mutants (TEN `PLYQ`, β17,
  Q794) that have weakened DNA binding.
* **O6** Processivity-defective IFD-TRAP mutants (K749E, R756E, F759E) are recruited to
  telomeres at WT levels; D788R keeps near-WT stimulated processivity while retaining only 24%
  activity.

---

## 3. Predictions

Notation: "∅" = the model makes no commitment; ✓ = prediction confirmed; ✗ = contradicted;
? = untested by existing data.

| Observation | Model A predicts | Model B predicts | Model C predicts |
|---|---|---|---|
| **O1** TPP1 far from DNA | must contact DNA or a DNA-binding element closely — a conduit has to end at the thing it changes | no commitment; consistent | ✓ **required**: the effect must be mediated by TERT's own DNA affinity |
| **O2** strand → turn → loop → helix, no continuous helix | ✗ a *sequential* conduit through one helix is impossible; only a junction mechanism survives | ✓ fine: any route through the body is allowed | ✓ fine: no route needed |
| **O3** only 766 of 760–766 contacts TPP1 | ✗ the named proximal segment is not at the interface | ? | ✓ fine |
| **O4** H762E near-WT, K760E severe | ✗ the middle of the conduit is dispensable; the model has to explain why the *distal* residue matters more than the interface-proximal ones, which is backwards | ? | ✓ fine: H762 is a surface residue not proximal to any partner; K760E can act on nucleic acid directly |
| **O5** TPT *rescues* anchor mutants | ✗ a missing signal cannot be restored by adding more of the same ligand; rescue contradicts a required pathway | ? rescue is not predicted by a pathway model but is not absolutely forbidden | ✓ **the signature prediction**: mass action compensating a weakened binding equilibrium |
| **O6** recruitment separable; D788R dissociates activity from stimulation | ? (A is about within-enzyme coupling, silent here) | ? | ✓ fine: the stimulation tracks DNA retention, not catalytic rate or telomere targeting |

Score: **A: 0 confirmed, 4 contradicted. B: 0 confirmed, 0 contradicted, 4 unaddressed.
C: 5 confirmed, 0 contradicted, 1 unaddressed.**

Model B's row of unaddressed entries is the important one. It is not supported; it is also not
refuted, because no existing experiment measures a TPP1-dependent conformational change at the
active site. Anyone who reports that "allostery has been excluded" is over-reading these data.

---

## 4. What the repository's own mechanical analysis can and cannot add

The allosteric-coupling calculation in `results/allostery/` (anisotropic network model on 7TRE,
BOUND vs APO vs a displaced DECOY built from the same coordinates) reported:

| quantity | value | reading |
|---|---|---|
| specific TEN↔IFD-TRAP coupling change | −0.041 (18 Å: −0.051) | small |
| z-score vs distance-matched null | +0.34 | indistinguishable from the null |
| z-score for the displaced DECOY | −1.78 | not a specific effect |
| short-range positive control | +0.041 | the same magnitude as the "effect" |
| catalytic site ↔ 3′ DNA | 0.103 vs 0.087 background | no separation from background |

Two conclusions follow, and only two.

1. **The calculation cannot adjudicate between the models.** Its "specific" effect is the same
   size as a short-range positive control and is within a third of a standard deviation of a
   distance-matched null. Reporting it as evidence for or against allostery would be a category
   error.
2. **The BOUND/APO axis does not exist.** WP1 §6 shows 7TRE and 7TRD differ by 0.36 Å over 934
   Cα — the same conformation. Building an "unbound" state by deleting TPP1 from a bound
   structure, and then reading mechanical differences off that deletion, measures the deletion,
   not the physics.

There is a further structural reason the calculation was doomed, which is new in WP1: it was
built to test propagation along a helix that is not there (O2). An ANM can only propagate motion
through the bonded network that actually exists, so the calculation inherited the error in the
hypothesis it was testing.

Consequence for the milestone: **the ANM result is recorded as inconclusive and is retired as
evidence.** It stays in the repository because removing it would hide a negative result, but no
model ranking in this document depends on it.

---

## 5. Verdict of WP3

1. **Model A, as stated (760–766 local conduit), is dead** — not by simulation but by two facts
   that were available in the literature: only 766 of the window is at the interface, and H762E,
   in the middle of the window, is the mildest mutant in its panel. A salvage version exists
   ("junction/hinge from helix α25 through the P771 loop to β17"), but it makes no sharp
   prediction, is not distinguished from Model C by any existing experiment, and should not be
   pursued without new data.
2. **Model C is the best-supported model** and is the only model that predicted the rescue
   result before seeing it. It is also the model the original publications themselves argue for.
3. **Model B is neither supported nor excluded.** Its only quantitative literature support is
   intra-TPP1 microsecond MD of the OB domain (Commun Chem 2023), which shows plasticity of the
   TPP1 surface, not coupling to the TERT active site. Any future claim for B needs a
   TPP1-dependent active-site measurement, not another network model.
4. **The ranking is not a discovery.** It is the published interpretation, recovered by refusing
   to let a structural gradient stand in for a mechanism. The original contribution here is the
   correction: the helix-conduit premise was wrong on secondary structure (O2/O3), and the
   mutation that was nominated to break the conduit is the one mutation that does not.
