# Work package 5 — mechanism assessment, ranked by confidence

This is the closing report for the mechanism-discovery audit. It answers one question and
deliberately does not answer a larger one.

> **Mission.** Determine whether TPP1 stimulates telomerase processivity through (A) a local
> IFD-TRAP structural conduit (760–766), (B) broader conformational allostery, or (C)
> recruitment/avidity/primer positioning requiring no mechanical propagation — without assuming
> that any of them is correct.

**Bottom line.** Model (C) is the best-supported explanation and Model (A) as stated is
contradicted by published experiments. Model (B) is untested rather than excluded. Separability
of TPP1-dependent processivity from basal catalysis **can** be established from existing data;
separability of a *mechanical* stimulation from DNA-affinity compensation **cannot**, because no
published perturbation changes one without the other. That gap is the deliverable: it converts an
open-ended modelling programme into a short list of two discriminating experiments
(`experimental_decision_matrix.md`).

---

## 1. Findings ranked by confidence

### 1.1 Certain — recomputed from coordinates, reproducible (`structural_validation.md`)

| # | Finding | Evidence |
|---|---|---|
| F1 | TPP1 does not touch DNA (28.7 Å) or hTR (14.2 Å). Its only route to the primer is through TERT. | min heavy-atom distances, 7TRE |
| F2 | Of the proposed conduit window 760–766, **only L766** is at the TPP1 interface (3.6 Å). Residues 760–765 are 12–18 Å away; five of them exceed 7 Å. | 4.5 Å contact census |
| F3 | **The anchor and the TPP1 footprint are not on one helix.** K757/F759 sit on an extended β-strand; 766–767 is a single turn; 768–771 is a loop containing P771; 772–780 is helix α25 carrying the closest contact R774 (2.3 Å). | backbone dihedrals + CA(i→i+4) distances |
| F4 | 7TRD is not an unbound state: TERT is the same conformation as 7TRE to 0.36 Å over 934 Cα. | Kabsch superposition |

F3 is a retraction. Earlier write-ups in this repository described 756–775 as one continuous
IFD-TRAP helix, and the "helix conduit" hypothesis was built on that description. The description
is wrong: the distance gradient along the sequence (26 Å at 756 → 2.3 Å at 774) is real, but it
runs across a strand→turn→loop→helix junction, so it cannot be a helix-propagation pathway.

### 1.2 High confidence — from published experiments (`experimental_evidence.md`)

| # | Finding | Evidence |
|---|---|---|
| F5 | **Recruitment is separable from processivity.** IFD-TRAP mutants with severe activity/processivity defects are recruited to telomeres at WT levels. | Tesmer 2023, Fig S2B/C (n.s.) |
| F6 | **H762E, inside the named conduit, is near-WT**: 86% of WT processivity with POT1–TPP1, 68% activity. It is the mildest mutant in its panel. | Tesmer 2023, Fig 3A–C |
| F7 | **Adding TPT rescues** the RAP of mutants whose DNA anchor is damaged (TEN `PLYQ`, β17, Q794) rather than failing to stimulate them. | Sekne 2022, Fig 4E–G |
| F8 | **Basal catalysis and stimulated processivity dissociate**: D788R has 24% activity but 64% WT-like stimulated processivity. | Tesmer 2023 |

### 1.3 Moderate — inference, stated as such

* **F9. TPP1 acts by stabilising the TEN/DNA interaction (Model C).** Supported by F1, F5, F7,
  F8 together, and it is the interpretation the structural paper itself gives. It is an
  inference, not a measurement: nobody has measured the DNA on-rate/off-rate with and without
  TPP1 in the same experiment at the single-molecule level.
* **F10. The 760–766 mutant phenotypes are not evidence about a conduit at all.** K760E is
  severe and H762E is mild, which is the opposite of the ordering a conduit requires; the authors
  attribute K760's defect to a possible transient nucleic-acid interaction and cannot assign
  H762's mildness to anything except that it is a surface residue with no partner.

### 1.4 Not established — explicitly open

* **F11. Broad conformational allostery to the active site (Model B) is neither supported nor
  excluded.** No published experiment measures a TPP1-dependent conformational change at the
  catalytic site. The only quantitative support often cited for B is a µs molecular-dynamics
  study of the **isolated TPP1 OB domain**, which shows plasticity of TPP1's own surface — not
  coupling to TERT. It should not be cited as evidence for TERT-wide allostery.
* **F12. Whether a mechanical component exists at all is undetermined**, and (see §3) is not
  decidable with the data in hand.

### 1.5 Retired — the repository's own computation adds nothing here

* **F13. The anisotropic-network-model coupling calculation is inconclusive and is withdrawn as
  evidence.** Its "specific" effect (−0.041) equals the short-range positive control (+0.041) and
  is +0.34 σ from a distance-matched null. The BOUND/APO axis it relies on does not exist (F4),
  and it was built to propagate along a helix that does not exist (F3).

---

## 2. Direct answer to the separability question

**Question:** can separability of TPP1-dependent processivity from basal catalysis and
recruitment be established from existing experimental data?

* **From recruitment: yes.** F5 is clean, replicated in one lab, and covers several mutants.
  Recruitment is a separate function; it is not the mechanism of the test-tube stimulation. Any
  model that leans on "TPP1 recruits telomerase, therefore processivity rises" is wrong.
* **From basal catalysis: partly.** F6 and F8 dissociate the two readouts in both directions of
  severity. What is missing is the **fold-stimulation (ΔRAP)** per mutant, which neither paper
  reports in the text. Without it, "kept 86% of processivity" cannot be decomposed into "stimulated
  normally" versus "was already processive".
* **From a mechanical signal: no.** Every published perturbation changes DNA affinity and TPP1
  binding simultaneously. There is no separation-of-function allele on either side of the
  interface.

So the answer is: **the separability that matters has not been established, and it is
establishable** — but by two cheap experiments (ΔRAP re-analysis; primer-affinity phenocopy), not
by more simulation.

---

## 3. Why more simulation cannot settle this

1. The two available deposits are the same conformation (F4), so any "bound vs unbound" contrast
   is a computational artefact of deleting a chain.
2. The structural premise the simulations were built on is wrong (F3).
3. The repository's trajectories are 3 ps. Allosteric and processivity-relevant motions in a
   megadalton RNP occur over µs–ms. A 3 ps trajectory cannot sample them, and no amount of
   analysis of a 3 ps trajectory recovers information it never contained.
4. The network model already returns an effect indistinguishable from a matched null (F13).

Continuing the simulation programme would produce a fifth document and no new knowledge.

---

## 4. The four things this repository must keep distinct

The earlier documentation mixed these together, which is how a distance gradient became a
mechanism. They are listed here once so that no future summary can re-merge them:

1. **The original bipartite TEN–TPP1–TRAP clamp hypothesis.** TPP1 and the TEN domain were
   proposed to hold the DNA between them, with a long-range allosteric switch coupling TPP1
   binding to the anchor site. *Status: refuted on its structural premise* — the TEN domain makes
   no primer contact at all (0 residues within 6 Å; Y176 is 11.1 Å away), and TPP1 is 28.7 Å from
   the DNA.
2. **The structural audit that invalidated that premise.** Independent recalculation from the
   coordinates (WP1; the earlier `results/verification/` audit reached the same conclusion for the
   interface). *Status: stands, and is now extended* — with the new secondary-structure correction
   (F3), the audit also invalidates the "one continuous helix" description that the follow-up
   hypothesis used.
3. **The conduit hypothesis built on top of the audit (760–766).** *Status: not supported* —
   only 766 is at the interface, and H762E inside the window is near-WT. Retired, with one
   defensible residue (L766) retained as a boundary marker, not as a pathway node.
4. **The mechanical (ANM) result.** *Status: inconclusive.* It is not evidence against allostery
   and must never be cited as such. It is also not evidence for it.

---

## 5. Claims rejected or overstated during this work

| Claim | Why it is rejected |
|---|---|
| "TPP1 sits on the same IFD-TRAP helix as the DNA anchor, so a helix conduit is plausible." | Wrong secondary structure (F3). The two elements are on opposite sides of a strand→turn→loop→helix junction. |
| "Mutations at 760–766 specifically break the TPP1 response." | H762E is near-WT; K760E has no assigned contact; no residue in the window has a TPP1-response-specific phenotype. |
| "The ANM shows allosteric coupling is absent." | Effect equals the positive control and is within noise of a matched null (F13). |
| "TPP1 fails to stimulate processivity when the anchor is mutated." | The opposite is published: TPT **rescues** those mutants (F7). |
| "Recruitment/avidity explains the stimulation." | Avidity in the cell is separable from the test-tube effect (F5). The best-supported version of C is *affinity compensation inside the enzyme*, not telomere targeting. |
| "Microsecond MD of the TPP1 OB domain supports TERT-wide allostery." | It shows plasticity within TPP1, not coupling to the TERT active site (F11). |
| Any claim supported by an "allosteric switch" distance of 25 Å. | That premise died with the first audit; R774 is 2.3 Å from TPP1. |

---

## 6. Standing hypothesis, its falsification criteria, and the next steps

**Leading hypothesis (Model C, restated narrowly).** TPP1 raises processivity by stabilising the
TEN domain in a DNA-binding-competent state, increasing the effective affinity of telomerase for
the primer and reducing dissociation between repeats. No propagating signal is required, and any
conformational change involved is a thermodynamic shift in the TEN ensemble rather than a
pathway to the active site.

**This hypothesis is falsified if:**
1. a TPP1 variant binds TERT with WT affinity yet does not stimulate RAP at any concentration
   (would require a signal), **or**
2. raising primer affinity (no TPP1) fails to phenocopy TPP1 for anchor-site mutants that TPT
   rescues (would break the compensation account), **or**
3. a TPP1-dependent conformational change is detected at the catalytic site by a method that
   resolves dynamics (HDX-MS, smFRET, NMR) and tracks ΔRAP quantitatively.

**Next steps, in order** (`experimental_decision_matrix.md`): re-quantify existing gels for
ΔRAP (E1); test the primer-affinity phenocopy (E2); keep H762E as a mandatory negative control in
any conduit panel (E3); only then engineer TPP1 separation-of-function variants (E4) or attempt a
conformational readout (E5). No further network modelling.

**What this audit does not claim.** It does not claim that TPP1 binding causes no conformational
change — stabilising a domain is a change. It does not claim that allostery is absent, only that
no existing experiment separates a mechanical component from an affinity component. It does not
claim that Model C is proven; it claims that C is the only model the data currently support and
that the two experiments which could overturn it are cheap and have not been done.
