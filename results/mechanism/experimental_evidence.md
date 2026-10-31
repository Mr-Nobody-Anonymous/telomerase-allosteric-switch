# Work package 2 — what the existing experiments already establish

**Question this document answers.** Is TPP1-dependent processivity stimulation separable, in
existing published data, from (a) basal catalysis and (b) recruitment of telomerase to telomeres?
Separability is the precondition for the whole mechanism question: if the three cannot be
separated experimentally, then no amount of simulation can identify the mechanism, and the next
experiment should be a separability experiment rather than another structural model.

**Provenance.** Every number below is taken from the published text of the cited papers. **No
figure was digitised for this document**, so nothing here needs the "digitised" caveat. Where a
number appears only in a figure and not in the text, it is marked *figure-only* and is not used
quantitatively. Independent recomputation of the structural facts quoted here is in
`structural_validation.md` (WP1).

**Statistical unit.** In both papers the unit is an independent transfection/primer-extension
experiment (n = 3), reported as mean ± SEM. Individual gel bands and individual cryo-EM frames
are not independent observations and are not counted here. Neither paper reports a formal power
calculation, and one of them states explicitly that it could not reach significance for activity
differences because of signal-to-noise (below). All conclusions in this document are therefore
"consistent with", not "demonstrated at a stated effect size".

---

## 1. The two datasets that matter

| | Sekne et al. 2022, *Science* 375:1173 | Tesmer et al. 2023, *Life Sci Alliance* 6:e202201727 |
|---|---|---|
| What was varied | TERT DNA-anchor residues: Y176A/Q177A (TEN `PLYQ`), K757A/F759A (IFD β17), Q794A (IFD) | Nine conserved IFD-TRAP residues, singly, as Glu (or D788R) |
| Assay | Primer extension with purified/reconstituted telomerase, ± the TPT complex (TERT–POT1–TPP1) | Primer extension in a *trans*-complementation system (FLAG-TERTΔTEN + Myc-TEN), ± 500 nM POT1–TPP1 |
| Readouts | Activity; RAP; effect of adding TPT | Activity; POT1–TPP1-stimulated processivity; recruitment by IF-FISH in full-length TERT |
| Recruitment measured? | No | Yes (full-length TERT, K749E/R756E/F759E and the IFD-TRAP panel) |
| Processivity metric | Repeat-addition processivity from the gel | Fraction of signal in products ≥16 repeats (main text); ≤12-repeat alternative metric in Fig S2E/F |
| Replication | Triplicate, SEM | Three independent experiments, SEM; some activity comparisons not significant |
| Structures | 7QXA/7QXB/7QXS; 7TRE is the related Liu et al. 2022 deposit | PDB 7QXA reused |

Both datasets use the same structural model of the interface, so they are not independent with
respect to *where* the residues are. They are independent with respect to the *phenotype*,
because they were done in different labs with different constructs.

---

## 2. The measurements, as published

| Variant | Basal activity vs WT | Processivity with POT1–TPP1 | Recruitment to telomeres | Interpretation offered by the authors |
|---|---|---|---|---|
| Y176A/Q177A (TEN `PLYQ`) | reduced | reduced; **RAP defect rescued by adding TPT** | not tested | anchor-site DNA contact lost; TPT compensates reduced DNA affinity |
| K757A/F759A (IFD β17) | reduced | reduced; **RAP defect rescued by adding TPT** | not tested | as above |
| Q794A (IFD) | reduced | reduced; **RAP defect rescued by adding TPT** | not tested | as above |
| K749E | severe defect | severe defect | **WT-like (n.s.)** | TR template contact (π–cation to C56) |
| R756E | severe defect | severe defect | **WT-like (n.s.)** | TR template contact (π–cation to A55) |
| K757E | severe defect | severe defect | WT-like (Fig S2B/C, n.s.) | DNA backbone contact |
| F759E | severe defect | severe defect | **WT-like (n.s.)** | DNA base stacking |
| K760E | severe defect | severe defect | not tested | **author cannot assign it**: "not close to TR or the primer … may be indirect or … a transient nucleic acid interaction" |
| **H762E** | **68% of WT** | **86% of WT with POT1–TPP1** | WT-like | polar surface residue, "not proximal to any of its binding partners" |
| F776E | ≈ dead | ≈ dead | WT-like | hydrophobic core of IFD-TRAP → misfolding |
| L786E | ≈ dead | ≈ dead | WT-like | hydrophobic core of IFD-TRAP → misfolding |
| D788R | **24% of WT** | **64% of WT-like stimulation** | WT-like | weakened TEN–(IFD-TRAP) interface |

Sources: Sekne et al. 2022 Fig 4E–G and text; Tesmer et al. 2023 Fig 3A–C, Fig S2B/C, text.
"Rescued by TPT" is the authors' wording: adding the TPT complex restored the RAP of the
anchor-site mutants, and they attribute this to compensation of reduced DNA affinity.

---

## 3. What is now established, and what is not

### 3.1 Recruitment is separable from catalysis and from processivity — **established**

Three IFD-TRAP mutants with severe activity and processivity defects (K749E, R756E, F759E) are
recruited to telomeres at WT levels in full-length TERT (Tesmer 2023, Fig S2B/C; n.s.), and the
wider IFD-TRAP panel likewise shows no significant recruitment defect. Recruitment therefore
does not require an intact IFD-TRAP, and a processivity-defective telomerase can still be
recruited normally.

This kills one of the models outright, and it is worth being blunt about it: **any model in which
TPP1 stimulates processivity by recruiting telomerase to the telomere is contradicted** by the
observation that processivity can be destroyed while recruitment is untouched. Recruitment
(avidity for telomeric DNA in the cell) is a *separate* function, not the mechanism of the
stimulation seen in the test tube.

### 3.2 Basal catalysis is separable from TPP1-stimulated processivity — **partially established**

D788R retains only 24% of WT activity but nearly WT-like POT1–TPP1-stimulated processivity
(64%). H762E retains 68% activity and 86% processivity. So the two readouts dissociate in both
directions of severity, and a telomerase that is bad at catalysis is not automatically bad at
being stimulated.

The caveat is that these are *normalised* processivity values measured in the presence of a
saturating concentration of POT1–TPP1. They show that stimulation still occurs; they do not
measure the fold-stimulation (ΔRAP = processivity with POT1–TPP1 ÷ processivity without), which
is the quantity a mechanistic model has to explain. Neither paper reports ΔRAP per mutant in the
text. **This is the single most important gap in the existing data.** See §5.

### 3.3 The proposed 760–766 conduit is not supported — **contradicted**

* Structurally, of residues 760–766 only **L766** is at the TPP1 interface (3.6 Å); 760–765 are
  12–18 Å from TPP1 and share no contact (WP1 §3).
* Genetically, the two residues in that window that have been mutated point in opposite
  directions from what a conduit requires: **K760E is severely defective** and **H762E is
  essentially WT** (86% of WT processivity with POT1–TPP1). If 760–766 were a conduit whose
  integrity is required to transmit the TPP1 signal, then H762E — a charge reversal in the middle
  of it — should have been the more damaging of the two mutations. It is the least damaging
  mutation in the entire panel.
* Tesmer et al. report that H762 is a polar surface residue not proximal to any binding partner,
  and cannot assign K760's defect to any observed contact.

The honest reading is that these two mutations do not probe a conduit at all: H762E is a
negative control showing that this stretch tolerates substitution, and K760E is a defect of
unknown origin that Tesmer explicitly declines to attribute to a structural pathway.

### 3.4 TPP1 stimulates by raising effective DNA affinity — **best-supported, but not proven**

The published mechanistic statement in Sekne et al. is that, in telomerase alone, the TEN domain
is conformationally flexible, which makes DNA binding unstable and increases dissociation during
repeat synthesis; TPP1 binding stabilises TEN, and POT1 stabilises the DNA–TEN interaction, so
RAP rises because DNA dissociates less often. Two independent observations fit this:

1. **The rescue result.** If TPP1's effect were a signal transmitted through TERT, then raising
   the concentration of the TPT complex should not repair a defect located in TERT's DNA anchor —
   the signal would still be absent. What is observed is the opposite: adding TPT repairs the RAP
   of mutations that damage DNA binding (TEN `PLYQ`, β17, Q794). That is the signature of a
   compensating *affinity* effect (mass action on a weakened binding equilibrium), not of a
   missing signal.
2. **The distance result.** TPP1 is 28.7 Å from the modelled primer and makes no contact with
   DNA or hTR (WP1 §7), so it cannot stabilise the hybrid directly; its only available route to
   DNA retention is through TERT's own affinity for DNA.

Note carefully what this does **not** say. It does not say that no conformational change happens
when TPP1 binds — stabilising the TEN domain *is* a change in the conformational ensemble. The
claim is narrower: the evidence supports a **thermodynamic stabilisation of a DNA-binding
competent state**, and there is no evidence for a **propagating signal** from the TPP1 interface
to the catalytic site.

---

## 4. Consequences for the three models

| Model | Status after WP2 | The specific observation that decides it |
|---|---|---|
| **A. Local structural conduit (760–766)** | **Not supported.** Needs to be restated or dropped. | H762E is near-WT; only 766 of 760–766 is even at the interface; the stretch is a β-strand, not a helix (WP1 §5). |
| **B. Broad conformational allostery** | **Weakened but not excluded.** | No experiment shows a TPP1-dependent conformational change at the active site. The strongest support in the literature is intra-TPP1 µs MD of the OB domain (Comm Chem 2023), which does not show TERT-wide coupling. The published TEN-stabilisation model is a *thermodynamic* version of B, not a kinetic-pathway version. |
| **C. Avidity / anchor stabilisation, no propagation required** | **Best supported.** | TPT rescues anchor-site mutants (compensation), D788R is catalytic-defective but still stimulated, recruitment is separable, TPP1 makes no DNA contact. |

Model C should not be treated as the winner of a modelling contest. It is the model that the
existing experiments *fail to distinguish from a mechanical one*, because every TERT mutation
tested changes DNA affinity and TPP1 binding at the same time. That is a limitation of the data,
not a demonstration.

---

## 5. The experiment that does not exist yet

Separability of a *mechanical* stimulation requires a perturbation that changes the TPP1→TERT
signal without changing the TPP1–TERT binding equilibrium or the DNA affinity. Nothing in the
two datasets does that. Concretely, what would settle it:

1. **A TPP1 separation-of-function mutant.** TPP1 variants that bind TERT normally (TEN α3/α5,
  IFD contacts; WP1 §2) but fail to stimulate RAP, tested *in trans* with WT TERT at saturating
  concentration. If such a mutant exists, stimulation requires a specific interface and Model A/B
  become live again. The published directions to look are the TPP1 OB-domain surface loops —
  but note that the µs MD evidence for those loops is intra-TPP1 motion, not TERT coupling.
2. **Measure ΔRAP, not processivity.** For every mutant already in hand, report processivity
  −POT1–TPP1 and +POT1–TPP1 (the fold-stimulation). If a mutant loses basal processivity but
  keeps the fold-stimulation, the stimulation is independent of the anchor; if a mutant keeps
  basal processivity but loses the fold-stimulation, the stimulation is a real signal. This is a
  re-analysis of existing gels, not a new experiment — it is the cheapest decisive test available.
3. **Vary DNA affinity independently.** The compensation interpretation predicts that raising
  primer concentration (or using a higher-affinity primer) should phenocopy TPP1 for the anchor
   mutants. If it does, TPP1's effect on those mutants is affinity, not signal. This is a
  one-experiment discriminator and, as far as these two papers show, it has not been done.

Falsification criteria for Model C: a TPP1 variant that binds TERT with WT affinity and does not
stimulate RAP at any concentration would falsify "avidity is sufficient". Falsification criteria
for Model A/B: a mutant in the interface that keeps TPP1 binding but loses ΔRAP, combined with a
detectable TPP1-dependent conformational change at the active site by HDX-MS, NMR, or
single-molecule FRET.

---

## 6. Verdict of WP2

**Can separability be established from existing experimental data? Partially, and not the part
that matters.**

* **Recruitment vs. processivity: yes, established** — and it removes recruitment/avidity-in-cell
  from the mechanism.
* **Basal catalysis vs. processivity: yes, partially** — D788R and H762E dissociate the two.
* **Mechanical stimulation vs. DNA-affinity compensation: no** — the two are confounded in every
  existing perturbation, and existing data lean towards compensation.

The mechanism question therefore does not become answerable by adding simulation. It becomes
answerable by re-analysing existing gels for ΔRAP (item 2 above) and by testing whether raising
DNA affinity phenocopies TPP1 (item 3). Both are within reach of a wet lab and neither requires
the computational machinery in this repository.
