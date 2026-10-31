# The Corrected Claim — TERT IFD-TRAP DNA-Anchor ↔ TPP1 Interface Coupling

**Version 2.0 — 2026-10-30** · supersedes the "TEN–TPP1–TRAP Bipartite Clamp / Lys757A
separation-of-function" claim in `README.md`, `report.md` and
`results/breakthrough_investigation/FINAL_RESEARCH_REPORT.md`
**Status:** pre-registered hypothesis, not a finding. No experiment in this repository tests it.

---

## 1. The claim, in one paragraph (abstract-ready)

> Human telomerase is recruited to telomeres by the TPP1 TEL patch, which docks onto a bipartite
> surface formed by the TERT TEN domain (Lys78) and the TERT IFD-TRAP (Arg774). TPP1 binding also
> stimulates repeat-addition processivity (RAP) 3–4 fold, but how a docking event 25 Å from the
> DNA exit channel is converted into primer retention is unknown. The IFD-TRAP residues
> Lys757, Phe759 and Gln794, together with the TEN `174PLYQ177` loop, form a DNA anchor surface at
> the primer exit channel and are required for processive synthesis. Here we test whether these two
> surfaces — the TPP1 docking surface and the DNA anchor surface — are **functionally coupled**, or
> whether they act as independent modules whose effects on processivity simply multiply. Using
> matched single and double mutations on both surfaces, reconstituted telomerase, and a
> pre-registered epistasis statistic on the TPP1 fold-stimulation of RAP, we distinguish allosteric
> coupling from independent-module behaviour. The answer determines whether "processivity
> stimulation" is a separate, tunable function of telomerase — the premise of any
> processivity-selective therapeutic — or an unavoidable by-product of recruitment plus anchor-site
> chemistry.

**One-sentence version for a talk slide:**

> We do not yet know whether TPP1 *switches on* a DNA-anchor mechanism or merely *recruits* an
> enzyme that already works; this is a two-experiment question, and we are going to answer it.

---

## 2. Why the old claim had to be replaced

| Old claim element | Why it fails |
|---|---|
| "TERT Lys757 is the critical mechanical hinge of an allosteric lever" | Lys757 is already assigned a function: it is part of the **DNA anchor site** (Sekne et al., *Science* 2022). It sits 27.9 Å from TPP1 and its closest primer contact is through the **backbone carbonyl O**, not the side chain. |
| "K757A uncouples TPP1-stimulated processivity while preserving basal catalysis" | **Already tested and contradicted.** Sekne 2022: TERT K757A/F759A *"significantly reduced telomerase activity and RAP"*. Tesmer 2023: K757E/F759E *"defective in telomerase activity and processivity"*. The repository's own criterion (k_cat ≥ 85 % of WT for K757A) is therefore predicted to fail. |
| "Tyr176–Gln794 is a state-dependent allosteric switch" | Tyr176/Gln794 are the paper's anchor-site residues (Sekne 2022, Y176A/Q177A and Q794A reduce activity and RAP). The 2.338 → 4.372 Å difference is real in the coordinates but sits in a region whose B-factors differ globally by 5× between the two depositions, so it cannot be read as a state-dependent switch. |
| "The 78.8 % B-factor drop at Lys757 proves TPP1 quenches the exit pore" | Not site-specific: Lys757's 7TRE/7TRD B-factor ratio (0.24) equals its own domain's (0.26) and the global ratio is 0.19. |
| "MD validates the mechanism" | 2.36 ps on a 97-residue fragment with no RNA and no DNA; the "K757A" system file is byte-identical to WT. |
| "TELIN-1 is a lead inhibitor" | Two disconnected fragments, 160 clashes to 0.17 Å, no measured affinity, no literature. |

What is **not** in dispute: the interface geometry (857 Å² per partner, R774–E168 at 2.32 Å,
K78–E215 at 3.58 Å), the pre-organisation of Leu183, the DE-loop disorder-to-order change, and the
coordinate invariance of the catalytic core. Those are measurements and they replicate.

---

## 3. Claim architecture (say it in this order)

**Tier 1 — Established (cite, do not claim).**
1. TPP1 TEL patch (E168, E169, E171, L183, L212) engages TERT TEN (K78) + IFD-TRAP (R774); TPP1
   binding damps TEN–TRAP dynamics (Liu 2022, *Nature* 604:578–583).
2. POT1–TPP1 stimulates RAP 3–4 fold by slowing primer dissociation and aiding translocation
   (Wang 2007, *Nature* 445:506–510; Latrick & Cech 2010, *EMBO J.* 29:924–933).
3. TPP1 K170Δ causes Hoyeraal–Hreidarsson syndrome and reduces both recruitment and processivity
   (Kocak 2014, *Genes Dev.* 28:2090–2102; Bisht 2016, *PNAS* 113:13021–13026).
4. TERT IFD-TRAP K757/F759/Q794 and TEN `174PLYQ177` form a DNA anchor site; K757A/F759A and Q794A
   reduce activity and RAP (Sekne 2022); K757E/F759E reduce activity and processivity but not
   telomere localisation (Tesmer 2023, *Life Sci. Alliance* 6:e202201727).

**Tier 2 — This work's own measurements (verified in `claim_recheck.csv`).**
5. 7TRE TERT–TPP1 interface: 857 Å² per partner, 5 residue pairs within 3.5 Å, R774:NH1–E168:OE1 =
   2.320 Å, K78:NZ–E215:OE2 = 3.584 Å, L769:CD2–L183:CD2 = 3.665 Å, P771:CG–L212:CD2 = 3.524 Å.
6. TPP1 OB-fold, bound vs unbound: global 1.966 Å RMSD; Leu183 0.512 Å; Lys170 4.293 Å;
   Glu169 3.809 Å; Glu171 3.320 Å; Glu168 1.253 Å.
7. Catalytic core invariance 7TRE vs 7TRD: 0.3377 Å over 770 Cα (below the 0.35–0.5 Å coordinate
   uncertainty; **not** evidence for or against allostery).
8. Anchor site is 25.4 Å (whole-residue) / 27.9 Å (NZ) from TPP1 and 5.75 Å (backbone O) / 6.31 Å
   (NZ) from the primer's 5′ phosphate — spatially segregated from the docking surface.

**Tier 3 — The new, testable claim (this is the breakthrough candidate).**
9. **Coupling claim.** TPP1's fold-stimulation of RAP is not a fixed multiple of whatever the
   anchor site allows. Mutating the anchor surface should therefore change χ_RAP *more* than it
   changes basal catalysis, and anchor × interface double mutants should be **non-additive** in
   log χ_RAP.
   *Null (independent modules):* anchor mutations scale basal activity and leave χ_RAP unchanged;
   interface mutations abolish binding and therefore abolish stimulation; the double mutant is the
   product of the two.
   *Alternative (recruitment-only):* TPP1 stimulates processivity purely by increasing the local
   enzyme concentration / residence time at the telomere; no anchor coupling exists; χ_RAP depends
   only on K_d.

**Tier 4 — Speculation (label as such).**
10. If Tier 3 is non-null, a small molecule or point mutant that lowers χ_RAP without lowering
    k_cat would be a processivity-selective tool. No such molecule exists; TELIN-1 is not one; and
    no cancer-over-stem-cell selectivity has been demonstrated for any telomerase modulator other
    than the observed myelosuppression of catalytic inhibitors (imetelstat).

---

## 4. Pre-registered design

### 4.1 Constructs (all in one reconstitution background, ≥ 2 independent protein preparations)

| Arm | Purpose |
|---|---|
| WT hTERT + hTR | reference |
| **K757A**, K757R, F759A, Q794A | anchor surface (Tier 1 function) |
| **K78A** (TEN), **R774A** (IFD-TRAP) | docking surface, TERT side |
| TPP1 E168A/E169A/E171A ("TEL-patch dead") | docking surface, TPP1 side; recruitment-negative control |
| **K757A + K78A**, **K757A + R774A**, **Q794A + R774A** | double mutants (the actual experiment) |
| D712A | catalytic-dead control |

### 4.2 Measurements

1. **Basal catalysis** — direct primer-extension, no TPP1: k_cat and product-length distribution.
2. **TPP1 responsiveness** — same assay ± saturating POT1–TPP1; χ_RAP = v(+TPP1)/v(−TPP1) from
   the band-intensity decay (the repository's `rap_gel_quantification.py` is adequate).
3. **Recruitment** — TPP1 binding by SPR (TERT construct) **and** telomere localisation by
   IF/FISH in cells (as Tesmer 2023 did).
4. **Primer retention** — pulse-chase dissociation, k_off of the bound primer ± TPP1 (the
   mechanistic endpoint that distinguishes "slower release" from "faster catalysis").

### 4.3 The statistic

For variant X: χ(X) = RAP(+TPP1)/RAP(−TPP1).
Epistasis on a multiplicative scale:

```
ε = ln χ(K757A+K78A) − ln χ(K757A) − ln χ(K78A) + ln χ(WT)
```

* ε ≈ 0 → independent modules (or a pure recruitment effect).
* ε < 0 (sub-multiplicative) → the anchor surface is **required** for TPP1 to stimulate →
  functional coupling exists.
* ε > 0 → the two surfaces are mutually antagonistic.

**Pre-registered thresholds.** n ≥ 4 independent reconstitutions per arm, triplicate gels; 95 % CI
by hierarchical bootstrap over reconstitutions; a claim of coupling requires **|ε| ≥ 0.30 log units
(≈ 35 % of a 3.5-fold stimulation)** with the CI excluding 0, *and* χ(WT) ≥ 2.0 in the same
experiment (assay sensitivity gate). Equivalence (null) requires the CI to exclude |ε| ≥ 0.30.

### 4.4 What would falsify the claim

* χ(K757A) ≈ χ(WT) once basal activity is accounted for (ε ≈ 0): anchor site is not part of the
  TPP1 conduit → **claim dead**, and the "allosteric switch" framing should be dropped entirely.
* χ(K757A) < 1 (TPP1 becomes inhibitory): the anchor site is permissive, not coupling → reframe as
  "TPP1 requires a functional exit channel" (a much weaker, already-published statement).
* K78A or R774A fails to reduce TPP1 binding in the same background: the docking surface is
  misassigned → the whole design collapses to a binding assay.

### 4.5 The second, cheaper experiment that should be done first

Before any new construct: **re-analyse the published Sekne 2022 Fig. 4E–G data as χ_RAP**.
If χ_RAP(K757A/F759A) is already visibly below χ_RAP(WT) at matched basal activity, the coupling
hypothesis is already supported by published data and the novelty claim reduces to "we quantified
an epistasis term that was already apparent". Know this before writing an abstract.

---

## 5. Text you can paste

**Abstract sentence.**
> We show that the TERT IFD-TRAP DNA-anchor surface and the TPP1 docking surface are functionally
> coupled, such that the processivity stimulation produced by POT1–TPP1 depends on an intact
> anchor site rather than on recruitment alone (epistasis ε = … log units, 95 % CI …).

**Results sentence (only after the experiment).**
> In four independent reconstitutions, POT1–TPP1 stimulated WT telomerase 3.4 ± 0.4 fold; the
> anchor mutant K757A retained 88 ± 6 % of WT basal catalysis but only 1.6 ± 0.3 fold stimulation,
> and the K757A/K78A double mutant was sub-multiplicative by ε = −0.62 (95 % CI −0.81 to −0.43),
> indicating that the two surfaces act on a common step.

**Limitations sentence.**
> Static cryo-EM models cannot establish coupling; the coordinate invariance of the catalytic core
> (0.34 Å) is within experimental error and is reported as such. Coupling is inferred from
> functional epistasis and remains to be localised to a structural transition.

**Do not write:** "breakthrough", "discovery", "validated", "proves", "allosteric switch
confirmed", "selective for cancer", "MD shows", "docking predicts", "TELIN-1 lead", or any number
from `data/benchmarks/*` presented as measurement.

---

## 6. Why this is the strongest honest version

| Criterion | Old claim | This claim |
|---|---|---|
| Novel? | No — K757A/F759A and K757E already published (Sekne 2022; Tesmer 2023) | Yes — the *coupling* between the two published surfaces has not been measured, and epistasis on χ_RAP has not been reported |
| Consistent with published data? | No — predicts preserved basal activity, published data show it is impaired | Yes — compatible with either outcome |
| Falsifiable in one experiment? | Its central prediction is already falsified | Yes — one double-mutant series with a pre-registered ε |
| Uses what the repo already has? | — | Yes: reconstitution SOP, RAP densitometry, SPR, and the verified interface geometry as the structural frame |
| Therapeutic overreach? | "Starve cancer cells, spare stem cells" | None; translational claim deferred until Tier 3 resolves |

---

## 7. Scope note on the therapeutic question

Nothing in this repository supports a cancer-versus-stem-cell therapeutic window. Telomerase
inhibition is limited by on-target marrow toxicity (imetelstat: grade 3/4 thrombocytopenia and
neutropenia), TEL-patch loss causes Hoyeraal–Hreidarsson syndrome, and the kinetic argument that
partial processivity loss would spare slowly dividing stem cells while killing tumour cells has
**not been measured** in any system. That argument is a hypothesis about telomere-length kinetics,
not a result, and should be stated as such — in one sentence, at the end, with no numbers.

---

## 8. Update 2026-10-30 — the bottleneck analysis changed this claim

`scripts/allosteric_coupling.py` tested the premise behind Tier 3 and produced a structural result
that supersedes the wording above. Full analysis:
[`../allostery/ALLOSTERIC_COUPLING_ANALYSIS.md`](../allostery/ALLOSTERIC_COUPLING_ANALYSIS.md).

**What changed**

1. The "25 Å gap" that motivated an allosteric framing is a mis-measurement of the problem.
   TPP1 sits on the **same IFD-TRAP helix** that carries the DNA anchor; the distance decays
   continuously along it (Arg756 26.0 → Lys757 25.4 → Phe759 20.0 → Ser761 15.2 → Thr765 7.2 →
   Leu766 3.6 → Arg774 2.3 Å). There is no void for a signal to cross.
2. An anisotropic network model with a distance-matched null finds **no** coupling beyond
   proximity (z = +0.34; the displaced-TPP1 decoy sits *below* its own null; the TEN↔IFD-TRAP
   staple metric changes by −0.04 once the non-specific effect is removed). Same model also fails
   a positive control (catalytic site ↔ primer, |corr| 0.10 vs 0.09 background), so this is
   "this model class cannot adjudicate", not "coupling is absent".
3. The two-hand clamp is not present in 7TRE: the TEN domain makes no contact with the modelled
   primer (Tyr176 is 11.1 Å away). Only the RT/IFD-TRAP side touches DNA.

**Superseding claim (H_conduit)**

> TPP1 stimulates processivity by binding one face of the IFD-TRAP helix that carries the DNA
> anchor, propagating a local change along residues 760–766 into the primer exit channel. It does
> not require long-range allosteric propagation.

Discriminating prediction: helix mutants at 760–766 (and L766P) should reduce χ_RAP with K_d
intact, whereas the already-published K757A/F759A anchor mutants lose basal catalysis as well
(Sekne et al. 2022). If both classes behave identically, H_conduit is falsified and the effect is
simply "anything that perturbs the exit channel".

The Tier 3 epistasis design above remains valid and complementary — but it should now be run on
the 760–766 panel rather than on the anchor-site mutants, because the anchor mutants cannot
distinguish recruitment from catalysis (their basal activity defect is already known).
