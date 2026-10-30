# FINAL EXPERIMENTAL PROTOCOL: HYPOTHESIS FALSIFICATION AND VALIDATION SUITE

**Document ID:** `BSV-DOC-004-REV-FINAL`  
**Date:** 2026-10-30  
**Classification:** Evidence Level 4 (Experimental Design) & Decision Framework  
**Target Repository:** `results/breakthrough_investigation/FINAL_EXPERIMENTAL_PROTOCOL.md`  

---

## 1. Executive Purpose & The Three Breakthrough Criteria

This protocol defines the definitive wet-lab experimental pipeline required to prove or falsify the **TEN–TPP1–TRAP Bipartite Clamp Mechanism**. 

To transition this discovery from a computational hypothesis to a verified biological breakthrough, an experimental candidate must satisfy **all three** of the following quantitative criteria:

```
+---------------------------------------------------------------------------------------------------------+
|                                    THE THREE BREAKTHROUGH CRITERIA                                      |
+---------------------------------------------------------------------------------------------------------+
| 1. STRUCTURAL CRITERION:                                                                                |
|    Physical TEN-TRAP coupling must be demonstrated experimentally.                                      |
|    - Threshold: HDX-MS demonstrates significant dynamic protection (ΔD >= 15%, p < 0.01) at the TEN     |
|      cleft (peptides 170–182) and TRAP N-terminus (peptides 752–765) upon TPP1 binding.                 |
+---------------------------------------------------------------------------------------------------------+
| 2. BIOCHEMICAL SEPARATION-OF-FUNCTION CRITERION:                                                        |
|    The mutation must preserve recruitment binding and basal catalysis, but abolish RAP stimulation.     |
|    - Recruitment binding: TPP1 OB domain affinity Kd <= 300 nM (SPR / BLI; WT Kd ~ 150 nM).             |
|    - Basal catalytic activity: k_cat >= 85% of Wild-Type in the absence of TPP1.                       |
|    - Selective RAP collapse: TPP1-stimulated processivity fold-enhancement chi_RAP <= 1.8               |
|      (versus WT chi_RAP >= 5.0 to 8.0).                                                                 |
+---------------------------------------------------------------------------------------------------------+
| 3. MECHANISTIC KINETIC CRITERION:                                                                       |
|    Processivity collapse must be directly driven by premature primer release during translocation.       |
|    - Primer dissociation rate during translocation (k_off^translocation) measured via pulse-chase or    |
|      smFRET remains elevated (>= 0.080 s^-1) in the mutant despite TPP1 presence (versus WT = 0.012 s^-1). |
+---------------------------------------------------------------------------------------------------------+
```

---

## 2. Experimental Candidate Panel

The test panel includes the primary uncoupling candidate, conservative stereochemical controls, domain bridge disruptors, and established reference standards:

```
========================================================================================================================
EXPERIMENTAL CANDIDATE MUTATION PANEL
========================================================================================================================
Identifier    Construct               Rationale                                     Predicted Outcome
------------------------------------------------------------------------------------------------------------------------
MUT-01        TERT Lys757A (Primary)  Pore exit charge removal; uncontacted by      Preserves Kd (<=300 nM) & basal k_cat;
                                      TPP1 (25.4 Å) or RNA template (8.3 Å).         collapses chi_RAP <= 1.8. BREAKTHROUGH HIT.

MUT-02        TERT Lys757R (Control)  Charge & packing preservation at pore rim;    Rescues chi_RAP >= 4.5. Proves basic charge
                                      maintains positive electrostatic exit guide.  and aliphatic packing are required.

MUT-03        TERT Tyr176F (Control)  Aromatic ring preserved; removes single       Disrupts H-bond to Gln794 (OE1); isolates
                                      phenolic hydroxyl forming H-bond to Gln794.    inter-domain hydrogen-bonding requirement.

MUT-04        TERT Tyr176A (Severe)   Complete aromatic ring ablation; disrupts     Severe structural collapse of TEN anchor
                                      Phe759 / Ile792 hydrophobic packing.           and TRAP coupling; loss of RAP.

MUT-05        TERT Gln794A (Partner)  TRAP-side partner of Tyr176 H-bond;           Reciprocal partner test; validates the
                                      removes OE1 acceptor atom.                    inter-domain hydrogen bond from TRAP side.

MUT-06        TERT Arg756A (Negative) Direct template contact (A55:C2 = 3.12 Å);     Confirmatory failure: basal k_cat < 30%
                                      template-alignment confounding control.       due to RNA displacement.

REF-01        TERT WT (Reference)     Native human TERT holoenzyme.                 Baseline: Kd ~ 150 nM, chi_RAP = 6.2 ± 0.8.

REF-02        TPP1 L104A (Reference)  Established TPP1-side separation-of-function  Trans-control: retains binding, fails to
                                      mutant (Nandakumar et al., 2012).             stimulate RAP in trans.
========================================================================================================================
```

---

## 3. Recombinant Expression, Purification, and Assembly

### 3.1 Expression System
- **Host**: Expi293F mammalian cells or Sf9 insect cells using the MultiBac baculovirus system.
- **Co-expression**:
  - Full-length human TERT (residues 1–1132) fused to an N-terminal Twin-Strep-tag and 3xFLAG-tag separated by a 3C protease cleavage site.
  - Full-length human telomerase RNA (hTR, residues 1–451) driven by the human U3 snoRNA promoter.
  - Human TPP1 OB domain (residues 89–249) fused to a C-terminal $6\times\text{His}$-tag.

### 3.2 Purification Protocol
1. **Cell Lysis**: Resuspend cell pellets in ice-cold lysis buffer ($50\text{ mM}$ HEPES-KOH pH 7.5, $300\text{ mM}$ KCl, $10\%\text{ glycerol}$, $1\text{ mM}$ EDTA, $1\text{ mM}$ TCEP, $0.1\%\text{ Triton X-100}$, supplemented with cOmplete EDTA-free protease inhibitors and RNasin Plus).
2. **Affinity Chromatography**:
   - Capture on Strep-Tactin XT Superflow resin.
   - Wash with 30 column volumes (CV) of high-salt buffer ($50\text{ mM}$ HEPES-KOH pH 7.5, $500\text{ mM}$ KCl, $1\text{ mM}$ TCEP).
   - On-column cleavage with 3C protease ($16\text{ h}$ at $4\text{ }^\circ\text{C}$) to liberate native holoenzyme.
3. **Size Exclusion Chromatography (SEC)**:
   - Inject concentrated eluate onto a Superose 6 Increase 10/300 GL column equilibrated in running buffer ($20\text{ mM}$ HEPES-KOH pH 7.5, $150\text{ mM}$ KCl, $2\text{ mM}$ $\text{MgCl}_2$, $1\text{ mM}$ TCEP).
   - Collect peak fractions corresponding to the monomeric holoenzyme complex (~450 kDa).
   - Assess purity and stoichiometric RNA incorporation via SDS-PAGE and denaturing urea-PAGE.

---

## 4. Multi-Tiered Quantitative Assay Pipeline

### 4.1 Assay 1: Structural Integrity & Folding (nanoDSF)
- **Objective**: Ensure mutations do not disrupt global tertiary folding.
- **Method**: Differential scanning fluorimetry using Prometheus NT.48 (NanoTemper).
- **Measurement**: Dual-UV intrinsic fluorescence ($330\text{ nm} / 350\text{ nm}$) across a thermal ramp ($20\text{ }^\circ\text{C}$ to $95\text{ }^\circ\text{C}$ at $1\text{ }^\circ\text{C/min}$).
- **Passing Threshold**: Apparent melting temperature $T_m \ge T_{m,\text{WT}} - 2.0\text{ }^\circ\text{C}$.

### 4.2 Assay 2: TPP1 Recruitment Binding Affinity (SPR / BLI)
- **Objective**: Measure equilibrium dissociation constant ($K_d$) for the TERT–TPP1 interaction.
- **Instrument**: Biacore 8K or Octet RH16.
- **Format**:
  - Sensor chip: Biotin CAP or SA chip coated with C-terminally biotinylated TPP1 OB-fold domain (residues 89–249).
  - Analyte: Purified WT or mutant TERT holoenzyme titrated at concentrations $0, 15, 30, 60, 120, 240, 480\text{ nM}$.
  - Running buffer: $20\text{ mM}$ HEPES-KOH pH 7.5, $150\text{ mM}$ KCl, $2\text{ mM}$ $\text{MgCl}_2$, $0.05\%\text{ Surfactant P20}$, $1\text{ mM}$ TCEP.
- **Passing Threshold**: Calculated $K_d \le 300\text{ nM}$ (WT baseline: $150 \pm 25\text{ nM}$).

### 4.3 Assay 3: Direct Primer Extension RAP Quantification
- **Objective**: Measure repeat addition processivity and fold-stimulation by TPP1.
- **Reaction Conditions**:
  - Substrate: $1\text{ }\mu\text{M}$ telomeric DNA primer 5'-$(\text{TTAGGG})_3$-3' or 5'-$(\text{dT})_{18}$-3'.
  - Nucleotides: $500\text{ }\mu\text{M}$ dATP, $500\text{ }\mu\text{M}$ dTTP, $5.0\text{ }\mu\text{M}$ dGTP, and $0.15\text{ }\mu\text{M}$ $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ ($3,000\text{ Ci/mmol}$).
  - Co-factor: Recombinant POT1–TPP1 heterodimer titrated at $0, 50, 100, 250, 500\text{ nM}$.
  - Incubation: $60\text{ min}$ at $30\text{ }^\circ\text{C}$.
- **Detection**: Quench with EDTA/formamide stop buffer, resolve on a $12\%$ denaturing polyacrylamide gel ($7\text{ M}$ urea), expose to phosphor storage screens, and scan on a Typhoon FLA 9500.
- **Mathematical Quantification**:
  - Activity: Total phosphor signal normalized to internal recovery control. Basal threshold: $k_{\text{cat}} \ge 85\%\text{ WT}$.
  - Repeat Addition Processivity ($R_p$): Calculated from band intensity ratios according to Latrick & Cech (2010):
    $$R_{p,i} = \frac{\sum_{j=i+1}^\infty I_j}{\sum_{j=i}^\infty I_j}$$
  - RAP Fold-Stimulation ($\chi_{\text{RAP}}$):
    $$\chi_{\text{RAP}} = \frac{R_p(+\text{TPP1})}{R_p(-\text{TPP1})}$$
- **Passing Threshold**: For uncoupling hit, $\chi_{\text{RAP}} \le 1.8$ (WT: $\chi_{\text{RAP}} = 6.2 \pm 0.8$).

### 4.4 Assay 4: Translocation Kinetic Chase (Pulse-Chase & smFRET)
- **Objective**: Measure the off-rate of primer dissociation during translocation ($k_{\text{off}}^{\text{translocation}}$).
- **Pulse-Chase Protocol**:
  - Initiate elongation complex on 5'-biotinylated $(\text{TTAGGG})_3$ captured on streptavidin beads.
  - Chase with a 200-fold excess of non-extendable 3'-dideoxy oligonucleotide competitor.
  - Take time points at $t = 0, 1, 2, 5, 10, 20, 45, 90\text{ min}$.
  - Fit complex decay to a single exponential: $I(t) = I_0 e^{-k_{\text{off}} t}$.
- **Passing Threshold**:
  - In WT + TPP1: $t_{1/2} \ge 20.0\text{ min}$ ($k_{\text{off}} \approx 0.012\text{ s}^{-1}$).
  - In `Lys757A` + TPP1: $t_{1/2} \le 3.0\text{ min}$ ($k_{\text{off}} \ge 0.080\text{ s}^{-1}$), confirming uncoupling of the retention gate.

### 4.5 Assay 5: Dynamic Structural Protection (HDX-MS)
- **Objective**: Verify physical TEN–TRAP inter-domain clamping in solution.
- **Labeling**: Dilute holoenzyme into $\text{D}_2\text{O}$ buffer ($90\%\text{ final D}_2\text{O}$) at $4\text{ }^\circ\text{C}$ for $10\text{ s}, 1\text{ min}, 10\text{ min}, 60\text{ min}, 240\text{ min}$.
- **Quench & Digestion**: Quench to $\text{pH } 2.5$ at $0\text{ }^\circ\text{C}$ with $2\text{ M}$ guanidine hydrochloride and $100\text{ mM}$ TCEP; pass through an online immobilized pepsin column.
- **LC-MS/MS**: Measure mass shifts on an Orbitrap Eclipse Tribrid.
- **Key Target Peptides**:
  - TEN loop (residues 170–182, contains `Tyr176`)
  - TRAP N-terminal rim (residues 752–765, contains `Lys757`)
  - TRAP TEL-patch interface (residues 766–780, contains `Leu769` and `Arg774`)
- **Passing Threshold**: Statistically significant reduction in deuteration ($\Delta D \ge 15\%$, $p < 0.01$) in both peptide 170–182 and peptide 752–765 upon TPP1 addition.

---

## 5. Experimental Decision Matrix & Falsification Criteria

```
========================================================================================================================
EXPERIMENTAL DECISION MATRIX (PRE-COMMITTED INTERPRETATION)
========================================================================================================================
Assay Result for TERT Lys757A                                         Scientific Verdict & Consequence
------------------------------------------------------------------------------------------------------------------------
Scenario A (BREAKTHROUGH HIT):
  - nanoDSF: Tm >= Tm,WT - 1.5 °C (folded)
  - SPR: Kd <= 300 nM (TPP1 binding intact)
  - Basal activity: k_cat >= 85% WT (catalysis intact)                HYPOTHESIS FULLY CONFIRMED.
  - TPP1-stimulated RAP: chi_RAP <= 1.8 (processivity collapsed)       First TERT-side molecular uncoupling switch
  - Pulse-chase: k_off^trans >= 0.080 s^-1 (retention failure)         proven. Eligible for Nature / Cell publication.
  - HDX-MS: Dynamic protection at TEN/TRAP junctions

Scenario B (RECRUITMENT ARTIFACT):
  - SPR: Kd > 1,500 nM (TPP1 binding severely disrupted)              HYPOTHESIS REFUTED.
                                                                      Lys757 is part of the binding interface; not a
                                                                      separation-of-function switch.

Scenario C (CATALYTIC DISRUPTION):
  - Basal activity: k_cat < 30% WT (catalytic core impaired)          HYPOTHESIS REFUTED.
                                                                      Lys757 directly affects basal active site geometry
                                                                      (similar to Arg756).

Scenario D (PHENOTYPIC RESCUE BY ARG):
  - Lys757R: chi_RAP >= 4.8 (rescues processivity)                    MECHANISM FULLY VALIDATED.
                                                                      Proves the basic electrostatic pore rim requirement.
========================================================================================================================
```
