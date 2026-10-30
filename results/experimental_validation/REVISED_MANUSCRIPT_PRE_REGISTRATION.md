# REVISED PRE-REGISTRATION STUDY: EXPERIMENTAL VALIDATION OF THE TERT LYS757 PROCESSIVITY UNCOUPLING SWITCH

**Document ID:** `EXP-REG-001-REV-FINAL`  
**Effective Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/REVISED_MANUSCRIPT_PRE_REGISTRATION.md`  
**Format:** Pre-Registered Experimental Protocol & Hypothesis Framework  

---

## 1. Abstract

Human telomerase reverse transcriptase (hTERT) maintains telomere length at eukaryotic chromosome ends, a catalytic function aberrantly reactivated in >85% of human cancers to enable indefinite cellular immortality. The shelterin component TPP1 recruits telomerase to chromosome ends via its TEL-patch and stimulates repeat addition processivity (RAP). However, clinical drug development targeting telomerase has been stymied by the cancer–aging dilemma: complete catalytic or recruitment blockade triggers fatal aplastic anemia and bone-marrow failure in patients. 

Through high-resolution coordinate analysis of human telomerase cryo-EM structures (PDB 7TRE, 7TRD, and 7TRF), we previously demonstrated that the TRAP domain does not actuate a rigid mechanical lever arm ($\Delta \text{C}\alpha \le 0.28\text{ \AA}$), but that TPP1 forms a bipartite clamp spanning the TEN and TRAP domains. In the TPP1-bound reconstruction (7TRE), `Tyr176:OH` (TEN domain) and `Gln794:OE1` (TRAP domain) are situated within 2.338 Å of each other (compared to 4.372 Å in unbound structures), and refined atomic displacement parameters ($B$-factors) at the single-stranded DNA exit pore show marked local reduction. 

Here, we present a pre-registered experimental validation study designed to test the central hypothesis: that `TERT Lys757` acts as an uncharacterized basic gatekeeper at the DNA exit pore ($5.75\text{ \AA}$ from the primer, $25.43\text{ \AA}$ from TPP1, $8.33\text{ \AA}$ from catalytic RNA) whose mutation to alanine (`Lys757A`) selectively abolishes TPP1-stimulated processivity while preserving nanomolar shelterin recruitment and basal catalytic elongation. We specify the complete recombinant expression protocol (Expi293F), Surface Plasmon Resonance (SPR) binding assays, direct radiometric primer extension RAP assays, statistical power calculations ($N = 3$, $>99\%$ power), and pre-committed falsification thresholds. This study establishes the empirical requirements to prove or refute the first TERT-side molecular uncoupling switch.

---

## 2. Epistemic Qualifications and Pre-Publication Demarcations

In strict adherence to rigorous scientific standards, this study formalizes the following critical distinctions:
1. **Structural Proximity vs. Biological Driving Force**: A static inter-atomic distance measured in a cryo-EM coordinate model (such as the 2.338 Å distance between `Tyr176:OH` and `Gln794:OE1`) identifies a candidate structural latch, but does not by itself prove that the hydrogen bond causally drives processivity in solution.
2. **Cryo-EM Displacement vs. Solution Thermal Variance**: Lower refined $B$-factors in cryo-EM models reflect local coordinate precision and particle alignment in vitrified vitreous ice; they motivate dynamic hypotheses but do not directly measure solution thermodynamic fluctuations.
3. **Harmonic Models vs. Empirical Dynamics**: Gaussian Network Model (GNM) normal mode calculations provide coarse-grained harmonic predictions of mode dampening; they do not establish biological solution dynamics.
4. **Pre-Registered Prediction vs. Experimental Fact**: All phenotypic values ascribed to `TERT Lys757A` ($K_d \le 300\text{ nM}$, $k_{\text{cat}} \ge 85\%$, $\chi_{\text{RAP}} \le 1.8$) are **pre-registered testable hypotheses**, not observed experimental facts.
5. **Kinetic Partition Models vs. Direct Measurements**: The apparent 12-fold difference in primer dissociation rate ($0.150 \to 0.012\text{ s}^{-1}$) is a model-derived partition parameter fitted to bulk pulse-chase half-lives ($1.8\text{ min} \to 22.5\text{ min}$; Latrick & Cech, 2010), rather than a directly observed single-molecule microscopic rate.

---

## 3. Pre-Registered Experimental Hypotheses

```
========================================================================================================================
PRE-REGISTERED HYPOTHESES AND BOUNDARY CRITERIA
========================================================================================================================
Hypothesis ID   Experimental Assay    Construct Tested        Predicted Outcome             Falsification Boundary
------------------------------------------------------------------------------------------------------------------------
H-REC           SPR Binding           TERT Lys757A            Kd <= 300 nM (Near WT)        Falsified if Kd > 1,500 nM
                (Biacore 8K)          (Recruitment test)      (Nanomolar affinity intact)   (Recruitment interface artifact)

H-CAT           Direct Extension      TERT Lys757A            k_cat >= 85% of WT            Falsified if k_cat < 30% of WT
                Radiometric Gel       (Basal catalysis test)  (Active site elongation intact)(Core catalytic disruption)

H-RAP           RAP Gel Titration     TERT Lys757A            chi_RAP <= 1.8 fold           Falsified if chi_RAP >= 4.0 fold
                (0-500 nM TPP1)       (Processivity test)     (Selective collapse of RAP)   (Processivity uncoupling absent)

H-RESCUE        RAP Gel Titration     TERT Lys757R            chi_RAP >= 4.5 fold           Falsified if chi_RAP <= 2.5 fold
                (0-500 nM TPP1)       (Stereochemical rescue) (Charge & geometry rescues)   (Non-specific charge artifact)
========================================================================================================================
```

---

## 4. Summary of Methodological Framework

The complete experimental protocols are specified across three accompanying documents:
- **Construct Design & Mutagenesis**: [`CONSTRUCT_DESIGN_WT_K757A_K757R.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md) (Q5 mutagenesis, BssHII screening, full-length sequencing primers).
- **Holoenzyme Reconstitution SOP**: [`RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md) (Expi293F co-transfection, Strep-Tactin XT capture, 500 mM KCl wash, PreScission 3C cleavage, Superose 6 SEC).
- **SPR Binding Data Package**: [`SPR_BINDING_DATA_PACKAGE.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/SPR_BINDING_DATA_PACKAGE.md) (Biacore 8K, 1:1 Langmuir kinetic fitting, concentration series 0–480 nM).
- **Quantitative Gel RAP Package**: [`QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md) (12% Urea-PAGE, $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$, Typhoon scanning, rolling-ball background subtraction, $3i$ normalization).
- **Statistical Analysis & Replication**: [`STATISTICAL_ANALYSIS_AND_REPLICATION_PLAN.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/STATISTICAL_ANALYSIS_AND_REPLICATION_PLAN.md) (Two-way repeated-measures ANOVA, Dunnett's post-hoc, blinded operator protocol, $N = 3$ independent biological replicates).

---

## 5. Decision Matrix & Scientific Consequences

```
========================================================================================================================
DECISIVE OUTCOME SCENARIOS
========================================================================================================================
Scenario                     Observed Phenotype (Lys757A)     Observed Control (Lys757R)   Definitive Conclusion
------------------------------------------------------------------------------------------------------------------------
Scenario 1: BREAKTHROUGH     Kd <= 300 nM; k_cat >= 85%;       Kd <= 300 nM; k_cat >= 90%;   HYPOTHESIS CONFIRMED.
                             chi_RAP <= 1.8                   chi_RAP >= 4.5                TERT-side uncoupling switch proven.
                                                                                            Eligible for Nature / Cell.

Scenario 2: CATALYTIC FAULT  k_cat < 30% WT                   k_cat < 30% WT                HYPOTHESIS REFUTED.
                                                                                            Lys757 affects active-site chemistry.

Scenario 3: RECRUITMENT LOSS Kd > 1,500 nM                    Kd > 1,500 nM                 HYPOTHESIS REFUTED.
                                                                                            Lys757 is part of recruitment surface.

Scenario 4: NO UNCOUPLING    Kd <= 300 nM; k_cat >= 85%;       Kd <= 300 nM; k_cat >= 90%;   HYPOTHESIS REFUTED.
                             chi_RAP >= 4.0                   chi_RAP >= 4.0                Lys757 does not control processivity.
========================================================================================================================
```
