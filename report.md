# Benchmark Audit & Structural Bioinformatics Evaluation of Human Telomerase Recruitment and the Cancer–Aging Dilemma

**Comparative Cryo-EM and Crystallographic Investigation of PDB Entries 7TRE, 7TRD, 7TRF, 7TRC, and 2I46**  
**Integrity Mode**: Benchmark (Epistemic Rigor, Multi-Agent Audit & Prior-Art Calibration)  
**Primary Authorship**: Multi-Agent Computational Structural Bioinformatics Research Team  
**Date**: October 2026  
**Status**: Final Verified Research Report & Retracted Claims Catalog  

---

## 1. Executive Summary & Forensic Audit Overview

Human telomerase is an essential specialized ribonucleoprotein (RNP) reverse transcriptase responsible for de novo synthesis of telomeric hexanucleotide repeats (`5'-TTAGGG-3'`) at chromosome termini. Reactivation or transcriptional upregulation of telomerase confers replicative immortality to over 85–90% of human malignancies, establishing it as a preeminent broad-spectrum oncology target. However, basal telomerase activity is simultaneously required for the lifelong proliferative capacity and self-renewal of adult somatic stem cell compartments—most critically hematopoietic stem and progenitor cells (HSPCs), intestinal crypt base columnar cells (`Lgr5+`), hair follicle stem cells, and pulmonary epithelial progenitors.

Direct catalytic inhibition of the human telomerase reverse transcriptase (hTERT) active site by the lipid-conjugated oligonucleotide competitive inhibitor **imetelstat (GRN163L)** demonstrates potent anti-neoplastic activity in clinical trials, but is severely limited by on-target myelosuppression, specifically Grade 3/4 thrombocytopenia and neutropenia.

To investigate whether targeting holoenzyme recruitment to the shelterin component TPP1 rather than the catalytic core can offer a viable therapeutic strategy, this project executed an automated, sequence-aware structural bioinformatics study across high-resolution experimental structures:
- `7TRE`: Primary human telomerase catalytic core complexed with the TPP1 OB-fold domain (cryo-EM at 3.50 Å).
- `7TRD`: Unbound human telomerase catalytic core reference (cryo-EM at 3.30 Å).
- `7TRF`: Telomerase catalytic core RNP complexed with histone H2A–H2B heterodimers (cryo-EM at 3.70 Å).
- `7TRC`: Human telomerase H/ACA RNP lobe reference (cryo-EM at 3.30 Å).
- `2I46`: High-resolution X-ray crystal structure of the unbound human TPP1 OB-fold domain (1.85 Å / 2.70 Å).

### Forensic Audit & De-Hyping Summary
A comprehensive multi-perspective audit was conducted to review earlier claims generated in this codebase. The audit established:
1. **No Breakthrough Discovery**: Earlier reports claiming a "breakthrough discovery" of a "two-step sequential induced-fit mechanism" and "proof of zero allosteric coupling" were scientifically overstated. The primary salt bridge (`Arg774–Glu168`) and recruitment contacts were already discovered and published by Liu et al. (*Nature* 2022).
2. **De-Hyped Biophysics**: An RMSD of 0.338 Å between static cryo-EM structures (`7TRE` vs `7TRD`) does **not** prove zero allosteric coupling; it falls within experimental coordinate uncertainty and cannot measure dynamic allostery or vibrational entropy shifts.
3. **Quarantine of Candidate Ligand `TELIN-1`**: The putative small-molecule inhibitor `TELIN-1` was found to be an unvalidated in silico construct comprising two disconnected chemical fragments with 160 severe steric clashes against TERT coordinates, with zero bibliographic backing or wet-lab validation. It has been placed under strict quarantine.
4. **Clinical Reality of Interface Disruption**: Total pharmacological shutdown of the TERT–TPP1 interface is clinically contraindicated; human germline genetics demonstrates that TPP1 TEL-patch loss (e.g. Lys170del) causes Hoyeraal–Hreidarsson syndrome and lethal bone marrow failure.

---

## 2. Repository Inventory & Evidence Ledger Summary

Every numerical and biological claim in the project was audited, indexed into `results/audit/evidence_ledger.csv`, and traced to its underlying scripts and raw coordinates.

### Table 2.1: Key Evidence Ledger Traceability Matrix
| Claim ID | Metric / Topic | Reported Value | Source Coordinates | Status | Epistemic Classification |
|:---|:---|:---|:---|:---|:---|
| `CLM-001` | Interface BSA (Bondi) | 856.82 Å² | `7TRE` (Chains A & C) | Established | Validated computational metric |
| `CLM-002` | Interface BSA (ProtOr) | 923.16 Å² | `7TRE` (Chains A & C) | Established | Validated (+7.74% radius shift) |
| `CLM-003` | Contact Count (3.5 Å) | 9 atomic contacts | `7TRE` (Chains A & C) | Established | Verified core persistent contacts |
| `CLM-004` | Contact Count (4.0 Å) | 46 atomic contacts | `7TRE` (Chains A & C) | Established | Verified intermediate shell |
| `CLM-005` | Contact Count (4.5 Å) | 123 atomic contacts | `7TRE` (Chains A & C) | Established | Verified outer contact shell |
| `CLM-006` | Catalytic Core RMSD | 0.338 Å (770 Cα) | `7TRE` vs `7TRD` | Established | Verified rigid-body coordinate fit |
| `CLM-007` | TEN Domain Excursion | 0.434 Å (162 Cα) | `7TRE` vs `7TRD` | Established | Verified rigid-body displacement |
| `CLM-008` | Primary Salt Bridge | 2.32 Å (`R774:NH1--E168:OE1`)| `7TRE` | Established | Compaction refinement artifact |
| `CLM-011` | TPP1 Leu183 Shift | 0.512 Å (Cα) | `7TRE` vs `2I46` (Chain A)| Established | Pre-organized hydrophobic anchor |
| `CLM-012` | TPP1 Lys170 Shift | 4.293 Å (Cα) | `7TRE` vs `2I46` (Chain A)| Established | State-dependent loop ordering |
| `CLM-021` | Two-Step Sequential Fit | Kinetic temporal claim | `7TRE` and `2I46` | **REFUTED** | Unvalidated kinetic hypothesis |
| `CLM-022` | Zero Allosteric Coupling| Zero dynamic allostery | `7TRE` vs `7TRD` | **REFUTED** | Biophysically invalid inference |
| `CLM-023` | TELIN-1 Lead Candidate | Putative nanomolar lead| `telin1_inhibitor.pdb` | **QUARANTINED**| Non-physical in silico artifact |

---

## 3. Prior-Art & Novelty Analysis (Liu et al. vs This Work)

To maintain absolute academic and scientific integrity, all findings were benchmarked against primary literature, specifically Liu et al. (*Nature* 2022; 604:578–585), Sekne et al. (*Science* 2022; 375:1173–1179), and Wang et al. (*Nature Struct. Mol. Biol.* 2007; 14:1011–1017).

### Table 3.1: Prior-Art Concordance Matrix
| Structural Observation | Liu et al. (Nature 2022) | This Work | Novelty Classification |
|:---|:---|:---|:---|
| **Arg774–Glu168 Salt Bridge** | Identified as central TEL-patch contact mediating recruitment. | Independently recalculated and verified at 2.32 Å. | **Positive Control (Replicated)** |
| **Lys78–Glu215 Contact** | Described as secondary N-terminal electrostatic interaction. | Independently recalculated at 2.89 Å (heavy atom) / 3.58 Å. | **Positive Control (Replicated)** |
| **Leu183 / Leu212 Hydrophobic Anchor**| Identified as core hydrophobic TEL-patch anchor. | Confirmed pre-organized with 0.512 Å Cα displacement. | **Positive Control (Replicated)** |
| **TEN–TRAP Plasticity vs Core Rigidity**| Noted that TPP1 binding dampens TEN-TRAP dynamics. | Rigorously quantified: 0.338 Å Core vs 0.434 Å TEN shift. | **Quantitative Clarification** |
| **DE-Loop Order-Disorder Transition** | Described bound conformation of DE-loop. | Analyzed B-factor transition (>100 Å² apo to 53 Å² holo) and lattice contacts. | **Refined Structural Insight** |
| **Two-Step Sequential Binding Claim** | Not claimed (static cryo-EM structure acknowledged). | Refuted earlier claims; categorized as open kinetic hypothesis. | **Epistemic Correction** |
| **Zero Allosteric Coupling Claim** | Not claimed (allosteric processivity stimulation acknowledged). | Refuted earlier claims; detailed dynamic and entropic allostery. | **Epistemic Correction** |

---

## 4. Falsifiable Hypotheses & Competing Models

Three competing biophysical hypotheses were formulated in `results/research/question_and_hypotheses.md` regarding the structural changes observed between apo and holo structures:

### Hypothesis 1 ($H_1$): Rigid Docking / Crystallographic Lattice Artifact
- **Premise**: The observed difference in the TPP1 DE-loop between `2I46` and `7TRE` is primarily an artifact of crystal packing in `2I46` space group $I 4_1 2 2$, rather than a physiologically relevant conformational change.
- **Evaluation**: Gemmi symmetry contact analysis (`check_lattice.py`) confirmed that Glu171 in `2I46` Chain A participates in crystal lattice contacts with symmetry-related molecules. Furthermore, Chain B lacks electron density entirely for residues 168–171. Thus, the crystal lattice stabilizes one arbitrary conformation out of a disordered ensemble.

### Hypothesis 2 ($H_2$): State-Dependent Induced-Fit Remodeling
- **Premise**: TPP1 binding induces a true transition from an intrinsically disordered apo loop ensemble to an ordered, receptor-bound conformation.
- **Evaluation**: Supported as a structural equilibrium phenomenon ($B$-factors drop from $>100\text{ \AA}^2$ to $53\text{ \AA}^2$). However, claiming that this occurs in a **temporal sequence** (hydrophobic docking first, then electrostatic zipping) is unsupported without stopped-flow kinetic data.

### Hypothesis 3 ($H_3$): Computational & Parameter Sensitivity
- **Premise**: Reported discrepancies in interface metrics (BSA, contact counts) arise from computational definitions rather than biophysical reality.
- **Evaluation**: Confirmed. Varying the probe radius convention from Bondi ($856.82\text{ \AA}^2$) to ProtOr ($923.16\text{ \AA}^2$) produces a $+7.74\%$ shift, whereas algorithm choice (FreeSASA vs Biopython ShrakeRupley) causes $<0.003\%$ variance.

---

## 5. Primary Structural Data & Quality Manifest

The structural foundation of this project comprises 5 experimental structures retrieved from the wwPDB and verified via SHA-256 cryptographic checksums.

### Table 5.1: Structural Data Manifest
| PDB ID | Experimental Method | Resolution | Modeled hTERT Range | Modeled TPP1 Range | Free Flag / Clashscore |
|:---|:---|:---|:---|:---|:---|
| **`7TRE`** | Cryo-EM (Single Particle) | 3.50 Å | 7–1132 (934 aa modeled) | 90–241 (152 aa modeled) | Clashscore 7.45; Ramachandran outliers 0.00% |
| **`7TRD`** | Cryo-EM (Single Particle) | 3.30 Å | 7–1132 (934 aa modeled) | Absent (Unbound Core) | Clashscore 7.59; Ramachandran outliers 0.00% |
| **`7TRF`** | Cryo-EM (Single Particle) | 3.70 Å | 7–1132 (934 aa modeled) | Absent (H2A-H2B Bound) | Clashscore 7.52; Ramachandran outliers 0.09% |
| **`7TRC`** | Cryo-EM (Single Particle) | 3.30 Å | Absent (H/ACA Lobe) | Absent | Clashscore 8.18; Ramachandran outliers 0.00% |
| **`2I46`** | X-Ray Diffraction | 1.85 / 2.70 Å | Absent (Apo TPP1) | Chain A: 90–241; Chain B: missing 168–171 | Clashscore 12.97; Ramachandran outliers 0.36% |

### Sequence Mapping & Offset Resolution:
Recombinant constructs contain non-canonical cloning sequences that shift author numbering:
- **TERT (`7TRE`, `7TRD`, `7TRF`)**: Contains an N-terminal 35-residue Twin-Strep tag ($\text{label\_seq\_id} = \text{auth\_seq\_id} + 35$). Canonical author numbering matches UniProt `O14746`.
- **TPP1 (`7TRE`)**: Carries an N-terminal cloning linker ($\text{auth\_seq\_id} = \text{label\_seq\_id} + 86$), matching UniProt `Q96AP0`.
- **TPP1 (`2I46`)**: Begins at residue Ser90 ($\text{auth\_seq\_id} = \text{label\_seq\_id} + 89$).

---

## 6. Quantitative Interface Geometry (Recalculated & Verified)

The hTERT–TPP1 interface was independently recalculated using orthogonal computational engines (`src/telotrade/independent/recalculate.py`) utilizing Biopython, Gemmi, and SciPy cKDTree architectures.

### Table 6.1: Master Recalculated Interface Benchmarks
| Metric | Original Reported Value | Independent Recalculation | Delta | Relative Error | Verification Status |
|:---|:---|:---|:---|:---|:---|
| **Interface BSA (Bondi)** | 856.82 Å² | 856.82 Å² | 0.00 Å² | 0.000% | VERIFIED |
| **Interface BSA (ProtOr)** | 923.16 Å² | 923.16 Å² | 0.00 Å² | 0.000% | VERIFIED |
| **Interface Contacts (3.5 Å)** | 9 contacts | 9 contacts | 0 | 0.000% | VERIFIED |
| **Interface Contacts (4.0 Å)** | 46 contacts | 46 contacts | 0 | 0.000% | VERIFIED |
| **Interface Contacts (4.5 Å)** | 123 contacts | 123 contacts | 0 | 0.000% | VERIFIED |
| **Arg774:NH1 – Glu168:OE1** | 2.32 Å | 2.322 Å | +0.002 Å | 0.086% | VERIFIED |
| **Lys78:NZ – Glu215:OE2** | 2.89 Å | 2.887 Å | -0.003 Å | 0.104% | VERIFIED |
| **Leu769:CD2 – Leu183:CD1** | 3.52 Å | 3.524 Å | +0.004 Å | 0.114% | VERIFIED |
| **Pro771:CG – Leu212:CD1** | 3.86 Å | 3.865 Å | +0.005 Å | 0.130% | VERIFIED |
| **Catalytic Core RMSD (770 Cα)**| 0.338 Å | 0.3377 Å | -0.0003 Å| 0.089% | VERIFIED |

### The Arg774–Glu168 Compaction Phenomenon:
The apparent $2.32\text{ \AA}$ distance between `TERT Arg774:NH1` and `TPP1 Glu168:OE1` is shorter than standard physiological hydrogen bonds ($2.7\text{--}3.1\text{ \AA}$). This reflects real-space refinement restraint compaction at $3.5\text{ \AA}$ cryo-EM map resolution rather than an anomalous low-barrier covalent interaction.

---

## 7. Parameter Sensitivity & Coordinate Uncertainty Stress-Testing

To establish the statistical robustness of the contact network, a multi-cutoff sensitivity analysis was conducted (`results/interface/sensitivity_analysis.csv`), partitioning 32 residue–residue interaction pairs into three distinct stability tiers:

```
[Core Persistent: 5 Pairs]  <-- Invariant across all cutoffs (Arg774, Glu168, Leu183)
         ▼
[Intermediate Boundary: 16 Pairs] <-- Formed between 3.5 Å and 4.0 Å (e.g. Lys78, Glu215)
         ▼
[Peripheral / Sensitive: 11 Pairs] <-- Present only at 4.5 Å shell (threshold-sensitive)
```

### Table 7.1: Contact Network Hierarchy
| Tier | Distance Cutoff | Contact Pairs | Key Residues Involved | Sensitivity / Biological Significance |
|:---|:---|:---|:---|:---|
| **Core Persistent** | $\le 3.5\text{ \AA}$ | 5 pairs (9 atomic contacts) | TERT Arg774, Leu769; TPP1 Glu168, Leu183 | Obligate anchoring core; intolerant to mutation |
| **Intermediate** | $3.5\text{ \AA} < d \le 4.0\text{ \AA}$ | 16 pairs (46 atomic contacts) | TERT Lys78, Gln775, Pro771; TPP1 Glu215, Leu212 | Secondary stabilizing clamp and peripheral hydrophobic rim |
| **Peripheral** | $4.0\text{ \AA} < d \le 4.5\text{ \AA}$ | 11 pairs (123 atomic contacts)| TERT Asp166, Glu778; TPP1 Lys170, Thr214 | Transient electrostatic halo; highly sensitive to thermal motion |

---

## 8. Cross-State Comparative Analysis & DE-Loop Remodeling

Superposition of bound TPP1 in `7TRE` against the unbound crystal structure `2I46` (Chain A) reveals a profound divergence in conformational plasticity across the TEL-patch:
1. **The Rigid Anchor (Leu183)**: Shifts by only **0.512 Å** in Cα position upon binding. Its sidechain rotamer is essentially pre-organized to dock into the hydrophobic cradle formed by TERT Leu769 and Pro771.
2. **The Dynamic DE-Loop (Glu168–Glu171)**: Undergoes extensive reorganization:
   - **Lys170**: Exhibits a peak Cα displacement of **4.293 Å**.
   - **Glu169**: Shifts by **3.809 Å**.
   - **Glu171**: Shifts by **3.320 Å**.
   - **Glu168**: Shifts by **2.748 Å** to engage Arg774 in the primary salt bridge.

### Structural Basis: Disorder-to-Order Transition
In apo `2I46`, the DE-loop exhibits extreme thermal motion ($B$-factors $>80\text{--}106\text{ \AA}^2$ in Chain A, completely disordered/unmodeled in Chain B). In holo `7TRE`, the loop is stabilized upon receptor engagement ($B$-factors drop to $\sim 50\text{--}61\text{ \AA}^2$). This demonstrates that TPP1 recruitment involves folding upon binding of an intrinsically disordered loop, rather than a classical rigid-body dock.

---

## 9. Mechanistic Evaluation of Allostery (Static vs Dynamic Coupling)

Earlier assertions that an RMSD of 0.338 Å proves "zero allosteric coupling" were critically audited and refuted:

### Why Static RMSD Cannot Disprove Allostery:
1. **Experimental Resolution Limit**: The cryo-EM maps of `7TRE` (3.5 Å) and `7TRD` (3.3 Å) have an estimated coordinate uncertainty (DPI / Cruickshank noise floor) of $0.35\text{--}0.50\text{ \AA}$. An RMSD of $0.338\text{ \AA}$ across 770 residues is within experimental noise.
2. **Dynamic & Entropic Allostery**: Conformational allostery can occur through altered low-frequency normal modes and vibrational entropy ($\Delta S_{\text{vib}}$) without altering the time-averaged mean atomic coordinates (Cooper & Dryden, *Nature* 1983).
3. **Biological Reality**: Published biochemical studies (Liu et al., 2022; Jansson et al., 2015) explicitly prove that TPP1 binding stimulates telomerase **Repeat Addition Processivity (RAP)** by dampening TEN–TRAP domain fluctuations and facilitating telomeric primer translocation. Thus, functional allosteric coupling is physically present despite mean-coordinate invariance.

---

## 10. Candidate Ligand Reality Check (TELIN-1 Provenance & Limits)

The small-molecule candidate `TELIN-1` (`telin1_inhibitor.pdb` and `7tre_telin1_complex.pdb`) was subjected to a thorough forensic audit:

### Audit Findings:
1. **Zero Prior Art**: 0 matches in PubChem, ChEMBL, CAS, PubMed, or patent records.
2. **Chemical Invalidity**: The coordinate file `telin1_inhibitor.pdb` is not a single compound: it consists of **two disconnected covalent fragments** (a 29-atom polycyclic core and an unbonded 7-atom pentanoic acid tail).
3. **Severe Steric Clashes**: Structural analysis of `7tre_telin1_complex.pdb` revealed **160 severe steric clashes** against TERT, with atom-atom interpenetrations down to **0.17 Å** (far below the van der Waals collision threshold).
4. **Zero Wet-Lab Evidence**: No synthesis protocol, no binding assay ($K_d$/$\text{IC}_{50}$), and no pharmacokinetic characterization.

**Disposition**: `TELIN-1` has been formally classified as an artifactual in silico construct and permanently quarantined.

---

## 11. Translation: The Cancer vs. Healthy Aging Dilemma

To evaluate clinical potential, the project maintains a strict separation between **Axis A (Molecular Mechanism)** and **Axis B (Therapeutic Selectivity)**.

### AXIS A: Molecular Mechanism
- The TERT–TPP1 interface is a nanomolar protein-protein interaction burying $\sim 857\text{ \AA}^2$ of surface area.
- The interface is mediated by an identical genetic sequence and identical biophysical contacts in both cancer cells and healthy human stem cells.

### AXIS B: Therapeutic Selectivity & Clinical Reality
- **The Kinetic Lag Paradox**: Cancer cells with $5\text{ kb}$ telomeres survive $40\text{--}50$ population doublings before reaching critical crisis length. In contrast, long-term hematopoietic stem cells (LT-HSCs) already possess short telomeres. Systemic telomerase shutdown triggers hematopoietic failure long before bulk tumor arrest.
- **The Hoyeraal–Hreidarsson Warning**: Human germline mutations in the TPP1 TEL-patch (e.g., deletion of Lys170) cause severe Hoyeraal–Hreidarsson syndrome and Dyskeratosis Congenita, characterized by fatal bone marrow failure, aplastic anemia, and early mortality.
- **Realistic Therapeutic Windows**:
  1. *Synthetic Lethality*: Targeting tumors with critically short telomeres using combinations with PARP or ATR/ATM inhibitors.
  2. *ALT Targeting*: Selectively targeting the Alternative Lengthening of Telomeres pathway (~10–15% of cancers), which is never utilized by healthy somatic or stem cells.
  3. *Transient mRNA Activation for Longevity*: Utilizing modified TERT mRNA with an intrinsic 48-hour half-life to expand stem cell telomeres transiently without conferring permanent oncogenic transformation.

---

## 12. Falsifiable Wet-Lab Experimental Protocols

Four concrete, reproducible experimental protocols were authored in `results/validation_experiments/falsifiable_predictions.md`:
1. **Stopped-Flow Fluorescence Kinetics**: Directly tests whether TPP1 binding follows single-step bimolecular kinetics or hyperbolic two-step induced fit ($H_2$).
2. **Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS)**: Measures conformational dynamics and entropic allostery in the TPP1 DE-loop and TERT catalytic cleft in apo vs. holo states.
3. **Single-Molecule FRET (smFRET)**: Monitors real-time opening/closing frequencies of the TEN–TRAP clamp upon TPP1 docking.
4. **Deep Mutational Scanning & RAP Assays**: Quantifies binding affinity changes ($\Delta\Delta G$) vs processivity stimulation ($\Delta\text{RAP}$) across TEL-patch variants (`L183A`, `E168A`, `K170A`, `R774A`).

---

## 13. Explicit Catalog of Overstated / Rejected Claims & Future Roadmap

To ensure permanent epistemic hygiene, this section catalogs all claims that were rejected, modified, or quarantined during this investigation:

### Table 13.1: Definitive Catalog of De-Hyped and Retracted Claims
| Original Claim | Initial Assertion | Corrected Finding / Audit Status | Rationale |
|:---|:---|:---|:---|
| **Claim 1** | "Two-step sequential induced-fit binding mechanism discovered." | **RETRACTED AS DEMONSTRATED FACT**. Reclassified as testable hypothesis ($H_2$). | Static endpoints (`7TRE` and `2I46`) cannot establish temporal reaction sequences. |
| **Claim 2** | "0.338 Å RMSD proves zero allosteric coupling between interface and core."| **REJECTED AS BIOPHYSICALLY INVALID**. | Ignores experimental coordinate uncertainty ($0.35\text{--}0.50\text{ \AA}$) and dynamic/entropic allostery. |
| **Claim 3** | "TELIN-1 is a validated nanomolar inhibitor lead for cancer." | **FORMALLY QUARANTINED**. | Non-physical construct with two disconnected fragments, 160 severe atomic clashes, and 0 wet-lab data. |
| **Claim 4** | "Complete TERT–TPP1 interface shutdown spares healthy stem cells." | **REFUTED BY CLINICAL GENETICS**. | Germline TPP1 TEL-patch deletion ($\Delta\text{K170}$) causes fatal Hoyeraal–Hreidarsson aplastic anemia. |
| **Claim 5** | "Computational pipeline discovered Arg774–Glu168 salt bridge." | **RECLASSIFIED AS POSITIVE CONTROL**. | Originally discovered and published by Liu et al. (*Nature* 2022). |

### Future Research Roadmap
1. Execute stopped-flow fluorescence kinetics to resolve the bimolecular vs induced-fit debate experimentally.
2. Perform HDX-MS profiling to map allosteric entropic pathways connecting the TEN–TRAP cradle to the RT catalytic pocket.
3. Advance transient LNP-delivered TERT mRNA therapies for ex vivo expansion of hematopoietic stem cells and acute regenerative rescue.
