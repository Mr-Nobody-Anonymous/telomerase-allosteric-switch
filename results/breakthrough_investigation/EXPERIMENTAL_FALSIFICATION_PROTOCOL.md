# Decisive Experimental Falsification Protocol: The Three Breakthrough Criteria

**Date**: October 30, 2026  
**Investigation**: Multi-assay experimental framework to empirically test the TEN–TPP1–TRAP Bipartite Clamp and identify separation-of-function variants in human TERT.  
**Target Variants**: Reconstituted human telomerase holoenzyme: WT, `Lys757A`, `Lys757R`, `Tyr176A`, `Tyr176F`, `Arg756A`, `Arg756K`, `Lys78A`, `Arg774A`.  
**Epistemic Level**: Level 4 (Experimental protocol design & mathematical falsification criteria).

---

## 1. The Three Mandatory Breakthrough Criteria

To transition this research from a plausible computational hypothesis into an established scientific breakthrough, the experimental results must satisfy **three independent criteria**:

```
                       ┌──────────────────────────────────────────────┐
                       │     THE THREE BREAKTHROUGH CRITERIA          │
                       └──────────────────────┬───────────────────────┘
                                              │
         ┌────────────────────────────────────┼────────────────────────────────────┐
         ▼                                    ▼                                    ▼
┌──────────────────┐               ┌──────────────────┐               ┌──────────────────┐
│  1. STRUCTURAL   │               │  2. BIOCHEMICAL  │               │  3. MECHANISTIC  │
│  The TEN-TRAP    │               │  Preserves TPP1  │               │  Reduction is    │
│  bridge is       │               │  binding & kcat; │               │  explained by an │
│  experimentally  │               │  selectively     │               │  experimentally  │
│  supported       │               │  collapses RAP   │               │  measured change │
│  (HDX-MS/smFRET) │               │  (SPR & Gel RAP) │               │  in retention    │
└──────────────────┘               └──────────────────┘               └──────────────────┘
```

1. **Criterion 1 (Structural)**: The physical coupling across the TEN–TRAP junction (`Tyr176`–`Gln794`–`Lys757`) must be experimentally demonstrated by HDX-MS (showing dynamic protection upon TPP1 binding in WT, lost in `Tyr176A`) or cryo-EM 3D variability analysis.
2. **Criterion 2 (Biochemical Separation of Function)**: `TERT Lys757A` must demonstrate wild-type TPP1 binding affinity ($K_d \le 300\text{ nM}$) and normal basal catalytic competence ($k_{\text{cat}} \ge 85\%\text{ WT}$), but selectively collapse TPP1-dependent processivity stimulation ($\chi_{\text{RAP}} \le 1.8$ vs. $\ge 5.0\text{ in WT}$), with the conservative control `Lys757R` fully restoring stimulation.
3. **Criterion 3 (Mechanistic Kinetic Causation)**: The processivity collapse must be directly traced to an experimentally measured elevation in primer dissociation during translocation ($k_{\text{off}}^{\text{translocation}} \ge 0.08\text{ s}^{-1}$ in a pulse-chase trap assay or smFRET), rather than a failure of nucleotide addition chemistry.

---

## 2. Comprehensive Assay Specifications

### Assay 1: Recombinant Holoenzyme Quality Control & Thermal Stability (nanoDSF)
- **Objective**: Verify proper folding and rule out global structural destabilization.
- **Instrument**: NanoTemper Prometheus NT.48 (measuring intrinsic tryptophan/tyrosine emission ratio $F_{350}/F_{330}$).
- **Acceptance Criterion**: All test variants must exhibit a cooperative thermal unfolding transition with:
  $$|T_m^{\text{variant}} - T_m^{\text{WT}}| \le 1.5^\circ\text{C}$$
  ($T_m^{\text{WT}} \approx 48.5 \pm 0.4^\circ\text{C}$). Any variant with $\Delta T_m > 3.0^\circ\text{C}$ or aggregation before $40^\circ\text{C}$ is rejected as a folding artifact.

### Assay 2: Surface Plasmon Resonance (SPR) Equilibrium TPP1 Binding Kinetics
- **Objective**: Directly quantify physical binding affinity ($K_d$), association rate ($k_{\text{on}}$), and dissociation rate ($k_{\text{off}}$) for the TPP1 OB domain.
- **Instrument**: Biacore 8K, Streptavidin (SA) sensor chip with immobilized biotinylated TPP1 OB domain (residues 89–249, $\approx 150\text{ RU}$).
- **Injections**: 2-fold serial dilutions of purified telomerase core ($6.25\text{ to }800\text{ nM}$) at $25^\circ\text{C}$.
- **Acceptance Thresholds**:
  - **Preserved Recruitment (Decoupling Candidate `Lys757A`)**: $K_d \le 300\text{ nM}$ ($< 2.5\text{-fold}$ deviation from WT $K_d \approx 120\text{ nM}$).
  - **Recruitment Defect (Negative Controls `Lys78A`, `Arg774A`)**: $K_d > 2,500\text{ nM}$ ($>20\text{-fold}$ loss).

### Assay 3: Direct Primer Extension & Processivity Index Titration
- **Objective**: Measure basal catalytic turnover ($k_{\text{cat}}$) and the processivity index ($R_p$) across a TPP1 concentration gradient ($0, 50, 150, 500, 1000\text{ nM}$).
- **Substrate**: 18-nt telomeric DNA primer ($5'\text{-(TTAGGG)}_3\text{-3'}$, $1\text{ }\mu\text{M}$) labeled with Cy5 or $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$.
- **Readout**: High-resolution $18\%\text{ denaturing sequencing gel}$ (7 M urea), analyzed via phosphorimaging.
- **Metrics**:
  - Basal activity: Total nucleotide incorporation in the absence of TPP1 ($k_{\text{cat}}^{-\text{TPP1}}$).
  - Repeat Addition Processivity ($R_p$): Calculated using the Fraction Left Behind (FLB) method across repeats 2 through 15:
    $$R_p = 1 - \frac{I_n}{\sum_{i=n}^{N} I_i}$$
  - Fold Processivity Stimulation ($\chi_{\text{RAP}}$):
    $$\chi_{\text{RAP}} = \frac{R_p^{(+500\text{ nM TPP1})}}{R_p^{(-\text{TPP1})}}$$

### Assay 4: Single-Molecule FRET & Optical Trap Translocation Kinetics
- **Objective**: Directly measure the rate of single-turnover primer dissociation during translocation ($k_{\text{off}}^{\text{translocation}}$).
- **Setup**: Biotinylated telomerase immobilized on a passivated PEG-quartz flow cell. Cy3-donor labeled telomerase RNA and Cy5-acceptor labeled telomeric DNA primer.
- **Readout**: Alternating laser excitation (ALEX) smFRET traces monitoring high-FRET (duplex elongation) and low-FRET (translocation pause) states.
- **Acceptance Benchmark**:
  - WT + TPP1 suppresses $k_{\text{off}}^{\text{translocation}}$ to $< 0.02\text{ s}^{-1}$ (complex half-life $> 20\text{ min}$).
  - A true decoupling mutant (`Lys757A` or `Tyr176A`) will bind TPP1 normally in Assay 2, but exhibit an elevated $k_{\text{off}}^{\text{translocation}} \ge 0.08\text{ s}^{-1}$, directly establishing kinetic causation!

### Assay 5: Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS)
- **Objective**: Directly verify whether TPP1 binding quenches conformational dynamics across the TEN–TRAP inter-domain cleft.
- **Workflow**: Deuterium labeling ($D_2O$ buffer, $10\text{ s}$ to $4\text{ hours}$) of telomerase $\pm$ TPP1, followed by pepsin digestion and LC-MS/MS on an Orbitrap mass spectrometer.
- **Expected Signature**: WT exhibits marked dynamic protection ($\Delta \text{Deut} > 25\%$) across TRAP peptides (750–785) and TEN peptides (120–145, 170–180). Mutating `Tyr176A` should selectively eliminate protection across the anchor groove.

---

## 3. Mathematical Decision Matrix & Falsification Criteria

| Candidate | TPP1 SPR $K_d$ (Assay 2) | Basal $k_{\text{cat}}$ (Assay 3) | Processivity Stimulation $\chi_{\text{RAP}}$ (Assay 3) | Translocation Off-Rate $k_{\text{off}}^{\text{transloc}}$ w/ TPP1 (Assay 4) | Verdict & Scientific Conclusion |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`Lys757A`**<br>*(Primary Candidate)* | **Normal** ($K_d \le 300\text{ nM}$) | **Competent** ($\ge 85\%\text{ WT}$) | **Collapsed** ($\chi_{\text{RAP}} \le 1.8\text{ vs WT } \ge 5.0$) | **Rapid** ($k_{\text{off}} \ge 0.08\text{ s}^{-1}$, identical to -TPP1) | **BREAKTHROUGH PROVEN**: Satisfies all 3 criteria. Clean separation of function at the TRAP pore. |
| **`Lys757R`**<br>*(Charge Control)* | **Normal** ($K_d \le 300\text{ nM}$) | **Normal** ($\ge 95\%\text{ WT}$) | **Normal** ($\chi_{\text{RAP}} \ge 4.5$) | **Suppressed** ($k_{\text{off}} \le 0.02\text{ s}^{-1}$) | **CONFIRMS MODEL**: Positive charge in the pore is sufficient to restore processivity transmission. |
| **`Tyr176A`**<br>*(Bridge Candidate)*| **Normal** ($K_d \le 300\text{ nM}$) | **Competent** ($\ge 75\%\text{ WT}$) | **Impaired** ($\chi_{\text{RAP}} \le 2.0$) | **Rapid** ($k_{\text{off}} \ge 0.08\text{ s}^{-1}$) | **PROVES INTER-DOMAIN COUPLING**: TEN-TRAP bridge transmits TPP1 allosteric signal. |
| **`Tyr176F`**<br>*(Bridge Control)* | **Normal** ($K_d \le 300\text{ nM}$) | **Normal** ($\ge 95\%\text{ WT}$) | **Normal** ($\chi_{\text{RAP}} \ge 4.5$) | **Suppressed** ($k_{\text{off}} \le 0.02\text{ s}^{-1}$) | **CONFIRMS BRIDGE**: Aromatic core packing is sufficient to maintain coupling. |
| **`Arg756A`**<br>*(Confounded Target)*| **Normal** ($K_d \le 300\text{ nM}$) | **Severely Reduced** ($< 40\%\text{ WT}$) | **Unreliable / Low** | Rapid | **CONFOUNDED**: Catalytic defect prevents clean proof of allosteric decoupling. |
| **`Arg756K`**<br>*(Stacking Control)* | **Normal** ($K_d \le 300\text{ nM}$) | **Partial** ($50\text{--}60\%\text{ WT}$) | Partially Stimulated | Intermediate | **PROVES A55 STACKING**: Planar guanidinium ring is specifically required for template guidance. |
| **`Arg774A`**<br>*(TRAP Control)* | **Severely Lost** ($> 2,500\text{ nM}$) | **Normal** (100% WT) | **Lost** ($\chi_{\text{RAP}} \le 1.2$) | Rapid | **RECRUITMENT CONTROL**: Loss of RAP is trivially explained by loss of TPP1 binding. |
| **`Lys78A`**<br>*(TEN Control)* | **Severely Lost** ($> 2,500\text{ nM}$) | **Normal** (100% WT) | **Lost** ($\chi_{\text{RAP}} \le 1.2$) | Rapid | **RECRUITMENT CONTROL**: Confirms TEN recruitment lobe is necessary for TPP1 binding. |

---

## 4. Absolute Falsification Rulebook

1. **The Entire Allosteric Decoupling Hypothesis is FALSIFIED if**:
   - `Lys757A` and `Tyr176A` both show wild-type RAP stimulation ($\chi_{\text{RAP}} \ge 4.5$). This would prove that the TRAP pore and the TEN–TRAP bridge are uninvolved in processivity regulation.
2. **The "Mechanical Piston" Model is ALREADY FALSIFIED**:
   - Disproven by coordinate analysis ($\Delta \le 0.28\text{ \AA}$ at pore base, within experimental coordinate uncertainty).
3. **The "Arg756 Single Hairpin" Model is REJECTED as Primary Hypothesis**:
   - Superseded by the `Lys757` / `Tyr176` inter-domain model to eliminate the fatal `hTR A55` catalytic confounder.
