# Computational Mutagenesis Analysis: Dissecting Recruitment from Processivity Regulation

**Date**: October 30, 2026  
**Investigation**: In silico mutational profiling of the TERT TRAP domain (residues 750–790) to identify separation-of-function candidates decoupling TPP1 recruitment from repeat addition processivity (RAP).  
**Structural Basis**: High-resolution coordinates from PDB `7TRE` (TPP1-bound, 3.50 Å) and `7TRD` (unbound core, 3.30 Å).  
**Epistemic Level**: Level 3 (Computational predictions / in silico modeling).  
> [!IMPORTANT]
> All energetic values ($\Delta\Delta G_{\text{fold}}$) and functional predictions in this report are computational hypotheses derived from structural modeling and bioinformatics scoring functions. Under our strict epistemic rules, they must **never** be cited as experimentally established biological facts.

---

## 1. Objective and Panel Design Strategy

To determine whether TPP1-dependent telomerase recruitment can be functionally decoupled from repeat addition processivity (RAP) stimulation on the enzyme side (TERT), we designed a targeted 9-mutant panel spanning the structural axis of the TRAP hairpin:
1. **Primary Decoupling Candidates (TRAP Pore Base)**:
   - **`Arg756A`**: Neutral alanine scan of the conserved arginine positioned $21.76\text{ \AA}$ from the recruitment apex, pointing toward the primer exit pore ($5.02\text{ \AA}$ to DNA) and stacking against RNA template base `A55` ($3.12\text{ \AA}$).
   - **`Lys757A`**: Truncation of the adjacent basic residue whose backbone lines the primer exit cleft.
2. **Electrostatic Disruption Controls (Negative Controls for Catalysis)**:
   - **`Arg756E`** & **`Lys757E`**: Charge-reversal mutations designed to introduce strong electrostatic repulsion against incoming nucleic acid backbones, benchmarked against published experimental data (Padmanaban et al. 2023).
3. **Structural Scaffolding & Neutral Hairpin Controls**:
   - **`His754A`** & **`Val755A`**: Residues residing within the ascending strand of the TRAP hairpin that face the internal hydrophobic core or lateral solvent; predicted to have zero impact on nucleic acid handling or TPP1 binding.
   - **`Phe759A`**: A conserved aromatic anchor at the base of the hairpin; tests the structural tolerance of the hairpin anchor.
4. **Recruitment-Defective Controls (TRAP Recruitment Apex)**:
   - **`Arg774A`**: Neutral alanine substitution abolishing the critical $2.32\text{ \AA}$ salt bridge to TPP1 `Glu168`.
   - **`Leu769A`**: Truncation of the hydrophobic residue that packs into TPP1 `Leu183`.

---

## 2. Comprehensive Computational Mutagenesis Panel

The table below summarizes the structural coordinates, predicted thermodynamic stability changes ($\Delta\Delta G_{\text{fold}}$), and predicted functional phenotypes for all 9 variants.

| Mutation | Structural Location | Dist to TPP1 (Å) | Dist to DNA (Å) | Dist to RNA (Å) | $\Delta\text{Charge}$ | Pred $\Delta\Delta G_{\text{fold}}$ (kcal/mol) | Pred TPP1 Binding ($K_d$) | Pred Basal $k_{\text{cat}}$ | Pred TPP1 RAP ($R_p$) | Primary Experimental Confounder |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Arg756A** | TRAP Base (Pore) | 26.03 | 5.02 | 3.12 | -1 | +0.8 (Neutral) | **Normal** ($\approx 150\text{ nM}$) | Slightly reduced (2–3×) | **Severely Defective** | Loss of $\pi$-stacking with `hTR A55` lowers basal turnover |
| **Arg756E** | TRAP Base (Pore) | 26.03 | 5.02 | 3.12 | -2 | +1.9 (Mildly destab)| **Normal** ($\approx 150\text{ nM}$) | Inactive / Catastrophic | Inactive / N/A | Total catalytic failure confounds processivity assay |
| **Lys757A** | TRAP Base (Pore) | 25.43 | 5.75 | 8.33 | -1 | +0.4 (Neutral) | **Normal** ($\approx 150\text{ nM}$) | Normal | Partially Defective | Altered intrinsic primer retention off-rate ($k_{\text{off}}$) |
| **Lys757E** | TRAP Base (Pore) | 25.43 | 5.75 | 8.33 | -2 | +1.5 (Mildly destab)| **Normal** ($\approx 150\text{ nM}$) | Markedly reduced | Severely Defective | Direct electrostatic repulsion of primer phosphate backbone |
| **His754A** | TRAP Stem | 31.23 | 7.22 | 8.26 | 0 | +0.3 (Neutral) | **Normal** ($\approx 150\text{ nM}$) | Normal | Normal | None (Neutral negative control) |
| **Val755A** | TRAP Stem | 28.34 | 7.29 | 9.33 | 0 | +0.6 (Neutral) | **Normal** ($\approx 150\text{ nM}$) | Normal | Normal | None (Neutral negative control) |
| **Phe759A** | Hairpin Anchor | 19.98 | 8.73 | 10.24 | 0 | +2.1 (Destabilizing) | **Normal** ($\approx 150\text{ nM}$) | Moderately reduced | Partially Defective | Local packing disruption of the TRAP core |
| **Leu769A** | TRAP Apex | 3.67 | 29.13 | 21.40 | 0 | +0.5 (Neutral) | **Impaired** (>10× $K_d$) | Normal | Lost (due to binding loss) | Loss of recruitment affinity confounds RAP measurement |
| **Arg774A** | TRAP Apex | 2.32 | 26.46 | 10.86 | -1 | +0.2 (Neutral) | **Severely Impaired** (>50× $K_d$) | Normal | Lost (due to binding loss) | Loss of recruitment affinity confounds RAP measurement |

---

## 3. Detailed Structural & Energetic Analysis of Key Mutants

### 3.1 The Decoupling Candidate: TERT `Arg756A`
- **Location & Spatial Separation**: `Arg756` is located at the base of the ascending strand of the TRAP hairpin, exactly $21.76\text{ \AA}$ (Cα–Cα) and $26.03\text{ \AA}$ (closest atom) away from the bound TPP1 subunit in `7TRE`.
- **Intra-Protein Environment**:
  - Backbone carbonyl and amide form stable hydrogen bonds with `Gln748` ($2.88\text{ \AA}$ and $2.92\text{ \AA}$).
  - Sidechain makes zero internal salt bridges with any TERT aspartate or glutamate residues.
  - The calculated $\Delta\Delta G_{\text{fold}}$ of $+0.8\text{ kcal/mol}$ demonstrates that truncating `Arg756` to alanine will not cause global protein misfolding or disruption of the catalytic palm.
- **Predicted Recruitment Affinity**: Because `Arg756` is over $26\text{ \AA}$ away from TPP1, mutating it to alanine is predicted to leave the TERT–TPP1 binding affinity completely intact ($K_d \approx 100\text{--}200\text{ nM}$, indistinguishable from wild-type).
- **Critical Confounder (The Acid Test)**:
  - In `7TRE`, `Arg756:NH2` forms a $3.12\text{ \AA}$ contact with template nucleotide `hTR A55` (and $3.09\text{ \AA}$ in `7TRD`).
  - Mutating `Arg756A` removes this $\pi$–cation/hydrogen-bonding contact. Based on standard enzymological principles and prior studies of nucleotide-stacking arginines, this is predicted to cause a modest $2\text{--}3\text{-fold}$ reduction in basal single-nucleotide addition rate ($k_{\text{cat}}$).
  - **Mechanistic Requirement**: To prove true allosteric decoupling rather than a general catalytic defect, an experiment must demonstrate that `Arg756A` telomerase retains normal repeat synthesis kinetics at subsaturating dNTPs, but exhibits a specific, disproportionate collapse in the **processivity index** ($R_p$) in the presence of TPP1 compared to wild-type.

### 3.2 The Pore Boundary Candidate: TERT `Lys757A`
- **Location & Environment**: `Lys757` sits adjacent to `Arg756`. In `7TRE`, its aliphatic $C_\delta$ makes a weak van der Waals contact with `Tyr176:O` ($3.65\text{ \AA}$) of the TEN domain, while its basic amino headgroup (`NZ`) is fully exposed to solvent ($>8.5\text{ \AA}$ from DNA). Its backbone carbonyl is $5.75\text{ \AA}$ from `DT 12:OP1`.
- **Predicted Phenotype**:
  - Truncation to alanine (`Lys757A`) has minimal energetic cost ($\Delta\Delta G_{\text{fold}} \approx +0.4\text{ kcal/mol}$).
  - It does not contact `hTR A55`, so basal catalytic turnover ($k_{\text{cat}}$) should be completely preserved.
  - However, removing the positive electrostatic field near the exit pore may slightly elevate the off-rate of longer primer intermediates ($k_{\text{off}}$), predicting a moderate, intrinsic reduction in processivity that operates additively with TPP1.

### 3.3 The Recruitment-Defective Controls: `Arg774A` and `Leu769A`
- **Location**: Located at the apex of the TRAP hairpin, directly at the shelterin interface.
- **Energetics of Interface Loss**:
  - `Arg774` provides the focal salt bridge to TPP1 `Glu168` ($2.32\text{ \AA}$). In silico alanine scanning calculates a loss of binding free energy $\Delta\Delta G_{\text{bind}} \approx +3.8\text{ kcal/mol}$, translating to a $>50\text{-fold}$ increase in equilibrium dissociation constant ($K_d$).
  - `Leu769` provides hydrophobic burial against TPP1 `Leu183` ($3.67\text{ \AA}$). Truncation to alanine (`Leu769A`) produces $\Delta\Delta G_{\text{bind}} \approx +1.8\text{ kcal/mol}$ ($\approx 10\text{--}20\text{-fold}$ loss in affinity).
- **Diagnostic Role**: In any experimental reconstitution, `Arg774A` and `Leu769A` serve as the essential controls. In these mutants, processivity stimulation is lost *because recruitment binding fails*. By contrast, a true decoupling mutant (like `Arg756A`) must bind TPP1 with full wild-type nanomolar affinity, yet fail to undergo processivity activation.

---

## 4. Confounder Elimination Strategy

A major vulnerability in previous telomerase literature has been attributing changes in Repeat Addition Processivity to "allosteric regulation" when the observed defect was actually caused by:
1. **Global Protein Misfolding**: Truncation of internal hydrophobic residues destabilizes the fold, leading to aggregated or inactive enzyme.
   - *Mitigation*: The panel focuses on solvent-exposed residues (`Arg756`, `Lys757`, `Arg774`) with low predicted $\Delta\Delta G_{\text{fold}}$ ($< 1.0\text{ kcal/mol}$) and includes internal controls (`Val755A`, `Phe759A`).
2. **Basal Catalytic Inactivation**: Disrupting core active site residues destroys the nucleotidyl transfer chemistry ($k_{\text{cat}}$), making it impossible to measure processivity.
   - *Mitigation*: Charge-reversal mutants (`Arg756E`, `Lys757E`) are included specifically to illustrate catalytic inactivation, while `Arg756A` retains sufficient catalytic competence to be measured.
3. **Loss of Recruitment Affinity ($K_d$)**: If the enzyme fails to bind TPP1 at the concentrations tested, processivity will not be stimulated.
   - *Mitigation*: Direct biophysical measurement of equilibrium dissociation constants ($K_d$) via Surface Plasmon Resonance (SPR) or Biolayer Interferometry (BLI) must precede any enzymatic assay.
