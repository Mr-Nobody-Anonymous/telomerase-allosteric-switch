# SPR BINDING DATA & ANALYSIS PACKAGE: TERT–TPP1 RECRUITMENT INTERACTION

**Document ID:** `EXP-DOC-002-REV-FINAL`  
**Effective Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/SPR_BINDING_DATA_PACKAGE.md`  
**Classification:** Experimental Validation Protocol & Data Capture Template  

---

## 1. Objective & Hypothesis Demarcation

This protocol evaluates the equilibrium dissociation constant ($K_d$) and binding kinetics between the human telomerase catalytic core and the shelterin factor TPP1 OB domain for wild-type (WT) and mutant variants (`Lys757A` and `Lys757R`).

> [!IMPORTANT]
> **Pre-Registered Epistemic Hypothesis (Not Yet Measured)**:  
> Because `TERT Lys757` is located **25.43 Å** away from the nearest TPP1 atom in cryo-EM structure PDB 7TRE, the mutational hypothesis predicts that `TERT Lys757A` will retain high-affinity, nanomolar binding to TPP1 comparable to wild-type ($K_d \le 300\text{ nM}$). If `Lys757A` disrupts binding ($K_d > 1,500\text{ nM}$), the processivity uncoupling hypothesis is formally refuted.

---

## 2. Experimental Instrumentation and Sensor Setup

- **Instrument**: Biacore 8K, Biacore T200, or Carterra LSA.
- **Sensor Chip**: Sensor Chip SA (Streptavidin, Cytiva BR100531) or Sensor Chip CAP (Biotin CAPture Kit, Cytiva 28920234).
- **Running & Dilution Buffer**:
  - $20\text{ mM}$ HEPES-KOH, $\text{pH } 7.5$
  - $150\text{ mM}$ $\text{KCl}$
  - $2.0\text{ mM}$ $\text{MgCl}_2$
  - $0.05\%\text{ (v/v) Surfactant P20}$ (Polysorbate 20)
  - $1.0\text{ mM}$ $\text{TCEP}$
  - $0.1\text{ mg/mL}$ Bovine Serum Albumin (BSA, fraction V, protease-free)
- **Temperature**: $25.0\text{ }^\circ\text{C} \pm 0.1\text{ }^\circ\text{C}$.

---

## 3. Immobilization & Concentration Series

### 3.1 Ligand Immobilization (TPP1 OB-Fold Domain)
- **Construct**: Recombinant human TPP1 OB domain (residues 89–249) containing a C-terminal AviTag chemically or enzymatically mono-biotinylated via BirA.
- **Surface Density Target**: Immobilize to a response level of **$50\text{--}80\text{ RU}$** on Channel 2 (active surface). Channel 1 is left blank or coated with biotinylated free peptide as a reference surface.
- **Rationale**: Low immobilization density prevents mass-transport limitation (MTL) and steric crowding artifacts.

### 3.2 Analyte Preparation (Holoenzyme Complex)
- **Analytes**:
  1. Purified Wild-Type TERT Holoenzyme (WT)
  2. Purified `TERT Lys757A` Holoenzyme
  3. Purified `TERT Lys757R` Holoenzyme
  4. Established negative control: `TPP1 Glu168A/Glu169A` (TEL-patch recruitment dead mutant)
- **Concentration Titration Series**:
  - Seven 2-fold serial dilutions: **$0, 15, 30, 60, 120, 240, 480\text{ nM}$**.
  - All dilutions performed in running buffer in low-binding 96-well microplates.
  - Replicate injections: Include duplicate injections for $60\text{ nM}$ and $120\text{ nM}$ to confirm sensor surface stability.

### 3.3 Cycle Parameters
- Flow Rate: $30\text{ }\mu\text{L/min}$.
- Association Phase: $180\text{ seconds}$.
- Dissociation Phase: $600\text{ seconds}$.
- Regeneration: $10\text{ mM}$ Glycine-HCl, $\text{pH } 2.5$ for $30\text{ seconds}$ at $30\text{ }\mu\text{L/min}$, followed by $120\text{ seconds}$ running buffer stabilization.

---

## 4. Kinetic Fitting and Data Analysis Pipeline

### 4.1 Double-Referencing Subtraction
1. Subtract the reference channel response ($Ch_1$) from the active channel response ($Ch_2$):
   $$R_{\text{ref-sub}}(t) = R_{\text{active}}(t) - R_{\text{ref}}(t)$$
2. Subtract the zero-concentration buffer blank injection ($0\text{ nM}$):
   $$R_{\text{double-sub}}(t) = R_{\text{ref-sub}}(t) - R_{\text{blank}}(t)$$

### 4.2 Mathematical Kinetic Model (1:1 Langmuir with MTL)
The double-subtracted sensorgrams are fitted to a 1:1 binding model with mass transport limitation:
$$\frac{d[AB]}{dt} = k_{\text{on}}[A]_{\text{surface}} [B]_{\text{free}} - k_{\text{off}}[AB]$$
$$\frac{d[A]_{\text{surface}}}{dt} = k_t([A]_{\text{bulk}} - [A]_{\text{surface}}) - \frac{d[AB]}{dt}$$
$$K_d = \frac{k_{\text{off}}}{k_{\text{on}}}$$

### 4.3 Acceptance Criteria for Kinetic Fits
- Residual standard deviation: $\chi^2 \le 0.10 \times R_{\text{max}}$.
- Uniqueness parameter: $U\text{-value} \le 15$ (verifies independent convergence of $k_{\text{on}}$ and $k_{\text{off}}$).
- Standard error of fitted parameters: $< 5\%$ of parameter value.

---

## 5. Pre-Registered Decision Thresholds

```
========================================================================================================================
SPR BINDING DECISION MATRIX
========================================================================================================================
Construct       Expected Kd Metric         Passing Threshold                  Scientific Interpretation
------------------------------------------------------------------------------------------------------------------------
TERT WT         Kd = 150 ± 25 nM           100 nM <= Kd <= 200 nM             Validates native shelterin recruitment baseline.

TERT Lys757A    Predicted Kd <= 300 nM     Kd <= 300 nM                       CRITERION 2A PASSED: TPP1 recruitment affinity
                                                                              is fully preserved.

                Falsified if:              Kd > 1,500 nM                      HYPOTHESIS REFUTED: Mutation disrupts recruitment;
                                                                              not a separation-of-function target.

TERT Lys757R    Predicted Kd <= 300 nM     Kd <= 300 nM                       RESCUE CONTROL CONFIRMED: Charge-preserving mutant
                                                                              maintains wild-type binding affinity.
========================================================================================================================
```

---

## 6. Raw Data Capture Template

Data collected from the Biacore 8K instrument must be tabulated in the following standard format for automated pipeline fitting via `scripts/spr_kinetic_analysis.py`:

```csv
Replicate_ID,Construct,Concentration_nM,k_on_M_s,k_off_s,Kd_nM,Rmax_RU,Chi2_RU2,U_value,QC_Pass
REP-01,TERT_WT,480,[ENTER_KON],[ENTER_KOFF],[ENTER_KD],[ENTER_RMAX],[ENTER_CHI2],[ENTER_U],[YES/NO]
REP-01,TERT_Lys757A,480,[ENTER_KON],[ENTER_KOFF],[ENTER_KD],[ENTER_RMAX],[ENTER_CHI2],[ENTER_U],[YES/NO]
REP-01,TERT_Lys757R,480,[ENTER_KON],[ENTER_KOFF],[ENTER_KD],[ENTER_RMAX],[ENTER_CHI2],[ENTER_U],[YES/NO]
```
