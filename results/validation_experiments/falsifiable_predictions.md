# Falsifiable Wet-Lab Experimental Protocols for Telomerase–TPP1 Recruitment and Allostery

## Objective and Epistemic Standard

This protocol suite translates computational findings from experimental structures (**PDB: 7TRE, 7TRD, 7TRF, 7TRC, 2I46**) into concrete, reproducible wet-lab experimental assays. 

Rather than presenting computational predictions as established facts, each protocol defines:
- **Hypothesis under test** (linking directly to $H_1$, $H_2$, or $H_3$ in `results/research/question_and_hypotheses.md`).
- **Construct specifications and controls** (positive, negative, and vehicle controls).
- **Assay parameters, instrumentation, and readouts**.
- **Quantitative thresholds for hypothesis confirmation vs. unequivocal falsification**.

---

## Protocol 1: Time-Resolved Binding Kinetics via Stopped-Flow Fluorescence & Surface Plasmon Resonance (SPR)

### 1.1 Objective
Test whether TPP1 binding to the TERT TEN–TRAP receptor cradle follows a **one-step bimolecular association** or a **two-step induced-fit mechanism** ($A + B \underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}} AB \underset{k_{-2}}{\overset{k_2}{\rightleftharpoons}} AB^*$), addressing the retracted claim `CLM-021` ("two-step sequential induced fit").

### 1.2 Construct Design & Labeling
- **Recombinant TPP1 OB-fold (residues 89–249)**:
  - Wild-type (WT) construct.
  - Fluorophore-labeled reporter: Site-specific cysteine insertion at position Ser165 (`S165C`) or Tyr180 (`Y180C`) labeled with maleimide-conjugated **Cy3** or **Alexa Fluor 488**.
- **Human Telomerase Core RNP (TERT + TERC)**:
  - Reconstituted in vitro via insect cell / rabbit reticulocyte lysate (RRL) systems or HEK293T overexpression, purified via tandem affinity tag.
  - Quencher construct: TERT labeled with **Cy5** or **Dabcyl** near the TEN domain cradle (e.g. `S770C`).
- **Negative Controls**:
  - TPP1 TEL-patch alanine triple mutant: `E168A / E169A / E171A`.
  - TPP1 hydrophobic anchor null mutant: `L183A`.
  - Mock reconstitution lacking TERT protein.
- **Positive Control**:
  - Full-length shelterin subcomplex (POT1–TPP1) binding to immobilized $(TTAGGG)_4$ primer.

### 1.3 Methodology
1. **Instrumentation**: Applied Photophysics SX20 stopped-flow spectrophotometer equipped with dual-syringe mixing ($1:1$ ratio, dead time $<1.2\text{ ms}$) at $25.0^\circ\text{C}$ and $37.0^\circ\text{C}$.
2. **Buffer Conditions**: $25\text{ mM}$ HEPES-KOH (pH 7.5), $150\text{ mM}$ KCl, $2\text{ mM}$ $\text{MgCl}_2$, $1\text{ mM}$ TCEP, $0.05\%$ Tween-20, $5\%$ glycerol.
3. **Concentration Series**: Constant labeled TPP1 ($10\text{ nM}$) mixed rapidly against excess pseudo-first-order concentrations of TERT RNP ($50\text{ nM}$ to $2.0\text{ }\mu\text{M}$, 8 concentration steps).
4. **Data Acquisition**: 2,000 logarithmic time points collected from $1.5\text{ ms}$ to $10.0\text{ s}$. 8 replicates averaged per concentration.

### 1.4 Kinetic Analysis & Falsification Criteria
- Fit observed relaxation rates ($k_{\text{obs}}$) as a function of [TERT]:
  - **One-Step Bimolecular Model**: $k_{\text{obs}} = k_{\text{on}}[\text{TERT}] + k_{\text{off}}$ (strictly linear dependence).
  - **Two-Step Induced-Fit Model**: Hyperbolic saturation curve:
    $$k_{\text{obs}} = k_{-2} + \frac{k_2 [\text{TERT}]}{K_{1} + [\text{TERT}]}$$
- **Falsification Thresholds**:
  - **Falsification of Induced-Fit ($H_2$)**: If $k_{\text{obs}}$ vs. [TERT] displays linear dependence ($R^2 > 0.98$) across the entire $50\text{ nM} - 2\text{ }\mu\text{M}$ range with no saturation, the two-step sequential induced-fit hypothesis is **FALSIFIED**; binding is an uncoupled single-step collision.
  - **Confirmation of Induced-Fit ($H_2$)**: Requires statistically significant hyperbolic curvature ($\Delta \text{AIC} > 10$ in favor of two-step model) with a measurable limiting isomerization rate constant $k_2 \ge 10\text{ s}^{-1}$.

---

## Protocol 2: Hydrogen-Deuterium Exchange Mass Spectrometry (HDX-MS) of TPP1 & TERT Core

### 1.1 Objective
Directly test whether TPP1 binding modulates conformational dynamics and vibrational entropy ($\Delta S_{\text{vib}}$) in the distant TERT catalytic core (fingers/palm/thumb), addressing `CLM-022` ("zero allosteric coupling").

### 1.2 Samples & State Comparisons
- **State 1**: Apo TPP1 OB-fold ($10\text{ }\mu\text{M}$).
- **State 2**: Apo Telomerase Catalytic Core ($10\text{ }\mu\text{M}$ TERT + TERC).
- **State 3**: Holo TERT–TERC–TPP1 complex ($10\text{ }\mu\text{M}$ equimolar complex, pre-equilibrated).
- **State 4**: Complex with TPP1 $\Delta\text{K170}$ mutant (Hoyeraal–Hreidarsson control).

### 1.3 Methodology
1. **On-Exchange Reaction**: Dilute protein $1:10$ into deuterated buffer ($25\text{ mM}$ HEPES-KOH, $\text{pD } 7.5$, $150\text{ mM}$ KCl in $99.9\%\text{ D}_2\text{O}$) at $4.0^\circ\text{C}$.
2. **Timepoints**: Quench at $10\text{ s}$, $30\text{ s}$, $1\text{ min}$, $5\text{ min}$, $30\text{ min}$, $2\text{ hr}$, and $8\text{ hr}$.
3. **Quench & Digestion**: Add ice-cold quench buffer ($100\text{ mM}$ potassium phosphate, $2\text{ M}$ guanidine hydrochloride, $200\text{ mM}$ TCEP, pH 2.3) bringing final pH to $2.5$ at $0.0^\circ\text{C}$. Online digestion through an immobilized pepsin/protease XIII column at $100\text{ }\mu\text{L/min}$.
4. **LC-MS/MS**: Desalt on a trap column and separate on a reverse-phase C18 column ($0^\circ\text{C}$ Peltier cooling) within $8\text{ minutes}$ connected to an Orbitrap high-resolution mass spectrometer.
5. **Data Processing**: Deuteration levels per peptic peptide extracted using DynamX 3.0. Back-exchange corrected using a fully deuterated control sample.

### 1.4 Falsification Criteria
- **Falsification of Dynamic Allostery**: If no peptide within the TERT catalytic core (fingers/palm residues 560–930) displays a statistically significant difference in deuterium uptake between State 2 and State 3 ($\Delta \text{D} < 0.20\text{ Da}$, $p > 0.05$ after FDR correction across all timepoints), the hypothesis of allosteric dynamic coupling across the interface is **FALSIFIED**.
- **Support for Dynamic Entropic Allostery**: Identification of protected or deprotected peptides ($\Delta \text{D} \ge 0.50\text{ Da}$, $p < 0.01$) in TERT regions $>25\text{ \AA}$ distant from the TPP1 binding site, accompanied by protection of the TPP1 DE-loop (residues 166–175) confirming rigidification from the disordered apo ensemble.

---

## Protocol 3: Single-Molecule FRET (smFRET) Analysis of TEN–TRAP Domain Fluctuations

### 1.1 Objective
Measure real-time conformational fluctuations and opening/closing dynamics between the TERT TEN domain and the TRAP cradle before and after TPP1 binding, testing whether TPP1 dampens or stabilizes specific hinge conformations.

### 1.2 Construct & Surface Immobilization
- **TERT Dual-Labeled Reporter**: Site-specific orthogonal amber-suppression labeling:
  - TEN domain: $p$-azido-L-phenylalanine at position Thr85 labeled with DBCO-Cy3 (donor).
  - TRAP cradle: Cysteine mutant at position Ser770 labeled with maleimide-Cy5 (acceptor).
  - Theoretical inter-dye distance in 7TRE: $R \approx 42\text{ \AA}$ (ideal Förster radius $R_0 = 54\text{ \AA}$ for Cy3–Cy5, giving predicted FRET efficiency $E_{\text{FRET}} \approx 0.72$).
- **Surface Tethering**: Biotinylated DNA primer $(TTAGGG)_3$ immobilized on a PEG/biotin-PEG-passivated quartz microfluidic flow cell via neutravidin.

### 1.3 Data Acquisition
- Prism-based Total Internal Reflection Fluorescence (TIRF) microscopy at 100 Hz frame rate ($10\text{ ms}$ temporal resolution) under oxygen scavenging conditions ($2\text{ mM}$ Trolox, $5\text{ mM}$ PCA, $50\text{ nM}$ PCD).
- Trace collection: $>500$ individual single-molecule trajectories per experimental condition.

### 1.4 Falsification Criteria
- **Falsification of Conformational Selection / Hinge Dampening**: If the FRET histogram of apo-TERT shows a single rigid Gaussian peak with full-width at half-maximum ($\text{FWHM} \le 0.10$) identical to holo-TERT ($E_{\text{apo}} = E_{\text{holo}}$), without any observable transitions in dwell-time analysis, the model of TPP1-mediated conformational dampening is **FALSIFIED**.
- **Support**: Bimodal or broad distribution in apo-TERT ($E \approx 0.35$ and $E \approx 0.72$) shifting to a narrow, high-FRET state ($E = 0.74 \pm 0.04$) upon addition of $100\text{ nM}$ TPP1, with transition frequencies $\tau_{\text{open}\to\text{closed}}$ dependent on TPP1 concentration.

---

## Protocol 4: Deep Mutational Scanning & Repeat Addition Processivity (RAP) Assays

### 1.1 Objective
Quantify the energetic contribution of each TEL-patch and TERT contact residue to overall binding affinity ($\Delta\Delta G_{\text{bind}}$) versus catalytic processivity stimulation ($\Delta\text{RAP}$).

### 1.2 Target Variants
- **TPP1 TEL-patch variants**: `E168A`, `E169A`, `K170A`, `E171A`, `L183A`, `L183F`, `L183V`, `L212A`, `E215A`.
- **TERT cradle variants**: `R774A`, `R774K`, `R774E`, `K78A`, `Q775A`, `L769A`, `P771A`.
- **Compound Controls**: TPP1 WT, TPP1 $\Delta\text{K170}$ (HHS positive control).

### 1.3 Primer Extension & TRAP Assay Methodology
1. **Reaction Setup**: $50\text{ fmol}$ purified telomerase RNP, $1.0\text{ }\mu\text{M}$ telomeric DNA primer `5'-(TTAGGG)3-3'`, $1.0\text{ mM}$ dATP, $1.0\text{ mM}$ dTTP, $10\text{ }\mu\text{M}$ dGTP, and $0.15\text{ }\mu\text{M}$ $[\alpha-^{32}\text{P}]$-dGTP ($3000\text{ Ci/mmol}$) in telomerase reaction buffer.
2. **TPP1 Titration**: Recombinant TPP1 variants titrated from $1\text{ nM}$ to $5\text{ }\mu\text{M}$.
3. **Electrophoresis**: Stop reactions at $60\text{ min}$ with phenol/chloroform, ethanol-precipitate DNA, and resolve on a $12\%$ polyacrylamide / $7\text{ M}$ urea denaturing sequencing gel.
4. **Phosphorimager Quantification**: Expose to phosphor-storage screen; quantify Repeat Addition Processivity (RAP) according to standard decay formula:
   $$\text{RAP} = \frac{\sum_{i=2}^N I_i \cdot i}{\sum_{i=1}^N I_i}$$
   where $I_i$ is the intensity of the $i$-th telomeric repeat band.

### 1.4 Expected Invariants & Falsification Criteria
- **Positive Control Concordance**:
  - `TPP1_L183A` must reduce binding affinity by $\Delta\Delta G \ge +2.0\text{ kcal/mol}$ ($>30$-fold increase in $K_d$).
  - `TPP1_K170A` and `TPP1_DeltaK170` must severely impair processivity ($\text{RAP} < 20\%$ of WT) while retaining measurable basal binding ($K_d \le 10\text{ }\mu\text{M}$).
- **Falsification of the Essential Salt Bridge**: If `TERT_R774A` or `TPP1_E168A` retains wild-type recruitment affinity and processivity stimulation ($>80\%$ of WT), the model placing `Arg774–Glu168` as an obligate electrostatic recruitment anchor is **FALSIFIED**.
