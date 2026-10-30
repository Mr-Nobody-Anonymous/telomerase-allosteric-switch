# Experimental Validation Protocol: Falsification Framework for TERT TRAP Decoupling

**Date**: October 30, 2026  
**Investigation**: Laboratory protocol to test whether TERT TRAP base mutation `Arg756A` decouples TPP1 recruitment binding from repeat addition processivity (RAP) stimulation.  
**Target Enzymes**: Reconstituted human telomerase holoenzyme variants (WT, `Arg756A`, `Arg756E`, `Lys757A`, `His754A`, `Leu769A`, `Arg774A`).  
**Epistemic Level**: Level 4 (Experimental protocol & falsification design).  
> [!IMPORTANT]
> This protocol specifies the definitive wet-lab biochemical assays required to empirically confirm or falsify the proposed allosteric decoupling mechanism. No computational model can substitute for these physical measurements.

---

## 1. Experimental Overview and Objectives

The central scientific question is:
> *Can human telomerase retain normal nanomolar binding affinity to shelterin TPP1 while selectively losing TPP1's ability to stimulate repeat addition processivity?*

To answer this conclusively without confounding by folding defects, catalytic inactivation, or binding failure, the experimental workflow executes four independent, orthogonal assays:
```
[Recombinant Holoenzyme Expression & Purification]
                      │
     ┌────────────────┼────────────────┐
     ▼                ▼                ▼
[Assay 1: nanoDSF/CD] [Assay 2: SPR]  [Assay 3: Direct RAP]
 (Folding & Stability) (TPP1 Affinity) (Processivity Index)
     │                │                │
     └────────────────┼────────────────┘
                      ▼
            [Assay 4: Primer Chase]
           (Dissociation Kinetics)
                      ▼
     [Rigorous Falsification Logic Matrix]
```

---

## 2. Detailed Assay Protocols

### Assay 1: Holoenzyme Reconstitution and Biophysical Stability (nanoDSF & Far-UV CD)
- **Objective**: Verify that test mutations do not cause global misfolding, thermal destabilization, or structural collapse of the TERT catalytic core.
- **Methodology**:
  - Reconstitute human telomerase in insect cells (Sf9/baculovirus) or human suspension cells (Expi293F) co-expressing full-length TERT, full-length hTR, and accessory proteins (dyskerin, NOP10, NHP2, GAR1).
  - Purify reconstituted complexes via anti-FLAG affinity chromatography followed by size-exclusion chromatography (Superose 6 Increase 10/300 GL) in $20\text{ mM HEPES-KOH (pH 7.9)}$, $150\text{ mM KCl}$, $1\text{ mM MgCl}_2$, $10\%\text{ glycerol}$, $1\text{ mM DTT}$.
  - **NanoDSF (Differential Scanning Fluorimetry)**: Measure intrinsic tryptophan/tyrosine fluorescence ($F_{350}/F_{330}$) across a thermal ramp ($20^\circ\text{C}$ to $95^\circ\text{C}$ at $1^\circ\text{C/min}$) on a Prometheus NT.48 instrument. Determine apparent melting temperature ($T_m$).
  - **Far-UV Circular Dichroism (CD)**: Record spectra between $200\text{--}260\text{ nm}$ on a Jasco J-815 spectropolarimeter at $20^\circ\text{C}$ in a $1\text{-mm}$ pathlength cuvette.
- **Acceptance Criteria**:
  - Recombinant mutant complexes must exhibit a single cooperative thermal unfolding transition with $\Delta T_m \le 1.5^\circ\text{C}$ relative to wild-type ($T_m^{\text{WT}} \approx 48.5^\circ\text{C}$).
  - Any mutant with $\Delta T_m > 3.0^\circ\text{C}$ or non-cooperative unfolding must be flagged as a folding artifact and excluded from mechanistic claims.

---

### Assay 2: Direct Equilibrium TERT–TPP1 Binding Affinity (Surface Plasmon Resonance / BLI)
- **Objective**: Directly measure the physical equilibrium dissociation constant ($K_d$) and kinetic rate constants ($k_{\text{on}}, k_{\text{off}}$) between TERT variants and the TPP1 OB-fold domain.
- **Methodology**:
  - Instrument: Biacore 8K (SPR) or Octet HTX (Biolayer Interferometry).
  - Immobilization: Chemically biotinylated wild-type human TPP1 OB domain (residues 89–249) immobilized on a Streptavidin (SA) biosensor chip to an immobilization level of $\approx 150\text{ RU}$.
  - Reference Surface: Sensor channel immobilized with inactive TPP1 mutant (`Glu168A / Leu183A`) or free biocytin to subtract bulk refractive index changes and non-specific binding.
  - Analyte Injections: Purified telomerase holoenzyme injected at concentrations ranging from $6.25\text{ nM}$ to $800\text{ nM}$ (2-fold serial dilutions, 6 concentrations in duplicate) in running buffer ($20\text{ mM HEPES-KOH pH 7.9}$, $150\text{ mM KCl}$, $1\text{ mM MgCl}_2$, $0.05\%\text{ Tween-20}$, $0.5\text{ mg/mL BSA}$).
  - Contact time: $180\text{ s}$; dissociation time: $600\text{ s}$ at $25^\circ\text{C}$.
  - Curve fitting: Global 1:1 Langmuir binding model with mass-transport limitation correction.
- **Quantitative Benchmark Expectations**:
  - **Wild-Type (WT)**: $K_d = 120 \pm 25\text{ nM}$; $k_{\text{on}} \approx 4.5 \times 10^4\text{ M}^{-1}\text{s}^{-1}$; $k_{\text{off}} \approx 5.5 \times 10^{-3}\text{ s}^{-1}$.
  - **Recruitment Control (`Arg774A`)**: $K_d > 5,000\text{ nM}$ ($>40\text{-fold}$ loss; no appreciable binding at $400\text{ nM}$).
  - **Decoupling Candidate (`Arg756A`)**: $K_d \le 300\text{ nM}$ ($<2.5\text{-fold}$ deviation from WT; normal recruitment binding).

---

### Assay 3: Direct Primer Extension & Repeat Addition Processivity (RAP) Assay
- **Objective**: Quantify basal single-repeat catalytic efficiency ($k_{\text{cat}}$) and TPP1-stimulated Repeat Addition Processivity ($R_p$) under steady-state and multi-repeat extension conditions.
- **Methodology**:
  - Substrate: 18-nt telomeric DNA primer ($5'\text{-(TTAGGG)}_3\text{-3'}$, $1\text{ }\mu\text{M}$) labeled at the 5' end with Cy5 or reconstituted with $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ incorporation.
  - Reaction Buffer: $50\text{ mM Tris-HCl (pH 8.0)}$, $50\text{ mM KCl}$, $1\text{ mM MgCl}_2$, $1\text{ mM spermidine}$, $5\text{ mM }\beta\text{-mercaptoethanol}$, $1\text{ mM dATP}$, $1\text{ mM dTTP}$, $10\text{ }\mu\text{M dGTP}$ ($+ [\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ $10\text{ }\mu\text{Ci}$).
  - TPP1 Titration: Reactions performed in parallel containing purified TPP1 OB domain at:
    $$[\text{TPP1}] = 0,\; 50,\; 150,\; 500,\; 1000\text{ nM}$$
  - Incubation: $30^\circ\text{C}$ for $60\text{ min}$. Terminated with stop buffer containing $20\text{ mM EDTA}$, $10\text{ mM Tris (pH 8.0)}$, $0.2\%\text{ SDS}$, and a $5'\text{-}^{32}\text{P}$-labeled 12-nt loading control.
  - Analysis: Products resolved on an $18\%\text{ polyacrylamide / 7 M urea}$ denaturing sequencing gel at $65\text{ W}$ for $2.5\text{ hours}$. Visualized by phosphorimaging (Typhoon FLA 9500) and quantified with ImageQuant TL.
- **Mathematical Formulations & Metrics**:
  1. **Basal Activity ($k_{\text{cat}}$)**:
     $$\text{Total Activity} = \sum_{i=1}^{N} \frac{I_i}{n_i}$$
     where $I_i$ is band intensity of repeat $i$, and $n_i$ is the number of incorporated dGTP residues in repeat $i$.
  2. **Processivity Index ($R_p$)**:
     Calculated using the decay ratio across repeats 2 through 15:
     $$R_p = 1 - \frac{I_n}{\sum_{i=n}^{N} I_i}$$
  3. **Fold-Stimulation of Processivity ($\chi_{\text{RAP}}$)**:
     $$\chi_{\text{RAP}} = \frac{R_p^{(+\text{TPP1})}}{R_p^{(-\text{TPP1})}}$$

---

### Assay 4: Primer Retention & Single-Turnover Trap Chase Assay
- **Objective**: Measure the physical dissociation rate of the DNA primer ($k_{\text{off}}^{\text{primer}}$) during repeat synthesis in the presence vs. absence of TPP1 to distinguish kinetic trap stabilization from true catalytic processivity.
- **Methodology**:
  - Pre-incubate telomerase ($20\text{ nM}$) with 5'-fluorophore-labeled telomeric primer ($20\text{ nM}$) in reaction buffer for $15\text{ min}$ at $25^\circ\text{C}$ to allow complete binary complex formation.
  - Initiate elongation by adding dNTPs and simultaneously adding a 100-fold molar excess ($2\text{ }\mu\text{M}$) of an unextendable "trap" oligonucleotide ($5'\text{-(TTAGGG)}_3\text{-ddC-3'}$).
  - Quench reaction aliquots at time points: $t = 0, 15, 30, 60, 120, 300, 600, 1200\text{ s}$.
  - Under these trap conditions, any primer that dissociates from telomerase is permanently sequestered by the trap, ensuring that only a single processive run is visualized per enzyme molecule.
  - Fit the fraction of elongated primer remaining over time to a single-exponential decay:
    $$F(t) = F_0 \cdot e^{-k_{\text{off}}^{\text{primer}} t}$$

---

## 3. Definitive Mathematical Falsification Matrix

The table below establishes the explicit, mutually exclusive decision rules for interpreting the combined experimental outcomes.

| Case | SPR TPP1 $K_d$ (Assay 2) | Basal Activity $k_{\text{cat}}$ (Assay 3) | TPP1 RAP Stimulation $\chi_{\text{RAP}}$ (Assay 3) | Primer Dissociation $k_{\text{off}}^{\text{primer}}$ w/ TPP1 (Assay 4) | Scientific Interpretation & Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Case A** | **Normal** ($K_d \le 300\text{ nM}$) | **Competent** ($\ge 40\%\text{ WT}$) | **Severely Impaired** ($\chi_{\text{RAP}} \le 1.5\text{ vs WT } \ge 5.0$) | **Rapid** (Identical to -TPP1) | **HYPOTHESIS CONFIRMED**: True separation of recruitment from processivity. Arg756 is required for TPP1-mediated allosteric processivity stimulation. |
| **Case B** | **Normal** ($K_d \le 300\text{ nM}$) | **Impaired** ($< 20\%\text{ WT}$) | **Unmeasurable** | **Very Rapid** | **ALTERNATIVE 1 (Catalytic Defect)**: Mutation disrupts fundamental template-translocation chemistry (A55 contact), preventing RAP measurement. Allostery not proven. |
| **Case C** | **Weakened** ($K_d > 1,000\text{ nM}$) | Normal | Severely Impaired | Rapid | **ALTERNATIVE 2 (Coupled Defect)**: TRAP base mutation propagates conformational change up to the apex, weakening TPP1 binding. Functions are NOT separable. |
| **Case D** | **Normal** ($K_d \le 300\text{ nM}$) | Normal | **Normal** ($\chi_{\text{RAP}} \ge 4.5$) | **Slowed** (Identical to WT +TPP1) | **HYPOTHESIS FALSIFIED**: Arg756 is uninvolved in TPP1-mediated processivity stimulation. The TRAP hairpin is not an allosteric conduit. |

---

## 4. Quantitative Thresholds for Scientific Discovery

To claim a **publication-quality breakthrough**, the experiment must satisfy all four of the following quantitative conditions:
1. **Binding Equivalence**: $K_d^{\text{Arg756A}} / K_d^{\text{WT}} \le 2.0$ (via SPR).
2. **Thermal Stability Equivalence**: $|T_m^{\text{Arg756A}} - T_m^{\text{WT}}| \le 1.5^\circ\text{C}$ (via nanoDSF).
3. **Selective Processivity Collapse**: $\chi_{\text{RAP}}^{\text{Arg756A}} \le 1.8$, while $\chi_{\text{RAP}}^{\text{WT}} \ge 5.0$ at $[TPP1] = 500\text{ nM}$.
4. **Recruitment-Defective Control Validation**: The `Arg774A` control must show $K_d > 5,000\text{ nM}$ and $\chi_{\text{RAP}} \le 1.2$, demonstrating that recruitment failure and allosteric failure produce distinguishable biophysical signatures.
