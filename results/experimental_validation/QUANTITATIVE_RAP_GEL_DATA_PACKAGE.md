# QUANTITATIVE REPEAT ADDITION PROCESSIVITY (RAP) GEL DATA PACKAGE

**Document ID:** `EXP-DOC-003-REV-FINAL`  
**Effective Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md`  
**Classification:** Experimental Assay SOP, Image Quantification Pipeline, & Decision Matrix  

---

## 1. Objective & Hypothesis Demarcation

This protocol provides the definitive standard for quantifying basal telomerase activity ($k_{\text{cat}}$) and TPP1-stimulated repeat addition processivity (RAP) across wild-type (WT), `Lys757A`, and `Lys757R` human telomerase holoenzymes using direct radiometric primer extension assays.

> [!IMPORTANT]
> **Pre-Registered Epistemic Hypothesis (Not Yet Measured)**:  
> - `TERT Lys757A` is predicted to maintain wild-type basal activity in the absence of TPP1 ($k_{\text{cat}} \ge 85\%\text{ WT}$), confirming that active site chemistry is intact.  
> - Upon titration of POT1–TPP1 ($0\text{--}500\text{ nM}$), `TERT Lys757A` is predicted to selectively fail to stimulate processivity ($\chi_{\text{RAP}} \le 1.8$ versus $\text{WT} \ge 5.0\text{--}8.0$).  
> - The conservative charge control `TERT Lys757R` is predicted to rescue processivity stimulation ($\chi_{\text{RAP}} \ge 4.5$).

---

## 2. Radiometric Direct Primer Extension Assay Protocol

### 2.1 Reagent Formulation & Master Mix Setup
- **Assay Buffer ($1\times$)**:
  - $50\text{ mM}$ Tris-HCl, $\text{pH } 8.0$
  - $50\text{ mM}$ Potassium acetate ($\text{CH}_3\text{COOK}$)
  - $1.0\text{ mM}$ Magnesium chloride ($\text{MgCl}_2$)
  - $1.0\text{ mM}$ Spermidine
  - $5.0\text{ mM}$ $\beta$-mercaptoethanol
  - $1.0\%\text{ (v/v) glycerol}$
- **Substrate & Nucleotides**:
  - DNA Primer: $1.0\text{ }\mu\text{M}$ telomeric primer 5'-$(\text{TTAGGG})_3$-3' (HPLC purified, 18 nt)
  - Deoxynucleotides: $500\text{ }\mu\text{M}$ dATP, $500\text{ }\mu\text{M}$ dTTP, $5.0\text{ }\mu\text{M}$ cold dGTP
  - Radiotracer: $0.15\text{ }\mu\text{M}$ $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ ($3,000\text{ Ci/mmol}$, $10\text{ mCi/mL}$, Revvity BLU008H)
- **Co-factor Titration**:
  - Recombinant human POT1–TPP1 heterodimer (full-length POT1 + TPP1 OB domain, purified via co-expression) titrated at concentrations: **$0, 50, 100, 250, 500\text{ nM}$**.
- **Internal Recovery & Loading Control**:
  - 5'-$^{32}\text{P}$-labeled 18-nt oligomer (5'-ATCGATCGATCGATCGAT-3', 2,000 cpm/reaction).

### 2.2 Reaction Incubation and Quenching
1. Assemble reactions on ice in $20\text{ }\mu\text{L}$ total volume:
   - $14\text{ }\mu\text{L}$ Master Mix (Buffer, primer, nucleotides, radiotracer).
   - $4\text{ }\mu\text{L}$ POT1–TPP1 heterodimer at appropriate titration stock.
   - $2\text{ }\mu\text{L}$ Purified holoenzyme ($50\text{ nM}$ final concentration).
2. Incubate at $30.0\text{ }^\circ\text{C} \pm 0.2\text{ }^\circ\text{C}$ for exactly **$60\text{ minutes}$** in a calibrated thermal cycler.
3. Terminate reaction by adding $80\text{ }\mu\text{L}$ ice-cold Stop Buffer ($20\text{ mM}$ EDTA, $10\text{ mM}$ Tris-HCl pH 7.5, $0.2\%\text{ SDS}$, $0.1\text{ mg/mL}$ Proteinase K, containing internal recovery control oligonucleotide).
4. Incubate at $45\text{ }^\circ\text{C}$ for $30\text{ min}$ to digest telomerase protein and prevent band-shifting.
5. Precipitate DNA by adding $10\text{ }\mu\text{L}$ $3\text{ M Sodium Acetate (pH 5.2)}$, $2.0\text{ }\mu\text{L}$ GlycoBlue ($15\text{ mg/mL}$), and $300\text{ }\mu\text{L}$ ice-cold 100% ethanol; chill at $-80\text{ }^\circ\text{C}$ for 1 hour.
6. Pellet at $21,000 \times g$ for 20 min at 4 °C, wash with 70% ethanol, air-dry for 5 min, and resuspend in $10\text{ }\mu\text{L}$ Gel Loading Buffer ($95\%\text{ formamide}$, $10\text{ mM EDTA}$, $0.05\%\text{ bromophenol blue}$, $0.05\%\text{ xylene cyanol}$).

---

## 3. High-Resolution Denaturing Gel Electrophoresis

1. Cast a $12\%\text{ denaturing polyacrylamide gel}$ ($38\text{ cm} \times 30\text{ cm} \times 0.4\text{ mm}$, 19:1 acrylamide:bisacrylamide, $7\text{ M urea}$, $1\times\text{ TBE}$).
2. Pre-run gel at $65\text{ W}$ constant power for 45 min until surface temperature reaches $50\text{ }^\circ\text{C}$.
3. Heat samples at $95\text{ }^\circ\text{C}$ for 5 min, snap-cool on ice, and load $5\text{ }\mu\text{L}$ per lane.
4. Electrophorese at $65\text{ W}$ until the bromophenol blue dye front migrates off the bottom of the gel (~2.5 hours).
5. Transfer gel onto Whatman 3MM chromatography paper, cover with plastic wrap, and dry under vacuum at $80\text{ }^\circ\text{C}$ for 60 min on a gel dryer.
6. Expose dried gel to a Storage Phosphor Screen (Cytiva BAS-IP MS 2040) inside an exposure cassette at room temperature for $18\text{--}24\text{ hours}$.
7. Scan screen on a Typhoon FLA 9500 or Amersham Typhoon biomolecular imager at $50\text{ }\mu\text{m}$ resolution, 16-bit linear dynamic range, and PMT voltage $650\text{ V}$.

---

## 4. Mathematical Image Densitometry Quantification Pipeline

Quantification is executed programmatically via `scripts/rap_gel_quantification.py` to prevent subjective bias:

```
+---------------------------------------------------------------------------------------------------------+
|                                  MATHEMATICAL FORMULATIONS FOR RAP                                      |
+---------------------------------------------------------------------------------------------------------+
| 1. GUANINE RESIDUE NORMALIZATION:                                                                       |
|    Each synthesized repeat i contains 3 incorporated dGTP residues:                                    |
|         n_G,i = 3 * i                                                                                   |
|    Normalized intensity of band i:                                                                      |
|         I_norm,i = (Raw_Band_Intensity_i - Background_i) / (3 * i * I_internal_control)                 |
|                                                                                                         |
| 2. TOTAL BASAL ACTIVITY (k_cat):                                                                        |
|         Total_Activity = sum(I_norm,i from i=1 to N)                                                    |
|                                                                                                         |
| 3. REPEAT ADDITION PROCESSIVITY FACTOR (Rp):                                                            |
|    According to Latrick & Cech (2010), the probability of synthesizing past repeat i is:               |
|         R_p,i = sum(I_norm,j from j=i+1 to N) / sum(I_norm,j from j=i to N)                            |
|                                                                                                         |
| 4. EXPONENTIAL DECAY SLOPE (Overall Processivity P_overall):                                            |
|    Fitting ln(I_norm,i) versus repeat number i yields slope -k:                                        |
|         ln(I_norm,i) = -k * i + C  ===>  P_overall = exp(-k)                                            |
|                                                                                                         |
| 5. FOLD-STIMULATION RATIO (chi_RAP):                                                                    |
|         chi_RAP = P_overall(+500 nM TPP1) / P_overall(0 nM TPP1)                                        |
+---------------------------------------------------------------------------------------------------------+
```

---

## 5. Pre-Registered Decision Thresholds & Falsification Matrix

```
========================================================================================================================
QUANTITATIVE RAP GEL ASSAY DECISION MATRIX
========================================================================================================================
Assay Metric               WT Reference           Lys757A Predicted Outcome          Lys757R Control Outcome
------------------------------------------------------------------------------------------------------------------------
Basal Elongation (0 nM)    100% (Reference)       >= 85% of WT                       >= 90% of WT
                                                  (Falsified if < 30%)               (Confirms catalytic integrity)

Processivity (-TPP1)       Rp = 0.38 ± 0.04       Rp = 0.35 ± 0.05                   Rp = 0.38 ± 0.04
                           P_overall ~ 0.35       P_overall ~ 0.35                   P_overall ~ 0.35

Processivity (+500 nM)     Rp = 0.89 ± 0.03       Rp <= 0.50                         Rp >= 0.82
                           P_overall ~ 0.88       P_overall <= 0.45                  P_overall >= 0.80

Fold-Stimulation (χ_RAP)   6.2 ± 0.8 fold         <= 1.8 fold                        >= 4.5 fold
                                                  (BREAKTHROUGH HIT)                 (RESCUE VALIDATION)
========================================================================================================================
```

### Pre-Committed Scientific Verdicts:
1. **Hypothesis Confirmed**: `Lys757A` demonstrates basal activity $\ge 85\%$ of WT, but $\chi_{\text{RAP}} \le 1.8$, while `Lys757R` demonstrates $\chi_{\text{RAP}} \ge 4.5$. The TERT-side molecular switch is proven.
2. **Hypothesis Refuted (Catalytic Failure)**: `Lys757A` exhibits basal activity $< 30\%$ of WT. The mutation disrupts fundamental reverse transcriptase elongation.
3. **Hypothesis Refuted (No Uncoupling)**: `Lys757A` exhibits $\chi_{\text{RAP}} \ge 4.0$. The residue is not required for processivity stimulation.
