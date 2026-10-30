# STATISTICAL ANALYSIS & INDEPENDENT REPLICATION PLAN

**Document ID:** `EXP-DOC-004-REV-FINAL`  
**Effective Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/STATISTICAL_ANALYSIS_AND_REPLICATION_PLAN.md`  
**Classification:** Pre-Registration Protocol, Power Analysis, & Statistical Framework  

---

## 1. Study Design & Epistemic Pre-Registration

This document establishes the pre-registered statistical analysis and replication framework for evaluating the TERT `Lys757A` processivity uncoupling hypothesis. In accordance with rigorous scientific standards (e.g., Nature Research reporting guidelines), all sample sizes, exclusion criteria, statistical tests, and significance thresholds are pre-committed prior to wet-lab data collection.

---

## 2. Statistical Power Analysis and Sample Size Justification

### 2.1 Hypothesized Effect Sizes
1. **Recruitment Affinity ($K_d$)**:
   - WT Baseline: $150 \pm 25\text{ nM}$.
   - Pre-registered null hypothesis ($H_0$): `Lys757A` binding is severely impaired ($K_d > 1,500\text{ nM}$).
   - Pre-registered alternative hypothesis ($H_1$): `Lys757A` retains nanomolar binding ($K_d \le 300\text{ nM}$).
2. **Processivity Stimulation ($\chi_{\text{RAP}}$)**:
   - WT Baseline: $\chi_{\text{RAP}} = 6.2 \pm 0.8\text{ fold}$.
   - `Lys757A` Predicted Outcome: $\chi_{\text{RAP}} \le 1.8 \pm 0.3\text{ fold}$.
   - Expected difference: $\Delta \chi_{\text{RAP}} = 4.4\text{ fold}$.
   - Pooled standard deviation: $s_{\text{pooled}} \approx 0.60$.
   - Standardized effect size (Cohen's $d$):
     $$d = \frac{|\mu_1 - \mu_2|}{s_{\text{pooled}}} = \frac{6.2 - 1.8}{0.60} = 7.33$$

### 2.2 Power Calculation ($G*\text{Power}$ Equivalent)
- Significance level: $\alpha = 0.05$ (two-tailed).
- Statistical power ($1 - \beta$): $\ge 0.95$ (95% power).
- For a conservative effect size $d = 3.0$ (substantially smaller than the predicted $d = 7.33$), the minimum required sample size to achieve $95\%$ power at $\alpha = 0.05$ is **$N = 3$ independent biological replicates**.
- Therefore, $N = 3$ independent biological replicates provides $>99.9\%$ statistical power to detect the hypothesized processivity collapse.

---

## 3. Biological vs. Technical Replication Hierarchy

```
========================================================================================================================
REPLICATION DEFINITIONS AND SAMPLING HIERARCHY
========================================================================================================================
Hierarchy Level        Definition                                                    Requirement
------------------------------------------------------------------------------------------------------------------------
Biological Replicate   Completely independent mammalian cell transfection,          Three distinct preparations
(N = 3)                distinct harvest day, and separate Strep-Tactin/SEC           (Batch A, Batch B, Batch C).
                       purification batch.

Technical Replicate    Independent assay reaction tubes, gel lanes, or sensor chip   Duplicate or triplicate measurements
(n = 2-3)              injections derived from the same purified protein batch.      averaged to produce 1 biological value.
========================================================================================================================
```

> **Mandatory Statistical Rule**: All formal hypothesis tests and error bars (mean $\pm$ standard error of the mean, SEM) must be computed across **$N = 3$ biological replicates**. Technical replicates ($n$) are averaged within each biological batch and cannot be treated as independent degrees of freedom (preventing pseudoreplication).

---

## 4. Pre-Registered Hypothesis Testing Pipeline

### 4.1 Assay 1: TPP1 Recruitment Binding Affinity (SPR)
- **Primary Endpoint**: Equilibrium dissociation constant ($K_d$ in nM).
- **Statistical Test**: One-way Analysis of Variance (ANOVA) comparing WT, `Lys757A`, and `Lys757R` across $N = 3$ biological replicates.
- **Post-Hoc Test**: Dunnett's multiple comparisons test against WT.
- **Equivalence Threshold**: Two one-sided tests (TOST) for equivalence testing to confirm `Lys757A` is statistically equivalent to WT within a bound of $\pm 150\text{ nM}$ ($\alpha = 0.05$).

### 4.2 Assay 2: Basal Telomerase Activity ($k_{\text{cat}}$)
- **Primary Endpoint**: Normalized total nucleotide incorporation in the absence of TPP1 ($0\text{ nM}$).
- **Statistical Test**: One-sample t-test testing whether `Lys757A` basal activity is significantly greater than the pre-set failure threshold of $30\%\text{ WT}$ ($H_0: \mu \le 30\%$ vs $H_1: \mu > 30\%$).
- **Equivalence Test**: Confirm basal activity of `Lys757A` does not differ significantly from WT ($p > 0.05$).

### 4.3 Assay 3: Repeat Addition Processivity Stimulation ($\chi_{\text{RAP}}$)
- **Primary Endpoint**: Fold-stimulation of processivity $\chi_{\text{RAP}} = \frac{P_{\text{overall}}(+500\text{ nM})}{P_{\text{overall}}(0\text{ nM})}$.
- **Statistical Model**: Two-way repeated-measures ANOVA:
  - Factor 1: Construct (WT, `Lys757A`, `Lys757R`).
  - Factor 2: POT1–TPP1 concentration ($0, 50, 100, 250, 500\text{ nM}$).
- **Key Interaction**: The interaction term $\text{Construct} \times \text{Concentration}$ must be statistically significant ($p < 0.001$).
- **Post-Hoc Pairwise Comparisons**:
  - `Lys757A` vs. WT at $500\text{ nM}$ TPP1: One-tailed t-test testing $\chi_{\text{RAP}}(\text{Lys757A}) < \chi_{\text{RAP}}(\text{WT})$ ($p < 0.001$).
  - `Lys757R` vs. `Lys757A` at $500\text{ nM}$ TPP1: One-tailed t-test testing rescue $\chi_{\text{RAP}}(\text{Lys757R}) > \chi_{\text{RAP}}(\text{Lys757A})$ ($p < 0.001$).

---

## 5. Quality Control, Outlier Policy, and Data Exclusion

1. **Precision & Reproducibility Metrics**:
   - Intra-assay Coefficient of Variation ($\text{CV}_{\text{intra}}$) between technical duplicates must be $\le 12.0\%$.
   - Inter-assay Coefficient of Variation ($\text{CV}_{\text{inter}}$) between biological triplicates must be $\le 18.0\%$.
2. **Outlier Identification**:
   - Outlier testing must be conducted using **Grubbs' test** ($\alpha = 0.01$, two-sided) or the **ROUT method** ($Q = 1\%$).
   - Data points may **never** be excluded based on subjective inspection or disagreement with the hypothesis. An outlier may only be excluded if:
     * A documented technical anomaly occurred (e.g., bubble in SPR flow cell, gel tear/bubble across lane, pipetting failure confirmed by internal recovery control).
     * Grubbs' test identifies the point at $p < 0.01$.
   - All excluded points and reasons must be reported in the manuscript supplementary data table.

---

## 6. Blinded Operator and Randomization Protocol

To prevent subconscious bias during gel densitometry and SPR sensorgram analysis:
1. **Sample Coding**: A non-participating laboratory member will assign randomized codes (e.g., Construct Alpha, Beta, Gamma) to purified protein aliquots prior to assay execution.
2. **Blinded Execution**: The experimentalist performing SPR injections, radioactive primer extensions, gel loading, and phosphorimager band integration will remain blinded to construct identity.
3. **Data Locking**: All band densitometry values and kinetic parameters must be committed and cryptographically timestamped in the data repository before unblinding.
