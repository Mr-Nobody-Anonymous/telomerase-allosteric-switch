# RECONSIDERED STATISTICAL POWER ANALYSIS: REALISTIC EFFECT SIZES & SAMPLE SIZE DETERMINATION

**Document ID:** `EXP-DOC-005-REV-FINAL`  
**Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/RECONSIDERED_POWER_ANALYSIS.md`  
**Classification:** Pre-Registration Statistical Power Analysis  

---

## 1. Executive Summary & Audit Mandate

The preliminary power calculation in this project assumed an idealized effect size of Cohen's $d = 7.33$, concluding that $N = 3$ biological replicates would provide $>99.9\%$ statistical power. 

In response to scientific audit, **we formally revise this calculation**. Assuming an effect size of $d > 7$ without prior pilot data is overly optimistic and carries a high risk of underpowering the study if the true biological uncoupling effect is more subtle or if experimental variance is elevated. 

This document computes exact statistical power across realistic biological effect sizes ($d \in [1.0, 4.0]$) using empirical literature variance estimates from published human telomerase direct extension assays (Latrick & Cech, 2010; Padmanaban et al., 2023). We demonstrate that:
1. Under a moderate-to-strong biological effect ($d = 2.0$), $N = 3$ biological replicates yields only **46.3% power** (a $>53\%$ risk of a Type II false negative).
2. To achieve the canonical benchmark of **$\ge 80\%$ statistical power** at $\alpha = 0.05$ across realistic effect sizes ($d = 2.0\text{--}2.5$), the laboratory must pre-commit to **$N = 4$ to $N = 5$ independent biological replicates**.

---

## 2. Empirical Literature Variance Modeling

From published radiometric direct primer extension assays measuring human telomerase processivity stimulation by POT1–TPP1:
- **Latrick & Cech (2010)**:
  - WT complex half-life: $t_{1/2} = 22.5 \pm 3.2\text{ min}$ ($\text{CV} = 14.2\%$).
  - Processivity factor $R_p$: $0.89 \pm 0.03$ ($\text{CV} = 3.4\%$).
  - Overall fold-stimulation: $\chi_{\text{RAP}} \approx 5.5\text{ to }6.8$, with standard deviation $\sigma_{\text{WT}} \approx 1.1$ ($\text{CV} \approx 18\%$).
- **Padmanaban et al. (2023)**:
  - Processivity reduction in TEN-domain mutants exhibited biological inter-preparation standard deviations of $\sigma \approx 0.6\text{ to }0.9$.
- **Realistic Scenario Modeling**:
  - Baseline WT: $\mu_1 = 6.0$, $\sigma_1 = 1.1$.
  - Candidate `Lys757A`:
    * *Ideal Uncoupling*: $\mu_2 = 1.8$, $\sigma_2 = 0.5 \implies \Delta = 4.2, s_{\text{pooled}} = 0.86, d = 4.88$.
    * *Moderate Uncoupling*: $\mu_2 = 2.8$, $\sigma_2 = 0.8 \implies \Delta = 3.2, s_{\text{pooled}} = 0.96, d = 3.33$.
    * *Conservative / Attenuated Uncoupling*: $\mu_2 = 3.6$, $\sigma_2 = 1.0 \implies \Delta = 2.4, s_{\text{pooled}} = 1.05, d = 2.29$.
    * *Subtle Uncoupling*: $\mu_2 = 4.2$, $\sigma_2 = 1.1 \implies \Delta = 1.8, s_{\text{pooled}} = 1.10, d = 1.64$.

---

## 3. Exact Statistical Power Calculation Matrix

Power was calculated for a two-sample two-tailed Student's t-test with $\alpha = 0.05$ across sample sizes $N \in [3, 8]$ using the non-central Student's t-distribution ($\text{df} = 2N - 2$, non-centrality parameter $\delta = d \sqrt{N/2}$):

```
========================================================================================================================
EXACT STATISTICAL POWER MATRIX (Two-Tailed, Alpha = 0.05)
========================================================================================================================
Standardized Effect Size (d)   N = 3         N = 4         N = 5         N = 6         N = 8
------------------------------------------------------------------------------------------------------------------------
d = 1.0 (Subtle biological effect)   15.9%         22.3%         28.6%         34.7%         46.1%
d = 1.5 (Moderate effect)            29.3%         43.0%         54.9%         65.0%         79.7%
d = 2.0 (Substantial effect)         46.3%         65.7%         79.1%         87.6%         96.0%
d = 2.5 (Strong biological effect)   63.6%         83.6%         93.2%         97.3%         99.6%
d = 3.0 (Very strong effect)         78.3%         93.9%         98.5%         99.6%        100.0%
d = 4.0 (Robust knockout-like effect)94.8%         99.6%        100.0%        100.0%        100.0%
========================================================================================================================
```

---

## 4. Key Takeaways & Pre-Committed Sample Size Decision

1. **The Risk of $N = 3$**: If the true biological effect of `Lys757A` is substantial ($d = 2.0$) rather than massive ($d > 4.0$), an $N = 3$ experiment has a **53.7% failure rate** to reach statistical significance ($p < 0.05$).
2. **The Optimal Sample Size**:
   - At **$N = 5$ biological replicates**, statistical power reaches **79.1%** for $d = 2.0$, **93.2%** for $d = 2.5$, and **98.5%** for $d = 3.0$.
3. **Mandatory Protocol Revision**:
   - The laboratory execution plan is hereby updated to require **$N = 4$ independent biological preparations** as the standard cohort, with a pre-registered sequential design permitting expansion to **$N = 6$** if the initial $p$-value falls into the ambiguous window $0.01 \le p \le 0.10$.
