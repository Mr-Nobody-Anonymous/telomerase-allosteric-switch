# BREAKTHROUGH EVIDENCE CRITERIA: EMPIRICAL REQUIREMENTS FOR MECHANISTIC DISCOVERY

**Document ID:** `EXP-EVI-001-REV-FINAL`  
**Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/BREAKTHROUGH_EVIDENCE_CRITERIA.md`  
**Classification:** Definitive Experimental Evidentiary Standards  

---

## 1. Principle & Epistemic Demarcation

A molecular mechanism cannot be claimed as a demonstrated scientific breakthrough based solely on computational docking, molecular dynamics trajectories, crystallographic B-factors, or mathematical predictions. 

### 1.1 Critical Computational Boundaries Established by Reproducibility Audit
1. **Simulation Duration Limit (3.0 ps)**: The all-atom molecular dynamics simulations executed in OpenMM (1,500 steps with $\Delta t = 2.0\text{ fs}$) sample only **$3.0\text{ ps}$** of physical time ($2.4\text{ ps}$ production post-burnin). This captures local coordinate relaxation but is orders of magnitude below biological turnover ($\approx 0.5\text{--}2\text{ s}$ per repeat addition).
2. **Cleft Compaction vs Latch Demarcation**: Independent replica analysis ($N = 3, df = 4$) establishes that TPP1 binding induces statistically significant cleft compaction ($3.82\text{ \AA}$ bound vs $5.62\text{ \AA}$ apo; $t = -4.71, p = 0.0093$). However, the hydrogen-bond occupancy difference is **not statistically significant** ($32.2\% \pm 28.7\%$ vs $0.0\% \pm 0.0\%$; $t = 1.95, p = 0.124$). Consequently, **a persistent Tyr176–Gln794 hydrogen-bond latch has not been demonstrated in silico**.
3. **Decisive Primacy of Wet-Lab Evidence**: Because computational models can only establish physical plausibility and generate testable hypotheses, only reproducible in vitro biochemical assays can prove or disprove whether `Lys757` acts as a biological processivity switch.

To claim that **TERT Lys757 defines a molecular switch uncoupling repeat addition processivity from telomere recruitment**, the research group must possess and deposit the following **three mandatory empirical primary datasets** and **one mechanistic validation dataset**.

---

## 2. The Three Mandatory Biochemical Observations

```
+---------------------------------------------------------------------------------------------------------+
|                                  MANDATORY PRIMARY DATASETS REQUIRED                                    |
+---------------------------------------------------------------------------------------------------------+
| 1. NORMAL RECRUITMENT AFFINITY (SPR / BLI):                                                             |
|    - Sensorgrams and fitted Langmuir equilibrium dissociation constants across N >= 4 biological        |
|      batches demonstrating:                                                                             |
|         Kd(TERT Lys757A) <= 300 nM  (WT reference: 150 ± 25 nM; p > 0.05 equivalence test).             |
|    - Raw data requirement: Multi-concentration sensorgrams (0–480 nM) showing identical association    |
|      and dissociation curves between WT and Lys757A.                                                    |
+---------------------------------------------------------------------------------------------------------+
| 2. PRESERVED BASAL CATALYTIC ELONGATION:                                                                |
|    - Direct primer extension gel bands in the absence of TPP1 (0 nM) showing total radiolabeled         |
|      nucleotide incorporation across N >= 4 biological batches:                                         |
|         k_cat(TERT Lys757A) >= 85% of WT  (p > 0.05 vs WT; p < 0.001 vs catalytic dead failure).         |
|    - Proves that the mutation does not distort catalytic core folding, template RNA base-pairing,       |
|      or nucleotide addition chemistry.                                                                  |
+---------------------------------------------------------------------------------------------------------+
| 3. SELECTIVELY IMPAIRED PROCESSIVITY STIMULATION WITH CONSERVATIVE RESCUE:                              |
|    - Phosphorimager lane profiles across POT1–TPP1 titration (0, 50, 100, 250, 500 nM) showing:        |
|         chi_RAP(Lys757A) <= 1.8 fold  (versus WT = 6.2 ± 0.8 fold; p < 0.001).                         |
|    - Stereochemical rescue control in the same gel:                                                     |
|         chi_RAP(Lys757R) >= 4.5 fold  (p < 0.001 vs Lys757A).                                           |
|    - Proves that the processivity defect is specific to exit-pore charge and sidechain geometry,        |
|      not a non-specific folding artifact.                                                               |
+---------------------------------------------------------------------------------------------------------+
```

---

## 3. Mechanistic Confirmation Dataset: Translocation Kinetics

To establish that the loss of processivity is directly caused by failure of primer retention during translocation (rather than altered dNTP binding or template alignment), one of the following dynamic datasets must be produced:

1. **Pulse-Chase Primer Retention Kinetics**:
   - Chase with 200-fold excess non-extendable competitor oligonucleotide over a 90-minute time course.
   - Measured complex half-life in the presence of TPP1:
     * In WT + TPP1: $t_{1/2} \ge 20.0\text{ min}$ (dissociation rate $k_{\text{off}}^{\text{translocation}} \approx 0.012\text{ s}^{-1}$).
     * In `Lys757A` + TPP1: $t_{1/2} \le 3.5\text{ min}$ ($k_{\text{off}}^{\text{translocation}} \ge 0.080\text{ s}^{-1}$).
2. **Single-Molecule FRET (smFRET) Translocation Dwell Times**:
   - Real-time TIRF imaging of Cy3-telomerase and Cy5-primer in microfluidic flow cells.
   - Demonstration that `Lys757A` exhibits premature low-FRET dissociation transitions during the translocation pause between successive repeats.

---

## 4. Minimum Reporting Package for Peer Review

Prior to submitting a manuscript claiming this discovery to a top-tier journal (*Nature*, *Science*, *Cell*), the experimental submission package must include:
1. **Raw Phosphorimager Gel Scans**: Full uncropped 16-bit `.gel` / `.tif` files for all replicates.
2. **SPR Surface Performance Files**: Biacore `.blr` / `.txt` raw sensorgrams including reference channel subtractions and blank cycles.
3. **Purity & Integrity Gels**: Uncropped SDS-PAGE (silver stain) and denaturing urea-PAGE showing stoichiometric 1:1 hTERT:hTR assembly.
4. **Independent Blinding Record**: Cryptographically signed blinding key unsealed only after data quantification.
5. **Pre-Registration Match**: Proof that the observed data was evaluated strictly against the pre-registered decision thresholds without post-hoc threshold adjustment.
