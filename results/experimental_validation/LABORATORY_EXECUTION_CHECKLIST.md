# FINAL LABORATORY EXECUTION CHECKLIST: TERT LYS757 VALIDATION

**Document ID:** `EXP-CHK-001-REV-FINAL`  
**Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/LABORATORY_EXECUTION_CHECKLIST.md`  
**Classification:** Operational Standard Operating Checklist  

---

## Phase 1: Recombinant Construct Engineering & Sequence Verification

- [ ] **1.1 Primer Procurement**: Synthesize high-purity HPLC-purified mutagenic primer pairs for `Lys757A` (K757A-F/R) and `Lys757R` (K757R-F/R) as specified in [`CONSTRUCT_DESIGN_WT_K757A_K757R.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md).
- [ ] **1.2 Mutagenesis PCR**: Perform Q5 High-Fidelity PCR on `pcDNA3.1(+)-TwinStrep-FLAG-3C-hTERT` template with 5x Q5 GC Enhancer (18 cycles).
- [ ] **1.3 DpnI Elimination**: Incubate PCR reaction with $1.0\text{ }\mu\text{L}$ DpnI at $37\text{ }^\circ\text{C}$ for 2 hours to eliminate parental methylated plasmid.
- [ ] **1.4 Bacterial Transformation**: Transform into NEB 5-alpha competent *E. coli*; plate on LB + Ampicillin ($100\text{ }\mu\text{g/mL}$).
- [ ] **1.5 Colony Restriction Screen**: Screen transformants using diagnostic digestion:
  - `Lys757A`: Cleave PCR amplicon with **BstUI** (`CG^CG`) or **HinP1I** (`G^CGC`).
  - `Lys757R`: Cleave PCR amplicon with **MspI / HpaII** (`C^CGG`).
- [ ] **1.6 Full-Length Sequence Verification**: Perform Sanger / Nanopore sequencing across the entire 3,399 bp open reading frame using the 8-primer panel to guarantee zero secondary mutations.

---

## Phase 2: Mammalian Cell Co-Expression (Expi293F)

- [ ] **2.1 Culture Conditioning**: Grow Expi293F suspension cells to $2.8\text{--}3.2 \times 10^6\text{ cells/mL}$ with $>95\%$ viability.
- [ ] **2.2 Lipid Transfection**: Co-transfect $0.25\text{ mg}$ hTERT plasmid + $0.75\text{ mg}$ hTR plasmid per liter using ExpiFectamine 293 in Opti-MEM I.
- [ ] **2.3 Enhancer Supplementation**: Add Transfection Enhancers 1 and 2 at 18–20 hours post-transfection.
- [ ] **2.4 Harvest**: Pellet cells at $48\text{ hours}$ post-transfection ($1,000 \times g$, 20 min, 4 °C); wash with ice-cold PBS.

---

## Phase 3: Purification and Holoenzyme Assembly

- [ ] **3.1 Lysis**: Lyse cells in Lysis Buffer with $300\text{ mM KCl}$, $0.1\%\text{ Triton X-100}$, protease inhibitors, and RNasin Plus; clarify by ultracentrifugation at $40,000 \times g$ for 45 min.
- [ ] **3.2 Affinity Capture**: Batch-bind clarified lysate onto Strep-Tactin XT Superflow resin ($4\text{ mL}$ resin/L) for 3 hours at 4 °C.
- [ ] **3.3 High-Salt Wash**: Wash resin with 30 column volumes of $500\text{ mM KCl}$ buffer to strip non-specific chaperone complexes.
- [ ] **3.4 Tag Removal**: Digest on-column with GST-PreScission 3C Protease ($100\text{ }\mu\text{g}$) for 16 hours at 4 °C.
- [ ] **3.5 Size-Exclusion Polish**: Polish eluate over a Superose 6 Increase 10/300 GL column; collect the monomeric holoenzyme peak ($V_e \approx 13.5\text{--}14.2\text{ mL}$).

---

## Phase 4: Biophysical Quality Control Gating

- [ ] **4.1 UV Stoichiometry Check**: Confirm $A_{260}/A_{280}$ ratio is within $1.15\text{--}1.35$ (indicates 1:1 TERT:hTR assembly).
- [ ] **4.2 RNA Integrity Check**: Denaturing 6% urea-PAGE with SYBR Gold staining; confirm clean 451-nt hTR band with zero degradation fragments.
- [ ] **4.3 Protein Purity Check**: SDS-PAGE (4–12% Bis-Tris) with Silver Stain; confirm $>90\%$ purity at ~125 kDa.
- [ ] **4.4 Folding QC (nanoDSF)**: Measure apparent melting temperature ($T_m$) on Prometheus NT.48; confirm $T_m \ge 48.5\text{ }^\circ\text{C}$ (within $1.5\text{ }^\circ\text{C}$ of WT).

---

## Phase 5: Surface Plasmon Resonance (SPR) Binding Assay

- [ ] **5.1 Sensor Preparation**: Coat Biacore Sensor Chip SA with biotinylated TPP1 OB domain (50–80 RU).
- [ ] **5.2 Holoenzyme Injections**: Titrate WT, `Lys757A`, and `Lys757R` at 0, 15, 30, 60, 120, 240, 480 nM ($30\text{ }\mu\text{L/min}$, 180 s association, 600 s dissociation).
- [ ] **5.3 Automated Kinetic Fitting**: Fit double-subtracted sensorgrams with `scripts/spr_kinetic_analysis.py` (1:1 Langmuir model with MTL).
- [ ] **5.4 Gating Check**: Confirm `Lys757A` exhibits $K_d \le 300\text{ nM}$. (If $K_d > 1,500\text{ nM}$, stop; hypothesis refuted).

---

## Phase 6: Quantitative Direct Primer Extension RAP Gel Assay

- [ ] **6.1 Reaction Assembly**: Incubate $50\text{ nM}$ holoenzyme, $1.0\text{ }\mu\text{M}$ $(\text{TTAGGG})_3$ primer, nucleotides, $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$, and POT1–TPP1 at 0, 50, 100, 250, 500 nM at 30 °C for 60 min.
- [ ] **6.2 Quenching & Recovery Control**: Add Stop Buffer containing $^{32}\text{P}$-labeled 18-nt internal recovery standard; digest with Proteinase K.
- [ ] **6.3 Electrophoresis**: Run 12% denaturing polyacrylamide gel ($7\text{ M urea}$, 65 W constant power, 50 °C surface temp).
- [ ] **6.4 Phosphorimaging**: Expose to storage phosphor screen for 18–24 h; scan on Typhoon FLA 9500 at $50\text{ }\mu\text{m}$ resolution.
- [ ] **6.5 Automated Quantification**: Run `scripts/rap_gel_quantification.py` to extract $R_{p,i}$, $P_{\text{overall}}$, basal $k_{\text{cat}}$, and $\chi_{\text{RAP}}$.

---

## Phase 7: Multi-Replicate Verification & Statistical Testing

- [ ] **7.1 Biological Triplicate / Quintuplicate**: Repeat Phases 2–6 across $N = 4\text{ to }5$ independent biological expressions.
- [ ] **7.2 Blinded Analysis**: Perform densitometry and fitting with randomized sample IDs.
- [ ] **7.3 Statistical Testing**: Run two-way repeated-measures ANOVA with Dunnett's post-hoc test.
- [ ] **7.4 Decision Lock**: Evaluate outcome against pre-registered breakthrough criteria.
