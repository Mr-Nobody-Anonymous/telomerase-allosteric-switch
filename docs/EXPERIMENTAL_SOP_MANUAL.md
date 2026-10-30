# Unified Experimental Standard Operating Procedures (SOP) Manual

**Document ID:** `SOP-MANUAL-001-REV-FINAL`  
**Classification:** Complete Wet-Lab Operational Protocol Suite  
**Target Repository:** `docs/EXPERIMENTAL_SOP_MANUAL.md`  

---

## 1. Recombinant Human Telomerase Holoenzyme Reconstitution

### 1.1 Co-Expression in Expi293F Cells
- **Host**: Human Expi293F suspension cells grown in Expi293 Expression Medium ($37^\circ\text{C}$, $8\%\text{ CO}_2$, $125\text{ rpm}$).
- **Transfection Plasmids (Molar Ratio 1:3:1:1:1:1)**:
  1. `pcDNA3.1-TwinStrep-PreScission-hTERT` (WT, `Lys757A`, or `Lys757R` — 3,399 bp CDS).
  2. `pBS-U1-hTR` (Full-length 451 nt human telomerase RNA driven by human U1 snRNA promoter).
  3. `pcDNA3.1-Flag-Dyskerin`.
  4. `pcDNA3.1-Myc-NOP10`.
  5. `pcDNA3.1-HA-NHP2`.
  6. `pcDNA3.1-V5-GAR1`.
- **Transfection**: ExpiFectamine 293 reagent according to manufacturer protocols at density $3.0 \times 10^6\text{ cells/mL}$. Harvest cells 48 hours post-transfection.

### 1.2 Affinity Purification & Stringency Wash
- **Lysis Buffer**: $50\text{ mM HEPES-KOH (pH 7.5)}$, $300\text{ mM KCl}$, $10\%\text{ glycerol}$, $1\text{ mM EDTA}$, $0.1\%\text{ Triton X-100}$, $1\text{ mM DTT}$, protease inhibitor cocktail (Roche cOmplete EDTA-free), $40\text{ U/mL}$ RNasin.
- **Clarification**: Centrifuge lysate at $20,000 \times g$ for 45 min at $4^\circ\text{C}$.
- **Affinity Capture**: Incubate supernatant with Strep-Tactin XT 4Flow resin (IBA Lifesciences) for 3 hours at $4^\circ\text{C}$.
- **High-Stringency Salt Wash (Crucial Step)**:
  - Wash resin with 30 column volumes of Wash Buffer containing **$500\text{ mM KCl}$** and $0.1\%\text{ Triton X-100}$ to strip endogenous shelterin, chaperones, and non-specific nucleic acids.
  - Wash with 10 column volumes of Isotonic Cleavage Buffer ($20\text{ mM HEPES-KOH pH 7.5}$, $150\text{ mM KCl}$, $1\text{ mM DTT}$).
- **Tag Removal**: Add PreScission 3C protease ($20\text{ U/mL}$) directly to resin; incubate overnight at $4^\circ\text{C}$. Collect cleaved native holoenzyme in flow-through.
- **Size-Exclusion Polishing**: Inject concentrated flow-through onto a Superose 6 Increase 10/300 GL column equilibrated in Assay Buffer ($20\text{ mM HEPES-KOH pH 7.5}$, $150\text{ mM KCl}$, $2\text{ mM MgCl}_2$, $1\text{ mM TCEP}$, $5\%\text{ glycerol}$). Collect the stoichiometric monomeric holoenzyme peak ($\approx 12.2\text{ mL}$).

---

## 2. Recombinant TPP1 OB-Fold Domain Preparation

- **Construct**: `pGEX-6P-1-TPP1(89–248)` encompassing the structured OB-fold and TEL-patch.
- **Expression**: *E. coli* BL21(DE3) induced with $0.4\text{ mM IPTG}$ at $18^\circ\text{C}$ for 16 hours.
- **Purification**: Glutathione Sepharose 4B capture, PreScission 3C protease cleavage to remove GST, followed by Superdex 75 Increase 10/300 GL SEC. Peak fractions pooled at purity $>98\%$ by SDS-PAGE.

---

## 3. Surface Plasmon Resonance (SPR) Kinetics (Biacore 8K)

- **Sensor Chip**: Series S Sensor Chip CM5 (Cytiva).
- **Immobilization**: Amine coupling via standard EDC/NHS chemistry. Anti-TERT monoclonal antibody or Strep-Tactin XT immobilized on Channel 2 (active); Channel 1 left unmodified (reference).
- **Running Buffer**: $20\text{ mM HEPES-KOH (pH 7.5)}$, $150\text{ mM KCl}$, $3\text{ mM EDTA}$, $0.05\%\text{ Surfactant P20}$, $1\text{ mM TCEP}$.
- **Multi-Cycle Kinetic Series**: Inject purified TPP1 OB-fold analyte at two-fold dilutions: $0, 30, 60, 120, 240, 480\text{ nM}$ at $30\text{ }\mu\text{L/min}$.
  - Association: $180\text{ s}$
  - Dissociation: $600\text{ s}$
  - Regeneration: $10\text{ mM Glycine-HCl (pH 2.0)}$ for $30\text{ s}$.
- **Data Analysis**: Subtract reference channel and blank cycles. Fit 1:1 Langmuir kinetic model using [`scripts/spr_kinetic_analysis.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/spr_kinetic_analysis.py).

---

## 4. Quantitative 12% Urea-PAGE RAP Assay

- **Primer Extension Assembly ($20\text{ }\mu\text{L}$ Total Volume)**:
  - $50\text{ mM Tris-HCl (pH 8.0)}$
  - $50\text{ mM KCl}$, $1\text{ mM MgCl}_2$, $5\text{ mM 2-mercaptoethanol}$, $1\text{ mM spermidine}$
  - $1.0\text{ }\mu\text{M}$ Telomeric DNA Primer: $5'\text{-(TTAGGG)}_3\text{-}3'$ (18 nt)
  - $1.0\text{ mM dATP}$, $1.0\text{ mM dTTP}$
  - $5.0\text{ }\mu\text{M}$ unlabeled dGTP supplemented with $10\text{ }\mu\text{Ci } [\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ ($3000\text{ Ci/mmol}$)
  - Purified telomerase holoenzyme ($10\text{ nM}$)
  - Purified TPP1 OB-fold ($0\text{ nM}$ for basal; $500\text{ nM}$ for stimulated).
- **Reaction**: Incubate at $30^\circ\text{C}$ for 60 min. Quench with $20\text{ }\mu\text{L}$ stop solution ($95\%\text{ formamide}$, $10\text{ mM EDTA}$, $0.05\%\text{ bromophenol blue}$, $0.05\%\text{ xylene cyanol}$) and 500 cpm of $5'\text{-}^{32}\text{P}$-labeled 12-nt oligonucleotide recovery control.
- **Denaturing Electrophoresis**: Resolve on a 40-cm 12% polyacrylamide / 7 M urea sequencing gel at $55\text{ W}$ constant power for 2 hours until bromophenol blue reaches the bottom.
- **Detection & Quantification**: Expose dried gel to a Storage Phosphor Screen for 16 hours. Scan on a Typhoon FLA 9500 scanner at $50\text{ }\mu\text{m}$ resolution. Quantify using [`scripts/rap_gel_quantification.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/rap_gel_quantification.py).

---

## 5. Double-Blind Protocol & Automated Pipeline Execution

1. **Blinding**: An independent researcher re-labels construct tubes as `Sample Alpha`, `Sample Beta`, `Sample Gamma`, and deposits the encryption key in a sealed envelope.
2. **Batch Processing**: The experimental operator runs SPR and RAP assays using blind identifiers.
3. **Execution**:
   ```bash
   python scripts/spr_kinetic_analysis.py --input raw_spr_data.csv --output results/spr_results.csv
   python scripts/rap_gel_quantification.py --input raw_rap_data.csv --output results/rap_results.csv
   ```
4. **Unblinding**: Only after script execution yields the cryptographic verification report are samples unblinded and matched against the pre-registered decision tree in [`results/experimental_validation/BREAKTHROUGH_EVIDENCE_CRITERIA.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/BREAKTHROUGH_EVIDENCE_CRITERIA.md).
