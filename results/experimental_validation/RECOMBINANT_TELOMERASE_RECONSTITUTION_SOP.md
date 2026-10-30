# STANDARD OPERATING PROCEDURE (SOP): RECOMBINANT HUMAN TELOMERASE EXPRESSION, PURIFICATION, AND RECONSTITUTION

**Document ID:** `SOP-RNP-001-REV-FINAL`  
**Effective Date:** 2026-10-30  
**Target Repository:** `results/experimental_validation/RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md`  
**Classification:** Laboratory Standard Operating Procedure  

---

## 1. Principle & Purpose

This SOP provides a validated, reproducible protocol for producing milligram-scale, monodisperse human telomerase ribonucleoprotein (RNP) complexes for wild-type (WT) and mutant variants (`Lys757A` and `Lys757R`). The protocol utilizes suspension Expi293F mammalian cells co-expressing TwinStrep-FLAG-3C-hTERT and full-length human telomerase RNA (hTR), followed by Strep-Tactin XT affinity chromatography, on-column proteolytic cleavage to remove affinity tags, and size-exclusion chromatography (SEC).

---

## 2. Reagents, Media, and Buffers

### 2.1 Media and Cell Lines
- **Host Cells**: Expi293F Suspension Cells (Thermo Fisher A14527).
- **Growth Medium**: Expi293 Expression Medium (Thermo Fisher A1435101).
- **Transfection Reagent**: ExpiFectamine 293 Transfection Kit (Thermo Fisher A14524).

### 2.2 Buffer Compositions (All filtered through 0.22 µm PES membranes)
- **Lysis Buffer**:
  - $50\text{ mM}$ HEPES-KOH, $\text{pH } 7.5$
  - $300\text{ mM}$ $\text{KCl}$
  - $10\%\text{ (v/v) glycerol}$
  - $1.0\text{ mM}$ $\text{EDTA}$
  - $1.0\text{ mM}$ $\text{TCEP}$ (freshly added)
  - $0.1\%\text{ (v/v) Triton X-100}$
  - $1\times$ cOmplete, EDTA-free Protease Inhibitor Cocktail (Roche 11873580001)
  - $40\text{ U/mL}$ RNasin Plus RNase Inhibitor (Promega N2615)
- **High-Salt Wash Buffer (Chaperone Removal)**:
  - $50\text{ mM}$ HEPES-KOH, $\text{pH } 7.5$
  - $500\text{ mM}$ $\text{KCl}$
  - $10\%\text{ (v/v) glycerol}$
  - $1.0\text{ mM}$ $\text{TCEP}$
  - $0.01\%\text{ (v/v) Triton X-100}$
- **On-Column Cleavage Buffer**:
  - $20\text{ mM}$ HEPES-KOH, $\text{pH } 7.5$
  - $150\text{ mM}$ $\text{KCl}$
  - $2.0\text{ mM}$ $\text{MgCl}_2$
  - $10\%\text{ (v/v) glycerol}$
  - $1.0\text{ mM}$ $\text{TCEP}$
- **Size-Exclusion Chromatography (SEC) Running Buffer**:
  - $20\text{ mM}$ HEPES-KOH, $\text{pH } 7.5$
  - $150\text{ mM}$ $\text{KCl}$
  - $2.0\text{ mM}$ $\text{MgCl}_2$
  - $1.0\text{ mM}$ $\text{TCEP}$
  - $0.02\%\text{ (w/v) sodium azide}$ (optional, omit for enzymatic assays)

---

## 3. Step-by-Step Experimental Procedure

### 3.1 Expi293F Co-Transfection (Day 0)
1. Maintain Expi293F cells in shaker flasks at $37\text{ }^\circ\text{C}$, $8\%\text{ CO}_2$, and $125\text{ rpm}$ until viability is $>95\%$ and cell density is $2.8\text{--}3.2 \times 10^6\text{ cells/mL}$.
2. Prepare plasmid DNA mixture per $1.0\text{ L}$ culture:
   - $0.25\text{ mg}$ `pcDNA3.1(+)-TwinStrep-FLAG-3C-hTERT` (WT, K757A, or K757R).
   - $0.75\text{ mg}$ `pcDNA3.1(+)-hTR` (expressing full-length 451-nt human telomerase RNA).
   - Dilute combined DNA ($1.0\text{ mg}$ total) into $50\text{ mL}$ Opti-MEM I Reduced Serum Medium.
3. Dilute $2.7\text{ mL}$ ExpiFectamine 293 Reagent into $50\text{ mL}$ Opti-MEM I. Incubate for 5 min at room temperature.
4. Combine DNA and ExpiFectamine mixtures, mix gently, and incubate for 15 min at room temperature to allow lipid-DNA complex formation.
5. Add the $100\text{ mL}$ transfection mixture slowly to the $1.0\text{ L}$ Expi293F cell suspension.
6. Return flasks to the $37\text{ }^\circ\text{C}$ incubator with $8\%\text{ CO}_2$ on an orbital shaker at $125\text{ rpm}$.
7. At 18–20 hours post-transfection, add $5.0\text{ mL}$ ExpiFectamine 293 Transfection Enhancer 1 and $50\text{ mL}$ Enhancer 2.

### 3.2 Cell Harvest and Lysis (Day 2 / 48 h Post-Transfection)
1. Harvest cells at $48\text{ h}$ post-transfection by centrifugation at $1,000 \times g$ for $20\text{ min}$ at $4\text{ }^\circ\text{C}$.
2. Wash cell pellets once with ice-cold $1\times\text{ PBS}$ ($150\text{ mM NaCl}, 10\text{ mM Na}_2\text{HPO}_4, \text{pH } 7.4$).
3. Resuspend pellet in $50\text{ mL}$ ice-cold Lysis Buffer per $1.0\text{ L}$ harvest.
4. Lyse cells by gentle rotation for $30\text{ min}$ at $4\text{ }^\circ\text{C}$, followed by 3 cycles of Dounce homogenization (tight pestle B, 15 strokes on ice).
5. Clarify lysate by ultracentrifugation at $40,000 \times g$ for $45\text{ min}$ at $4\text{ }^\circ\text{C}$ in a Beckman Ti45 rotor.
6. Collect the clear supernatant, filtering through a $0.45\text{ }\mu\text{m}$ syringe filter to remove residual lipids.

### 3.3 Affinity Capture & On-Column Proteolysis (Day 2–3)
1. Equilibrate $4.0\text{ mL}$ Strep-Tactin XT Superflow high-capacity resin (IBA Lifesciences 2-4010-010) with 10 column volumes (CV) of Lysis Buffer.
2. Incubate clarified supernatant with equilibrated resin in batch mode with gentle tumbling at $4\text{ }^\circ\text{C}$ for $3\text{ hours}$.
3. Pack resin into a gravity column. Collect flow-through for diagnostic Western blotting.
4. Wash column with $20\text{ CV}$ of Lysis Buffer.
5. Wash column with $30\text{ CV}$ of High-Salt Wash Buffer ($500\text{ mM KCl}$) to strip non-specifically bound chaperone complexes (e.g., Hsp90/p23).
6. Equilibrate column with $15\text{ CV}$ of On-Column Cleavage Buffer ($150\text{ mM KCl}$).
7. Resuspend resin in $4.0\text{ mL}$ Cleavage Buffer containing $100\text{ }\mu\text{g}$ recombinant GST-tagged PreScission 3C Protease (Cytiva 27084301).
8. Seal column and incubate on a tube rotator at $4\text{ }^\circ\text{C}$ for $16\text{ hours}$ (overnight).

### 3.4 Elution and Polish via Size-Exclusion Chromatography (Day 3)
1. Collect the flow-through containing cleaved, untagged hTERT–hTR holoenzyme complex (the PreScission protease and cleaved tags remain bound to resin and glutathione beads).
2. Wash resin with an additional $4.0\text{ mL}$ Cleavage Buffer; combine eluates ($8.0\text{ mL}$ total).
3. Pass eluate through $0.5\text{ mL}$ Glutathione Sepharose 4B resin to remove any trace GST-3C protease.
4. Concentrate the flow-through to $<0.5\text{ mL}$ using an Amicon Ultra-15 $100\text{ kDa}$ MWCO centrifugal filter unit ($3,000 \times g$ at $4\text{ }^\circ\text{C}$).
5. Inject concentrated sample onto a Superose 6 Increase 10/300 GL column pre-equilibrated in SEC Running Buffer ($20\text{ mM}$ HEPES pH 7.5, $150\text{ mM}$ KCl, $2\text{ mM}$ $\text{MgCl}_2$, $1\text{ mM}$ TCEP) on an ÄKTA pure system at $0.4\text{ mL/min}$.
6. Collect $0.5\text{ mL}$ fractions while monitoring UV absorbance at $260\text{ nm}$ and $280\text{ nm}$.
7. The intact holoenzyme elutes as a sharp, monodisperse peak at elution volume $V_e \approx 13.5\text{--}14.2\text{ mL}$ (apparent molecular weight $\approx 450\text{ kDa}$).

---

## 4. Quality Control & Acceptance Criteria

```
========================================================================================================================
QUALITY CONTROL GATES AND PASS/FAIL THRESHOLDS
========================================================================================================================
Quality Metric          Method                           Passing Criteria                           Action if Failed
------------------------------------------------------------------------------------------------------------------------
Purity (Protein)        SDS-PAGE (4-12% Bis-Tris)        > 90% single band at ~125 kDa (TERT)       Re-chromatograph SEC
                        Coomassie / Silver Stain

RNA Integrity           Denaturing 6% Urea-PAGE          Clean single band at 451 nt;               Discard; re-transfect
                        SYBR Gold RNA Staining           zero low-MW degradation fragments          fresh hTR plasmid

Stoichiometry           UV A260 / A280 Ratio             1.15 to 1.35                               Investigate RNA excess;
                                                         (indicates 1:1 TERT:hTR complex)           adjust plasmid ratio

Dispersity / Folding    nanoDSF (Prometheus NT.48)       Monodisperse; single thermal transition;   Discard aggregated
                                                         Tm >= 48.5 °C (WT baseline)                preparations

Concentration           BCA Assay / NanoDrop             >= 1.0 µM (>= 0.45 mg/mL)                  Concentrate with 100 kDa
                                                                                                    Amicon filter
========================================================================================================================
```

---

## 5. Storage and Sample Handling

1. Pool peak SEC fractions meeting all QC criteria.
2. Determine final concentration via $A_{280}$ ($\epsilon_{280} = 145,000\text{ M}^{-1}\text{cm}^{-1}$) and BCA protein assay.
3. Flash-freeze $10\text{ }\mu\text{L}$ single-use aliquots in pre-chilled thin-wall PCR tubes in liquid nitrogen.
4. Store at $-80\text{ }^\circ\text{C}$. Avoid repeated freeze-thaw cycles (maximum 1 freeze-thaw allowed).
