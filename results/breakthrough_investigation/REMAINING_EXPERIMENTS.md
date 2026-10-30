# REMAINING EXPERIMENTS: THE EXPERIMENTAL VALIDATION ROADMAP

**Document ID:** `BSV-DOC-005-REV-FINAL`  
**Date:** 2026-10-30  
**Classification:** Experimental Requirements & Validation Roadmap  
**Target Repository:** `results/breakthrough_investigation/REMAINING_EXPERIMENTS.md`  

---

## 1. Executive Summary

This catalog outlines all outstanding laboratory and high-performance computational procedures required to convert the **TEN–TPP1–TRAP Bipartite Clamp Discovery** from a validated computational hypothesis into a demonstrated biological breakthrough.

Each experiment is defined with its specific objective, biological system, technical protocol, instrumentation, quantitative deliverables, estimated timeline, and decisive falsification criteria.

---

## 2. Master Table of Outstanding Experiments

```
========================================================================================================================
MASTER EXPERIMENTAL ROADMAP
========================================================================================================================
ID    Category       Experiment Name                         System / Reagent          Key Readout          Timeline
------------------------------------------------------------------------------------------------------------------------
EX-01 Biochemistry   Recombinant Expression & Purification   Expi293F / MultiBac Sf9   Yield, Stoichiometry 3 weeks
EX-02 Biophysics     Thermal Stability Screening (nanoDSF)   Prometheus NT.48          Melting Temp (Tm)    1 week
EX-03 Biophysics     TPP1 Binding Affinity (SPR / BLI)       Biacore 8K / Octet RH16   Kd, k_on, k_off      1 week
EX-04 Enzymology     Direct Primer Extension Gel RAP Assay   12% Urea-PAGE / Typhoon   chi_RAP, Rp          2 weeks
EX-05 Kinetics       Pulse-Chase & smFRET Primer Retention   TIRF / Sepharose Chase    k_off^translocation  3 weeks
EX-06 Structural Bio Solution HDX-MS Dynamic Footprinting    Orbitrap Eclipse Tribrid  Deuteration (ΔD %)   4 weeks
EX-07 Structural Bio Cryo-EM 3D Variability Analysis (3DVA)  Titan Krios / CryoSPARC   Continuous modes     4 weeks
EX-08 Computation    HPC All-Atom Molecular Dynamics (18 µs) GPU Cluster (H100)        RMSF, PMF Barrier    3 weeks
EX-09 Cell Biology   Cellular Telomere Dynamics & Viability  HeLa / CD34+ HSCs         PDL, Senescence, TRF 8 weeks
========================================================================================================================
```

---

## 3. Detailed Experimental Specifications

### Experiment EX-01: Recombinant Holoenzyme Expression & Purification
- **Objective**: Produce milligram-scale, monodisperse human telomerase holoenzymes for wild-type and the 6 candidate mutants (`K757A`, `K757R`, `Y176A`, `Y176F`, `Q794A`, `R756A`).
- **Protocol**:
  - Transfect Expi293F cells at $3.0 \times 10^6\text{ cells/mL}$ with pCDNA3.1-TwinStrep-FLAG-TERT and pBluescript-hTR (1:3 plasmid ratio).
  - Harvest cells after 48 hours; lyse in hypotonic buffer supplemented with 0.1% Triton X-100, protease inhibitors, and RNasin Plus.
  - Affinity capture on Strep-Tactin XT Superflow resin, wash with 500 mM KCl buffer, and cleave with 3C protease at 4 °C for 16 hours.
  - Polish via Superose 6 Increase 10/300 GL size-exclusion chromatography.
- **Deliverables**: Purified holoenzymes (>90% purity by SDS-PAGE, stoichiometric hTR RNA incorporation verified by urea-PAGE, concentration $\ge 1.0\text{ }\mu\text{M}$).

### Experiment EX-02: Thermal Stability Screening (nanoDSF)
- **Objective**: Confirm that point mutations do not cause global misfolding or destabilization of the catalytic core.
- **Protocol**:
  - Load $10\text{ }\mu\text{L}$ of $0.5\text{ }\mu\text{M}$ holoenzyme into Prometheus high-sensitivity capillaries.
  - Monitor intrinsic tryptophan/tyrosine fluorescence ($F_{350}/F_{330}$) across a thermal gradient from $20\text{ }^\circ\text{C}$ to $95\text{ }^\circ\text{C}$ at $1.0\text{ }^\circ\text{C/min}$.
- **Falsification Threshold**: If $T_m(\text{mutant}) < T_m(\text{WT}) - 2.5\text{ }^\circ\text{C}$, the mutant is classified as a global folding defect and disqualified.

### Experiment EX-03: TPP1 Recruitment Binding Affinity (SPR / BLI)
- **Objective**: Quantify equilibrium dissociation constants ($K_d$) and binding kinetics between the TERT catalytic core and TPP1 OB domain.
- **Protocol**:
  - Immobilize chemically biotinylated TPP1 OB domain (residues 89–249) onto Biacore Sensor Chip SA.
  - Inject TERT holoenzymes at concentrations of 0, 15, 30, 60, 120, 240, and 480 nM at $30\text{ }\mu\text{L/min}$ for $180\text{ s}$ association and $600\text{ s}$ dissociation.
  - Fit sensorgrams using a 1:1 Langmuir binding model with mass-transport limitation correction.
- **Breakthrough Threshold**: Primary candidate `K757A` must exhibit $K_d \le 300\text{ nM}$ (WT baseline: $150 \pm 25\text{ nM}$). If $K_d > 1,500\text{ nM}$, the uncoupling hypothesis is refuted.

### Experiment EX-04: Direct Primer Extension Radiometric Gel RAP Assay
- **Objective**: Measure basal nucleotide addition rate ($k_{\text{cat}}$) and TPP1-stimulated processivity enhancement factor ($\chi_{\text{RAP}}$).
- **Protocol**:
  - Incubate $50\text{ nM}$ TERT holoenzyme, $1.0\text{ }\mu\text{M}$ $(\text{TTAGGG})_3$ primer, $500\text{ }\mu\text{M}$ dATP, $500\text{ }\mu\text{M}$ dTTP, $5.0\text{ }\mu\text{M}$ dGTP, and $0.15\text{ }\mu\text{M}$ $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ in reaction buffer at $30\text{ }^\circ\text{C}$ for 60 min.
  - Titrate recombinant POT1–TPP1 at 0, 50, 100, 250, and 500 nM.
  - Resolve on 12% denaturing polyacrylamide gels (7 M urea), expose to phosphor screens, and quantify band intensities using ImageQuant TL.
- **Breakthrough Threshold**:
  - Basal activity: $k_{\text{cat}}(\text{mutant}) \ge 85\%\text{ of } k_{\text{cat}}(\text{WT})$.
  - RAP enhancement: $\chi_{\text{RAP}}(\text{K757A}) \le 1.8$ (versus $\chi_{\text{RAP}}(\text{WT}) \ge 5.0$).
  - Rescued control: $\chi_{\text{RAP}}(\text{K757R}) \ge 4.5$.

### Experiment EX-05: Translocation Chase & Single-Molecule FRET (smFRET)
- **Objective**: Prove that RAP collapse is driven specifically by premature primer dissociation during translocation.
- **Protocol**:
  - Surface-immobilize Cy3-labeled telomerase holoenzyme in a microfluidic flow cell.
  - Deliver Cy5-labeled $(\text{TTAGGG})_3$ primer in the presence of nucleotides and recombinant TPP1.
  - Monitor donor/acceptor fluorescence via total internal reflection fluorescence (TIRF) microscopy at 50 ms time resolution.
  - Extract translocation dwell times and determine dissociation rate $k_{\text{off}}^{\text{translocation}}$.
- **Breakthrough Threshold**: In `K757A` + TPP1, $k_{\text{off}}^{\text{translocation}}$ must remain $\ge 0.080\text{ s}^{-1}$ (versus WT + TPP1 = $0.012\text{ s}^{-1}$), confirming failure of primer retention during translocation.

### Experiment EX-06: Solution HDX-MS Dynamic Footprinting
- **Objective**: Demonstrate physical dynamic coupling between the TEN domain and TRAP domain upon TPP1 binding in solution.
- **Protocol**:
  - Expose purified holoenzyme (± TPP1) to $\text{D}_2\text{O}$ buffer for time points $10\text{ s}, 60\text{ s}, 600\text{ s}, 3,600\text{ s}$.
  - Quench with 2 M guanidine-HCl, pH 2.5, digest on-column with immobilized pepsin, and analyze deuteration levels via LC-MS on an Orbitrap Eclipse.
- **Breakthrough Threshold**: Peptides covering residues 170–182 (TEN loop) and 752–765 (TRAP pore rim) must show statistically significant reduction in deuterium uptake ($\Delta D \ge 15\%$, $p < 0.01$) upon TPP1 addition.

### Experiment EX-07: Cryo-EM 3D Variability Analysis (3DVA)
- **Objective**: Visualize continuous conformational dynamics and verify that TPP1 binding quenches inter-domain rocking.
- **Protocol**:
  - Collect 5,000 micrographs of apo-holoenzyme and TPP1-bound holoenzyme on a Titan Krios (300 kV, Gatan K3 detector).
  - Reconstruct consensus maps in CryoSPARC v4.x.
  - Perform 3D Variability Analysis (3DVA) solving for 3 continuous principal eigenmodes.
- **Expected Result**: Apo-holoenzyme shows continuous domain breathing of TEN relative to TRAP (amplitude $\sim 3\text{--}4\text{ \AA}$); TPP1-bound complex shows complete suppression of breathing modes.

### Experiment EX-08: GPU All-Atom Explicit Solvent Molecular Dynamics
- **Objective**: Quantify thermodynamic free-energy barriers of primer retention via high-throughput simulation.
- **Protocol**:
  - Run $3 \times 1\text{ }\mu\text{s}$ per system (WT apo, WT bound, K757A apo, K757A bound, Y176F apo, Y176F bound; aggregate $18\text{ }\mu\text{s}$) under CHARMM36m at 303 K.
  - Perform Steered MD and Umbrella Sampling across 35 reaction coordinate windows to compute Potential of Mean Force (PMF).
- **Expected Result**: PMF unbinding barrier $\Delta G^\ddagger$ drops by $\ge 2.5\text{ kcal/mol}$ in `K757A` relative to WT bound.

### Experiment EX-09: Cellular Telomere Dynamics in Cancer vs Stem Cells
- **Objective**: Evaluate the therapeutic window hypothesis: does uncoupling RAP suppress cancer cell immortalization while sparing normal stem-cell viability?
- **Protocol**:
  - Generate CRISPR/Cas9 homozygous knock-in cell lines for `TERT-K757A` in HeLa (cancer) and primary CD34+ human hematopoietic stem cells (HSCs).
  - Monitor population doubling levels (PDL), senescence markers ($\text{p}16^{\text{INK4a}}$, $\beta$-galactosidase), and telomere length via Telomere Restriction Fragment (TRF) Southern blotting over 120 days.
- **Expected Result**: HeLa-K757A cells undergo progressive telomere shortening and enter crisis/apoptosis after 40–50 PDL. CD34+ HSCs maintain sufficient short-term telomere capping to support multilineage reconstitution.

---

## 4. Resource Requirements and Execution Timeline

```
========================================================================================================
BUDGET & RESOURCE REQUIREMENTS
========================================================================================================
Resource Type               Specification                                Estimated Cost    Source
--------------------------------------------------------------------------------------------------------
Reagents & Synthesis        Genes, plasmids, primers, radiochemicals     $18,000 USD       Institutional Grant
Instrument Time             Biacore 8K (40 h), Prometheus (20 h)         $8,500 USD        Core Facility
Instrument Time             Orbitrap Eclipse HDX-MS (60 h)               $12,000 USD       Core Facility
Instrument Time             Titan Krios Cryo-EM (48 h)                   $15,000 USD       Cryo-EM Consortium
HPC Compute Hours           100,000 GPU-hours (NVIDIA H100)              $25,000 USD       National Supercomputing
Cell Culture & In Vivo      Expi293F, CD34+ HSCs, consumables            $22,000 USD       Operating Grant
--------------------------------------------------------------------------------------------------------
TOTAL ESTIMATED BUDGET                                                   $100,500 USD      18 Months Total
========================================================================================================
```
