# Research Framing & Falsifiable Hypotheses: Human Telomerase TERT–TPP1 Recruitment Interface

> **Superseded — 2026-11-01.** These questions were narrowed during the causality phase. The
> current statement of what can and cannot be asked of a processivity measurement is in
> [`paper/MANUSCRIPT.md`](../../paper/MANUSCRIPT.md) and
> [`results/causality/MECHANISM_CLASSIFICATION.md`](../causality/MECHANISM_CLASSIFICATION.md).

**Document Role**: Primary Research Question & Falsifiable Hypotheses Formulation  
**Milestone**: Benchmark Milestone B2 (Falsifiable Hypotheses Formulation)  
**Deliverable Path**: `results/research/question_and_hypotheses.md`  
**Date**: October 2026  
**Status**: Authoritative Milestone B2 Deliverable  

---

## 1. Executive Summary & Problem Formulation

Human telomerase recruitment to the telomere is mediated by a specialized macromolecular interaction between the catalytic subunit **TERT** (telomerase reverse transcriptase) and the shelterin accessory protein **TPP1** (specifically its N-terminal oligosaccharide/oligonucleotide-binding [OB]-fold domain). Cryo-EM reconstructions—including the TPP1-bound active catalytic core (PDB: `7TRE`, 3.5 Å), the unbound catalytic core (PDB: `7TRD`, 3.3 Å), and the histone-bound holoenzyme (PDB: `7TRF`, 3.7 Å)—together with high-resolution X-ray crystallography of the isolated TPP1 OB-fold (PDB: `2I46`, 2.7 Å), provide the structural foundation for understanding recruitment mechanics.

However, an acute scientific dilemma emerges when interpreting these coordinate sets:
1. **The Rigid Docking Paradigm**: The TERT catalytic core (residues 326–1132) aligns between TPP1-bound (`7TRE`) and free (`7TRD`) states with an extraordinarily low Cα root-mean-square deviation (**RMSD = 0.338 Å**), below typical cryo-EM coordinate precision (~0.4–0.5 Å), suggesting a rigid, pre-organized evolutionary joint.
2. **The Dynamic Remodeling Paradigm**: Comparison of the bound TPP1 subunit in `7TRE` with the unbound crystal structure `2I46` reveals extensive localized conformational remodeling, concentrated sharply within the acidic **DE-loop** (residues 168–171), peaking at **Lys170 with a 4.293 Å Cα displacement**, alongside subtle rigid-body clamping of the TERT TEN domain (0.434 Å excursion) and TRAP/IFD hairpin (0.427 Å excursion).
3. **The Parameter Sensitivity Paradigm**: Across standard computational threshold boundaries (3.5 Å, 4.0 Å, 4.5 Å), the apparent interface contact network undergoes an explosive 13.67-fold expansion (from **9 contacts** across 5 residue pairs to **123 contacts** across 32 residue pairs), while buried surface area (BSA) varies substantially with probe radius and atomic van der Waals conventions.

This document establishes the **Primary Research Question** and formulates **Three Competing, Mutually Falsifiable Scientific Hypotheses (H1, H2, H3)** with explicit, quantitative criteria for empirical support and rejection.

---

## 2. The Primary Research Question

> **Primary Research Question**:  
> **Does human telomerase recruitment via the TERT–TPP1 interface operate as an evolutionarily conserved, rigid structural joint whose observed conformational variations are experimental and crystallographic artifacts (H1); or does it function as a genuine state-dependent dynamic interface driven by induced-fit DE-loop remodeling and bipartite domain clamping that mediates functional regulation (H2); or are the reported structural differences and contact architectures predominantly artifacts of computational cutoff and parameter sensitivities (H3)?**

### Sub-Questions
1. **Mechanistic Sub-Question (Axis A)**: To what extent does TPP1 binding alter the structural or dynamic landscape of the TERT bipartite cradle (TEN and TRAP/IFD) and the catalytic core, and can static coordinate identity (0.338 Å RMSD) legitimately rule out dynamic allostery?
2. **Conformational Plasticity Sub-Question**: Does the 4.293 Å peak displacement of TPP1 Lys170 reflect a physiological induced-fit binding transition, or is it an artifact of crystal packing lattice forces in `2I46` and cryo-EM map coordinate uncertainty in `7TRE`?
3. **Computational Methodological Sub-Question**: How much of the reported interface complexity (BSA = 856.82 Å², 9 vs 46 vs 123 contacts) is physically invariant versus computationally tunable by radius conventions, probe radii, and distance cutoffs?
4. **Translational Sub-Question (Axis B)**: If the interface is an invariant evolutionary joint (H1), does any realistic therapeutic window exist for recruitment-targeted small-molecule or peptide inhibitors (e.g., TELIN-1 candidates) to selectively suppress cancer cell proliferation without triggering catastrophic on-target stem cell failure (Dyskeratosis Congenita / bone marrow aplasia)?

---

## 3. Competing Scientific Hypotheses

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        COMPETING SCIENTIFIC HYPOTHESES                                 │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                        │
│   ┌────────────────────────┐  ┌────────────────────────┐  ┌────────────────────────┐   │
│   │     HYPOTHESIS 1       │  │     HYPOTHESIS 2       │  │     HYPOTHESIS 3       │   │
│   │  Rigid Conservation &  │  │ Genuine State-Dependent│  │     Computational      │   │
│   │ Experimental Artifact  │  │  Remodeling / Induced  │  │ Parameter Sensitivity  │   │
│   └───────────┬────────────┘  └───────────┬────────────┘  └───────────┬────────────┘   │
│               │                           │                           │                │
│               ▼                           ▼                           ▼                │
│   • Interface is rigid joint  • Genuine induced-fit dock  • Network is computationally │
│   • RMSD <= coord error       • DE-loop remodeling real   │  tunable artifact          │
│   • 2I46 loop shift is        • TEN/TRAP dynamic clamping │ • 9 -> 46 -> 123 contacts  │
│     crystal packing artifact  • Dynamic allostery in RAP  │ • Cutoff & radius driven   │
│   • Static key-in-lock dock   • Multi-state energy barrier│ • Pose / affinity unproven │
│                                                                                        │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

### Hypothesis 1 (H1): Evolutionary Conservation & Rigid Architectural Joint (Experimental / Crystal Artifact)

#### 1. Theoretical Premise & Physical Model
The TERT–TPP1 interface is an evolutionarily conserved, structurally pre-organized rigid macromolecular joint. Telomerase recruitment relies on an invariant "key-in-lock" geometric fit between pre-formed hydrophobic anchors on TPP1 (Leu183, Leu212) and the pre-formed TERT receptor cradle (Leu769, Pro771, Arg774, Lys78).

Observed coordinate discrepancies across PDB entries (`7TRE`, `7TRD`, `7TRF`, `2I46`) do **not** reflect functional biological plasticity, but are entirely explained by:
1. **Experimental Coordinate Uncertainty**: At 3.3–3.5 Å cryo-EM resolution, coordinate precision is bounded by Cruickshank/DPI limits of ~0.35–0.50 Å. The global Cα RMSD of 0.356 Å and catalytic core RMSD of 0.338 Å fall strictly within this experimental noise floor.
2. **Crystallographic Lattice Packing Artifacts**: In `2I46` (2.7 Å), Chain A participates in extensive crystal lattice contacts that artificially constrain the DE-loop, while Chain B exhibits complete loop disorder (unmodeled residues 168–171) and high temperature factors ($B > 80\ \text{Å}^2$). The 4.293 Å displacement of Lys170 is a lattice contact artifact rather than a physiological binding pathway.
3. **Recombinant Construct Tags**: Truncation of the N-terminus and engineered affinity tags introduce minor peripheral coordinate shifts that have no bearing on physiological assembly.

#### 2. Quantitative Support Criteria (Empirical Verification)
To establish empirical support for **H1**, all of the following criteria must be met:
1. **Core Coordinate Invariance**: The catalytic core Cα RMSD between `7TRE` and `7TRD` satisfies:
   $$\text{RMSD}_{\text{core}} \le 0.40\ \text{Å}$$
   and the local coordinate error $\sigma_{\text{coord}}$ estimated by Cruickshank DPI satisfies:
   $$\text{RMSD}_{\text{core}} \le \sigma_{\text{coord}} \quad (\text{where } \sigma_{\text{coord}} \approx 0.35\text{--}0.50\ \text{Å})$$
2. **Hydrophobic Anchor Pre-Organization**: The structural anchor residues responsible for initial engagement (TPP1 Leu183 and Leu212) show negligible Cα displacement between unbound (`2I46` Chain A) and bound (`7TRE` Chain C):
   $$\Delta r_{\text{C}\alpha}(\text{Leu183}) \le 0.60\ \text{Å} \quad (\text{Observed: } 0.512\ \text{Å})$$
   $$\Delta r_{\text{C}\alpha}(\text{Leu212}) \le 1.10\ \text{Å} \quad (\text{Observed: } 1.040\ \text{Å})$$
3. **Crystallographic Lattice Correlation**: Across all modeled residues of `2I46` Chain A, per-residue Cα displacement relative to `7TRE` correlates positively and significantly with crystallographic B-factors and proximity to neighboring crystal symmetry mates:
   $$\text{Pearson } r(\Delta r_{\text{C}\alpha}, B_{\text{2I46}}) \ge 0.55 \quad (p < 0.001)$$
   $$\text{Distance to symmetry mate in 2I46 for DE-loop (168--171)} < 4.0\ \text{Å}$$
4. **Solution Ensemble Equivalence**: Molecular dynamics (MD) simulations of isolated apo-TPP1 in explicit solvent sample conformations with DE-loop RMSD $\le 1.5\ \text{Å}$ relative to `7TRE` without crossing a major free energy barrier:
   $$\Delta G^{\ddagger}_{\text{remodel}} \le 1.0\ \text{kcal/mol}$$

#### 3. Quantitative Falsification Criteria (Empirical Rejection)
**H1** is decisively **falsified** if any of the following occur:
1. **Excursion Beyond Error Envelope**: Cryo-EM local map coordinate uncertainty analysis (e.g., via 3D flex, cryoSPARC local resolution, or ResMap) demonstrates that the TEN domain excursion (0.434 Å) or TRAP/IFD shift (0.427 Å) exceeds the local directional coordinate variance with statistical significance:
   $$p < 0.01 \quad (\text{paired directional displacement test})$$
2. **Crystal Artifact Disproof**: High-resolution solution NMR spectroscopy or cryo-EM of isolated apo-TPP1 (free of lattice contacts) confirms that the ground-state solution conformation of the DE-loop deviates from the `7TRE` bound state by:
   $$\text{C}\alpha\text{ RMSD}(\text{DE-loop}_{168\text{--}171}) > 2.50\ \text{Å}$$
   demonstrating that the open or displaced state is an intrinsic solution property, not a crystal artifact.
3. **Conformational Free Energy Barrier**: Free energy profiling (umbrella sampling or metadynamics) shows that transitioning the DE-loop from the apo-like state to the `7TRE` bound state requires overcoming a substantial thermodynamic barrier:
   $$\Delta G^{\ddagger} > 3.0\ \text{kcal/mol}$$
   precluding passive, rigid encounter docking.

---

### Hypothesis 2 (H2): Genuine State-Dependent Remodeling & Induced-Fit Clamping

#### 1. Theoretical Premise & Physical Model
The TERT–TPP1 recruitment interface is a state-dependent, functionally plastic macromolecular junction. Telomerase recruitment proceeds via a multi-step **induced-fit conformational remodeling** mechanism:
1. **Encounter Complex Formation**: Pre-organized hydrophobic anchors (Leu183, Leu212) make initial, low-affinity encounter contacts with the TERT TRAP/IFD hydrophobic groove (Leu769, Pro771).
2. **Induced-Fit DE-Loop Clamping**: Upon encounter, the intrinsically mobile and negatively charged TPP1 DE-loop (Glu168, Glu169, Lys170, Glu171) undergoes an extensive structural rearrangement (peaking at Lys170 with 4.293 Å Cα displacement) to collapse into the basic cradle of TERT.
3. **Bipartite Cradle Accommodation**: The TERT TEN domain undergoes a coordinated rigid-body hinge motion (0.434 Å excursion) and the TRAP/IFD hairpin flexes (0.427 Å excursion) to lock the secondary salt latch (Lys78–Glu215, 3.58 Å) and primary salt bridge (Arg774–Glu168, 2.32 Å).
4. **Functional Allosteric Regulation**: This conformational transition alters low-frequency normal modes and internal dynamics across the catalytic core, stabilizing the primer-template duplex and stimulating repeat addition processivity (RAP).

#### 2. Quantitative Support Criteria (Empirical Verification)
To establish empirical support for **H2**, all of the following criteria must be met:
1. **Statistically Significant Loop Remodeling**: The Cα displacement of the DE-loop (residues 168–172) exceeds the combined experimental coordinate standard deviation $\sigma_{\text{combined}} = \sqrt{\sigma_{\text{7TRE}}^2 + \sigma_{\text{2I46}}^2} \approx 0.65\ \text{Å}$ by $> 4\sigma$:
   $$\Delta r_{\text{C}\alpha}(\text{Lys170}) \ge 3.0\ \text{Å} \quad (\text{Observed: } 4.293\ \text{Å} > 6.6\sigma)$$
   $$\text{Mean }\Delta r_{\text{C}\alpha}(\text{DE-loop}_{168\text{--}171}) \ge 2.5\ \text{Å} \quad (\text{Observed: } 3.169\ \text{Å} > 4.8\sigma)$$
2. **Directionally Correlated Cradle Adaptation**: The rigid-body excursion vectors of the TEN domain ($\Delta \vec{r}_{\text{TEN}} = 0.434\ \text{Å}$) and TRAP/IFD hairpin ($\Delta \vec{r}_{\text{TRAP}} = 0.427\ \text{Å}$) are oriented vectorially toward the center-of-mass of the docked TPP1 OB-fold:
   $$\cos \theta(\Delta \vec{r}_{\text{domain}}, \vec{r}_{\text{dock}}) > 0.70 \quad (p < 0.05 \text{ vs isotropic random motion})$$
3. **Dynamic Allostery & Normal Mode Remodeling**: Elastic network models (ENM) or atomistic molecular dynamics reveal that TPP1 binding induces a significant shift in low-frequency vibrational eigenvalues ($\Delta \omega_1 / \omega_1 > 15\%$) and alters root-mean-square fluctuations (RMSF) in the RNA-template duplex cleft by $> 0.5\ \text{Å}$, demonstrating dynamic allosteric coupling despite static Cα invariance.
4. **Kinetic Multi-Step Binding**: Stopped-flow fluorescence or surface plasmon resonance (SPR) biophysical kinetics fit significantly better to a two-step induced-fit binding model ($A + B \rightleftharpoons AB^* \rightleftharpoons AB$) over a simple one-step bimolecular model ($A + B \rightleftharpoons AB$):
   $$\Delta\text{AIC} > 10.0 \quad (\text{Akaike Information Criterion preference for 2-step model})$$
   $$k_{\text{isomerization}} > 0 \quad (\text{statistically significant forward rate constant } p < 0.01)$$

#### 3. Quantitative Falsification Criteria (Empirical Rejection)
**H2** is decisively **falsified** if any of the following occur:
1. **Zero Dynamic Coupling**: Solution-state HDX-MS (hydrogen-deuterium exchange mass spectrometry) or NMR relaxation dispersion shows zero change in backbone amide exchange rates or conformational entropy in the TERT catalytic core upon TPP1 binding:
   $$\Delta\text{HDX}_{\text{core}} < 2\% \quad (\text{across all peptides in RT/TRBD domains, } p > 0.20)$$
2. **Pre-Existing Rigid Solution State**: Solution NMR or cryo-EM reconstruction of isolated apo-TPP1 demonstrates that the DE-loop adopts the `7TRE` conformation in the absence of TERT:
   $$\text{RMSD}_{\text{NMR}}(\text{DE-loop vs 7TRE}) \le 1.0\ \text{Å}$$
3. **Strict One-Step Association Kinetics**: Real-time kinetic measurements (stopped-flow / BLI) demonstrate strict pseudo-first-order single-exponential binding kinetics across an escalating concentration series ($k_{\text{obs}}$ linear with $[L]$):
   $$R^2(\text{linear } k_{\text{obs}} \text{ vs } [L]) \ge 0.995 \quad (\text{no hyperbolic saturation indicative of induced fit})$$
4. **Conformational Rigidity Invariance**: Engineering rigid crosslinks (e.g., disulfide staples) locking the DE-loop in its `7TRE` bound conformation yields binding affinity ($K_d$) and telomerase repeat addition processivity (RAP) indistinguishable from wild-type:
   $$\Delta K_d \le 1.2\text{-fold} \quad \text{and} \quad \Delta\text{RAP} \le 10\%$$

---

### Hypothesis 3 (H3): Computational & Parameter Sensitivity (Methodological Artifact)

#### 1. Theoretical Premise & Physical Model
Reported quantitative metrics describing the TERT–TPP1 interface—including contact counts, polar network complexity, buried surface area (BSA), and domain excursion metrics—are predominantly driven by arbitrary computational threshold choices, geometric cutoff parameters, and algorithmic conventions rather than invariant physical ground truth.

Specifically:
1. **Non-Linear Contact Horizon**: Adjusting the heavy-atom distance threshold from 3.5 Å to 4.5 Å produces a non-linear, >13-fold explosion in detected contacts ($9 \rightarrow 46 \rightarrow 123$), shifting the qualitative interpretation from a discrete salt-bridge latch into an extensive, diffuse surface.
2. **SASA / BSA Convention Sensitivity**: Buried surface area calculations are heavily dependent on probe sphere radius (1.2 Å vs 1.4 Å vs 1.6 Å), van der Waals radius sets (Bondi vs ProtOr vs Chothia), and numerical surface tessellation algorithms (Shrake-Rupley vs Lee-Richards vs analytical), producing variations up to 250 Å² (~30%).
3. **Alignment Subset Bias**: Structural superposition RMSD is highly sensitive to the choice of alignment subset (Global TERT Cα 934 residues vs Catalytic Core 770 residues vs secondary structure elements) and whether outlier rejection (trimmed Kabsch / PyMOL `align`) is permitted.

#### 2. Quantitative Support Criteria (Empirical Verification)
To establish empirical support for **H3**, all of the following criteria must be met:
1. **Contact Network Multiplicative Volatility**: Contact counts exhibit extreme threshold dependence, where $> 60\%$ of total contacts at 4.5 Å fall into the "Threshold Sensitive" category (absent at 4.0 Å):
   $$\frac{N_{\text{contacts}}(4.5\ \text{Å})}{N_{\text{contacts}}(3.5\ \text{Å})} \ge 10.0 \quad (\text{Observed: } 123 / 9 = 13.67)$$
   $$\frac{N_{\text{Threshold Sensitive}}}{N_{\text{total}}(4.5\ \text{Å})} \ge 60\% \quad (\text{Observed: } 77 / 123 = 62.6\%)$$
2. **Buried Surface Area Dispersion**: Systematically varying the probe radius from 1.2 Å to 1.6 Å and swapping vdW definitions between Bondi (1964) and ProtOr (Tsai et al., 1999) yields a coefficient of variation ($\text{CV}$) in calculated BSA:
   $$\text{CV}(\text{BSA}) \ge 15.0\% \quad \text{or} \quad \Delta\text{BSA} \ge 200\ \text{Å}^2$$
3. **Alignment Selection Instability**: Altering the reference alignment selection between Global TERT (934 Cα) and Catalytic Core (770 Cα) alters reported domain excursions by $> 25\%$:
   $$\frac{|\text{Displacement}_{\text{global}} - \text{Displacement}_{\text{core\_fit}}|}{\text{Displacement}_{\text{core\_fit}}} \ge 25.0\%$$
4. **Cross-Software Divergence**: Running independent structural biology packages (e.g., Biopython vs Gemmi vs FreeSASA vs PyMOL) with default settings on identical raw mmCIF coordinates yields discrepancies in reported interface areas:
   $$\Delta\text{BSA}_{\text{packages}} > 100\ \text{Å}^2$$

#### 3. Quantitative Falsification Criteria (Empirical Rejection)
**H3** is decisively **falsified** if any of the following occur:
1. **Thermodynamic Invariance of Core Contacts**: Electrostatic binding energy calculations (e.g., Poisson-Boltzmann / MMPBSA) demonstrate that the 5 Core Robust residue pairs (present at 3.5 Å) account for $> 80\%$ of total electrostatic binding free energy:
   $$\frac{\sum \Delta G_{\text{elec}}(\text{Core Robust})}{\Delta G_{\text{elec}}(\text{Total Interface})} \ge 0.80$$
   proving that the threshold-sensitive contacts ($4.0\text{--}4.5\ \text{Å}$) represent negligible energetic noise rather than meaningful structural interactions.
2. **BSA Algorithmic Convergence**: Under strictly harmonized atomic radii and probe radii ($r = 1.40\ \text{Å}$), BSA calculated across independent software implementations (FreeSASA, Shrake-Rupley, Gemmi) converges within:
   $$\Delta\text{BSA} \le 2.0\% \quad (\le 17\ \text{Å}^2 \text{ out of } 856.8\ \text{Å}^2)$$
3. **Rank-Order Invariance of Domain Mobility**: The relative mobility ranking of TERT domains ($\text{TEN} > \text{TRAP} > \text{CTE} > \text{RT} > \text{TRBD}$) remains 100% invariant (Spearman $\rho = 1.0, p < 0.001$) regardless of whether global, core, or outlier-trimmed superposition is employed.

---

## 4. Multi-Hypothesis Comparative Decision Matrix

The following decision matrix contrasts the three competing hypotheses across physical drivers, quantitative predictions, translational implications, and decisive empirical discriminators:

| Evaluation Dimension | Hypothesis 1 (H1): Rigid Conservation & Lattice Artifact | Hypothesis 2 (H2): State-Dependent Induced Fit | Hypothesis 3 (H3): Computational Parameter Sensitivity |
| :--- | :--- | :--- | :--- |
| **Primary Physical / Computational Driver** | Evolutionary conservation, static thermodynamic minimum, experimental noise floor | Multi-state conformational landscape, dynamic allostery, induced-fit binding | Arbitrary threshold selection, algorithmic conventions, coordinate discretization |
| **Catalytic Core Cα RMSD (7TRE vs 7TRD)** | Expected $\le 0.40\ \text{Å}$ (within cryo-EM map coordinate error) | Compatible with low static RMSD ($\le 0.40\ \text{Å}$) via dynamic allostery | Artifactually sensitive to reference selection ($\Delta\text{RMSD} > 0.15\ \text{Å}$) |
| **TEN Domain Adaptation (7TRE vs 7TRD)** | Rigid-body drift within experimental noise ($\le 0.50\ \text{Å}$) | Genuine clamping excursion ($\ge 0.40\ \text{Å}$) correlated with docking axis | Metric dependent on alignment anchor residues |
| **TPP1 DE-Loop Behavior (7TRE vs 2I46)** | 4.293 Å shift is an artifact of 2I46 crystal lattice contacts & high B-factors | Genuine physiological induced-fit remodeling upon TERT engagement | Loop shift is real, but contact consequences depend on cutoff choice |
| **Interface Contacts (3.5 Å $\rightarrow$ 4.5 Å)** | Only the 9 core contacts are biophysically relevant; peripheral contacts are noise | Intermediate contacts (37) represent essential outer energetic stabilization | Contact count expansion ($9 \rightarrow 46 \rightarrow 123$) is a pure mathematical artifact |
| **Buried Surface Area (BSA)** | Standard non-obligate interface (~850 Å²) pre-organized for binding | Interface area expands dynamically during induced-fit clamp closure | BSA varies by $> 20\%$ depending on probe radius and atomic radius set |
| **Axis B: Cancer vs Aging Therapeutic Window** | **Extremely Narrow / Non-Existent**: Small molecules targeting static interface inevitably phenocopy germline DC mutations | **Potentially Viable**: Conformation-selective inhibitors could trap intermediate or open states, achieving partial/allosteric inhibition | **Unquantifiable from Docking**: Docked poses (e.g. TELIN-1) reflect parameter tuning rather than true binding free energy |
| **Decisive Empirical Test** | Solution NMR / cryo-EM of apo-TPP1 showing DE-loop pre-organized in bound state | Stopped-flow 2-step kinetics ($\Delta\text{AIC} > 10$) and HDX-MS dynamic coupling | Harmonized algorithmic benchmark proving metric volatility across parameter sweeps |

---

## 5. Multi-Phase Empirical Discrimination Roadmap

To systematically discriminate between H1, H2, and H3 across the remainder of this benchmark evaluation, the following experimental and computational phases are established:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        EMPIRICAL DISCRIMINATION PHASES                                 │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                        │
│   PHASE A (Computational / Milestone B2 & B3):                                         │
│   ├─ Independent structural recalculation engine (src/telotrade/independent/)          │
│   ├─ Parameter sensitivity sweep (cutoffs 3.0-5.0 Å; Bondi vs ProtOr; probe 1.2-1.6 Å)  │
│   ├─ Domain excursion vector alignment against TPP1 center-of-mass                     │
│   └─ Lattice contact analysis of 2I46 crystal symmetry mates                           │
│                                                                                        │
│   PHASE B (Allostery & Ligand Audit / Milestone B3):                                   │
│   ├─ DE-loop B-factor vs coordinate displacement correlation analysis                  │
│   ├─ Dynamic allostery evaluation (decoupling 0.338 Å RMSD from entropic coupling)     │
│   └─ TELIN-1 reality check (proving zero experimental binding data / pure pose artifact)│
│                                                                                        │
│   PHASE C (Translational Modeling & Wet-Lab Design / Milestone B4):                    │
│   ├─ Cancer vs aging evidence matrix (ClinVar DC mutations vs cancer dependencies)     │
│   ├─ Stopped-flow / SPR kinetic protocol design (falsifying 1-step vs 2-step binding)  │
│   └─ HDX-MS and cryo-EM apo-TPP1 validation assay protocols                            │
│                                                                                        │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### Phase A: In Silico Parameter & Uncertainty Stress-Testing (Milestones B2 & B3)
1. **Independent Metric Recalculation (`independent_metrics.csv`)**:
   - Re-implement SVD Kabsch alignment, Shrake-Rupley BSA, and contact distance algorithms using independent code paths (e.g. Biopython and pure NumPy vs Gemmi).
   - Test whether discrepancies exceed $\pm 1\%$ under identical parameter configurations.
2. **Multi-Cutoff Sensitivity Analysis (`sensitivity_analysis.csv`)**:
   - Sweep heavy-atom cutoffs across a fine grid (3.0 Å to 5.0 Å in 0.1 Å increments).
   - Evaluate whether the 5 Core Robust residue pairs dominate total contact surface area across the entire range, supporting H1/H2 over H3.
3. **Crystallographic Packing Audit of 2I46**:
   - Calculate inter-molecular distances between the DE-loop in Chain A and all crystallographic symmetry mates generated from space group $P2_12_12_1$.
   - If symmetry mate atoms are found within $< 3.5\ \text{Å}$ of Glu169–Lys170, lattice artifact contribution (H1) is quantitatively substantiated.

### Phase B: Mechanistic & Allostery Rigor Assessment (Milestone B3)
1. **DE-Loop Remodeling vs B-Factor Analysis (`de_loop_mechanism_report.md`)**:
   - Statistically evaluate whether Lys170 peak displacement (4.293 Å) is accompanied by anomalous B-factor inflation.
   - Deconstruct any claim of sequential "Leu183 first, DE-loop second" binding, reclassifying it as an unvalidated kinetic hypothesis.
2. **Allosteric Coupling Rigor Analysis (`allostery_evidence_assessment.md`)**:
   - Formalize the biophysical decoupling of static coordinate invariance (0.338 Å) from dynamic/energetic coupling.
   - Document how entropic allostery and DNA exit path stabilization by TPP1-POT1 stimulate repeat addition processivity (RAP).
3. **TELIN-1 Candidate Ligand Reality Check (`TELIN1_evidence_audit.md`)**:
   - Confirm complete absence of wet-lab provenance (0 literature citations, 0 patent filings).
   - Demonstrate that TELIN-1's predicted affinity is a parameter-sensitive docking artifact (supporting H3).

### Phase C: Translational & Wet-Lab Experimental Protocols (Milestone B4)
1. **Cancer–Aging Evidence Matrix (`cancer_aging_evidence_matrix.csv`)**:
   - Synthesize human genetics data (germline TEL-patch mutations causing Dyskeratosis Congenita, IPF, bone marrow failure) against cancer dependency datasets (DepMap CRISPR Achilles screens).
   - Establish whether any therapeutic index separates complete cancer arrest from hematologic toxicity.
2. **Falsifiable Wet-Lab Validation Protocols (`falsifiable_predictions.md`)**:
   - Specify concrete experimental designs (stopped-flow kinetic assays, HDX-MS, Cryo-EM of isolated full-length TPP1, and site-directed mutagenesis of core vs peripheral contact pairs) to decisively arbitrate between H1, H2, and H3 in future empirical investigations.

---

## 6. Traceability to Evidence Ledger & Prior-Art Audits

The hypotheses and criteria defined herein directly map to the audited findings recorded in `results/audit/evidence_ledger.csv` and `results/literature/novelty_gap_report.md`:

| Hypothesis Focus | Mapped Evidence Ledger Claims | Literature Prior-Art Status | Milestone Target |
| :--- | :--- | :--- | :--- |
| **H1 (Core Invariance & Lattice Artifact)** | CLM-010 (0.338 Å RMSD, `uncertain`), CLM-016 (1.966 Å OB-fold RMSD, `established`), CLM-017 (4.293 Å Lys170 displacement, `uncertain`) | Established baseline (Liu et al. 2022; Wang et al. 2007) | Milestone B2 / B3 |
| **H2 (Induced-Fit Remodeling & Clamping)** | CLM-013 (0.434 Å TEN displacement, `established`), CLM-015 (0.427 Å TRAP displacement, `established`), CLM-021 (Sequential binding, `untested`), CLM-022 (Zero allostery claim, `contradicted`) | Dynamic coupling established in MD literature (*Comms Chem* 2023); sequential kinetics unproven | Milestone B3 |
| **H3 (Parameter & Threshold Sensitivity)** | CLM-001 (856.82 Å² BSA, `established`), CLM-002 (9 contacts at 3.5 Å, `established`), CLM-003 (46 contacts at 4.0 Å, `established`), CLM-004 (123 contacts at 4.5 Å, `established`), CLM-005 (2.32 Å compaction, `uncertain`), CLM-023 (TELIN-1 ligand, `untested`) | Multi-cutoff robustness is authentic novel contribution of this study; TELIN-1 quarantined as unvalidated | Milestone B2 / B3 |

---

## 7. Conclusion

By establishing explicit, mathematically rigorous criteria for the support and falsification of **H1 (Rigid Conservation)**, **H2 (Induced-Fit Plasticity)**, and **H3 (Parameter Sensitivity)**, this research framework moves the analysis of human telomerase recruitment from speculative narrative interpretation to rigorous, falsifiable computational biophysics. 

All subsequent deliverables in Benchmark Milestones B2, B3, B4, and B5 shall be evaluated against the operational criteria defined in this document.
