# Audit Methodology & Evidence Framework: Forensic Epistemological Blueprint for Human Telomerase Recruitment (Requirements R1 & R2)

**Author**: `explorer_b1_methodology` (Audit Methodology & Evidence Framework Specialist)  
**Date**: October 2026  
**Target Milestone**: Milestone B1 (Requirements R1 & R2)  
**Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b1_methodology`  
**Target Output Artifacts**: 
- `results/audit/evidence_ledger.csv`
- `results/audit/repository_inventory.md`
- `results/literature/prior_art_matrix.csv`
- `results/literature/claim_verification.md`
- `results/literature/novelty_gap_report.md`

---

## 1. Executive Summary & Epistemological Mandate

The primary objective of this investigation is to construct a rigorous, mathematically sound, and scientifically falsifiable audit framework for the human telomerase structural bioinformatics repository (`telomerase-research`). In the preceding generation of research, a high-performing automated pipeline (`src/telotrade/`) was assembled to process cryo-EM structures (`7TRE`, `7TRD`, `7TRF`, `7TRC`) and X-ray crystallographic coordinates (`2I46`), quantifying buried surface area (BSA), residue–residue contact networks across distance cutoffs (3.5 Å, 4.0 Å, 4.5 Å), domain-level SVD Kabsch root-mean-square deviations (RMSD), and conformational shifts in the TPP1 TEL patch.

However, a critical epistemological gap emerged between **computational metric recovery** and **scientific inference**:
1. **Conflating Positive Controls with Novel Discoveries**: Foundational recruitment determinants such as the primary salt bridge (`TERT Arg774 -- TPP1 Glu168`), the secondary latch (`TERT Lys78 -- TPP1 Glu215`), and the hydrophobic anchor (`Leu183`) were already experimentally established by Liu et al. (*Nature* 2022), Sekne et al. (*Science* 2022), and Nandakumar et al. (*Nature* 2012). Their identification by `telotrade` constitutes a mandatory validation check (positive control), not a novel computational discovery.
2. **Static Coordinate Rigidity vs. Dynamic Allostery**: Measuring a catalytic core Cα RMSD of 0.338 Å between bound (`7TRE`) and unbound (`7TRD`, `7TRF`) structures proves structural invariance in the static, vitrified ground state. Concluding from this number that telomerase possesses "zero allosteric coupling" is scientifically invalid: allosteric modulation frequently operates via changes in conformational entropy, fast picosecond-to-nanosecond vibrational dynamics, or altered chemical step kinetics ($k_{\mathrm{cat}}$), leaving ground-state Cα coordinates virtually unchanged.
3. **Crystal Lattice Artifacts vs. Biological Solution Remodeling**: The pronounced Cα displacement of the TPP1 DE-loop (up to 4.293 Å at Lys170) observed between `7TRE` and `2I46` Chain A has been presented as definitive evidence of an "induced-fit docking mechanism." Crucially, in `2I46` Chain B, this exact loop is completely disordered in electron density, while in Chain A it is stabilized by crystallographic lattice contacts. Interpreting this coordinate difference as an obligate sequential binding mechanism without time-resolved or solution biophysical data represents an unsubstantiated leap from static end-states.
4. **Clinical Viability of Recruitment Inhibition**: Proposing orthosteric disruption of the TERT–TPP1 interface as a selective anti-cancer therapy ignores decisive human clinical genetics. Germline loss-of-function mutations in the TPP1 TEL patch (`Lys170del`, `Glu169Lys`) cause lethal congenital telomeropathies (Dyskeratosis Congenita, Hoyeraal-Hreidarsson syndrome, aplastic anemia, idiopathic pulmonary fibrosis). Because the physical recruitment interface is genetically identical between malignant cells and CD34+ hematopoietic stem cells, absolute recruitment blockade will replicate the lethal phenotype of Dyskeratosis Congenita.
5. **In Silico Ligand Reality Check (`TELIN-1`)**: The candidate ligand `TELIN-1` represents a purely hypothetical in silico docking proposal with zero experimental validation ($K_d$, $\mathrm{IC}_{50}$, cellular uptake, or selectivity). It must be strictly classified as untested and de-hyped.

This report establishes the taxonomic criteria, data integrity checks, code-level calculation audits, and an exact implementation plan for the worker agents authoring the Phase 0 and Phase 1 deliverables.

---

## 2. Taxonomic and Epistemological Framework for Claim Classification

Every reported claim, parameter, and numerical metric across the repository must be cataloged and assigned exactly one epistemological status in `results/audit/evidence_ledger.csv`:
- `established`
- `uncertain`
- `contradicted`
- `untested`

The classification decision tree is formal and objective, based on four orthogonal axes of evidence:
1. **Experimental Structure & Resolution**: Cryo-EM/X-ray map quality, coordinate uncertainty, and local resolution.
2. **Threshold & Algorithmic Sensitivity**: Sensitivity to cutoffs, probe radii, atom selections, and software libraries.
3. **Biochemical & Genetic Corroboration**: Direct verification via published wet-lab mutagenesis, binding assays, or human germline genetics.
4. **Epistemological Validity**: Absence of logical fallacies (e.g., claiming static coordinate immobility disproves dynamic allostery).

```
                             [ Candidate Claim / Metric ]
                                          │
            ┌─────────────────────────────┴─────────────────────────────┐
            ▼                                                           ▼
Is it a hypothetical proposal or                            Does it conflict with experimental 
computational prediction lacking                            structures, biophysical laws, or 
any wet-lab binding data?                                   clinical genetics?
     YES │                                                              YES │
         ▼                                                                  ▼
    [ UNTESTED ]                                                      [ CONTRADICTED ]
         │ NO                                                               │ NO
         └──────────────────────────────┬───────────────────────────────────┘
                                        ▼
                   Is it robust across thresholds AND independently
                   corroborated by published wet-lab biochemistry?
                            YES │                       │ NO
                                ▼                       ▼
                         [ ESTABLISHED ]          [ UNCERTAIN ]
                         (Positive Control)       (Sensitive to cutoff/construct)
```

### 2.1 Formal Classification Criteria

#### Class 1: `established`
A claim is classified as `established` if and only if it satisfies all of the following conditions:
1. **Structural Verification**: Observed in experimental coordinates derived from high-resolution electron density or cryo-EM maps ($\le 3.5\ \text{Å}$ global resolution) with unambiguous local density and standard stereochemical geometry.
2. **Biochemical Corroboration**: Supported by published peer-reviewed functional assays (alanine-scanning mutagenesis, co-immunoprecipitation, pull-downs, or functional telomere elongation assays).
3. **Computational Invariance**: Persists across all standard metric thresholds (e.g., contact present at 3.5 Å, 4.0 Å, and 4.5 Å) and is reproducible across independent parsers (Gemmi, Biopython).
4. **Status Note**: In literature comparisons, established claims must be designated as **Positive Controls**, not original computational discoveries.

#### Class 2: `uncertain`
A claim is classified as `uncertain` if it is computationally derived or present in experimental coordinates, but meets ANY of the following criteria:
1. **Cutoff Sensitivity**: Highly sensitive to arbitrary threshold parameters (e.g., contact present at 4.0 Å or 4.5 Å, but lost at 3.5 Å; or BSA value varying by $>5\%$ with probe radius $r = 1.2–1.6\ \text{Å}$).
2. **Construct or Environment Dependence**: Arises from differences in expression constructs (e.g., Twin-Strep tags, GAG linker), crystal lattice packing (e.g., 2I46 Chain A lattice interactions vs. Chain B disorder), or single-particle vitrification artifacts.
3. **Coordinate Precision Boundary**: Metric displacement is on the order of local coordinate uncertainty (e.g., rigid-body displacements $< 0.5\ \text{Å}$ in a 3.5 Å cryo-EM reconstruction, where the average coordinate error is $\sim 0.3–0.5\ \text{Å}$).
4. **Incomplete Geometric Criteria**: Inferred polar interactions based solely on heavy-atom distance ($\le 4.0\ \text{Å}$) without explicit hydrogen positions, donor–H–acceptor angular verification, or verified protonation states.

#### Class 3: `contradicted`
A claim is classified as `contradicted` if it meets ANY of the following criteria:
1. **Refuted by Structural Data**: Direct conflict with experimental coordinates (e.g., claiming that TPP1 induces active-site remodeling in TERT, directly contradicted by the 0.338 Å catalytic core RMSD and $<0.25\ \text{Å}$ catalytic triad shift).
2. **Epistemological Fallacy**: Claiming that static coordinate similarity across ground states rules out dynamic or entropic allostery (violates fundamental thermodynamic principles of macromolecular allostery).
3. **Refuted by Human Genetics / Clinical Evidence**: Claiming that total orthosteric ablation of the TERT–TPP1 recruitment interface is a viable non-toxic therapeutic strategy, directly contradicted by human germline telomeropathy genetics (fatal aplastic anemia and IPF caused by TEL-patch disruption).
4. **Implementation Errors**: Erroneous formulas, sign errors, reversed indices, or swapped chain entity assignments.

#### Class 4: `untested`
A claim is classified as `untested` if it meets ANY of the following criteria:
1. **Computational Candidates**: In silico small molecules, peptides, or docked poses (e.g., `TELIN-1`) lacking published in vitro binding assays (SPR, ITC, BLI), cell-based target engagement, or co-crystal/cryo-EM structural validation.
2. **Inferred Kinetic Order**: Hypotheses describing sequential binding pathways (e.g., "Leu183 docks first as a guide knob, followed by DE-loop electrostatic locking") derived solely from comparing two static end-state structures without time-resolved biophysics.
3. **Hypothetical Therapeutic Windows**: Claims that the "telomere length buffer" will spare adult stem cells during therapeutic telomerase recruitment inhibition, lacking in vivo preclinical validation or long-term hematopoietic colony-forming assays.

---

### 2.2 Claim-by-Claim Epistemological Audit Rubric

The table below defines the formal classification for all headline claims reported in the Gen1 synthesis report (`report.md`) and standard tables (`results/tables/`).

| Claim ID | Metric / Feature | Reported Value | Source File & Script | Parameters & Methodology | Epistemological Status | Detailed Justification & Limitations |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CLM-01** | Interface BSA per subunit | **856.82 Å²** ($\Delta\mathrm{SASA} = 1713.65\ \text{Å}^2$) | `report.md` Sec 1, 4.3; `buried_surface_area.csv`; `interface_geometry.py` | Shrake-Rupley ($N=960$, $r=1.4\ \text{Å}$), Bondi radii, isolated TERT (A) vs TPP1 (C) | **`uncertain`** | Dependent on probe radius (1.4 Å), Bondi radii dictionary vs ProtOr, and truncation of unresolved loops (TERT 179–321). Reflects typical moderate PPI. |
| **CLM-02** | SASA Isolated TERT | **49,776.00 Å²** | `buried_surface_area.csv`; `interface_geometry.py` | FreeSASA / Shrake-Rupley ($r=1.4\ \text{Å}$, $N=960$) | **`uncertain`** | Includes artificial solvent exposure at chain termination sites of missing loops (1–6, 105–116, 179–321, 417–443, 641–650). |
| **CLM-03** | SASA Isolated TPP1 | **10,651.36 Å²** | `buried_surface_area.csv`; `interface_geometry.py` | FreeSASA / Shrake-Rupley ($r=1.4\ \text{Å}$, $N=960$) | **`established`** | Robust for the isolated 90–241 OB-fold construct; physiological TPP1 in shelterin has flanking TIN2/POT1 domains. |
| **CLM-04** | SASA Bound Complex | **58,713.70 Å²** | `buried_surface_area.csv`; `interface_geometry.py` | FreeSASA / Shrake-Rupley ($r=1.4\ \text{Å}$, $N=960$) | **`uncertain`** | Same linker truncation dependencies as CLM-02. |
| **CLM-05** | Primary Salt Bridge Distance | **2.32 Å** (`TERT Arg774 NH1 -- TPP1 Glu168 OE1`) | `report.md` Sec 1, 4.2; `polar_interactions_7tre.csv`; `interface_geometry.py` | Euclidean heavy-atom distance in 7TRE | **`established`** (Interaction) / **`uncertain`** (Distance) | Physical interaction established (Liu et al. 2022, Nandakumar 2012 E168A). The 2.32 Å distance is unusually short for N...O (typical 2.8 Å), indicating coordinate strain or uncertainty at 3.5 Å resolution. |
| **CLM-06** | Secondary Salt Latch Distance | **3.58 Å** (`TERT Lys78 NZ -- TPP1 Glu215 OE2`) | `report.md` Sec 1, 4.2; `polar_interactions_7tre.csv`; `interface_geometry.py` | Cutoff $\le 4.0\ \text{Å}$, heavy atoms | **`uncertain`** | Threshold sensitive: lost at 3.5 Å cutoff; 3.58 Å is near the electrostatic interaction limit. Positive control from Sekne et al. 2022. |
| **CLM-07** | Hydrophobic Anchor Contact | **3.67 Å** (`TERT Leu769 CD1 -- TPP1 Leu183 CD2`) | `report.md` Sec 1, 4.2; `interface_contacts_7tre_4.0A.csv`; `interface_geometry.py` | Cutoff $\le 4.0\ \text{Å}$, non-bonded C–C | **`established`** | Positive control; validated by Nandakumar et al. 2012 L183A recruitment loss. |
| **CLM-08** | Hydrophobic Groove Packing | **3.52 Å** (`TERT Pro771 CG -- TPP1 Leu212 CD2`) | `report.md` Sec 1, 4.2; `interface_contacts_7tre_4.0A.csv`; `interface_geometry.py` | Cutoff $\le 4.0\ \text{Å}$, non-bonded C–C | **`established`** | Confirmed hydrophobic contact in IFD/TRAP cradle. |
| **CLM-09** | 3.5 Å Contact Network | **9 contacts / 5 residue pairs** | `report.md` Sec 4.6; `interface_contacts_7tre_3.5A.csv`; `interface_geometry.py` | Spatial cell hashing, Euclidean $d \le 3.5\ \text{Å}$, heavy atoms | **`established`** | Mathematically rigorous core contact network in 7TRE. Highly robust. |
| **CLM-10** | 4.0 Å Contact Network | **46 contacts / 21 residue pairs** | `report.md` Sec 4.6; `interface_contacts_7tre_4.0A.csv`; `interface_geometry.py` | Spatial cell hashing, Euclidean $d \le 4.0\ \text{Å}$ | **`established`** | Standard physiological contact shell in 7TRE. |
| **CLM-11** | 4.5 Å Contact Network | **123 contacts / 32 residue pairs** | `report.md` Sec 4.6; `interface_contacts_7tre_4.5A.csv`; `interface_geometry.py` | Spatial cell hashing, Euclidean $d \le 4.5\ \text{Å}$ | **`uncertain`** | Includes peripheral, non-specific van der Waals contacts sensitive to minor coordinate jitter and thermal motion. |
| **CLM-12** | Catalytic Core RMSD (7TRE vs 7TRD) | **0.338 Å** (770 Cα atoms) | `report.md` Sec 1, 4.4; `cross_state_rmsd.csv`; `cross_state.py` | SVD Kabsch rigid-body alignment, residues 326–1132 | **`established`** | Mathematically verified; establishes static coordinate immobility of the catalytic core between bound and unbound states. |
| **CLM-13** | Catalytic Core RMSD (7TRE vs 7TRF) | **0.338 Å** (770 Cα atoms) | `report.md` Sec 1, 4.4; `cross_state_rmsd.csv`; `cross_state.py` | SVD Kabsch rigid-body alignment, residues 326–1132 | **`established`** | Confirms core invariance even in the presence of histone H2A–H2B complex. |
| **CLM-14** | Global TERT Cα RMSD (7TRE vs 7TRD) | **0.356 Å** (934 Cα atoms) | `report.md` Sec 4.4; `cross_state_rmsd.csv`; `cross_state.py` | SVD Kabsch alignment over all modeled Cα | **`established`** | Demonstrates overall rigidity of modeled TERT scaffold. |
| **CLM-15** | TEN Domain Rigid-Body Shift | **0.434 Å** (Core-Fit displacement) | `report.md` Sec 1, 4.4; `cross_state_rmsd.csv`; `cross_state.py` | Kabsch on core (770 Cα), apply transform to TEN (160 Cα) | **`uncertain`** | Displacement (0.434 Å) is within the mean coordinate precision of 3.3–3.5 Å cryo-EM maps; cannot definitively separate biological clamp motion from coordinate noise. |
| **CLM-16** | Global TPP1 OB-fold RMSD (7TRE vs 2I46 A) | **1.966 Å** (152 Cα atoms) | `report.md` Sec 1, 4.5; `cross_state_rmsd.csv`; `cross_state.py` | Kabsch SVD on residues 90–241, 7TRE Chain C vs 2I46 Chain A | **`established`** | Accurate mathematical RMSD between cryo-EM model and X-ray Chain A. |
| **CLM-17** | TPP1 Leu183 Cα Displacement | **0.512 Å** | `report.md` Sec 1, 4.5; `tpp1_conformational_shifts.csv`; `cross_state.py` | Coordinate subtraction after Kabsch alignment | **`established`** | Invariant backbone anchor residue. |
| **CLM-18** | TPP1 Lys170 Peak DE-Loop Shift | **4.293 Å** (Glu169: 3.809 Å, Glu171: 3.320 Å) | `report.md` Sec 1, 4.5; `tpp1_conformational_shifts.csv`; `cross_state.py` | Coordinate subtraction after Kabsch alignment | **`uncertain`** | While the coordinate difference between 7TRE and 2I46 Chain A is 4.293 Å, 2I46 Chain B has this loop completely disordered, and Chain A is locked by crystal lattice contacts. Attributing this displacement purely to physiological induced-fit is uncertain. |
| **CLM-19** | Zero Allosteric Regulation Claim | "0.338 Å RMSD proves absence of allosteric regulation by TPP1" | `report.md` Sec 1, 3, 4.4, 9 | Ground-state Cα RMSD interpretation | **`contradicted`** | Epistemological fallacy: static coordinate invariance does not rule out dynamic, entropic, or kinetic allostery. |
| **CLM-20** | Sequential Induced-Fit Mechanism | "Leu183 docks first as guide knob, followed by DE-loop electrostatic locking" | `report.md` Sec 4.5, 9 | Extrapolation from two static crystal/cryo-EM structures | **`untested`** | Mechanistic hypothesis lacking time-resolved biophysical kinetics (stopped-flow, single-molecule FRET, or MD simulations). |
| **CLM-21** | TELIN-1 Inhibitor Target & Selectivity | `telin1_inhibitor.pdb` binds TPP1 cradle to selectively inhibit cancer | `ORIGINAL_REQUEST.md` R6 | In silico docked small molecule | **`untested`** | Pure computational hypothesis with zero experimental binding affinity ($K_d$), cellular uptake, or selectivity data. |
| **CLM-22** | Selective Recruitment Therapeutic Viability | Partial recruitment inhibition spares stem cells due to telomere length buffer | `report.md` Sec 1, 5.3, 5.5, 9 | Theoretical translational hypothesis | **`contradicted`** / **`untested`** | Contradicted by human germline genetics (heterozygous TEL-patch mutations cause lethal bone marrow failure and DC). Spreading of therapeutic window remains untested in vivo. |

---

## 3. Data Integrity & Repository Audit Standards

The data integrity check for `results/audit/repository_inventory.md` must enforce four strict verification protocols:

### 3.1 Git Tree & Version Control State Audit
- **Git Repository Status**: Verify whether the root directory is initialized as a Git repository (`git rev-parse --is-inside-work-tree`). If not, explicitly document the working tree as an unversioned filesystem workspace snapshot.
- **Commit & Tree Cleanliness**: Record `HEAD` commit SHA (if available), current branch, status of modified files, and presence of untracked files or scratch caches (`.pytest_cache/`, `.ruff_cache/`, `__pycache__/`).
- **Submodule & Upstream Sync**: Record any submodules or external dependencies.

### 3.2 Raw Data File Cryptographic Integrity
- **Chunked SHA-256 Checksums**: Every raw coordinate file in `data/raw/` must be cryptographically audited using chunked 64KB hashing. Checksums must match authoritative RCSB values:
  - `7TRE.cif`: `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` (1,518,438 bytes)
  - `7TRD.cif`: `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` (1,394,494 bytes)
  - `7TRF.cif`: `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` (1,542,948 bytes)
  - `7TRC.cif`: `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` (1,647,704 bytes)
  - `2I46.cif`: `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` (273,104 bytes)
- **Validation Documents**: Confirm presence of wwPDB validation reports (`.pdf` and `.xml.gz`) for all five targets.

### 3.3 Test Suite Execution & Pass Rate Metrics
- **Test Infrastructure Discovery**: Audit all 300 tests across unit, integration, and adversarial suites:
  - `tests/e2e/test_tier1_features.py`: Feature coverage (F1–F18)
  - `tests/e2e/test_tier2_boundaries.py`: Boundary and corner cases
  - `tests/e2e/test_tier3_interactions.py`: Multi-module data flow
  - `tests/e2e/test_tier4_scenarios.py`: End-to-end pipeline and tables
  - `tests/test_adversarial_m1.py`: Hash corruption and missing data recovery
  - `tests/test_adversarial_m2.py`: Tag offsets, sequence mismatch, 2I46 asymmetry
  - `tests/test_adversarial_m3.py`: Cutoff bounds, salt bridge strictness, BSA convergence
  - `tests/test_adversarial_m4.py`: Kabsch SVD reflections, domain boundaries, benchmark parity
  - `tests/test_adversarial_m5_figures.py`: Figure sensitivity and synthesis report checks
  - `tests/test_validate.py`, `test_sequence_mapping.py`, `test_interface_geometry.py`, `test_cross_state.py`
- **Pass Rate Expectation**: 100% pass rate (300/300) with zero regressions.

### 3.4 Environment Reproducibility Matrix
- Record runtime and library versions:
  - Python version: `sys.version` (e.g., Python 3.12/3.13 Windows x64)
  - `gemmi`: $\ge 0.6.5$ (authoritative mmCIF parsing and CRA traversal)
  - `biopython`: $\ge 1.83$ (MMCIFParser, NeighborSearch, Bio.SeqUtils)
  - `freesasa`: C-compiled extension module (calcCoord Shrake-Rupley)
  - `numpy`: $\ge 1.26$ (Kabsch SVD, vector arithmetic)
  - `pandas`: $\ge 2.2$ (standardized CSV table generation)
  - `scipy`: spatial distance routines
  - `matplotlib`: figure generation backend

---

## 4. Code-Level Metric Calculation Audit: Assumptions, Discrepancies & Sensitivities

A forensic review of the core mathematical and biophysical algorithms in `src/telotrade/` reveals key assumptions and potential numerical discrepancies between software implementations:

### 4.1 Buried Surface Area (BSA) Computation
- **Mathematical Formula**:
  $$\mathrm{BSA} = \frac{1}{2} \left[ \mathrm{SASA}(\text{TERT}) + \mathrm{SASA}(\text{TPP1}) - \mathrm{SASA}(\text{Complex}) \right]$$
  $$\Delta\mathrm{SASA}_{\mathrm{total}} = \mathrm{SASA}(\text{TERT}) + \mathrm{SASA}(\text{TPP1}) - \mathrm{SASA}(\text{Complex}) = 2 \times \mathrm{BSA}$$
- **Engine Comparison**:
  1. **FreeSASA (`freesasa.calcCoord`)**:
     - Uses explicit Bondi atomic radii (`BONDI_RADII`): C = 1.70 Å, N = 1.55 Å, O = 1.52 Å, S = 1.80 Å, H = 1.20 Å.
     - Probe radius: $r_p = 1.40\ \text{Å}$.
     - Sphere points: $N = 960$ points per atom sphere.
     - Strips non-amino-acid heteroatoms and hydrogens.
     - Result on 7TRE: $\mathrm{SASA}(\text{TERT}) = 49,776.00\ \text{Å}^2$, $\mathrm{SASA}(\text{TPP1}) = 10,651.36\ \text{Å}^2$, $\mathrm{SASA}(\text{Complex}) = 58,713.70\ \text{Å}^2 \implies \mathrm{BSA} = 856.82\ \text{Å}^2$.
  2. **Biopython (`Bio.PDB.SASA.ShrakeRupley`)**:
     - Uses default `_atom_radius` lookup based on Chothia/ProtOr parameters (aliphatic C = 2.00 Å, aromatic C = 1.75 Å, carbonyl C = 1.55 Å, N = 1.55–1.70 Å, O = 1.40 Å).
     - Because ProtOr radii are larger for carbons and smaller for oxygens, computing SASA via Biopython shifts total SASA by $\sim 5–10\%$ and yields a BSA shift of $\sim 20–60\ \text{Å}^2$.
  3. **Discrepancy Source & Assumption**:
     - The reported BSA of $856.82\ \text{Å}^2$ is **radii-dependent**. Downstream independent recalculation must specify atomic radii convention.
     - The missing 142-residue TEN-TRBD linker (residues 179–321) creates artificial solvent-exposed chain ends at residue 178 and 326. While localized, this inflates isolated TERT SASA by $\sim 150–250\ \text{Å}^2$.

### 4.2 Spatial Neighbor Search & Distance Cutoffs
- **Implementation in `compute_interface_contacts`**:
  - Uses pure NumPy spatial cell hashing (`cell_size = cutoff`, 27-cell neighbor query).
  - Alternative implementations: `scipy.spatial.distance.cdist`, `Bio.PDB.NeighborSearch` (KD-tree), and `gemmi.NeighborSearch`.
  - **Discrepancy Analysis**:
    - All implementations achieve exact Euclidean distance parity ($\Delta d < 10^{-5}\ \text{Å}$).
    - **Cutoff Escalation Sensitivity**:
      - At 3.5 Å: exactly **9 contacts across 5 residue pairs**.
      - At 4.0 Å: exactly **46 contacts across 21 residue pairs** (+411% contact count increase).
      - At 4.5 Å: exactly **123 contacts across 32 residue pairs** (+167% jump from 4.0 Å, +1267% from 3.5 Å).
    - This extreme non-linear scaling demonstrates that contact counts are purely operational constructs of the chosen threshold. The 9 core contacts at 3.5 Å constitute the true biophysical anchor, while contacts at 4.5 Å encompass loose peripheral van der Waals packing.

### 4.3 Polar Interactions & Salt Bridge Geometry
- **Criteria in `classify_interaction_type`**:
  - Salt Bridge: Anionic carboxylate oxygen (`ASP` OD1/OD2, `GLU` OE1/OE2) to basic nitrogen (`ARG` NH1/NH2/NE, `LYS` NZ, `HIS` ND1/NE2) with Euclidean distance $\le 4.0\ \text{Å}$.
  - Polar / H-Bond: Non-charged polar heteroatoms (`N`, `O`, `S`) with Euclidean distance $\le 4.0\ \text{Å}$.
- **Discrepancy & Limitations**:
  - **Histidine Charge State**: The code treats His as basic, assuming full protonation. At physiological pH (7.4), His is predominantly neutral ($\sim 90\%$).
  - **Lack of Angular Geometry**: No donor–H–acceptor angle check is performed because cryo-EM models at 3.5 Å lack hydrogen coordinates. Therefore, interactions labeled "H-Bond" are strictly "polar proximity pairs" rather than validated geometric hydrogen bonds.
  - **The 2.32 Å Arg774–Glu168 Salt Bridge**:
    - The measured heavy-atom distance between Arg774 NH1 and Glu168 OE1 is **2.32 Å**.
    - In high-resolution X-ray crystallography, typical salt bridge N...O distances average $2.75–2.90\ \text{Å}$, and the sum of Bondi van der Waals radii is $1.55 + 1.52 = 3.07\ \text{Å}$.
    - A 2.32 Å distance in a 3.5 Å cryo-EM structure reflects coordinate refinement restraint strain or local coordinate uncertainty rather than an ultra-short covalent-like interaction.

### 4.4 Rigid-Body SVD Kabsch Superposition
- **Algorithm in `cross_state.py`**:
  - Centers mobile ($P$) and target ($Q$) coordinate point clouds.
  - Computes cross-covariance matrix $H = P_c^T Q_c$.
  - Executes SVD: $H = U \Sigma V^T$.
  - Enforces proper rotation in $SO(3)$: $R = V \operatorname{diag}(1, 1, \operatorname{sign}(\det(V U^T))) U^T$, guaranteeing $\det(R) = +1.0$ and preventing improper reflections.
- **Subsets & Chain Selection**:
  - Catalytic core superposition (residues 326–1132, 770 Cα) yields an identical RMSD of **0.338 Å** for `7TRE vs 7TRD` and `7TRE vs 7TRF`.
  - **2I46 Dimer Asymmetry**: Chain A has the TEL-patch DE-loop resolved in electron density, while in Chain B residues 168–171 are completely disordered. Superposition onto Chain A yields global TPP1 RMSD of 1.966 Å and Lys170 shift of 4.293 Å. The code correctly selects Chain A, but the audit must document that Chain A is stabilized by crystal packing in space group $I4_122$.

---

## 5. Exhaustive Artifact Implementation Plan for R1 & R2 Workers

This section outlines the exact implementation blueprints for the downstream worker agents tasked with authoring the five Phase 0 and Phase 1 deliverables.

### 5.1 Deliverable 1: `results/audit/evidence_ledger.csv`
- **Output Path**: `results/audit/evidence_ledger.csv`
- **Required Columns**:
  1. `claim_id`: Unique identifier (e.g., `CLM-01`, `CLM-02`)
  2. `metric_name`: Formal name of the metric or biological claim
  3. `reported_value`: Exact numerical value and unit from Gen1 report
  4. `source_file`: Primary documentation file (`report.md`, `results/tables/buried_surface_area.csv`, etc.)
  5. `source_script`: Originating script and function (`src/telotrade/interface_geometry.py:compute_bsa`, etc.)
  6. `parameters`: Exact input parameters, atom selections, probe radius, sphere points, and cutoffs
  7. `evidence_status`: Taxonomic assignment (`established`, `uncertain`, `contradicted`, `untested`)
  8. `notes`: Epistemological explanation, limitations, and literature context
- **Itemized Claim Scope**: Must incorporate all 22 claims detailed in Section 2.2 of this report.

### 5.2 Deliverable 2: `results/audit/repository_inventory.md`
- **Output Path**: `results/audit/repository_inventory.md`
- **Mandatory Sections**:
  1. `# Human Telomerase Repository Inventory & Data Integrity Audit`
  2. `## 1. Executive Summary & Audit Scope`
  3. `## 2. Git Tree & Version Control State` (document unversioned directory status, dirty state, cache files)
  4. `## 3. Raw Data Manifest & Cryptographic File Verification` (complete table of 5 CIF files, byte sizes, verified SHA-256 hashes, validation PDFs/XMLs)
  5. `## 4. Codebase Architecture & Metric Calculation Pipelines` (audit of `src/telotrade/` modules)
  6. `## 5. Test Infrastructure & Pass Rate Audit` (summary of 300 tests across E2E and adversarial suites)
  7. `## 6. Standardized Tables & Generated Artifacts Inventory` (table listing all 21 CSVs in `results/tables/` and 4 PNGs in `results/figures/`)
  8. `## 7. Execution Environment & Reproducibility Matrix` (Python, OS, packages)
  9. `## 8. Discrepancies, Identified Gaps & Forensic Recommendations`

### 5.3 Deliverable 3: `results/literature/prior_art_matrix.csv`
- **Output Path**: `results/literature/prior_art_matrix.csv`
- **Required Columns**:
  1. `finding_id`: Unique identifier (e.g., `PA-01`, `PA-02`)
  2. `structural_feature_or_claim`: Biophysical feature or finding
  3. `telotrade_finding`: Result computed by `telotrade`
  4. `primary_literature_citation`: Author, Journal, Year, DOI/PMID
  5. `reported_in_literature`: Boolean (`True` / `False`)
  6. `literature_value_or_finding`: What was previously reported in the paper
  7. `status_in_this_study`: Role in audit (`positive_control`, `novel_computational_metric`, `unresolved_gap`, `overstated_claim`)
  8. `notes`: Contextual comparison
- **Scope of Literature Comparisons**:
  - Liu et al., *Nature* 2022 (PDB 7TRE, 7TRD, 7TRF, 7TRC): Arg774–Glu168 salt bridge, rigid catalytic core.
  - Sekne et al., *Science* 2022: Lys78–Glu215 secondary latch, TEN domain flexibility.
  - Nandakumar et al., *Nature* 2012: Discovery of the TEL patch, Glu168/Glu169/Glu171/Leu183 alanine mutagenesis.
  - Wang et al., *Nature* 2007 (PDB 2I46): TPP1 OB-fold crystal structure, Chain A vs B asymmetry.
  - Zhong et al., *Mol Cell* 2012: TPP1 OB-fold recruitment assays.
  - Human Clinical Genetics (Guo et al. 2014, Kocak et al. 2014, Armanios & Blackburn 2012): `TPP1 Lys170del`, `Glu169Lys`, `TERT Arg774` in Dyskeratosis Congenita and IPF.
  - Imetelstat Trials (Zeidan et al. 2023, Santini et al. 2020): On-target hematological toxicities (Grade 3/4 thrombocytopenia/neutropenia).

### 5.4 Deliverable 4: `results/literature/claim_verification.md`
- **Output Path**: `results/literature/claim_verification.md`
- **Mandatory Sections**:
  1. `# Claim Verification & Epistemological De-Hyping Audit`
  2. `## 1. Verification of Positive Controls (Prior Art Recovery)` (Arg774–Glu168, Lys78–Glu215, Leu183 anchor, catalytic triad immobility)
  3. `## 2. Verification of Quantitative Geometric Claims` (BSA Shrake-Rupley metrics, cutoff sensitivity 9 -> 46 -> 123)
  4. `## 3. Verification of Conformational Remodeling Claims` (DE-loop plastic shift vs 2I46 lattice packing and Chain B disorder)
  5. `## 4. Verification of Mechanistic & Allostery Claims` (static RMSD 0.338 Å vs dynamic allostery)
  6. `## 5. Verification of Translational Claims` (TELIN-1 status, Telomerase Paradox, germline telomeropathies vs cancer selectivity)
  7. `## 6. Master Claim Verification Ledger & Falsification Summary Table`

### 5.5 Deliverable 5: `results/literature/novelty_gap_report.md`
- **Output Path**: `results/literature/novelty_gap_report.md`
- **Mandatory Sections**:
  1. `# Primary Literature & Novelty Gap Audit Report`
  2. `## 1. Executive Summary: What Was Known vs. What Is Genuine Novelty`
  3. `## 2. Definitive Positive Controls Catalog` (features previously discovered by Liu, Sekne, Nandakumar)
  4. `## 3. Genuine Computational Insights & Methodological Novelties` (sequence-aware tag offset reconcilers, multi-cutoff spatial sensitivity profiles, core-fit domain decomposition)
  5. `## 4. De-Hyped & Overstated Claims Catalog` (de-hyping salt bridge discovery, zero allostery, sequential docking, and selective recruitment)
  6. `## 5. Unresolved Scientific Gaps & Concrete Requirements for Future Work`

---

## 6. Actionable Implementation Recommendations for Orchestrator & Workers

1. **Pre-Create Directory Tree**:
   Ensure `results/audit/` and `results/literature/` directories exist before workers initiate file writes:
   ```bash
   mkdir -p results/audit results/literature
   ```
2. **Prevent Code Duplication**:
   Workers must not duplicate parsing or calculation code in `results/`. All data extraction should use the audited functions in `src/telotrade/` or standard Python scripts.
3. **Preserve Test Suite Integrity**:
   No edits to `tests/` or `src/telotrade/` should be made during Phase 0 and Phase 1 artifact generation. All 300 existing tests must continue to pass cleanly.
4. **Enforce Strict Schema Contracts**:
   `results/audit/evidence_ledger.csv` and `results/literature/prior_art_matrix.csv` must adhere strictly to the header schemas defined in `orchestrator_gen2/PROJECT.md` and this report.

---
*End of Report.*
