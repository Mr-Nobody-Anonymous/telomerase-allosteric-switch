# Repository Inventory & Baseline Forensic Audit

**Document Designation**: R1 / Milestone B1 Deliverable  
**Investigation Lead**: `explorer_b1_repo` (Codebase & Artifact Auditor)  
**Date**: October 2026  
**Repository Path**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research`  
**Execution Context**: Accelerated Single-Run Pipeline (`telotrade`)  

---

## 1. Executive Summary & Inventory Scope

This repository inventory establishes a comprehensive, immutable forensic baseline for the computational structural bioinformatics codebase `telotrade`, designed to investigate human telomerase reverse transcriptase (TERT) recruitment by shelterin subunit TPP1 across comparative cryo-EM and X-ray crystallographic datasets (PDB IDs: `7TRE`, `7TRD`, `7TRF`, `7TRC`, and `2I46`).

This audit inspects every source script, unit/adversarial/E2E test, raw macromolecular data file, standardized output table in `results/tables/`, publication figure in `results/figures/`, and headline claim in `report.md`. Every reported numerical claim is traced to its exact source code, computational parameters, mathematical formulation, and analytical boundaries.

### Key Audit Findings:
1. **Cryptographic Data Integrity**: 100% of raw coordinate mmCIF files in `data/raw/` match the authoritative SHA-256 hashes recorded in `data_manifest.csv`.
2. **Environment & Git Status**: The repository workspace is currently uninitialized under Git (no `.git` directory present). A total of 8 root files, 8 Python source files in `src/telotrade/`, 12 test modules in `tests/`, 21 CSV tables in `results/tables/`, and 4 high-resolution figures in `results/figures/` were inventoried.
3. **Execution Baseline**: The test infrastructure consists of 300 automated test cases across unit, adversarial stress, and E2E tiers.
4. **Audited Code Anomaly (Defect Identified)**: In `src/telotrade/sequence_mapping.py::generate_tpp1_correspondence_table()`, an off-by-one indexing error (`ob_seq[resnum - 89]`) shifts the canonical amino acid sequence (`resname_canonical`) in `results/tables/residue_correspondence_tpp1.csv` by 1 position (e.g. TPP1 residue 166 is mislabeled as `TRP` instead of canonical `ASP`). Interface geometry and RMSD modules read directly from coordinate records and are unaffected.
5. **Physical Coordinate Anomaly**: The primary salt bridge distance (`TERT Arg774:NH1 -- TPP1 Glu168:OE1`) is measured at **2.32 Å** in `7TRE.cif`. While mathematically verified from deposited coordinates, a 2.32 Å heavy-atom distance indicates severe unphysical coordinate compaction/clash in the 3.50 Å cryo-EM Phenix real-space refinement model.
6. **Interpretive Over-Claims in Existing Report**: Claims that static Cα RMSD (0.338 Å) rules out allosteric regulation, that TPP1 DE-loop displacement (4.293 Å) demonstrates a two-step kinetic induced-fit pathway, and that TELIN-1 represents a viable lead candidate are unsubstantiated by static coordinates or peer-reviewed literature.

---

## 2. Environment & Version Control Audit

### 2.1 Git Repository Status
- **Status**: `UNINITIALIZED`
- **Error Command**: `git status` -> `fatal: not a git repository (or any of the parent directories): .git`
- **Forensic Assessment**: The workspace is operated as a standalone directory tree without Git history tracking. To establish full audit provenance and change-tracking across multi-agent milestones, Git initialization (`git init`) or automated manifest hashing is required.

### 2.2 System & Runtime Environment
- **Operating System**: Microsoft Windows (win32)
- **Python Version**: Python 3.10.11 (`C:\Users\hp\AppData\Local\Programs\Python\Python310\python.exe`)
- **Key Dependencies & Tooling**:
  - `gemmi`: 0.6.x (Primary high-performance mmCIF parser & coordinate extraction engine)
  - `biopython`: 1.8x (Secondary validation parser `Bio.PDB.MMCIFParser` and `Bio.PDB.SASA`)
  - `freesasa`: C-accelerated Shrake-Rupley / Lee-Richards surface area engine
  - `numpy`: 1.24+ (SVD Kabsch linear algebra, grid cell spatial hashing)
  - `pandas`: 2.0+ (Standardized tabular data generation and alignment manipulation)
  - `matplotlib`: 3.7+ (300 DPI publication figure rasterization)
  - `pytest`: 9.1.1 (Testing framework)
  - `pytest-cov`: 7.1.0 (Coverage analysis)

---

## 3. Raw Experimental Data Inventory & Cryptographic Validation

All five raw structures were downloaded directly from the Worldwide Protein Data Bank (wwPDB / RCSB PDB) alongside official wwPDB validation XML reports and full validation PDFs.

### Table 1: Primary Structural Dataset Manifest & Hash Verification
| PDB ID | Experimental Method | Resolution | File Size (Bytes) | SHA-256 Checksum (Computed) | Manifest Checksum (Expected) | Integrity Verdict |
| :--- | :--- | :---: | :---: | :--- | :--- | :---: |
| **7TRE** | Cryo-EM (Single Particle) | 3.50 Å | 1,518,438 | `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` | `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` | **MATCH (100%)** |
| **7TRD** | Cryo-EM (Single Particle) | 3.30 Å | 1,394,494 | `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` | `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` | **MATCH (100%)** |
| **7TRF** | Cryo-EM (Single Particle) | 3.70 Å | 1,542,948 | `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` | `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` | **MATCH (100%)** |
| **7TRC** | Cryo-EM (Single Particle) | 3.30 Å | 1,647,704 | `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` | `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` | **MATCH (100%)** |
| **2I46** | X-Ray Diffraction | 2.70 Å | 273,104 | `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` | `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` | **MATCH (100%)** |

### 3.1 Macromolecular Entity & Chain Disambiguation
- **7TRE**:
  - `Chain A` (Author) / `Chain B` (Label): Human TERT catalytic subunit (934 modeled residues, 7,488 heavy atoms, 7–1132). Recombinant Twin-Strep tag creates offset `auth = label - 35`.
  - `Chain C` (Author) / `Chain D` (Label): Human TPP1 OB-fold domain (152 modeled residues, 1,171 heavy atoms, 90–241). Recombinant N-terminal GAG cloning tag creates offset `auth = label + 86`.
  - `Chain B` (Author) / `Chain A` (Label): Human Telomerase RNA component (hTR / TER; 217 modeled nucleotides).
  - `Chain D` (Author) / `Chain C` (Label): Telomeric DNA substrate (7 modeled nucleotides, `5'-TTAGGGT-3'`).
- **7TRD**:
  - Reference unbound catalytic core. TERT (Chain A/B), hTR (Chain B/A), DNA (Chain N/C). TPP1 absent.
- **7TRF**:
  - Telomerase catalytic core RNP complexed with human histone H2A–H2B heterodimers. TERT (Chain A/D), H2B (Chain F/A), H2A (Chain E/B), hTR (Chain B/C), DNA (Chain D/E).
- **7TRC**:
  - Human telomerase H/ACA RNP lobe. Contains TCAB1, DKC1, NOP10, NHP2, GAR1, and 3' hTR lobe. TERT and TPP1 absent.
- **2I46**:
  - Crystal structure of human TPP1 OB-fold domain (homodimer in asymmetric unit).
  - **Chain A**: Fully modeled residues 90–241 (152 aa). Entire TEL-patch DE-loop (Glu168–Glu171) is fully resolved.
  - **Chain B**: DE-loop residues 168–171 completely disordered (unmodeled). Crucial: any comparative alignment must strictly target Chain A.

---

## 4. Source Code Architecture & Module Audit

The `src/telotrade/` package comprises 8 Python source files totaling 5,036 lines of code.

### Table 2: Source Code Module Inventory
| Module Path | Lines | Size (Bytes) | SHA-256 (Prefix) | Functional Responsibility | Audited Architectural Notes |
| :--- | :---: | :---: | :---: | :--- | :--- |
| `src/telotrade/__init__.py` | 151 | 4,153 | `d9ac060d4b98` | Package entrypoint & public API exports | Clean namespace exports for all pipeline functions. |
| `src/telotrade/__main__.py` | 6 | 97 | `976e3e64de11` | CLI invocation proxy | Forwards directly to `telotrade.pipeline:main`. |
| `src/telotrade/validate.py` | 595 | 22,406 | `b404c1d4d02b` | Baseline integrity & wwPDB validation parser | Zero-crash handling of missing residues, 0-byte files, and altlocs. |
| `src/telotrade/sequence_mapping.py` | 1,331 | 46,376 | `873a8e98096d` | Sequence extraction, tag offsets & homology mapping | **Defect identified**: off-by-one indexing in `ob_seq[resnum - 89]`. |
| `src/telotrade/interface_geometry.py` | 880 | 30,991 | `8692d076d29e` | Spatial contact hashing, BSA & polar networks | Employs Bondi radii and 960-point Shrake-Rupley sphere discretization. |
| `src/telotrade/cross_state.py` | 951 | 35,419 | `6d8c82348a18` | Kabsch SVD superposition & Cα RMSD profiling | Enforces $SO(3)$ proper rotation ($\det(R) = +1.0$); strictly selects 2I46 Chain A. |
| `src/telotrade/figures.py` | 793 | 33,543 | `3bf1fc3c7620` | Publication-grade Matplotlib figure rendering | 300 DPI, atomic tmp-file saving, proper `plt.close()` cleanup. |
| `src/telotrade/pipeline.py` | 327 | 11,990 | `23bbef5e2848` | End-to-end automated runner & orchestrator | Manages directories, executes Stages 1–5, verifies deliverables. |

---

## 5. Standardized Tabular Results Audit (`results/tables/`)

A total of 21 CSV files reside in `results/tables/`, containing the quantitative structural bioinformatics metrics generated by the pipeline.

### Table 3: Standardized Tables Inventory & Schema Verification
| File Name | Rows | Columns | Size (Bytes) | Core Measured Metrics & Verified Values | Schema Compliance |
| :--- | :---: | :---: | :---: | :--- | :---: |
| `buried_surface_area.csv` | 1 | 17 | 427 | Interface BSA = **856.82 Å²**; $\Delta\mathrm{SASA} = 1713.65\ \text{Å}^2$; isolated TERT = 49,776.00 Å²; isolated TPP1 = 10,651.36 Å² | PASSED |
| `chain_entity_mapping.csv` | 19 | 11 | 2,648 | Complete entity and chain mappings across all 5 structures | PASSED |
| `cross_state_rmsd.csv` | 11 | 9 | 887 | Catalytic Core RMSD = **0.338 Å** (7TRE vs 7TRD / 7TRF); Global TERT = 0.356 Å; TPP1 OB-fold = **1.966 Å** (7TRE vs 2I46 A) | PASSED |
| `data_manifest.csv` | 5 | 8 | 2,461 | Master manifest of PDB IDs, resolutions, SHA-256 hashes | PASSED |
| `domain_rmsd.csv` | 12 | 10 | 985 | TEN independent = 0.411 Å; TEN core-fit = **0.434 Å**; TRAP/IFD core-fit = **0.427 Å**; TRBD = 0.270 Å; RT = 0.295 Å; CTE = 0.238 Å | PASSED |
| `domain_rmsds.csv` | 12 | 10 | 985 | Duplicate alias of `domain_rmsd.csv` | PASSED |
| `interface_contacts_3.5A.csv` | 9 | 14 | 754 | Exactly **9 contacts** across 5 residue pairs (Core Robust) | PASSED |
| `interface_contacts_4.0A.csv` | 46 | 14 | 3,038 | Exactly **46 contacts** across 21 residue pairs (37 Intermediate) | PASSED |
| `interface_contacts_4.5A.csv` | 123 | 14 | 7,772 | Exactly **123 contacts** across 32 residue pairs (77 Sensitive) | PASSED |
| `interface_contacts_7tre_3.5A.csv` | 9 | 14 | 754 | Canonical naming for 3.5 Å contacts | PASSED |
| `interface_contacts_7tre_4.0A.csv` | 46 | 14 | 3,038 | Canonical naming for 4.0 Å contacts | PASSED |
| `interface_contacts_7tre_4.5A.csv` | 123 | 14 | 7,772 | Canonical naming for 4.5 Å contacts | PASSED |
| `interface_contacts_robustness.csv` | 123 | 11 | 12,323 | Robustness classification: 9 Core Robust, 37 Intermediate, 77 Sensitive | PASSED |
| `per_residue_displacement_tert.csv` | 934 | 6 | 31,561 | Per-residue Cα displacements across 934 modeled TERT residues | PASSED |
| `polar_interactions.csv` | 11 | 10 | 802 | Primary salt bridge (Arg774-Glu168: **2.32 Å**); secondary clamp (Lys78-Glu215: **3.58 Å**) | PASSED |
| `polar_interactions_7tre.csv` | 11 | 10 | 802 | Canonical naming for polar interactions table | PASSED |
| `residue_correspondence_tert.csv` | 1,132 | 12 | 56,958 | 1,132-row canonical UniProt O14746 alignment across 7TRE, 7TRD, 7TRF | PASSED |
| `residue_correspondence_tpp1.csv` | 161 | 10 | 7,301 | 161-row canonical UniProt Q96AP0 alignment (contains 1-aa offset anomaly in resname column) | AUDIT FLAG |
| `tel_patch_contacts.csv` | 53 | 15 | 4,216 | Contacts involving TEL patch residues (DE-loop, Leu183, Leu212) at 4.5 Å | PASSED |
| `tpp1_bound_vs_unbound_rmsd.csv` | 152 | 7 | 7,578 | Per-residue shifts: Lys170 peak = **4.293 Å**; Leu183 anchor = **0.512 Å**; Leu212 = 1.040 Å | PASSED |
| `tpp1_conformational_shifts.csv` | 152 | 7 | 7,578 | Duplicate alias of `tpp1_bound_vs_unbound_rmsd.csv` | PASSED |

---

## 6. Publication Figures Audit (`results/figures/`)

### Table 4: Figures Inventory
| File Name | Dimensions (Pixels) | DPI | File Size (Bytes) | Biological Visualization |
| :--- | :---: | :---: | :---: | :--- |
| `figure1_interface_contacts.png` | 3300 x 2400 | 300 | 548,365 | 2D contact matrix highlighting TEN (Lys78) and TRAP/IFD (Arg774) cradle docking onto TPP1 TEL-patch. |
| `figure2_cross_state_rmsd.png` | 3300 x 2400 | 300 | 404,756 | Multi-panel bar and displacement plot illustrating catalytic core rigidity (0.338 Å) vs TEN excursion (0.434 Å). |
| `figure3_tpp1_conformational_shifts.png` | 3300 x 2400 | 300 | 912,897 | Per-residue Cα displacement profile highlighting the DE-loop peak at Lys170 (4.293 Å) vs invariant Leu183 (0.512 Å). |
| `figure4_cutoff_robustness.png` | 3300 x 2400 | 300 | 400,073 | Threshold escalation staircase (9 -> 46 -> 123 contacts) and persistence partition. |

---

## 7. Testing Infrastructure & Verification Status

The test suite covers unit logic, boundary limits, adversarial stress, and end-to-end integration across 12 test files.

### Table 5: Test Suite Inventory
| Test File | Test Suite Tier | Test Focus |
| :--- | :--- | :--- |
| `tests/e2e/test_tier1_features.py` | E2E Tier 1 | Features F1 through F18 functional verification |
| `tests/e2e/test_tier2_boundaries.py` | E2E Tier 2 | Numerical boundaries, zero cutoffs, missing residue edge cases |
| `tests/e2e/test_tier3_interactions.py` | E2E Tier 3 | Cross-module pipeline dataflow and subset monotonicity |
| `tests/e2e/test_tier4_scenarios.py` | E2E Tier 4 | Full pipeline execution and report deliverables |
| `tests/test_adversarial_m1.py` | Adversarial M1 | Checksum tampering, corrupted CIFs, XML schema errors |
| `tests/test_adversarial_m2.py` | Adversarial M2 | Truncated sequences, tag insertion offsets, dimer asymmetry |
| `tests/test_adversarial_m3.py` | Adversarial M3 | Spatial neighbor grid jitter, Shrake-Rupley sphere point density |
| `tests/test_adversarial_m4.py` | Adversarial M4 | Kabsch SVD reflection detection ($\det(R)=-1.0$), non-matching point sets |
| `tests/test_adversarial_m5_figures.py` | Adversarial M5 | Figure perturbation sensitivity, resource leak checks |
| `tests/test_cross_state.py` | Integration M4 | Comprehensive cross-state RMSD and domain evaluations |
| `tests/test_interface_geometry.py` | Integration M3 | Interface contacts, BSA formulas, polar networks |
| `tests/test_sequence_mapping.py` | Integration M2 | Residue correspondence and chain entity mappings |
| `tests/test_validate.py` | Integration M1 | Data validation and parsing engine tests |

---

## 8. Headline Claims Evidence Ledger Summary

Across `report.md` and `results/tables/`, exactly 24 discrete headline metrics and mechanistic claims were cataloged and audited.

### Evidence Status Distribution:
- **Established (17 claims, 70.8%)**: Computationally and mathematically confirmed against primary structural coordinates in `7TRE.cif` and `2I46.cif`. Includes BSA (856.82 Å²), contact counts (9, 46, 123), Leu183 rigidity (0.512 Å), Cα displacements across domain subsets, and sequence alignments.
- **Uncertain (4 claims, 16.7%)**:
  1. `CLM-005`: Primary salt bridge distance Arg774-Glu168 (2.32 Å) — physically uncharacteristic, likely a refinement artifact of 3.5 Å cryo-EM map.
  2. `CLM-010`: Catalytic core RMSD (0.338 Å) — within experimental coordinate error (~0.35–0.5 Å) and shared refinement bias.
  3. `CLM-014`: Catalytic core RMSD vs 7TRF (0.338 Å) — same coordinate uncertainty boundary.
  4. `CLM-017`: TPP1 DE-loop peak shift (4.293 Å) — affected by 2I46 crystallographic packing contacts and high B-factors (>80 Å²).
- **Contradicted (2 claims, 8.3%)**:
  1. `CLM-022`: Claim that 0.338 Å static RMSD proves complete absence of allosteric regulation — contradicted by biophysical principles of dynamic/entropic allostery.
  2. `CLM-024`: TPP1 canonical residue name column in `residue_correspondence_tpp1.csv` — contradicted by actual UniProt sequence due to off-by-one indexing in `sequence_mapping.py`.
- **Untested (2 claims, 8.3%)**:
  1. `CLM-021`: Two-step sequential induced-fit recruitment mechanism (Leu183 first, DE-loop second) — cannot be established from static crystal and cryo-EM endpoints.
  2. `CLM-023`: Candidate inhibitor TELIN-1 efficacy, binding affinity, and selectivity — zero wet-lab, literature, or patent evidence exists.

---

## 9. Identified Discrepancies & Recommendations for Phase 1–5

1. **Fix Sequence Correspondence Indexing Bug**:
   In `src/telotrade/sequence_mapping.py::generate_tpp1_correspondence_table`, correct the index offset:
   ```python
   # Current (defective):
   idx = resnum - 89
   # Correct:
   idx = resnum - 90  # Since ob_seq begins at residue 90 (Gly90)
   ```
2. **Explicitly Contextualize the 2.32 Å Salt Bridge**:
   In `report.md` Section 4.2 and Section 6, explicitly note that 2.32 Å is an overcompacted coordinate distance from 3.5 Å cryo-EM refinement, not a physically validated ground-truth bond length.
3. **Decouple Static Coordinate RMSD from Dynamic Allostery**:
   Replace the claim that "0.338 Å proves absence of allostery" with a nuanced discussion separating static rigid-body invariance from dynamic/entropic allosteric coupling and repeat addition processivity (RAP).
4. **Quarantine TELIN-1 as a Purely Hypothetical In Silico Construct**:
   Ensure TELIN-1 is clearly designated as an unvalidated computational candidate with zero demonstrated affinity or cellular selectivity.
