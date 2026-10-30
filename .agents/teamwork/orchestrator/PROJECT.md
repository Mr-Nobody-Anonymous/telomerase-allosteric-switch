# Project: Human Telomerase TERT–TPP1 Structural Bioinformatics Study

## Architecture
- **Language & Runtime**: Python 3.10.11 AMD64
- **Core Dependencies**: Biopython 1.88, Gemmi 0.7.5, FreeSASA 2.2.1, NumPy 1.26.4, SciPy 1.15.3, Pandas 2.3.3, Matplotlib 3.10.7, Requests 2.32.5, Pytest 9.1.1
- **Target Dataset**:
  - `7TRE`: Human telomerase catalytic core complexed with TPP1 (cryo-EM 3.5 Å, EMD-26087). Primary recruitment interface structure.
  - `7TRD`: Human telomerase catalytic core reference (cryo-EM 3.3 Å, EMD-26086).
  - `7TRF`: Telomerase catalytic core RNP with H2A–H2B histones (cryo-EM 3.5 Å, EMD-26088).
  - `7TRC`: Human telomerase H/ACA RNP lobe (cryo-EM 3.4 Å, EMD-26085).
  - `2I46`: Unbound TPP1 OB-fold reference crystal structure (X-ray 2.70 Å, Space Group I 41 2 2).
- **Directory Layout**:
  - `data/raw/`: Immutable raw mmCIF files and official validation reports.
  - `data/manifest/`: Master `data_manifest.csv` with SHA-256 hashes, source URLs, resolutions, chain mappings.
  - `src/telotrade/`: Modular, production-grade structural analysis package:
    - `validate.py`: mmCIF integrity, checksum verification, validation metrics, exception-free parsing.
    - `sequence_mapping.py`: Sequence-aware residue correspondence, UniProt mapping, expression tag offset handling.
    - `interface_geometry.py`: Contact network at 3.5 Å, 4.0 Å, 4.5 Å, Shrake-Rupley / FreeSASA BSA, polar/salt bridge networks, hydrophobic packing.
    - `cross_state.py`: SVD Kabsch superposition, domain-specific Cα RMSDs, TEN-TRAP rearrangements, 7TRE vs 2I46 Chain A TEL-patch conformational shift.
    - `pipeline.py`: End-to-end execution engine running all analyses and generating artifacts.
  - `results/tables/`: Standardized CSV tables (11 tables covering all metrics).
  - `results/figures/`: Publication-quality structural figures and comparative plots.
  - `report.md`: Publication-grade synthesis report separating Axis A (Biophysics) from Axis B (Therapeutic Selectivity & Stem Cell Trade-Off).
  - `tests/`: Comprehensive test suite (Tiers 1–5).

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| 1 | Raw Data Acquisition | Download 7TRE, 7TRD, 7TRF, 7TRC, 2I46 mmCIF and validation data to `data/raw/` | M1 | ORIGINAL_REQUEST R1 |
| 2 | Immutable Data Manifest | Compile `data_manifest.csv` with SHA-256 hashes, URLs, resolutions, chain mappings | M1 | ORIGINAL_REQUEST R1 |
| 3 | Baseline Validation Script | Implement `src/telotrade/validate.py` with zero-crash handling of missing residues/altlocs | M1 | ORIGINAL_REQUEST R1 |
| 4 | Sequence-Aware Chain Mapping | Automated Gemmi/Biopython parsers establishing sequence correspondence across homologous chains | M2 | ORIGINAL_REQUEST R2 |
| 5 | Construct Tag Offset Handling | Reconcile TERT +35 Twin-Strep tag and TPP1 +86 GAG tag to canonical UniProt | M2 | ORIGINAL_REQUEST R2 |
| 6 | 2I46 Dimer Asymmetry Handling | Select Chain A (resolved DE-loop 168-171) and handle Chain B disorder | M2 | ORIGINAL_REQUEST R2 |
| 7 | Multi-Cutoff Contact Networks | Compute 7TRE TERT-TPP1 contacts at 3.5 Å, 4.0 Å, 4.5 Å with atom selections | M3 | ORIGINAL_REQUEST R3 |
| 8 | Buried Surface Area (BSA) | Compute BSA using Shrake-Rupley / FreeSASA with formula $(\mathrm{SASA}_A + \mathrm{SASA}_B - \mathrm{SASA}_{AB})/2$ | M3 | ORIGINAL_REQUEST R3 |
| 9 | Polar Networks & Salt Bridges | Map TERT Arg774 -- TPP1 Glu168 (2.32 Å), Lys78 -- Glu215 (3.58 Å), and H-bonds | M3 | ORIGINAL_REQUEST R3 |
| 10 | Hydrophobic Packing & TEL-Patch | Map TERT Leu769, Pro771 with TPP1 Leu183, Leu212; map TEL patch (168-171, 183, 212) | M3 | ORIGINAL_REQUEST R3 |
| 11 | Cross-State SVD Kabsch RMSD | Compute global and domain-specific Cα RMSD across 7TRE, 7TRD, 7TRF | M4 | ORIGINAL_REQUEST R4 |
| 12 | TEN-TRAP Domain Rearrangements | Quantify conformational shifts in TEN and TRAP/IFD upon TPP1 binding | M4 | ORIGINAL_REQUEST R4 |
| 13 | Interface Cutoff Robustness Table | Analyze contact persistence vs threshold sensitivity across 3.0 Å to 5.0 Å | M4 | ORIGINAL_REQUEST R4 |
| 14 | TPP1 Bound vs Unbound Shift | Quantify TEL-patch DE-loop remodeling (7TRE bound vs 2I46 Chain A unbound) | M4 | ORIGINAL_REQUEST R4 |
| 15 | Standardized Results CSV Tables | Populate `results/tables/` with 11 CSV outputs for all calculated metrics | M5 | ORIGINAL_REQUEST R5 |
| 16 | Publication-Grade Synthesis Report | Author `report.md` with executive summary, methodology, results, Axis A vs B | M5 | ORIGINAL_REQUEST R5 |
| 17 | Publication Figures Generation | Generate publication-grade figures in `results/figures/` | M5 | ORIGINAL_REQUEST R5 |
| 18 | Automated Single-Run Pipeline | End-to-end `run_pipeline.py` script reproducing all tables, figures, and report | M5 | ORIGINAL_REQUEST R5 |
| 19 | E2E Opaque-Box Test Suite | Comprehensive 4-tier test suite covering feature, boundary, pairwise, and application tests | M6 | ORIGINAL_REQUEST Acceptance Criteria |
| 20 | Adversarial Coverage & Forensic Audit | White-box adversarial testing (Tier 5) and forensic integrity audit | M6 | Orchestration Pattern & Hard Constraints |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| M1 | Data Acquisition & Baseline Validation | Features 1, 2, 3 (`data/raw/`, `data_manifest.csv`, `src/telotrade/validate.py`) | none | DONE |
| M2 | Sequence-Aware Chain & Residue Correspondence | Features 4, 5, 6 (`src/telotrade/sequence_mapping.py`, UniProt alignment, offset handling) | M1 | DONE |
| M3 | Quantitative Interface Geometry & Interaction Profiling | Features 7, 8, 9, 10 (`src/telotrade/interface_geometry.py`, contact networks, BSA, polar, TEL-patch) | M2 | DONE |
| M4 | Cross-State Structural Analysis & Robustness | Features 11, 12, 13, 14 (`src/telotrade/cross_state.py`, SVD Kabsch, RMSDs, 7TRE vs 2I46, robustness) | M2, M3 | DONE |
| M5 | Evidence Matrix, Report, Figures & Pipeline Runner | Features 15, 16, 17, 18 (`report.md`, `results/tables/`, `results/figures/`, `run_pipeline.py`) | M3, M4 | DONE |
| M6 | E2E Testing, Adversarial Hardening & Forensic Audit | Features 19, 20 (`tests/`, Tiers 1-5 pass 100%, Forensic Audit CLEAN) | M1–M5 | IN_PROGRESS |

## Interface Contracts

### `src/telotrade/validate.py` ↔ Data Pipeline
- `download_structures(target_ids, output_dir) -> Dict[str, Path]`
- `verify_checksums(manifest_path, raw_dir) -> bool`
- `parse_validation_metrics(entry_id) -> Dict[str, Any]`
- CLI: `python -m telotrade.validate [--raw-dir DIR] [--manifest PATH]`

### `src/telotrade/sequence_mapping.py` ↔ Interface & Cross-State Analysis
- `get_chain_sequence(cif_path, chain_id, use_auth=True) -> Tuple[str, List[int]]`
- `map_residue_correspondence(cif_path_1, chain_1, cif_path_2, chain_2) -> Dict[int, int]`
- `get_uniprot_mapping(cif_path, chain_id) -> Dict[int, int]`
- Canonical offsets:
  - TERT (7TRE, 7TRD, 7TRF): `auth_seq_num = label_seq_id - 35` (UniProt O14746)
  - TPP1 (7TRE): `auth_seq_num = label_seq_id + 86` (UniProt Q96AP0)
  - TPP1 (2I46): `auth_seq_num = label_seq_id + 89` (Chain A strictly)

### `src/telotrade/interface_geometry.py` ↔ Output Tables
- `compute_interface_contacts(cif_path, chain_tert, chain_tpp1, cutoffs=[3.5, 4.0, 4.5]) -> Dict[float, pd.DataFrame]`
- `compute_bsa(cif_path, chain_tert, chain_tpp1, probe_radius=1.4) -> Dict[str, float]`
  - Formula: $\mathrm{BSA} = (\mathrm{SASA}_A + \mathrm{SASA}_B - \mathrm{SASA}_{AB}) / 2$
- `compute_polar_interactions(contacts_df) -> pd.DataFrame`
- `map_tel_patch_contacts(contacts_df) -> pd.DataFrame`

### `src/telotrade/cross_state.py` ↔ Output Tables
- `svd_kabsch_align(coords_target, coords_mobile) -> Tuple[np.ndarray, float]`
- `compute_domain_rmsds(cif_ref, cif_mobile, domain_defs) -> pd.DataFrame`
- `compute_tpp1_remodeling(cif_7tre, cif_2i46_a) -> pd.DataFrame`
- `compute_cutoff_robustness(cif_path, chain_a, chain_b, cutoffs_range) -> pd.DataFrame`

## Code Layout
```
C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\
├── data\
│   ├── raw\                  # 7TRE.cif, 7TRD.cif, 7TRF.cif, 7TRC.cif, 2I46.cif, validation reports
│   └── manifest\             # data_manifest.csv
├── src\
│   └── telotrade\
│       ├── __init__.py
│       ├── validate.py
│       ├── sequence_mapping.py
│       ├── interface_geometry.py
│       ├── cross_state.py
│       └── pipeline.py
├── results\
│   ├── tables\               # 11 standardized CSV files
│   └── figures\              # publication-grade figures
├── tests\                    # unit and E2E test suites (Tiers 1-5)
├── report.md                 # synthesis report with Axis A vs Axis B
├── run_pipeline.py           # single-entry automated pipeline runner
└── pyproject.toml / requirements.txt
```
