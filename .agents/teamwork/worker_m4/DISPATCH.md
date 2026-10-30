# DISPATCH — worker_m4

## Mission
Implement Milestone 4: Cross-State Comparative Structural Analysis & Robustness Testing (R4) in `src/telotrade/cross_state.py` and unit tests in `tests/test_cross_state.py`.

## Working Directory
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m4`

## Mandatory Reference Documents
Read these files thoroughly before starting:
1. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md` (architecture, milestones, interface contracts)
3. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md` (Section 3: SVD Kabsch, domain definitions, benchmarks; Section 4: CSV table schemas)
4. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md` (offsets, 2I46 asymmetry: Chain A vs Chain B DE-loop disorder)
5. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\sequence_mapping.py` (sequence extraction, offsets, pairwise alignment)
6. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\interface_geometry.py` (contact functions, robustness classification)

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Detailed Scope & Technical Requirements
1. **Module Implementation in `src/telotrade/cross_state.py`**:
   - Rigid-Body Superposition Engine:
     - `svd_kabsch_align(coords_target: np.ndarray, coords_mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]`:
       - Implements SVD Kabsch superposition via NumPy:
         - Centroid centering: $P_c, Q_c$.
         - Cross-covariance matrix $H = P_c^T Q_c$.
         - SVD: $H = U \Sigma V^T$.
         - Proper rotation enforcement: $d = \det(V U^T)$; $R = V \mathrm{diag}(1, 1, d) U^T$.
         - Translation vector: $\mathbf{t} = \bar{\mathbf{q}} - R \bar{\mathbf{p}}$.
         - Returns aligned coordinates, rotation matrix $R$, translation vector $\mathbf{t}$, and Cα RMSD.
   - Cross-State Superposition & RMSD Calculations:
     - `compute_cross_state_rmsd(raw_dir='data/raw') -> pd.DataFrame`:
       - Evaluates:
         - 7TRE vs 7TRD (Global Cα, Catalytic Core, TEN domain).
         - 7TRE vs 7TRF (Global Cα, Catalytic Core, TEN domain).
         - 7TRE vs 2I46 Chain A (Global Cα TPP1, TEL-patch).
       - Benchmarks:
         - 7TRE vs 7TRD catalytic core ($n=770$ matched Cα) is $\sim 0.338\ \text{Å}$.
         - 7TRE vs 7TRD global Cα is $\sim 0.356\ \text{Å}$.
         - 7TRE vs 2I46 Chain A TPP1 global Cα is $\sim 1.966\ \text{Å}$.
       - Exports/formats matching `cross_state_rmsd.csv` schema in `survey_methodology/report.md` Section 4.9.
   - Domain-Specific Structural Rearrangements:
     - `compute_domain_rmsds(raw_dir='data/raw') -> pd.DataFrame`:
       - Dissects TERT domains using canonical definitions:
         - `TEN`: residues 7–196 (UniProt 1–199)
         - `TRBD`: residues 326–540
         - `RT`: residues 541–935
         - `IFD_TRAP`: residues 755–785
         - `CTE`: residues 936–1132
         - `Catalytic Core`: residues 326–1132
       - Computes RMSD for each domain independently and under Catalytic Core fit.
       - Generates `domain_rmsd.csv` (or `domain_rmsds.csv`).
   - Per-Residue Cα Displacements on Core Fit:
     - `compute_per_residue_displacement(raw_dir='data/raw') -> pd.DataFrame`:
       - Superimposes 7TRD and 7TRF onto 7TRE based on the catalytic core (residues 326–1132).
       - Computes Euclidean distance for each residue of TERT (7–1132) between 7TRE and 7TRD/7TRF.
       - Proves that the catalytic core is a rigid monolith while TEN domain exhibits hinge displacement (~0.434 Å).
       - Generates `per_residue_displacement_tert.csv`.
   - Bound TPP1 (7TRE) vs Unbound TPP1 (2I46) Plasticity:
     - `compute_tpp1_conformational_shifts(raw_dir='data/raw') -> pd.DataFrame`:
       - Superimposes 7TRE TPP1 and 2I46 Chain A TPP1 (strictly Chain A!).
       - Computes per-residue Cα displacement for residues 90–241.
       - Verifies that:
         - Leu183 is structurally pre-organized in the unbound state (displacement $\sim 0.512\ \text{Å}$, `Rigid Anchor`).
         - Acidic DE-loop (Glu168–Glu171) undergoes extensive plastic remodeling (displacements up to $4.293\ \text{Å}$ at Lys170, `Extensive Remodeling`).
         - Classifies residues: `Rigid Anchor` (<1.0 Å), `Moderate Shift` (1.0–2.5 Å), `Extensive Remodeling` (>2.5 Å).
         - Extracts local B-factor / Q-score if available.
       - Generates `tpp1_bound_vs_unbound_rmsd.csv` (and `tpp1_conformational_shifts.csv`).
   - Cutoff Robustness Analysis:
     - `compute_cutoff_robustness(cif_path='data/raw/7TRE.cif') -> pd.DataFrame`:
       - Evaluates persistence across cutoffs (3.5 Å, 4.0 Å, 4.5 Å), classifying contacts into Core Robust, Intermediate, and Threshold Sensitive.
       - Generates `interface_contacts_robustness.csv` (and `contact_robustness.csv`).
   - Table Generation:
     - `generate_cross_state_tables(output_dir='results/tables', raw_dir='data/raw')`:
       - Saves all generated tables to `results/tables/`.
2. **Expose Core APIs in `src/telotrade/__init__.py`**:
   - Export `svd_kabsch_align`, `compute_cross_state_rmsd`, `compute_domain_rmsds`, `compute_per_residue_displacement`, `compute_tpp1_conformational_shifts`, `generate_cross_state_tables`.
3. **Unit Tests in `tests/test_cross_state.py`**:
   - Comprehensive tests validating Kabsch alignment on known rotation/translation test sets, exact RMSD benchmark values (7TRE vs 7TRD catalytic core 0.338 Å, 7TRE vs 2I46 TPP1 1.966 Å), DE-loop plasticity vs Leu183 rigidity, domain partition calculations, CSV schemas, and handling missing/corrupted files gracefully.
4. **Verification**:
   - Run `pytest tests/test_cross_state.py -v`.
   - Run full test suite `pytest tests/ -v`.
   - Run `ruff check src/ tests/`.
5. **Handoff**:
   - Write standard 5-component `handoff.md` (Observation, Logic Chain, Caveats, Conclusion, Verification Method) in your working directory and notify orchestrator.
