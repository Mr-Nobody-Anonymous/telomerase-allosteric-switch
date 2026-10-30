# Forensic Audit Report — Milestone 4 (R4)

**Work Product**: `src/telotrade/cross_state.py`, `tests/test_cross_state.py`, and `results/tables/`  
**Profile**: General Project (Integrity Mode: `development`, per `ORIGINAL_REQUEST.md` line 13)  
**Verdict**: **CLEAN**  
**Auditor**: `auditor_m4`  
**Parent Agent**: `eb46fecb-7562-4ce7-aa98-998216cee2c4`  
**Date**: 2026-10-16  

---

## 1. Observation

### 1.1 Static Code & Integrity Analysis
1. **No Hardcoded Output Stubs**:
   - Inspected `src/telotrade/cross_state.py` (952 lines).
   - In lines 48–58, `BENCHMARK_RMSDS` is declared as an informational constant dictionary. Grep analysis across the entire project (`grep_search Query="BENCHMARK_RMSDS"`) revealed that it is referenced on line 48 only. None of the computational functions (`svd_kabsch_align`, `compute_cross_state_rmsd`, `compute_domain_rmsds`, `compute_per_residue_displacement`, `compute_tpp1_conformational_shifts`) access or return values from `BENCHMARK_RMSDS`.
   - All return values are computed dynamically from `gemmi` atomic structures and NumPy linear algebra routines.

2. **No Mocking or Test Suppression**:
   - Inspected `tests/test_cross_state.py` (609 lines, 29 test methods).
   - Grep searches for `mock`, `unittest`, `patch`, and `monkeypatch` returned zero matches. All tests execute live calculations against real mmCIF files (`7TRE.cif`, `7TRD.cif`, `7TRF.cif`, `2I46.cif`) or dynamically generated synthetic matrices.

3. **Pre-Populated Artifact Integrity**:
   - Re-executed `generate_cross_state_tables(output_dir=tmpdir)` and compared outputs with files committed in `results/tables/`.
   - Results:
     * `cross_state_rmsd.csv`: EXACT MATCH (11 rows, 9 cols)
     * `domain_rmsd.csv`: EXACT MATCH (12 rows, 10 cols)
     * `domain_rmsds.csv`: EXACT MATCH (12 rows, 10 cols)
     * `per_residue_displacement_tert.csv`: EXACT MATCH (934 rows, 6 cols)
     * `tpp1_bound_vs_unbound_rmsd.csv`: EXACT MATCH (152 rows, 7 cols)
     * `tpp1_conformational_shifts.csv`: EXACT MATCH (152 rows, 7 cols)
     * `interface_contacts_robustness.csv`: EXACT MATCH (123 rows, 11 cols)
   - Zero discrepancies between pre-existing artifacts and fresh execution.

### 1.2 Mathematical & Runtime Behavioral Verification
1. **Rigid-Body Superposition Ground-Truth Parity**:
   - Tested `svd_kabsch_align` with a synthetic 50-point cloud subjected to a known random 3D rotation $Q_{rot} \in SO(3)$ and translation $\mathbf{t} = [5.2, -3.1, 7.8]$:
     * Returned RMSD: $< 1 \times 10^{-12}$ Å.
     * Recovered rotation matrix: $\|\mathbf{R} - Q_{rot}\|_\infty < 1 \times 10^{-10}$.
     * Recovered translation vector: $\|\mathbf{t}_{calc} - \mathbf{t}_{true}\|_\infty < 1 \times 10^{-10}$.
     * $\det(\mathbf{R}) = 1.0000000000$.

2. **Chiral Reflection & Proper Rotation Enforcement**:
   - Inverted coordinates along the Z-axis ($z \to -z$, improper reflection matrix with determinant $-1$).
   - `svd_kabsch_align` applied reflection factor $d = \mathrm{sign}(\det(VU^T))$, yielding proper rotation with $\det(\mathbf{R}) = +1.000000$, strictly preventing illegal chirality inversion.

3. **Coordinate Sensitivity & Noise Injection**:
   - Tested actual 7TRE vs 7TRD C-alpha coordinates with additive Gaussian noise $\sigma \in \{0.0, 0.1, 0.5, 1.0, 5.0\}$ Å:
     * $\sigma = 0.0$ Å: $\mathrm{RMSD} = 0.355718$ Å
     * $\sigma = 0.1$ Å: $\mathrm{RMSD} = 0.390712$ Å
     * $\sigma = 0.5$ Å: $\mathrm{RMSD} = 0.914564$ Å
     * $\sigma = 1.0$ Å: $\mathrm{RMSD} = 1.730474$ Å
     * $\sigma = 5.0$ Å: $\mathrm{RMSD} = 8.508548$ Å
   - Confirms that calculations dynamically depend on 3D spatial coordinates.

4. **Empirical Biological Metric Parity**:
   - 7TRE vs 7TRD Catalytic Core (residues 326–1132, 770 atoms): $\mathrm{RMSD} = 0.3377$ Å (Benchmark: ~0.338 Å).
   - 7TRE vs 7TRD Global (934 atoms): $\mathrm{RMSD} = 0.3557$ Å (Benchmark: ~0.356 Å).
   - 7TRE vs 7TRD TEN domain alone (160 atoms): $\mathrm{RMSD} = 0.4112$ Å (Benchmark: ~0.411 Å).
   - 7TRE vs 7TRD TEN domain on Core fit: $\mathrm{RMSD} = 0.4341$ Å (Benchmark: ~0.434 Å).
   - Bound TPP1 (7TRE Chain C) vs Unbound TPP1 (2I46 Chain A, 152 atoms): $\mathrm{RMSD} = 1.9661$ Å (Benchmark: ~1.966 Å).
   - TEL-patch per-residue displacements:
     * Leu183: 0.5119 Å (`Rigid Anchor`, $<1.0$ Å)
     * Glu168: 1.2534 Å (`Moderate Shift`, 1.0–2.5 Å)
     * Glu169: 3.8093 Å (`Extensive Remodeling`, $>2.5$ Å)
     * Lys170: 4.2929 Å (`Extensive Remodeling`, peak loop displacement)
     * Glu171: 3.3204 Å (`Extensive Remodeling`)
   - 2I46 Dimer Asymmetry: Confirmed Chain A has all 6 TEL-patch residues modeled, while Chain B is missing residues 168–171. Code strictly selects Chain A.

### 1.3 Test Suite & Quality Execution
1. **Milestone 4 Tests (`tests/test_cross_state.py`)**:
   - Command: `pytest tests/test_cross_state.py -v`
   - Result: `29 passed in 11.67s` (exit code 0).
2. **Full Repository Test Suite (`tests/`)**:
   - Command: `pytest tests/ -v`
   - Result: `257 passed, 6 skipped in 183.06s` (exit code 0).
   - The 6 skipped tests are Tier 4 end-to-end scenario placeholders reserved for Milestone 5 (`report.md` presence and `run_pipeline.py` script runner).
3. **Lint Check**:
   - Command: `ruff check src/ tests/`
   - Result: `All checks passed!`.

---

## 2. Logic Chain

1. **Integrity Mode Specification**:
   - Observation: `ORIGINAL_REQUEST.md` line 13 specifies `Integrity mode: development`.
   - Implication: Focus is to catch fabricated outputs, facade implementations, hardcoded test values, and mocked test runs.

2. **Absence of Hardcoded or Facade Logic**:
   - Observation: Static search confirms `BENCHMARK_RMSDS` is an unreferenced informational dict. All functions execute genuine mathematical routines using `gemmi` and NumPy SVD.
   - Deduction: The implementation contains no facade stubs or hardcoded lookup shortcuts.

3. **Verification of Mathematical Authenticity**:
   - Observation: Independent testing proved `svd_kabsch_align` computes exact rotations and translations minimizing RMSD, enforces proper right-handed rotations in $SO(3)$, and responds dynamically to coordinate perturbations and noise.
   - Deduction: The SVD Kabsch engine and coordinate extraction pipeline are fully authentic and functional.

4. **Consistency of Deliverables**:
   - Observation: Re-generating the standardized CSV tables yields exact frame-for-frame equality with the tables in `results/tables/`.
   - Deduction: Deliverable tables are authentic products of the code and not pre-fabricated artifacts.

---

## 3. Caveats

- **Scope Boundary**: This audit evaluated Milestone 4 (R4) artifacts and regression stability across Milestones 1–3. Tier 4 end-to-end integration tests (6 skipped tests) remain pending until Milestone 5 creates `report.md` and `run_pipeline.py`.
- **Static Coordinate Context**: Displacements and RMSDs reflect static cryo-EM models (3.3–3.5 Å) and X-ray coordinates (2.70 Å).

---

## 4. Conclusion

**Verdict: CLEAN**.

Milestone 4 (Cross-State Comparative Structural Analysis & Robustness Testing) satisfies all technical, architectural, and integrity constraints without violations:
- SVD Kabsch algorithm is mathematically sound, enforces $SO(3)$ chirality preservation, and handles boundary conditions.
- Cross-state RMSDs and per-residue displacement profiles match empirical benchmarks authentically.
- 2I46 Chain A is safely disambiguated from disordered Chain B.
- All 29 Milestone 4 tests and all 257 project-wide tests pass cleanly.
- Code style and linting are 100% compliant.

Milestone 4 is approved for progression to Milestone 5.

---

## 5. Verification Method

To independently reproduce this audit:

1. **Run Milestone 4 Unit Test Suite**:
   ```bash
   pytest tests/test_cross_state.py -v
   ```
   *Expected outcome*: `29 passed` with 0 failures.

2. **Run Full Test Suite**:
   ```bash
   pytest tests/ -v
   ```
   *Expected outcome*: `257 passed, 6 skipped`.

3. **Run Code Quality Check**:
   ```bash
   ruff check src/ tests/
   ```
   *Expected outcome*: `All checks passed!`.

4. **Run Dynamic Perturbation & Table Parity Verification**:
   ```bash
   python -c "
   import tempfile, pandas as pd
   from pathlib import Path
   from telotrade.cross_state import generate_cross_state_tables, svd_kabsch_align, extract_ca_coordinates, resolve_cif_path
   import numpy as np

   # Verify dynamic coordinate response
   p = np.random.randn(50, 3)
   r_base = svd_kabsch_align(p, p)
   assert r_base.rmsd == 0.0
   p_noisy = p + np.random.randn(*p.shape) * 1.0
   assert svd_kabsch_align(p, p_noisy).rmsd > 0.5

   # Verify table generation parity
   with tempfile.TemporaryDirectory() as tmpdir:
       gen = generate_cross_state_tables(output_dir=tmpdir)
       for name, pth in gen.items():
           ref = Path('results/tables') / name
           pd.testing.assert_frame_equal(pd.read_csv(pth), pd.read_csv(ref))
   print('Verification complete: 100% CLEAN')
   "
   ```
