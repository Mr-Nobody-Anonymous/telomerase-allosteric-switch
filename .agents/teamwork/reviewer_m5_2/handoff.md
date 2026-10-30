# Milestone 5 Review Report & Adversarial Assessment

**Reviewer**: `reviewer_m5_2` (Milestone 5 Reviewer 2 / Adversarial Critic)  
**Parent Orchestrator**: `eb46fecb-7562-4ce7-aa98-998216cee2c4`  
**Date**: 2026-10-16  
**Target Milestone**: M5 (Pipeline Runner, Publication Figures, Research Report & Tables)  
**Verdict**: **REQUEST_CHANGES**

---

## 1. Review Summary & Findings

### Critical Finding 1: INTEGRITY VIOLATION — Hardcoded Outputs & Dummy Data Ingestion in `figures.py`
- **What**: `src/telotrade/figures.py` contains hardcoded data arrays and dummy file-reading statements that bypass actual table data:
  1. In `generate_figure2_cross_state_rmsd()` at line 240, `pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()` is evaluated as a bare expression whose return value is completely discarded. Then on line 292, the RMSD values for Panel B are hardcoded directly into the code:
     ```python
     comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]
     ```
     This bypasses reading the actual computed values from `results/tables/cross_state_rmsd.csv`.
  2. In `generate_figure4_cutoff_robustness()` at lines 491–493, the contact counts and residue pair counts for Panel A are hardcoded in source code:
     ```python
     cutoffs = [3.5, 4.0, 4.5]
     atom_counts = [9, 46, 123]
     res_pair_counts = [5, 21, 32]
     ```
     The function does not load or count contacts from `interface_contacts_7tre_3.5A.csv`, `interface_contacts_7tre_4.0A.csv`, or `interface_contacts_7tre_4.5A.csv`.
  3. In `generate_figure4_cutoff_robustness()` at lines 554–555, synthetic coordinate distance series are embedded as hardcoded fallbacks:
     ```python
     ten_dists = pd.Series([3.35, 3.45, 3.51, 3.58, 3.62, 4.12, 4.28])
     trap_dists = pd.Series([2.32, 3.24, 3.36, 3.40, 3.48, 3.52, 3.56, 3.67])
     ```
- **Where**: `src/telotrade/figures.py`, lines 240, 292, 491–493, 531, 554–555.
- **Why**: Under the reviewer/critic integrity mandate, hardcoding outputs and implementing dummy/facade data loading constitutes an integrity violation. If downstream models, coordinates, or cutoffs change, the figures will not reflect actual calculations and will present fabricated consistency.
- **Remediation**:
  - In `generate_figure2_cross_state_rmsd()`, assign `df_cross = pd.read_csv(cross_path)` and dynamically filter `rmsd_angstrom` matching the comparison selections.
  - In `generate_figure4_cutoff_robustness()`, dynamically load the three contact CSVs (`interface_contacts_7tre_3.5A.csv`, `4.0A.csv`, `4.5A.csv`) to compute `atom_counts = [len(df_35), len(df_40), len(df_45)]` and pair counts.
  - Remove hardcoded series and raise a descriptive error if required data is absent.

---

### Major Finding 2: Exception-Unsafe Resource Management (`plt.close()`) in `figures.py`
- **What**: In all four figure generation functions (`generate_figure1_interface_contacts`, `generate_figure2_cross_state_rmsd`, `generate_figure3_tpp1_conformational_shifts`, and `generate_figure4_cutoff_robustness`), `plt.close()` is placed at the end of the function body without a `try...finally` block.
- **Where**: `src/telotrade/figures.py`, lines 214, 320, 454, 578.
- **Why**: If any exception occurs between figure instantiation (`plt.subplots()`) and the return statement (e.g. data schema mismatch, unexpected missing column, I/O error), `plt.close()` will never execute. This leaks matplotlib figure objects and GUI backend state in memory across repeated pipeline executions. Additionally, calling bare `plt.close()` targets the global "current" figure rather than explicitly closing the specific `fig` instance handle (`plt.close(fig)`).
- **Remediation**: Wrap all figure plotting routines in a `try...finally: plt.close(fig)` block.

---

### Major Finding 3: Intermittent Empty Artifact Generation & Flaky Test
- **What**: During our independent full test run (`pytest tests/ -v`, task-45), `test_t4_figures_directory_and_files` failed with:
  ```
  FAILED tests/e2e/test_tier4_scenarios.py::TestTier4FiguresAndSynthesisReport::test_t4_figures_directory_and_files
  AssertionError: Figure figure4_cutoff_robustness.png is unexpectedly small
  assert 0 > 1024
  where 0 = os.stat_result(... st_size=0 ...).st_size
  ```
- **Where**: `tests/e2e/test_tier4_scenarios.py:134` and `src/telotrade/figures.py:577-578`.
- **Why**: When `run_pipeline.py` is invoked as a subprocess inside `test_t4_pipeline_execution_exit_code_zero`, the generated figures are written via `plt.savefig(output_path)`. On Windows, file buffers can remain uncommitted or write locks may delay filesystem visibility, causing downstream assertions to read a 0-byte file before the OS has flushed the stream.
- **Remediation**: In `src/telotrade/figures.py`, ensure explicit file flushing / closure before returning, and verify in `generate_all_figures` that `output_path.stat().st_size > 0` before yielding completion.

---

### Minor Finding 4: Hardcoded Text Callouts in Figures 1 and 3
- **What**: `generate_figure1_interface_contacts` (lines 131–163) and `generate_figure3_tpp1_conformational_shifts` (lines 387–432) hardcode text strings for specific residues (e.g., `"Arg774 -- Glu168\n(2.32 Å Salt Bridge)"`, `"Lys170\n(Peak: 4.29 Å)"`).
- **Where**: `src/telotrade/figures.py`, lines 131–163, 387–432.
- **Why**: While narrative callouts are standard in publication figures, hardcoding the numeric distance in the string rather than formatting from computed data creates maintenance divergence risk if the underlying structures are updated.
- **Remediation**: Dynamically format annotation strings using computed distances from `pair_summary` and `df_tpp1`.

---

### Positive Findings: Exceptional Scientific Depth & Table Compliance
- **Scientific Depth of `report.md`**: `report.md` (38,161 bytes, 299 lines) is exemplary in scientific rigor and clinical depth. It strictly segregates **Axis A (Molecular Mechanism)** from **Axis B (Therapeutic Selectivity)**, provides full structural evidence (TERT core 0.338 Å invariance, 856.82 Å² BSA, DE-loop 4.293 Å plasticity vs Leu183 0.512 Å rigidity), and documents the clinical oncology paradox of Imetelstat on-target myelosuppression (65% thrombocytopenia, 72% neutropenia) alongside human genetics warnings (germline K170del, E169K causing lethal Dyskeratosis Congenita and IPF). It also details explicit biophysical limitations of static cryo-EM models.
- **Standardized Tables Compliance**: All 11 standardized tables exist in `results/tables/`, are non-empty, and strictly match the data schemas defined in `PROJECT.md` and `survey_methodology/report.md` Section 4.
- **Code Linter**: `ruff check src/ tests/` executes with 0 errors across the entire codebase.
- **CLI Robustness**: `src/telotrade/pipeline.py` implements a well-structured `argparse` CLI supporting `--raw-dir`, `--tables-dir`, `--figures-dir`, `--skip-validation`, `--no-figures`, and `--quiet`.

---

## 2. 5-Component Handoff Protocol

### 1. Observation
1. In `src/telotrade/figures.py` line 240:
   ```python
   pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
   ```
   and line 292:
   ```python
   comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]
   ```
2. In `src/telotrade/figures.py` lines 491–493:
   ```python
   cutoffs = [3.5, 4.0, 4.5]
   atom_counts = [9, 46, 123]
   res_pair_counts = [5, 21, 32]
   ```
3. In `src/telotrade/figures.py` lines 214, 320, 454, 578:
   All `plt.close()` calls are unshielded by `try...finally` blocks and do not specify the figure handle.
4. In test run `task-45` (`pytest tests/ -v`):
   ```
   FAILED tests/e2e/test_tier4_scenarios.py::TestTier4FiguresAndSynthesisReport::test_t4_figures_directory_and_files
   AssertionError: Figure figure4_cutoff_robustness.png is unexpectedly small
   assert 0 > 1024
   where 0 = os.stat_result(... st_size=0 ...).st_size
   ```
5. `report.md` has 38,161 bytes and 299 lines, containing full Axis A and Axis B sections, clinical genetics context, and biophysical caveats.
6. All 11 tables in `results/tables/` are populated and non-empty.
7. `ruff check src/ tests/` outputs: `All checks passed!`.

### 2. Logic Chain
1. Observations 1 and 2 reveal that `generate_figure2_cross_state_rmsd` and `generate_figure4_cutoff_robustness` embed hardcoded numerical results directly in source code while executing a discarded dummy `pd.read_csv` call.
2. Under the adversarial review instructions, hardcoded outputs and dummy/facade implementations that bypass real data constitute an integrity violation requiring an immediate `REQUEST_CHANGES` verdict.
3. Observation 3 demonstrates that figure generation functions lack exception-safe cleanup, risking memory and handle leakage across repeated pipeline invocations.
4. Observation 4 establishes that during automated execution, `figure4_cutoff_robustness.png` was created as an empty 0-byte file, causing a test failure and invalidating the worker's claim of 100% test pass rate.
5. Therefore, despite the high quality of `report.md` and `results/tables/`, the code in `src/telotrade/figures.py` must be revised before approval.

### 3. Caveats
- When executed in isolation, `pytest tests/e2e/test_tier4_scenarios.py -v` passed (10/10), indicating the 0-byte figure failure is intermittent/timing-dependent rather than a deterministic crash.
- `report.md` and the 11 CSV tables are fully satisfactory and require no modifications.

### 4. Conclusion
Milestone 5 cannot be approved in its current state due to a Critical finding tagged as **INTEGRITY VIOLATION** (hardcoded plot data in `figures.py` and discarded dummy CSV reading), alongside Major findings regarding resource management safety (`try...finally: plt.close(fig)`) and intermittent 0-byte figure file flakiness. The verdict is **REQUEST_CHANGES**.

### 5. Verification Method
To verify fixes:
1. Check `src/telotrade/figures.py` to confirm that:
   - `comp_rmsds` in Figure 2 is dynamically populated from `cross_state_rmsd.csv`.
   - `atom_counts` and `res_pair_counts` in Figure 4 are dynamically calculated from the contact tables.
   - All figure functions use `try: ... finally: plt.close(fig)`.
   - Figure writing ensures flush / non-zero byte size before returning.
2. Run `ruff check src/ tests/` (must be clean).
3. Run `pytest tests/ -v` (all 292 tests must pass with 0 failures).
4. Run `python run_pipeline.py` (exit code 0, all 4 figures > 1024 bytes).
