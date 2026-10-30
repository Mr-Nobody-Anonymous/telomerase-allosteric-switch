# Dispatch: Milestone 5 Remediation Worker (worker_m5_fix)

**Target Milestone**: M5 Remediation (Iteration 2)  
**Parent Orchestrator**: eb46fecb-7562-4ce7-aa98-998216cee2c4  
**Date**: 2026-10-16  

---

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

---

## Remediation Objectives
Apply the complete dynamic data query patches developed by `explorer_m5_fix_1`, `explorer_m5_fix_2`, and `explorer_m5_fix_3` to `src/telotrade/figures.py`:

1. **Figure 2 (`generate_figure2_cross_state_rmsd`)**:
   - Assign `df_cross = pd.read_csv(cross_path)` if it exists.
   - Dynamically lookup RMSD values from `df_cross` and `df_domain` for each comparison bar.
   - Eradicate hardcoded `comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]`.
   - Eradicate hardcoded fallback arrays on lines 262–263.
   - Wrap generation in `try...finally: plt.close(fig)`.

2. **Figure 4 (`generate_figure4_cutoff_robustness`)**:
   - Dynamically compute `atom_counts` and `res_pair_counts` by reading `interface_contacts_7tre_3.5A.csv`, `interface_contacts_7tre_4.0A.csv`, `interface_contacts_7tre_4.5A.csv` or querying `interface_contacts_robustness.csv`.
   - Eradicate hardcoded lists `atom_counts = [9, 46, 123]` and `res_pair_counts = [5, 21, 32]`.
   - Derive persistence category counts dynamically from `df_rob["robustness_class"].value_counts()`. Eradicate hardcoded `.get()` defaults and fallback list `values = [9, 37, 77]`.
   - Derive distance distributions dynamically from contact tables; eradicate synthetic fallback Series `ten_dists` and `trap_dists`.
   - Wrap generation in `try...finally: plt.close(fig)`.

3. **Figures 1 & 3 Dynamic Annotations & Resource Cleanup**:
   - In Figure 1, replace hardcoded annotation strings (e.g. 2.32 Å, 3.58 Å) with dynamic lookups from contact and polar interaction DataFrames.
   - In Figure 3, dynamically query coordinates for TEL-patch residues (Glu168, Glu169, Lys170, Glu171, Leu183, Leu212) from `df_tpp1` instead of hardcoded coordinates.
   - Wrap Figure 1 and Figure 3 in `try...finally: plt.close(fig)`.
   - Implement safe atomic saving (flush buffer before close, verify file size > 1024 bytes) to eliminate any intermittent 0-byte file generation.

4. **Pipeline & Test Verification**:
   - Re-run `python run_pipeline.py` and verify all 4 figures render dynamically from CSV tables with exit code 0.
   - Run `pytest tests/e2e/test_tier4_scenarios.py -v` (assert 10 passed, 0 skipped, 0 failed).
   - Run `pytest tests/ -v` (assert all 292 tests pass, 0 failed, 0 skipped).
   - Run `ruff check src/ tests/` (assert 0 errors).
   - Author a 5-component `handoff.md` and report to the parent orchestrator via `send_message`.


## 2026-10-16T19:17:25Z
You are Milestone 5 Remediation Worker (worker_m5_fix).
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m5_fix
Scope:
1. In src/telotrade/figures.py:
   - Apply dynamic table loading to generate_figure2_cross_state_rmsd: assign df_cross = pd.read_csv(cross_path), lookup RMSDs dynamically for all comparison bars, eradicate hardcoded comp_rmsds list and fallback arrays.
   - Apply dynamic table loading to generate_figure4_cutoff_robustness: dynamically compute atom_counts and res_pair_counts from contact CSVs, derive robustness category counts dynamically from df_rob, eliminate hardcoded lists and mock fallback arrays.
   - Update Figure 1 & 3 annotations to query DataFrames dynamically rather than embedding static strings/coordinates.
   - Wrap all 4 figure generators in try...finally: plt.close(fig) for zero memory leaks.
   - Implement safe atomic figure saving to prevent intermittent 0-byte file generation.
2. Run python run_pipeline.py and verify all 4 figures render dynamically with exit code 0.
3. Run tests: pytest tests/e2e/test_tier4_scenarios.py -v and pytest tests/ -v. Run ruff check src/ tests/.
4. Author 5-component handoff.md and send completion message to parent.
