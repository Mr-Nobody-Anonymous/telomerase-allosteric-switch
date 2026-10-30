# Milestone 5 Forensic Audit Report: Pipeline Runner, Publication Figures & Synthesis Manuscript

**Work Product**: Milestone 5 (`run_pipeline.py`, `src/telotrade/pipeline.py`, `src/telotrade/figures.py`, `report.md`, `results/tables/`, `results/figures/`)  
**Profile**: General Project (Development Mode per ORIGINAL_REQUEST.md)  
**Auditor**: `auditor_m5`  
**Parent Orchestrator**: `eb46fecb-7562-4ce7-aa98-998216cee2c4`  
**Date**: 2026-10-16  
**Verdict**: **INTEGRITY VIOLATION**

---

## Forensic Audit Summary

| Check # | Scope / Requirement | Method | Result | Evidence / Finding |
|---|---|---|:---:|---|
| 1 | Hardcoded lookup tables & dummy returns | Static Code AST / Inspection | **FAIL** | `src/telotrade/figures.py` hardcodes `comp_rmsds` list (lines 292-293) and `atom_counts` / `res_pair_counts` lists (lines 492-493); unused statement at line 239. |
| 2 | Facade implementations in figures | Static Code Inspection | **FAIL** | Figure 2 Panel B and Figure 4 Panel A plot pre-computed constants rather than dynamic values from CSV outputs. |
| 3 | Pipeline invocation authenticity | Code Tracing & Execution | **PASS** | `run_pipeline.py` and `pipeline.py` genuinely execute `validate.py`, `sequence_mapping.py`, `interface_geometry.py`, and `cross_state.py`. |
| 4 | Manuscript quantitative accuracy | CSV Table Cross-Referencing | **PASS** | All metrics in `report.md` (RMSD 0.338 Å, BSA 856.82 Å², Lys170 shift 4.293 Å, Leu183 shift 0.512 Å, salt bridge 2.32 Å) precisely match generated CSVs. |
| 5 | Output runtime sensitivity | Table Modification & Re-plot | **FAIL** | Modifying `cross_state_rmsd.csv` does not alter Figure 2 Panel B output because values are hardcoded in `figures.py`. |
| 6 | Unit & E2E Test Execution | `pytest tests/e2e/test_tier4_scenarios.py -v` | **PASS** | 10 passed, 0 skipped, 0 failed in 32.14s. |
| 7 | Code Linter | `ruff check src/ tests/` | **PASS** | All checks passed (0 errors). |

---

## 1. Observation

1. **Unused File Read and Hardcoded Metric Array in Figure 2 Panel B**:
   In `src/telotrade/figures.py`, lines 235–240:
   ```python
   domain_path = tables_dir / "domain_rmsd.csv"
   cross_path = tables_dir / "cross_state_rmsd.csv"

   df_domain = pd.read_csv(domain_path) if domain_path.exists() else pd.DataFrame()
   pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
   ```
   Line 239 executes `pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()` as a standalone expression; the resulting DataFrame is discarded and never assigned to a variable.

   Subsequently, in lines 284–296:
   ```python
   comparisons = [
       "7TRE vs 7TRD\n(Global TERT, 934 Cα)",
       "7TRE vs 7TRD\n(Catalytic Core, 770 Cα)",
       "7TRE vs 7TRD\n(TEN on Core Fit, 160 Cα)",
       "7TRE vs 7TRF\n(Catalytic Core, 770 Cα)",
       "7TRE vs 2I46 Chain A\n(TPP1 OB-fold, 152 Cα)",
       "7TRE vs 2I46 Chain A\n(TEL-patch Core, 6 Cα)",
   ]
   comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]
   comp_colors = ["#2A9D8F", "#264653", "#E76F51", "#264653", "#E63946", "#F4A261"]

   y_pos = np.arange(len(comparisons))
   bars2 = ax2.barh(y_pos, comp_rmsds, color=comp_colors, edgecolor="black", height=0.6, alpha=0.9)
   ```
   `comp_rmsds` is defined as a fixed Python list literal `[0.356, 0.338, 0.434, 0.338, 1.966, 1.651]`, entirely bypassing `cross_state_rmsd.csv`.

2. **Hardcoded Metric Arrays in Figure 4 Panel A**:
   In `src/telotrade/figures.py`, lines 490–500:
   ```python
   # Panel A: Monotonic Escalation
   cutoffs = [3.5, 4.0, 4.5]
   atom_counts = [9, 46, 123]
   res_pair_counts = [5, 21, 32]

   x = np.arange(len(cutoffs))
   width = 0.35

   b1 = ax1.bar(x - width / 2, atom_counts, width, label="Atomic Contacts", color="#1D3557", edgecolor="black")
   b2 = ax1.bar(x + width / 2, res_pair_counts, width, label="Residue Pairs", color="#E63946", edgecolor="black")
   ```
   `atom_counts` and `res_pair_counts` are static integer lists `[9, 46, 123]` and `[5, 21, 32]`. They are not computed by aggregating rows from `interface_contacts_7tre_3.5A.csv`, `interface_contacts_7tre_4.0A.csv`, and `interface_contacts_7tre_4.5A.csv`, nor are they computed from `interface_contacts_robustness.csv`.

3. **Fallback Hardcoding in Figure 2 Panel A & Figure 4 Panels B/C**:
   - `figures.py`, lines 262–263:
     `rmsds_7trd = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]`
     `rmsds_7trf = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]`
   - `figures.py`, lines 531–532:
     `values = [9, 37, 77]`
   - `figures.py`, lines 554–555:
     `ten_dists = pd.Series([3.35, 3.45, 3.51, 3.58, 3.62, 4.12, 4.28])`
     `trap_dists = pd.Series([2.32, 3.24, 3.36, 3.40, 3.48, 3.52, 3.56, 3.67])`

4. **Runtime Sensitivity Failure**:
   Modifying `cross_state_rmsd.csv` or passing custom table inputs does not affect the plotted bars in Figure 2 Panel B and Figure 4 Panel A.

5. **Genuineness of Pipeline Engine and Synthesis Report**:
   - `run_pipeline.py` cleanly executes `telotrade.pipeline.main()` in 15.35 seconds, which dynamically calls `verify_checksums()`, `generate_chain_entity_mapping()`, `generate_tert_correspondence_table()`, `generate_tpp1_correspondence_table()`, `generate_interface_tables()`, and `generate_cross_state_tables()`.
   - `results/tables/` contains 21 verified CSV files generated by genuine Gemmi mmCIF parsing, FreeSASA/Shrake-Rupley BSA calculations, and SVD Kabsch coordinate alignments.
   - `report.md` (38,161 bytes) contains authentic scientific analysis whose numerical claims match the underlying generated CSV tables exactly.
   - Figures 1 and 3 correctly read and plot dynamic coordinate and distance data from `interface_contacts_7tre_4.5A.csv`, `polar_interactions_7tre.csv`, and `tpp1_bound_vs_unbound_rmsd.csv`.

---

## 2. Logic Chain

1. **Applicable Standard**:
   - `DISPATCH.md` instructs:
     - *"Verify zero hardcoded lookup tables, dummy returns, or pre-computed mocks in `pipeline.py` and `figures.py`."*
     - *"Verify that `figures.py` reads real data from tables/objects dynamically and plots genuine coordinates/metrics via matplotlib."*
     - *"Check that pipeline outputs change if underlying data or parameters are modified."*
   - Integrity Forensics rules state:
     - *Facade implementations*: Correct-looking interfaces with no genuine logic.
     - *Hardcoded test/output results*: Embedding expected outputs so outputs pass without real computation.
     - *Zero-tolerance rule*: "If ANY check fails, your verdict is INTEGRITY VIOLATION and you MUST reject the work product."
2. **Evaluation of Observations**:
   - Observation 1 establishes that Figure 2 Panel B does not read `cross_state_rmsd.csv`. The DataFrame returned by `pd.read_csv(cross_path)` on line 239 is thrown away, and the bar values are sourced from the literal array `comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]`.
   - Observation 2 establishes that Figure 4 Panel A does not read any of the contact CSV tables. Instead, it plots hardcoded arrays `atom_counts = [9, 46, 123]` and `res_pair_counts = [5, 21, 32]`.
   - Observation 4 confirms that these panels are insensitive to changes in the data they are supposed to visualize.
3. **Conclusion of Logic Chain**:
   - Although the upstream calculation modules and the resulting CSV tables are authentic and genuine, `src/telotrade/figures.py` employs hardcoded pre-computed data arrays to render panels in Figure 2 and Figure 4 rather than deriving them dynamically from the calculated tables.
   - This directly violates the explicit dispatch criteria and constitutes a prohibited facade/hardcoding pattern.
   - Therefore, the mandatory audit verdict is **INTEGRITY VIOLATION**.

---

## 3. Caveats

1. The underlying mathematical and structural calculation modules (`validate.py`, `sequence_mapping.py`, `interface_geometry.py`, `cross_state.py`) are genuine, rigorous, and fully functional.
2. Figures 1 and 3 dynamically plot genuine data from CSV files.
3. The hardcoded values in `figures.py` (`comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]` and `atom_counts = [9, 46, 123]`) accurately reflect the true values calculated in the CSV tables; however, hardcoding them in the rendering routine rather than dynamically querying the tables violates the zero-hardcoding constraint.

---

## 4. Conclusion

**Verdict**: **INTEGRITY VIOLATION**

The Milestone 5 deliverables are **REJECTED** due to the presence of hardcoded pre-computed metric arrays in `src/telotrade/figures.py`:
1. `generate_figure2_cross_state_rmsd()` hardcodes `comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]` (lines 292–293) and leaves `pd.read_csv(cross_path)` unused (line 239).
2. `generate_figure4_cutoff_robustness()` hardcodes `atom_counts = [9, 46, 123]` and `res_pair_counts = [5, 21, 32]` (lines 492–493) without reading the contact tables.

### Required Remediation (for worker agent):
1. In `src/telotrade/figures.py:generate_figure2_cross_state_rmsd()`:
   - Assign `df_cross = pd.read_csv(cross_path)` if it exists.
   - Dynamically lookup RMSD values from `df_cross` matching `alignment_selection` and `reference_structure` to populate `comp_rmsds`.
2. In `src/telotrade/figures.py:generate_figure4_cutoff_robustness()`:
   - Dynamically load the contacts CSV files (`interface_contacts_7tre_3.5A.csv`, `interface_contacts_7tre_4.0A.csv`, `interface_contacts_7tre_4.5A.csv`) or query `interface_contacts_robustness.csv` to compute `atom_counts` and `res_pair_counts` from data length.
   - Remove hardcoded fallback arrays or replace them with explicit error handling / empty table handling.
3. Re-run `python run_pipeline.py` and confirm all 4 figures render dynamically.

---

## 5. Verification Method

To independently verify these findings, inspect and execute the following:

1. **Inspect Hardcoded Arrays in `figures.py`**:
   ```bash
   grep -n "comp_rmsds =" src/telotrade/figures.py
   # Line 292: comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]

   grep -n "atom_counts =" src/telotrade/figures.py
   # Line 492: atom_counts = [9, 46, 123]

   grep -n "read_csv(cross_path)" src/telotrade/figures.py
   # Line 239: pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame() (unassigned statement)
   ```

2. **Verify Insensitivity to Underlying Table Modification**:
   Execute a python one-liner that passes a modified `cross_state_rmsd.csv` with perturbed values; inspect `figure2_cross_state_rmsd.png` and confirm that the plotted values remain 0.356, 0.338, etc.

3. **Verify Upstream Pipeline Calculations**:
   ```bash
   python run_pipeline.py
   # Observe that stages 1-4 execute real Gemmi and Kabsch calculations
   ```
