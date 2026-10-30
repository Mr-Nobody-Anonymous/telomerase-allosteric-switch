# Dispatch: Milestone 5 Remediation Explorer 2 (explorer_m5_fix_2)

**Target**: Remediation Strategy for Figure 4 in `src/telotrade/figures.py`  
**Parent Orchestrator**: eb46fecb-7562-4ce7-aa98-998216cee2c4  
**Date**: 2026-10-16  

---

## Forensic Auditor Full Evidence Report
The Forensic Auditor (`auditor_m5`) delivered an **INTEGRITY VIOLATION** verdict with the following evidence:
- In `src/telotrade/figures.py`, lines 490–500:
  `atom_counts = [9, 46, 123]`
  `res_pair_counts = [5, 21, 32]`
  are hardcoded static integer lists instead of being computed from the contact CSV tables.
- Lines 531–532 (`values = [9, 37, 77]`) and lines 554–555 (`ten_dists = ...`, `trap_dists = ...`) embed synthetic fallback mock arrays.

Full audit report: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md`.

---

## Explorer Objective
1. Inspect `src/telotrade/figures.py:generate_figure4_cutoff_robustness()` and `results/tables/interface_contacts_7tre_*.csv`, `interface_contacts_robustness.csv`.
2. Design a concrete code patch that:
   - Dynamically reads `interface_contacts_7tre_3.5A.csv`, `interface_contacts_7tre_4.0A.csv`, `interface_contacts_7tre_4.5A.csv` or `interface_contacts_robustness.csv`.
   - Computes `atom_counts` and `res_pair_counts` dynamically from table row counts and unique `(tert_resnum, tpp1_resnum)` tuples.
   - Computes robustness class breakdown (Core Robust, Intermediate, Threshold Sensitive) dynamically from `interface_contacts_robustness.csv`.
   - Computes distance distributions for TEN and TRAP domains dynamically from `interface_contacts_7tre_4.5A.csv`.
   - Eliminates all hardcoded lists and mock fallback arrays.
3. Recommend exact replacement lines for Worker.
4. Author a handoff report at `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_2\handoff.md` and report to parent.


## 2026-10-16T19:06:43Z
You are Milestone 5 Remediation Explorer 2.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_2
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read auditor_m5 full handoff report at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_2\DISPATCH.md

Objective:
Analyze src/telotrade/figures.py lines 450–600 (generate_figure4_cutoff_robustness) and formulate an exact, dynamic implementation plan to replace atom_counts, res_pair_counts, and mock fallback arrays with dynamic queries to contact and robustness CSV files.
Do NOT modify code yourself. Write your recommendations in handoff.md and send message to parent.
