# Dispatch: Milestone 5 Remediation Explorer 1 (explorer_m5_fix_1)

**Target**: Remediation Strategy for Figure 2 in `src/telotrade/figures.py`  
**Parent Orchestrator**: eb46fecb-7562-4ce7-aa98-998216cee2c4  
**Date**: 2026-10-16  

---

## Forensic Auditor Full Evidence Report
The Forensic Auditor (`auditor_m5`) delivered an **INTEGRITY VIOLATION** verdict with the following evidence:
- In `src/telotrade/figures.py`, lines 235–240:
  `cross_path = tables_dir / "cross_state_rmsd.csv"`
  `pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()` is executed as an unused expression statement without assignment.
- In lines 292–293:
  `comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]` is hardcoded as a literal list.
- Fallback arrays are embedded on lines 262–263.

Full audit report: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md`.

---

## Explorer Objective
1. Inspect `src/telotrade/figures.py:generate_figure2_cross_state_rmsd()` and `results/tables/cross_state_rmsd.csv`, `domain_rmsd.csv`.
2. Design a concrete code patch that:
   - Assigns `df_cross = pd.read_csv(cross_path)` if it exists.
   - Dynamically extracts RMSD values from `df_cross` and `df_domain` for each comparison bar.
   - Eliminates hardcoded `comp_rmsds = [...]`.
   - Eliminates hardcoded fallback arrays on lines 262–263.
3. Recommend exact replacement lines for Worker.
4. Author a handoff report at `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_1\handoff.md` and report to parent.


## 2026-10-16T19:06:43Z
You are Milestone 5 Remediation Explorer 1.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_1
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read auditor_m5 full handoff report at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_1\DISPATCH.md

Objective:
Analyze src/telotrade/figures.py lines 200–350 (generate_figure2_cross_state_rmsd) and formulate an exact, dynamic implementation plan to replace the hardcoded comp_rmsds list and unassigned pd.read_csv statement with genuine DataFrame queries.
Do NOT modify code yourself. Write your recommendations in handoff.md and send message to parent.
