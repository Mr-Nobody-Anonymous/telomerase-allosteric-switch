# Dispatch: Milestone 5 Remediation Explorer 3 (explorer_m5_fix_3)

**Target**: Remediation Strategy for Resource Cleanup, Atomicity, and Figure 1/3 Audits in `src/telotrade/figures.py`  
**Parent Orchestrator**: eb46fecb-7562-4ce7-aa98-998216cee2c4  
**Date**: 2026-10-16  

---

## Forensic Auditor & Reviewer Evidence Reports
- Forensic Auditor (`auditor_m5`): Verdict **INTEGRITY VIOLATION** (`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md`).
- Reviewer 2 (`reviewer_m5_2`): Verdict **REQUEST_CHANGES** (`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\reviewer_m5_2\handoff.md`):
  * Figure generation functions lack exception-safe resource cleanup (`try...finally: plt.close(fig)`), risking memory/handle leaks.
  * In test suite runs, `figure4_cutoff_robustness.png` was created as an empty 0-byte file due to unhandled file saving or buffering issues.

---

## Explorer Objective
1. Inspect all 4 figure generator functions in `src/telotrade/figures.py` (`generate_figure1_interface_contacts`, `generate_figure2_cross_state_rmsd`, `generate_figure3_tpp1_conformational_shifts`, `generate_figure4_cutoff_robustness`).
2. Audit Figure 1 and Figure 3 for any covert hardcoding or synthetic fallback lists.
3. Design a pattern using `try...finally: plt.close(fig)` across all 4 functions to ensure zero memory or file-handle leaks.
4. Ensure atomic saving of image files (e.g., flush buffers before close) to guarantee non-zero file sizes (>1024 bytes).
5. Recommend exact code patterns for Worker.
6. Author a handoff report at `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_3\handoff.md` and report to parent.


## 2026-10-16T19:06:43Z
You are Milestone 5 Remediation Explorer 3.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_3
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read reviewer_m5_2 handoff report at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\reviewer_m5_2\handoff.md
Read auditor_m5 handoff report at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m5\handoff.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_3\DISPATCH.md

Objective:
Analyze all figure generators in src/telotrade/figures.py for resource leaks and buffering issues (0-byte file generation). Formulate a clean try...finally: plt.close(fig) design and audit Figures 1 and 3 for any covert hardcoded data.
Do NOT modify code yourself. Write your recommendations in handoff.md and send message to parent.
