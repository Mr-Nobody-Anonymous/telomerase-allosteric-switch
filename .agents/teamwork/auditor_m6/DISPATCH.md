# Dispatch: auditor_m6

## 2026-10-16T19:55:00Z
Role: Milestone 6 Forensic Integrity Auditor (Project-Wide Anti-Cheat & Forensic Certification)
Agent Type: teamwork_preview_auditor
Working Directory: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\auditor_m6

Target Scope:
Conduct exhaustive project-wide static AST and dynamic integrity forensics across all files in:
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\ (validate.py, sequence_mapping.py, interface_geometry.py, cross_state.py, pipeline.py, figures.py)
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\run_pipeline.py
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\tests\ (unit, adversarial, e2e)

Forensic Integrity Checks (Zero-Tolerance):
1. Hardcoded outputs: Search AST and source code for static dictionaries, lists, or returns hardcoding benchmark outputs (RMSD, BSA, contact counts, residue names, Q-scores, figure arrays).
2. Facade/Dummy implementations: Verify every function performs genuine calculations (Biopython/Gemmi coordinate extraction, SVD Kabsch in SO(3), FreeSASA surface calculation, Euclidean distance matrices, dynamic pandas queries).
3. Test authenticity: Verify tests assert on genuine outputs and do not mock out production logic or assert against hardcoded dummy facades.
4. Dynamic sensitivity: Verify that perturbing input coordinates or input data alters computed results and generated figures proportionally.
5. Resource hygiene: Confirm no leaked matplotlib figure handles, open file descriptors, or 0-byte output files.

Deliver an authoritative binary verdict: CLEAN or INTEGRITY VIOLATION in handoff.md, with full supporting evidence, and report to parent orchestrator via send_message.
