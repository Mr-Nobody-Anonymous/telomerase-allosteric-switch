## 2026-10-17T10:32:34Z
You are worker_b1 (role: Implementation & Verification Worker).
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_b1

MANDATORY FIRST STEP:
Read C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\ORIGINAL_REQUEST.md
Pay special attention to the follow-up request from 2026-10-17T10:11:04Z regarding Requirements R1 & R2:
"R1. Repository Inventory & Evidence Ledger (Phase 0)
Audit existing repository assets, code, tests, and tables. Generate results/audit/evidence_ledger.csv and results/audit/repository_inventory.md tracing every reported numerical claim (BSA, contact counts, salt bridge distances, RMSDs) to exact inputs, scripts, parameters, and limitations. Classify each claim as established, uncertain, contradicted, or untested."
"R2. Primary Literature & Novelty Gap Audit (Phase 1)
Conduct an exhaustive prior-art audit comparing claims against the original 7TRE publication (Liu et al., Nature 2022), Sekne et al. (Science 2022), and subsequent conformational dynamics studies. Treat recovery of Arg774–Glu168 and Lys78–Glu215 as positive controls, not discoveries. Generate results/literature/prior_art_matrix.csv, results/literature/claim_verification.md, and results/literature/novelty_gap_report.md."

Read C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator_gen2\PROJECT.md.
Also read the handoff reports from:
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b1_repo\handoff.md
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\spec_miner_b1_lit\handoff.md
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b1_methodology\handoff.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

EXCLUSIVE WRITE OWNERSHIP:
- `src/telotrade/sequence_mapping.py`
- `results/tables/residue_correspondence_tpp1.csv`
- `results/audit/evidence_ledger.csv`
- `results/audit/repository_inventory.md`
- `results/literature/prior_art_matrix.csv`
- `results/literature/claim_verification.md`
- `results/literature/novelty_gap_report.md`

YOUR TASKS:
1. Fix the audited sequence mapping defect in `src/telotrade/sequence_mapping.py` in `generate_tpp1_correspondence_table()`:
   Line ~1154: `ob_seq[resnum - 89]` assumes `ob_seq` begins at 89, whereas 7TRE Entity 4 starts at Gly90 (after stripping the 3-aa GAG tag). Fix the index calculation so `resnum_canonical` and `resname_canonical` match UniProt sequence without any offset phase shift (e.g., residue 166 should be ASP, not TRP).
2. Execute sequence mapping regeneration to update `results/tables/residue_correspondence_tpp1.csv`.
3. Verify that the 5 deliverables in `results/audit/` and `results/literature/`:
   - `results/audit/evidence_ledger.csv`
   - `results/audit/repository_inventory.md`
   - `results/literature/prior_art_matrix.csv`
   - `results/literature/claim_verification.md`
   - `results/literature/novelty_gap_report.md`
   are completely aligned with this fix and with all requirements in R1 and R2.
4. Run the full test suite (`python -m pytest tests/ -v`) and linter (`ruff check src/ tests/`) to ensure 100% pass (300+ tests pass, 0 fail).
5. Document all commands, test outputs, and changes in your handoff report at `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_b1\handoff.md`.
6. Send a message to your orchestrator when finished.
