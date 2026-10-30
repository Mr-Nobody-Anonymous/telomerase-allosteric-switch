## 2026-10-17T11:04:29Z
You are explorer_b2_recalc (role: Independent Recalculation Architect).
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b2_recalc

MANDATORY FIRST STEP:
Read C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\ORIGINAL_REQUEST.md
Pay special attention to the follow-up request from 2026-10-17T10:11:04Z regarding Requirement R4:
"R4. Independent Structural Recalculation & Stress-Testing (Phases 3 & 4)
Using an independent implementation (separate code paths), recalculate Cα RMSDs, buried surface area (BSA), and heavy-atom contacts across distance cutoffs (3.5 Å, 4.0 Å, 4.5 Å) and alternative alignment subsets. Deliver results/validation/independent_metrics.csv and results/interface/sensitivity_analysis.csv."

Also read:
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator_gen2\PROJECT.md
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\interface_geometry.py
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\cross_state.py
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\results\audit\evidence_ledger.csv

YOUR TASK:
1. Design an independent recalculation architecture using alternative libraries/algorithms that do NOT share code paths with the primary `telotrade` modules:
   - Alternative Cα RMSD: Use Gemmi superposition or Biopython SVD Superimposer with explicit coordinate extraction from raw mmCIF files (`data/raw/7TRE.cif`, `7TRD.cif`, `7TRF.cif`, `2I46.cif`).
   - Alternative BSA: Compare FreeSASA calculation against Biopython ShrakeRupley (or exact numerical integration), noting differences between Bondi and ProtOr radii.
   - Alternative Contact Networks: Direct Euclidean distance matrix or SciPy `cKDTree` / Biopython `NeighborSearch` at 3.5 Å, 4.0 Å, 4.5 Å.
2. Outline the exact implementation script (`src/telotrade/independent/recalculate.py` or similar) that will compute:
   - Primary vs Independent metrics comparison table: `results/validation/independent_metrics.csv` (Columns: `structure_id,metric_name,original_value,independent_value,delta,percent_diff,methodology_notes`).
3. Detail how to handle sequence offsets and chain mappings cleanly without duplicating previous errors.
4. Write your architecture plan to `report.md` and `handoff.md` in your directory, and send a message when finished.
