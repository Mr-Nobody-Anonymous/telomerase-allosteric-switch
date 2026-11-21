# Repository Inventory & Baseline Forensic Audit

**Document Designation**: R1 / Milestone B1 Deliverable  
**Investigation Lead**: `explorer_b1_repo` (Codebase & Artifact Auditor)  
**Date**: October 2026  
**Repository Path**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research`  
**Execution Context**: Accelerated Single-Run Pipeline (`telotrade`)  

### Status distribution after the no-lab audit (all 29 claims)

- **Established (19)**: the original 17, plus `CLM-027` (one of nine panel mutants survives Holm
  correction on processivity) and `CLM-029` (the estimator's window bias is real but bounded to
  phi >= 2%, established within the model only).
- **Uncertain (5)**: the original 4, plus `CLM-026` — the 2.34 A Tyr176-Gln794 contact is below
  hydrogen-bond distance at 3.5 A resolution, so it is a property of the deposited model rather than
  a measurement of the protein.
- **Contradicted (2)**: the original 1, plus `CLM-025` — Lys757 is 9.9-11.9 A from both Tyr176 OH and
  Gln794 OE1 in 7TRE/7TRD/7TRF, so the allosteric-lever premise is not supported by the coordinates.
- **Untested (3)**: the original 2, plus `CLM-028` — an equivalence claim on activity needs 8-93
  independent preparations and cannot be made from n = 3.

---

## 9. Identified Discrepancies & Recommendations for Phase 1–5

1. **Fix Sequence Correspondence Indexing Bug (RESOLVED IN MILESTONE B1)**:
   In `src/telotrade/sequence_mapping.py::generate_tpp1_correspondence_table`, the index offset was corrected from `ob_seq[resnum - 89]` to `ob_seq[resnum - 90]` (with residue 89 mapped from `full_poly[2]` Gly89), restoring 100% residue name alignment (res 166 is ASP, res 168 is GLU) and regenerating `results/tables/residue_correspondence_tpp1.csv`.
2. **Explicitly Contextualize the 2.32 Å Salt Bridge**:
   In `report.md` Section 4.2 and Section 6, explicitly note that 2.32 Å is an overcompacted coordinate distance from 3.5 Å cryo-EM refinement, not a physically validated ground-truth bond length.
3. **Decouple Static Coordinate RMSD from Dynamic Allostery**:
   Replace the claim that "0.338 Å proves absence of allostery" with a nuanced discussion separating static rigid-body invariance from dynamic/entropic allosteric coupling and repeat addition processivity (RAP).
4. **Quarantine TELIN-1 as a Purely Hypothetical In Silico Construct**:
   Ensure TELIN-1 is clearly designated as an unvalidated computational candidate with zero demonstrated affinity or cellular selectivity.
