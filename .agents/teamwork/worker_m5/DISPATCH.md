# Dispatch: Milestone 5 Implementation Worker (worker_m5)

**Target Milestone**: M5 (Evidence Matrix, Publication Figures, Report & Automated Pipeline Runner)  
**Parent Orchestrator**: eb46fecb-7562-4ce7-aa98-998216cee2c4  
**Date**: 2026-10-16  

---

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

---

## Mission Objectives
1. Implement `src/telotrade/pipeline.py` and top-level CLI `run_pipeline.py`.
   - Single-entry script executable via `python run_pipeline.py [--raw-dir DIR] [--output-dir DIR]`.
   - Re-executes the full pipeline cleanly and deterministically: validation, sequence mapping, interface geometry, cross-state RMSD, table generation, and figure generation.
   - Ensures exit code 0 under standard execution.
2. Generate all 4 publication-quality figures in `results/figures/` (PNG and/or SVG, >= 1024 bytes each):
   - `figure1_interface_contacts.png`: Interface contact map and polar/hydrophobic interaction network for 7TRE TERT-TPP1.
   - `figure2_cross_state_rmsd.png`: Domain RMSD comparison and catalytic core rigidity across 7TRE, 7TRD, 7TRF, and 2I46.
   - `figure3_tpp1_conformational_shifts.png`: TPP1 Cα displacement profile highlighting DE-loop plasticity (Glu168–Glu171 peaking at Lys170 at 4.29 Å) vs Leu183 rigidity (0.51 Å).
   - `figure4_cutoff_robustness.png`: Contact count and persistence distribution across 3.5 Å, 4.0 Å, and 4.5 Å distance cutoffs.
3. Ensure all 11 standardized CSV tables in `results/tables/` are populated and validated:
   - `data_manifest.csv`
   - `chain_entity_mapping.csv`
   - `residue_correspondence_tert.csv`
   - `residue_correspondence_tpp1.csv`
   - `interface_contacts_7tre_3.5A.csv`
   - `interface_contacts_7tre_4.0A.csv`
   - `interface_contacts_7tre_4.5A.csv`
   - `interface_contacts_robustness.csv`
   - `buried_surface_area.csv`
   - `polar_interactions_7tre.csv`
   - `cross_state_rmsd.csv`
   (Plus auxiliary tables: `domain_rmsd.csv`, `per_residue_displacement_tert.csv`, `tpp1_bound_vs_unbound_rmsd.csv`).
4. Author comprehensive publication-grade research manuscript `report.md` (>5,000 bytes, ideally 15,000–25,000 bytes) at the project root:
   - Executive Summary & Background.
   - Structural Evidence Matrix.
   - Strict bifurcation between **Axis A (Molecular Mechanism of Recruitment)** and **Axis B (Therapeutic Selectivity & Stem Cell Trade-Off)**.
   - Axis A: High-resolution recruitment anatomy, 856.82 Å² BSA, bidentate Arg774-Glu168 salt bridge, Lys78-Glu215 clamp, Leu769/Pro771-Leu183/Leu212 hydrophobic anchor, rigid catalytic core (0.338 Å RMSD), plastic DE-loop remodeling (up to 4.29 Å).
   - Axis B: The telomerase paradox (cancer addiction vs adult stem cell exhaustion), active site inhibition (Imetelstat on-target hematological toxicity: thrombocytopenia, neutropenia), allosteric recruitment inhibition hypothesis, human genetic validation (germline TEL-patch mutations cause Dyskeratosis Congenita, Hoyeraal-Hreidarsson syndrome, aplastic anemia, and Idiopathic Pulmonary Fibrosis).
   - Biophysical caveats: static cryo-EM limitations, lack of explicit hydrogens, missing flexible linkers (179–321), absence of POT1/TIN2 shelterin context.
   - Translational validation roadmap: TRAP processivity assays, human CD34+ hematopoietic stem cell colony-forming unit (CFU) safety screens, humanized knock-in mouse models.
5. Verification:
   - Run `pytest tests/e2e/test_tier4_scenarios.py -v` — all 6 scenarios must PASS (0 skipped, 0 failed).
   - Run `pytest tests/ -v` — all 292 tests across unit, adversarial, and E2E suites must PASS 100%.
   - Run `ruff check src/ tests/` — 0 violations.
   - Author a 5-component `handoff.md` (Observation, Logic Chain, Caveats, Conclusion, Verification Method) in your working directory and notify the parent orchestrator.

## 2026-10-16T18:22:43Z
From: eb46fecb-7562-4ce7-aa98-998216cee2c4
Content: You are Milestone 5 Implementation Worker (worker_m5).
Working directory: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m5
Scope:
1. Implement src/telotrade/pipeline.py and top-level run_pipeline.py
2. Generate all 4 publication figures in results/figures/
3. Author comprehensive research manuscript report.md (>5,000 bytes) at project root
4. Run verification tests: pytest tests/e2e/test_tier4_scenarios.py -v, pytest tests/ -v, ruff check src/ tests/
5. Write 5-component handoff.md and notify parent orchestrator via send_message.
