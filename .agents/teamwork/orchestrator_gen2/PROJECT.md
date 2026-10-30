# Project: Human Telomerase Benchmark Audit & Falsifiable Mechanistic Evaluation

## Architecture
- **Audit & Evidence Infrastructure**: `src/telotrade/` and independent verification scripts in `src/telotrade/independent/`
- **Output Repositories**:
  - `results/audit/`: Evidence ledger (`evidence_ledger.csv`), repository inventory (`repository_inventory.md`)
  - `results/literature/`: Prior-art matrix (`prior_art_matrix.csv`), claim verification (`claim_verification.md`), novelty gap report (`novelty_gap_report.md`)
  - `results/research/`: Falsifiable research question and competing hypotheses (`question_and_hypotheses.md`)
  - `results/validation/`: Independently recalculated metrics (`independent_metrics.csv`)
  - `results/interface/`: Multi-cutoff parameter sensitivity analysis (`sensitivity_analysis.csv`)
  - `results/mechanism/`: DE-loop conformational analysis (`de_loop_mechanism_report.md`), allosteric coupling evidence assessment (`allostery_evidence_assessment.md`)
  - `results/ligand/`: Provenance and reality check of candidate ligand (`TELIN1_evidence_audit.md`)
  - `results/translation/`: Cancer-aging trade-off evidence matrix (`cancer_aging_evidence_matrix.csv`), therapeutic window analysis (`therapeutic_window_limits.md`)
  - `results/validation_experiments/`: Concrete wet-lab validation designs (`falsifiable_predictions.md`)
  - `results/`: Adversarial review report (`final_adversarial_review.md`), comprehensive publication-grade report (`report.md`)

## Feature Inventory
| # | Feature | Description | Milestone | Source |
|---|---------|-------------|-----------|--------|
| F1 | Repository Inventory & Git Audit | Detailed inventory of code, tests, structures, and git state | B1 | ORIGINAL_REQUEST R1 |
| F2 | Evidence Ledger | Tracing every reported numerical claim to inputs, scripts, parameters, and status | B1 | ORIGINAL_REQUEST R1 |
| F3 | Prior-Art Matrix | Systematic comparison vs Liu et al. (Nature 2022), Sekne et al. (Science 2022), etc. | B1 | ORIGINAL_REQUEST R2 |
| F4 | Claim Verification & Novelty Gap | Distinguishing positive controls from genuine computational insights | B1 | ORIGINAL_REQUEST R2 |
| F5 | Falsifiable Hypotheses Formulation | Competing hypotheses H1, H2, H3 with concrete falsification criteria | B2 | ORIGINAL_REQUEST R3 |
| F6 | Independent Recalculation Engine | Independent code path recalculating RMSDs, BSA, and contact networks | B2 | ORIGINAL_REQUEST R4 |
| F7 | Parameter Sensitivity Analysis | Multi-cutoff (3.5 Å, 4.0 Å, 4.5 Å) and alignment subset stress testing | B2 | ORIGINAL_REQUEST R4 |
| F8 | DE-Loop Mechanistic Evaluation | Testing state-dependent remodeling vs crystallographic/construct differences | B3 | ORIGINAL_REQUEST R5 |
| F9 | Allosteric Coupling Rigor Analysis | Decoupling static RMSD (0.338 Å) from dynamic/energetic allostery | B3 | ORIGINAL_REQUEST R5 |
| F10 | TELIN-1 Provenance & Reality Check | Critical audit of docked ligand, stereochemistry, affinity, and selectivity limits | B3 | ORIGINAL_REQUEST R6 |
| F11 | Cancer-Aging Evidence Matrix | Systematic matrix separating biological facts from clinical speculation | B4 | ORIGINAL_REQUEST R7 |
| F12 | Therapeutic Window Limits | Quantitative requirements for cancer-selective disruption vs stem cell maintenance | B4 | ORIGINAL_REQUEST R7 |
| F13 | Falsifiable Wet-Lab Predictions | Experimental assays, controls, and readouts for empirical hypothesis testing | B4 | ORIGINAL_REQUEST R7 |
| F14 | Independent Adversarial Review | Multi-perspective challenge of findings, methodology, and limitations | B4 | ORIGINAL_REQUEST R7 |
| F15 | 13-Section Report Rewrite | Complete report rewrite with rejected/overstated claims catalog | B5 | ORIGINAL_REQUEST R7 |

## Milestones
| # | Name | Scope | Dependencies | Status |
|---|------|-------|-------------|--------|
| B1 | Repository Inventory, Evidence Ledger & Prior-Art Audit | F1, F2, F3, F4 (`evidence_ledger.csv`, `repository_inventory.md`, `prior_art_matrix.csv`, `claim_verification.md`, `novelty_gap_report.md`) | none | DONE |
| B2 | Falsifiable Hypotheses & Independent Structural Recalculation | F5, F6, F7 (`question_and_hypotheses.md`, `independent_metrics.csv`, `sensitivity_analysis.csv`) | B1 | DONE |
| B3 | Mechanistic DE-Loop, Allostery & TELIN-1 Ligand Audit | F8, F9, F10 (`de_loop_mechanism_report.md`, `allostery_evidence_assessment.md`, `TELIN1_evidence_audit.md`) | B1, B2 | PLANNED |
| B4 | Translation, Wet-Lab Predictions & Adversarial Review | F11, F12, F13, F14 (`cancer_aging_evidence_matrix.csv`, `therapeutic_window_limits.md`, `falsifiable_predictions.md`, `final_adversarial_review.md`) | B1, B2, B3 | PLANNED |
| B5 | 13-Section Report Rewrite & Final Forensic Verification | F15 (`report.md` complete rewrite, full test pass, final gate & audit) | B1, B2, B3, B4 | PLANNED |

## Interface Contracts
### Evidence Audit ↔ Independent Recalculation
- `results/audit/evidence_ledger.csv`: Columns: `claim_id`, `metric_name`, `reported_value`, `source_file`, `source_script`, `parameters`, `evidence_status` (`established`, `uncertain`, `contradicted`, `untested`), `notes`
- `results/validation/independent_metrics.csv`: Columns: `structure_id`, `metric_name`, `original_value`, `independent_value`, `delta`, `percent_diff`, `methodology_notes`
- `results/interface/sensitivity_analysis.csv`: Columns: `contact_pair`, `chain_tert`, `res_tert`, `chain_tpp1`, `res_tpp1`, `distance_3_5A`, `distance_4_0A`, `distance_4_5A`, `persistent`, `sensitivity_class`

### Translation ↔ Final Report
- `results/translation/cancer_aging_evidence_matrix.csv`: Columns: `target_or_interface`, `cancer_dependency`, `stem_cell_dependency`, `therapeutic_window_evidence`, `experimental_status`, `citation`
- `report.md`: 13 mandatory sections:
  1. Executive Summary & Forensic Audit Overview
  2. Repository Inventory & Evidence Ledger Summary
  3. Prior-Art & Novelty Analysis (Liu et al. vs This Work)
  4. Falsifiable Hypotheses & Competing Models
  5. Primary Structural Data & Quality Manifest
  6. Quantitative Interface Geometry (Recalculated & Verified)
  7. Parameter Sensitivity & Coordinate Uncertainty Stress-Testing
  8. Cross-State Comparative Analysis & DE-Loop Remodeling
  9. Mechanistic Evaluation of Allostery (Static vs Dynamic Coupling)
  10. Candidate Ligand Reality Check (TELIN-1 Provenance & Limits)
  11. Translation: The Cancer vs. Healthy Aging Dilemma
  12. Falsifiable Wet-Lab Experimental Protocols
  13. Explicit Catalog of Overstated / Rejected Claims & Future Roadmap

## Code Layout
- `src/telotrade/`: Core pipeline modules (validated from Gen1)
- `src/telotrade/independent/`: Independent recalculation scripts using alternative libraries/parsing paths
- `results/audit/`: Audit ledgers and inventories
- `results/literature/`: Literature matrix and novelty reports
- `results/research/`: Hypotheses and research framing
- `results/validation/`: Recalculation validation tables
- `results/interface/`: Sensitivity tables
- `results/mechanism/`: Mechanistic and allostery reports
- `results/ligand/`: TELIN-1 audit reports
- `results/translation/`: Cancer-aging matrix and therapeutic window reports
- `results/validation_experiments/`: Wet-lab prediction protocols
- `results/tables/`: Master numerical tables
- `results/figures/`: Publication figures
- `tests/`: Automated unit, adversarial, and E2E test suites
