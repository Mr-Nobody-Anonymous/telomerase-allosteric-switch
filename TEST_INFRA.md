# Test Infrastructure & Methodology (`TEST_INFRA.md`)

## 1. Testing Philosophy & Guiding Principles

This test infrastructure implements a rigorous, requirement-driven, opaque-box testing framework for the **Human Telomerase TERT–TPP1 Structural Bioinformatics Study**. The testing regime is derived strictly from `ORIGINAL_REQUEST.md` (Requirements R1 through R5) and `PROJECT.md` (Features 1 through 18 and Interface Contracts).

### Core Principles
1. **Opaque-Box Contract Verification**: Tests validate behavioral specifications and interface contracts rather than private implementation details. Every module is verified against its declared function signatures, input schemas, and expected outputs.
2. **Empirical Ground-Truth Derivation**: Expected outputs are anchored to authoritative references:
   - Structural ground truths from official wwPDB / RCSB mmCIF coordinate files (`7TRE`, `7TRD`, `7TRF`, `7TRC`, `2I46`).
   - Canonical sequence databases (UniProt `O14746` for human TERT, `Q96AP0` for human TPP1).
   - Validated biophysical and mathematical formulas (Shrake-Rupley sphere discretization for SASA/BSA, SVD Kabsch algorithm with proper rotation determinant check).
3. **No Facades / Anti-Cheat Guarantee**: Every test executes genuine computational logic (coordinate arithmetic, KD-tree/cell-hashing spatial search, Kabsch SVD decomposition, SASA numerical integration, or strict CSV schema and statistical distribution checks). Facades, dummy mocks of core algorithms, or trivial assertions are strictly prohibited.
4. **Progressive Testability & Isolation**: Tests are self-contained and execute cleanly. Where implementation modules are loaded dynamically, clear diagnostics guide the pipeline engineer toward compliance without obscuring functional regressions.

---

## 2. 4-Tier Test Architecture

The E2E test suite is organized under `tests/e2e/` across four distinct tiers:

```
tests/e2e/
├── conftest.py                   # Shared fixtures, structure loaders, reference data
├── test_tier1_features.py        # Tier 1: Core Feature Coverage (>=5 tests per feature group)
├── test_tier2_boundaries.py      # Tier 2: Boundary, Corner Cases & Adversarial Inputs
├── test_tier3_interactions.py    # Tier 3: Cross-Feature Interactions & Pairwise Workflows
└── test_tier4_scenarios.py       # Tier 4: Real-World Scenarios, Pipeline Execution & Artifacts
```

### Tier Breakdown

| Tier | Focus | Description | Marker |
| :--- | :--- | :--- | :--- |
| **Tier 1** | **Feature Coverage** | Systematic verification of all 18 features from `PROJECT.md`. Validates mmCIF acquisition, parsing, sequence mapping, contact networks, BSA calculation, RMSD superpositions, and table outputs. | `@pytest.mark.tier1` |
| **Tier 2** | **Boundary & Corner Cases** | Edge condition stress tests: unmodeled loops (TERT 197–325), 2I46 Chain B TEL-patch disorder, alternate locations (`altloc`), distance threshold boundaries (3.500 Å, 4.000 Å, 4.500 Å), zero-distance handling, and malformed inputs. | `@pytest.mark.tier2` |
| **Tier 3** | **Cross-Feature Interactions** | Multi-module integration: sequence-aware alignment output feeding into interface distance calculations; UniProt index harmonizer feeding into Kabsch SVD superpositions; contact networks feeding into robustness classification and polar interaction mapping. | `@pytest.mark.tier3` |
| **Tier 4** | **Real-World Scenarios** | Complete pipeline reproduction: end-to-end execution of `run_pipeline.py`, full validation of all 11 standardized CSV tables in `results/tables/`, figure generation in `results/figures/`, and synthesis report validation in `report.md`. | `@pytest.mark.tier4` |

---

## 3. Feature Inventory Coverage Matrix

| Feature ID | Feature Name | Tier 1 Tests | Tier 2 Tests | Tier 3 Tests | Tier 4 Tests |
| :---: | :--- | :--- | :--- | :--- | :--- |
| **F1** | Raw Data Acquisition | `test_f1_cif_files_exist`, `test_f1_sha256_checksums` | `test_t2_missing_cif_handling` | `test_t3_raw_data_to_parser_flow` | `test_t4_full_pipeline_data_audit` |
| **F2** | Immutable Data Manifest | `test_f2_manifest_exists`, `test_f2_manifest_columns` | `test_t2_manifest_corrupted_hash` | `test_t3_manifest_to_validation_sync` | `test_t4_manifest_artifact_integrity` |
| **F3** | Baseline Validation Script | `test_f3_validate_module_api`, `test_f3_validate_cli_execution` | `test_t2_validate_altloc_and_disorder` | `test_t3_validate_feeds_downstream` | `test_t4_pipeline_validation_step` |
| **F4** | Sequence-Aware Chain Mapping | `test_f4_chain_extraction`, `test_f4_auth_vs_label_chain` | `test_t2_unmodeled_terminal_residues` | `test_t3_seq_mapping_to_contacts` | `test_t4_all_structures_mapped` |
| **F5** | Construct Tag Offset Handling | `test_f5_tert_tag_offset_35`, `test_f5_tpp1_tag_offset_86` | `test_t2_negative_or_oob_residues` | `test_t3_offset_to_uniprot_reconciliation`| `test_t4_uniprot_tables_populated` |
| **F6** | 2I46 Dimer Asymmetry Handling | `test_f6_2i46_chain_a_tel_patch`, `test_f6_2i46_chain_b_disorder` | `test_t2_chain_b_tel_patch_nan_safety`| `test_t3_2i46_chain_a_to_cross_state` | `test_t4_tpp1_bound_unbound_table` |
| **F7** | Multi-Cutoff Contact Networks | `test_f7_contacts_3_5A`, `test_f7_contacts_4_0A`, `test_f7_contacts_4_5A` | `test_t2_cutoff_exact_boundary_floating_point` | `test_t3_contacts_monotonic_progression` | `test_t4_contact_tables_generated` |
| **F8** | Buried Surface Area (BSA) | `test_f8_bsa_positive_value`, `test_f8_bsa_formula_convention` | `test_t2_bsa_probe_radius_sensitivity` | `test_t3_bsa_isolated_vs_complex_sum` | `test_t4_bsa_table_verification` |
| **F9** | Polar Networks & Salt Bridges | `test_f9_arg774_glu168_salt_bridge`, `test_f9_lys78_glu215` | `test_t2_donor_acceptor_angle_cutoff` | `test_t3_polar_subset_of_contacts` | `test_t4_polar_table_schema` |
| **F10** | Hydrophobic Packing & TEL-Patch | `test_f10_leu769_leu183_contact`, `test_f10_pro771_leu212` | `test_t2_tel_patch_mutation_simulation`| `test_t3_tel_patch_to_b_factor_profile` | `test_t4_tel_patch_representation` |
| **F11** | Cross-State SVD Kabsch RMSD | `test_f11_kabsch_svd_exact_rotation`, `test_f11_core_rmsd_7tre_7trd` | `test_t2_kabsch_reflection_determinant` | `test_t3_kabsch_with_filtered_subsets` | `test_t4_cross_state_rmsd_table` |
| **F12** | TEN-TRAP Rearrangements | `test_f12_ten_domain_displacement`, `test_f12_trap_stability` | `test_t2_unmodeled_linker_197_325_gap` | `test_t3_domain_decomposition_rmsd` | `test_t4_per_residue_disp_table` |
| **F13** | Contact Robustness Table | `test_f13_robustness_classification`, `test_f13_persistence_score` | `test_t2_robustness_fine_grained_sweep`| `test_t3_robustness_to_polar_crossref` | `test_t4_robustness_table_schema` |
| **F14** | TPP1 Bound vs Unbound Shift | `test_f14_tpp1_tel_patch_displacement`, `test_f14_leu183_rigidity` | `test_t2_loop_gap_handling_2i46` | `test_t3_tpp1_rmsd_with_local_bfactor` | `test_t4_tpp1_shift_report_metrics` |
| **F15** | Standardized CSV Tables | `test_f15_tables_exist`, `test_f15_tables_schema_validation` | `test_t2_empty_or_nan_in_tables` | `test_t3_relational_foreign_keys` | `test_t4_all_11_tables_present` |
| **F16** | Synthesis Report Authoring | `test_f16_report_exists`, `test_f16_report_axis_a_and_b` | `test_t2_report_missing_sections` | `test_t3_report_numbers_match_csvs` | `test_t4_report_markdown_linting` |
| **F17** | Publication Figures Generation | `test_f17_figures_directory_populated`, `test_f17_figure_resolutions` | `test_t2_corrupted_figure_images` | `test_t3_figures_referenced_in_report`| `test_t4_visual_artifacts_present` |
| **F18** | Single-Run Pipeline Runner | `test_f18_run_pipeline_cli_syntax`, `test_f18_pipeline_exit_code` | `test_t2_pipeline_partial_failure_recovery` | `test_t3_pipeline_idempotency` | `test_t4_pipeline_reproducibility` |

---

## 4. Test Execution & Runner Invocation

### Full Test Suite Run
```powershell
python -m pytest tests/e2e/ -v
```

### Running by Tier
```powershell
# Tier 1: Core Feature Coverage
python -m pytest tests/e2e/test_tier1_features.py -v -m tier1

# Tier 2: Boundary & Corner Cases
python -m pytest tests/e2e/test_tier2_boundaries.py -v -m tier2

# Tier 3: Cross-Feature Interactions
python -m pytest tests/e2e/test_tier3_interactions.py -v -m tier3

# Tier 4: Real-World Scenarios & Full Pipeline
python -m pytest tests/e2e/test_tier4_scenarios.py -v -m tier4
```

### Pytest Configuration (`pyproject.toml` / `pytest.ini`)
The project utilizes strict markers:
- `tier1`: Feature Coverage tests
- `tier2`: Boundary & Corner Case tests
- `tier3`: Cross-Feature Interaction tests
- `tier4`: Application Scenario & Pipeline tests
- `slow`: Long-running numerical integrations (e.g., FreeSASA, full SVD on 10k atoms)

---

## 5. Authoritative Derivation & Ground Truth Benchmarks

All test assertions compare against analytically derived and empirically verified benchmarks:

1. **Interface Distance Contacts (7TRE TERT Chain A vs TPP1 Chain C)**:
   - At 3.5 Å: Exactly 9 atomic contacts across 5 residue pairs (`TERT Arg774 -- TPP1 Glu168` at 2.32 Å; `TERT Ala46 -- TPP1 Glu169` at 3.35 Å; `TERT Ala45 -- TPP1 Glu171` at 3.45 Å).
   - At 4.0 Å: Exactly 46 atomic contacts across 21 residue pairs (`TERT Leu769 -- TPP1 Leu183` at 3.67 Å, `TERT Pro771 -- TPP1 Leu212` at 3.52 Å; `TERT Lys78 -- TPP1 Glu215` at 3.58 Å).
   - At 4.5 Å: Exactly 123 atomic contacts across 32 residue pairs.
2. **Buried Surface Area (BSA)**:
   - $\mathrm{SASA}_{\text{TERT}} = 49,776.00\ \text{Å}^2$
   - $\mathrm{SASA}_{\text{TPP1}} = 10,651.36\ \text{Å}^2$
   - $\mathrm{SASA}_{\text{Complex}} = 58,713.70\ \text{Å}^2$
   - $\Delta\mathrm{SASA}_{\text{total}} = 1,713.65\ \text{Å}^2$
   - $\mathrm{BSA}_{\text{interface}} = \frac{1}{2}\Delta\mathrm{SASA} = 856.82\ \text{Å}^2$ ($\pm 5.0\ \text{Å}^2$ numerical probe discretization allowance).
3. **Cross-State SVD Kabsch Cα RMSD**:
   - 7TRE vs 7TRD Global Cα RMSD: $0.356\ \text{Å}$ (934 Cα pairs).
   - 7TRE vs 7TRD Catalytic Core (residues $\ge 326$): $0.338\ \text{Å}$ (770 Cα pairs).
   - 7TRE vs 2I46 Chain A TPP1 Global Cα RMSD: $1.966\ \text{Å}$ (152 Cα pairs).
   - TEL-patch Cα Displacements:
     - Glu169: $3.809\ \text{Å}$
     - Lys170: $4.293\ \text{Å}$
     - Glu171: $3.320\ \text{Å}$
     - Leu183: $0.512\ \text{Å}$ (strictly invariant anchor).


---

## Manuscript and analysis guards (added 2026-11-01)

Two modules exist for a reason the rest of this suite does not cover: they guard **prose**. A
document can drift from its own computation without any code failing, so the numbers are asserted
rather than trusted.

| module | guards | method |
|---|---|---|
| `tests/test_mechanism_classification.py` | 25 | re-runs `scripts/mechanism_classification_checks.py`, then asserts each quoted value equals the recomputed one; asserts the SHA-256 freeze is current; asserts no retracted claim has reappeared in the prose |
| `tests/test_paper_manuscript.py` | 28 | the same discipline applied to `paper/MANUSCRIPT.md`, plus: every reference carries a DOI, every figure exists in PNG and PDF, every figure has a legend, and the abstract still says no mechanism is demonstrated |

The design principle: **a negative assertion is worth as much as a positive one.** Both modules
contain assertions of the form "this claim must *not* appear", because the failures in this
repository's history have been overstatements that survived review, not errors in arithmetic.

They run in about 3 seconds and are part of the default suite.
