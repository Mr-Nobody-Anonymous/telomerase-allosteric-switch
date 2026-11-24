<!-- generated, do not edit by hand: python scripts/solve/wp49_master_record.py -->

# Artefact index — every file, what it is, and whether it is generated

    python scripts/solve/wp49_master_record.py

`generated` means the file is rewritten by a script and must not be hand-edited.

## Scripts

| script | produces | self-test |
|---|---|---|
| `wp10_model_discriminator` | — | no |
| `wp11_audit_identity` | — | no |
| `wp12_recovery` | — | no |
| `wp13_identifiability` | — | no |
| `wp14_robustness` | — | no |
| `wp15_design_validation` | — | no |
| `wp18_candidate_model` | — | no |
| `wp19_literature_closure` | — | no |
| `wp1_anchor_geometry` | — | no |
| `wp20_exhaustive_feasibility` | — | no |
| `wp22_rap_model_v2` | — | no |
| `wp23_tpp1_free_series` | — | no |
| `wp24_gate3_harness` | — | no |
| `wp25_resolving_tests` | — | no |
| `wp26_simulation_readiness` | — | no |
| `wp27_transcription_layer` | — | no |
| `wp28_operating_point` | — | yes |
| `wp29_maintenance` | — | yes |
| `wp2_model_audit` | — | no |
| `wp30_ligand_completion` | — | yes |
| `wp31_epitope` | — | yes |
| `wp32_physics_scoring` | — | yes |
| `wp33_valid_scoring` | — | yes |
| `wp34_status` | — | yes |
| `wp35_p2_identifiability` | — | yes |
| `wp36_circularity` | — | yes |
| `wp37_missing` | — | yes |
| `wp38_archive` | RETRACTIONS.md, retractions.json/.csv | yes |
| `wp39_executable` | — | yes |
| `wp3_step_resolved_model` | — | no |
| `wp40_design_to_fifty` | — | yes |
| `wp41_status_ledger` | STATUS_LEDGER.md, status_ledger.json/.csv | yes |
| `wp42_designs` | wp42_designs.json/.csv, WP42_DESIGNS.md | yes |
| `wp43_test_battery` | wp43_test_battery.json/.csv, WP43_BATTERY.md | yes |
| `wp44_modular_designs` | wp44_modular_designs.json/.csv, WP44_MODULAR.md | yes |
| `wp45_graft` | wp45_graft.json/.csv, WP45_GRAFT.md | yes |
| `wp46_closure` | wp46_closure.json/.csv, CLOSURE.md | yes |
| `wp47_close_open` | wp47_close_open.json/.csv, WP47_CLOSE_OPEN.md | yes |
| `wp48_new_designs` | wp48_new_designs.json/.csv, WP48_NEW_DESIGNS.md | yes |
| `wp49_master_record` | MASTER_RECORD.md, METHODS.md, ARTEFACT_INDEX.md, NEXT_STEPS.md | yes |
| `wp4_coupling_paths` | — | no |
| `wp50_partials` | wp50_partials.json/.csv, WP50_PARTIALS.md | yes |
| `wp51_ordering_md` | wp51_ordering_md.json/.csv, WP51_ORDERING_MD.md | yes |
| `wp52_interface` | wp52_interface.json/.csv, WP52_INTERFACE.md | yes |
| `wp53_conformational_change` | wp53_conformational_change.json/.csv, WP53_CONFORMATIONAL_CHANGE.md | yes |
| `wp54_remaining` | wp54_remaining.json/.csv, WP54_REMAINING.md | yes |
| `wp5_decisive_experiment` | — | no |
| `wp6_mutant_panel` | — | no |
| `wp7_entry_independence` | — | no |
| `wp8_dna_polarity` | — | no |
| `wp9_ten_face` | — | no |

## Result documents

| file | contents | generated |
|---|---|---|
| `MASTER_RECORD.md` | the whole programme, all numbers | yes |
| `METHODS.md` | how every number was produced | yes |
| `ARTEFACT_INDEX.md` | this file | yes |
| `NEXT_STEPS.md` | materials, protocols, decisions | yes |
| `CLOSURE.md` | status of every sub-problem | yes |
| `RETRACTIONS.md` | every retracted/corrected claim | yes |
| `STATUS_LEDGER.md` | mechanism-probability ledger | yes |
| `WP42_DESIGNS.md` | interface geometry and the 5 ligand + 5 protein designs | by its script |
| `WP43_BATTERY.md` | the filter battery and its failure | by its script |
| `WP44_MODULAR.md` | per-site modules and merge stability (carries CORRECTED BY WP45) | by its script |
| `WP45_GRAFT.md` | the S3 graft (carries CORRECTED BY WP47) | by its script |
| `WP47_CLOSE_OPEN.md` | crosslink screen, filter ceiling, assay protocols | by its script |
| `WP48_NEW_DESIGNS.md` | four agents, four new designs, the critic | by its script |
| `WP50_PARTIALS.md` | the three PARTIAL items attacked by other means | by its script |

## Machine-readable

| file | contents |
|---|---|
| `retractions.json` | machine-readable output of its generator |
| `status_ledger.json` | machine-readable output of its generator |
| `wp1_chains.json` | machine-readable output of its generator |
| `wp36_circularity.json` | machine-readable output of its generator |
| `wp39_executable.json` | machine-readable output of its generator |
| `wp3_predictions.json` | machine-readable output of its generator |
| `wp40_design_to_fifty.json` | machine-readable output of its generator |
| `wp42_designs.json` | machine-readable output of its generator |
| `wp43_test_battery.json` | machine-readable output of its generator |
| `wp44_modular_designs.json` | machine-readable output of its generator |
| `wp45_graft.json` | machine-readable output of its generator |
| `wp46_closure.json` | machine-readable output of its generator |
| `wp47_close_open.json` | machine-readable output of its generator |
| `wp48_new_designs.json` | machine-readable output of its generator |
| `wp49_master_record.json` | machine-readable output of its generator |
| `wp50_partials.json` | machine-readable output of its generator |
| `wp51_ordering_md.json` | machine-readable output of its generator |
| `wp52_interface.json` | machine-readable output of its generator |
| `wp53_conformational_change.json` | machine-readable output of its generator |
| `wp54_remaining.json` | machine-readable output of its generator |

## Infrastructure

| file | role |
|---|---|
| `scripts/setup_environment.sh` | rebuild the venv and restore the branch tip |
| `tests/` | the suite; 727 tests |
| `data/raw/*.cif` | the five deposits, immutable |
