# Dispatch for E2E Test Writer (`e2e_test_writer`)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Read PROJECT.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\e2e_test_writer`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## E2E Testing Track Scope & Deliverables
Design a comprehensive, requirement-driven, opaque-box E2E test suite derived strictly from `ORIGINAL_REQUEST.md` and `PROJECT.md § Feature Inventory`.

1. **Create `TEST_INFRA.md`** at project root (`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\TEST_INFRA.md`):
   - Document test philosophy, feature inventory coverage, test architecture, and runner invocation (`python -m pytest tests/e2e/ -v`).

2. **Implement 4-Tier Test Suite** under `tests/e2e/`:
   - **Tier 1 - Feature Coverage**: Minimum 5 test cases per feature for basic data existence, parsing, contacts, RMSD, table outputs.
   - **Tier 2 - Boundary & Corner Cases**: Edge cases (missing residues, unmodeled loops, altlocs, distance cutoff boundaries 3.5/4.0/4.5 Å, chain B disorder in 2I46).
   - **Tier 3 - Cross-Feature Combinations**: Pairwise interactions (e.g. sequence mapping feeding into interface contact calculation, SVD Kabsch superposition using mapped sequence indices).
   - **Tier 4 - Real-World Application Scenarios**: End-to-end pipeline execution reproducing all 11 tables in `results/tables/`, figure generation in `results/figures/`, and `report.md` presence.

3. **Publish `TEST_READY.md`** at project root:
   - Provide test runner command, tier breakdown table, feature checklist, and expected exit code 0.

Write your 5-component `handoff.md` in your directory and send a completion message to parent.


## 2026-10-16T14:39:25Z
From: eb46fecb-7562-4ce7-aa98-998216cee2c4
Content:
You are the E2E Test Writer.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\e2e_test_writer
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\e2e_test_writer\DISPATCH.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope:
1. Create TEST_INFRA.md at project root with testing architecture and methodology.
2. Implement 4-tier E2E test suite in tests/e2e/ covering Feature Coverage (Tier 1), Boundary & Corner Cases (Tier 2), Cross-Feature Interactions (Tier 3), and Real-World Application Scenarios (Tier 4).
3. Publish TEST_READY.md at project root with runner instructions and coverage checklist.
4. Write handoff.md and send message to parent when complete.
