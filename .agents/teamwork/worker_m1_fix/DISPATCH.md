# Dispatch for Milestone 1 Remediation Worker (`worker_m1_fix`)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Read PROJECT.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1_fix`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Remediation Scope for Milestone 1 (Iteration 2)
The Milestone 1 Gate reported 3 unanimous, concrete findings that must be remediated:

1. **Remove Static Q-score Fallback Dictionary in `src/telotrade/validate.py`**:
   - In `parse_validation_metrics()` (around lines 318–325), remove the static dictionary `q_scores = {"7TRE": 0.392, "7TRD": 0.447, ...}`.
   -wwPDB validation XML files (`data/raw/*_validation.xml.gz`) natively supply `Q_score` in the root `<Entry>` element (e.g. `entry_elem.attrib.get("Q_score")`). Parse this value dynamically and cast to float if present.
   - Do NOT use hardcoded dictionary lookups for validation metrics.

2. **Widen Exception Handling in `src/telotrade/validate.py` for Zero-Crash Guarantee**:
   - In `parse_structure_summary()`:
     - Add an explicit 0-byte check: `if path.stat().st_size == 0: ...` setting `gemmi_status = "EMPTY_FILE (0 bytes)"` and `biopython_status = "EMPTY_FILE (0 bytes)"`.
     - In Gemmi block (around line 424): catch `(RuntimeError, ValueError, OSError, AttributeError, IndexError, KeyError, Exception)`.
     - In Biopython block (around line 442): catch `(Bio.PDB.PDBExceptions.PDBConstructionException, ValueError, OSError, AttributeError, KeyError, IndexError, Exception)`.
   - In `parse_validation_metrics()` (around line 326): catch `(ValueError, KeyError, RuntimeError, OSError, IndexError, Exception)`.
   - In `main()` (around line 553): catch `(FileNotFoundError, ValueError, RuntimeError, OSError, IndexError, KeyError, Exception)`.

3. **Correct `data_manifest.csv` and `data/manifest/data_manifest.csv`**:
   - For 2I46 (line 6), change `Space Group: P 21 21 21` to `Space Group: I 41 2 2` (Space Group #98, Tetragonal).
   - In 2I46 Chain B unmodeled segments, ensure `126-139` and `242-250` are documented: `"Chain_A:aa 242-250;Chain_B:aa 126-139,168-171(DE-loop),242-250"`.

4. **Verify All Tests Pass 100%**:
   - Run `python -m pytest tests/test_validate.py tests/test_adversarial_m1.py -v`.
   - Run `python -m pytest tests/e2e/ -v`.
   - Run `python -m telotrade.validate`.
   - Run `ruff check src tests/`.
   - Ensure all 90 tests pass (27 validate + 23 adversarial + 40 e2e) with 0 failures.

Write your 5-component `handoff.md` and report to parent.


## 2026-10-16T15:13:20Z
Received from eb46fecb-7562-4ce7-aa98-998216cee2c4:
You are Milestone 1 Remediation Worker (Iteration 2).
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1_fix
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1_fix\DISPATCH.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope:
1. In src/telotrade/validate.py: Remove static q_scores dictionary in parse_validation_metrics(); dynamically parse attrib.get('Q_score') from XML.
2. In src/telotrade/validate.py: Widen exception handling in parse_structure_summary, parse_validation_metrics, and main() to catch IndexError, KeyError, and general Exception; add explicit 0-byte file check.
3. In data_manifest.csv and data/manifest/data_manifest.csv: Fix 2I46 space group to 'I 41 2 2' and complete unmodeled boundaries for Chain B ('126-139,168-171(DE-loop),242-250').
4. Run all tests: pytest tests/test_validate.py tests/test_adversarial_m1.py tests/e2e/ -v. Verify 90 tests pass with 0 failures.
5. Write your 5-component handoff.md and send completion message to parent.
