# Dispatch for Milestone 2 Remediation Worker (`worker_m2_fix`)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Read PROJECT.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2_fix`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Remediation Scope for Milestone 2 (Iteration 2)
The Milestone 2 Gate reported 4 concrete findings across Reviewer 2 and Challenger 1:

1. **Fix `_get_chain_sequence_biopython()` to Honor `use_auth=False`**:
   - In `src/telotrade/sequence_mapping.py`: When `use_auth=False`, match chains and extract residues using CIF `label_asym_id` and `label_seq_id` (from mmCIF dict `_atom_site.label_asym_id` and `_atom_site.label_seq_id`) so that CIF label assignments strictly match the Gemmi parser branch.

2. **Entity Classification Guard in `get_uniprot_mapping()`**:
   - Only map protein chains (TERT or TPP1) to UniProt.
   - If called on nucleic acid chains (e.g. 7TRE hTR RNA `label=A` / `auth=B` or DNA `label=C` / `auth=D`), do NOT apply protein offsets (which produce negative numbers). Return an empty dictionary or raise `ValueError` for non-protein entities.

3. **Empty Sequence Guard in `align_sequences()`**:
   - In `align_sequences(seq1, seq2, ...)`:
     ```python
     if not seq1 or not seq2:
         return AlignmentResult(
             query_seq=seq1,
             target_seq=seq2,
             score=0.0,
             identity=0.0,
             query_to_target={},
             target_to_query={},
         )
     ```
   - This prevents Biopython's PairwiseAligner from raising `ValueError: sequence has zero length`.

4. **Dynamic Entity Mapping in `generate_chain_entity_mapping()`**:
   - Ensure `generate_chain_entity_mapping` dynamically parses entity records from `data/raw/*.cif` or `data_manifest.csv` rather than relying on any static dictionaries for 7TRC or 7TRF.

5. **Test Verification**:
   - Run `python -m pytest tests/test_sequence_mapping.py tests/test_adversarial_m2.py -v`.
   - Run full suite: `python -m pytest tests/ -v`.
   - Run `ruff check src tests/`.
   - Ensure all unit, adversarial, and E2E tests pass with 0 failures and clean linting.

Write your 5-component `handoff.md` and report to parent.


## 2026-10-16T16:21:28Z
You are Milestone 2 Remediation Worker (Iteration 2).
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2_fix
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2_fix\DISPATCH.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope:
1. In src/telotrade/sequence_mapping.py: Make _get_chain_sequence_biopython honor use_auth=False by extracting CIF label_asym_id and label_seq_id.
2. In get_uniprot_mapping: Guard against non-protein chains (return empty dict for RNA/DNA).
3. In align_sequences: Add early check for empty string inputs ('').
4. In generate_chain_entity_mapping: Dynamically extract entity counts from data_manifest.csv or mmCIF files.
5. Run pytest tests/test_sequence_mapping.py tests/test_adversarial_m2.py tests/ -v and ruff check. Verify 100% pass.
6. Write 5-component handoff.md and send message to parent.
