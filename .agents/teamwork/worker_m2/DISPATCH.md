# Dispatch for Milestone 2 Worker (`worker_m2`)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Read PROJECT.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md`

Read survey reports:
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md`
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Milestone 2 Scope & Deliverables
Implement `src/telotrade/sequence_mapping.py` and comprehensive unit tests:

1. **Automated mmCIF Sequence Extraction**:
   - Extract polymer sequences and modeled coordinate sequences using Gemmi and Biopython.
   - Disambiguate `label_asym_id` vs `auth_asym_id` according to `data_manifest.csv`.

2. **Construct Tag Offset Reconcilers**:
   - TERT in 7TRE, 7TRD, 7TRF: Construct contains N-terminal 35-aa Twin-Strep tag. `auth_seq_num = label_seq_id - 35`, exactly matching canonical UniProt O14746. Modeled residues span 7–178 and 322–1132 (residues 179–321 unmodeled).
   - TPP1 in 7TRE: Construct contains 3-aa GAG tag. `auth_seq_num = label_seq_id + 86`, matching canonical UniProt Q96AP0 (residues 90–241 modeled).
   - TPP1 in 2I46: `auth_seq_num = label_seq_id + 89`.

3. **2I46 Dimer Asymmetry Handling**:
   - In 2I46, Chain A has TEL-patch DE-loop (168–171) modeled. In Chain B, the DE-loop is disordered (`auth_seq_num = '?'`).
   - Implement explicit selection of 2I46 Chain A and provide functions to detect/warn when Chain B is queried for TEL-patch residues.

4. **Sequence-Aware Homology Alignment**:
   - Pairwise sequence alignment (using Biopython `Bio.Align.PairwiseAligner` or dynamic programming) to establish residue-level correspondence between homologous chains across structures before any geometric superposition.
   - Map equivalent residues between:
     - 7TRE TERT (chain B/A) ↔ 7TRD TERT (chain B/A)
     - 7TRE TERT (chain B/A) ↔ 7TRF TERT (chain D/A)
     - 7TRE TPP1 (chain D/C) ↔ 2I46 TPP1 (chain A)

5. **Expose in `src/telotrade/__init__.py`**:
   - Expose core mapping functions in `__init__.py`.

6. **Unit Tests & Quality Verification**:
   - Implement `tests/test_sequence_mapping.py` verifying all mappings, offsets, alignments, and edge cases.
   - Run `python -m pytest tests/test_sequence_mapping.py -v`.
   - Run `python -m pytest tests/ -v`.
   - Run `ruff check src tests/`.

Write your 5-component `handoff.md` and report to parent.

## 2026-10-16T15:29:39Z
You are Milestone 2 Implementation Worker.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m2\DISPATCH.md
Read survey reports:
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope:
1. Implement src/telotrade/sequence_mapping.py:
   - Extract polymer sequences from mmCIF via Gemmi and Biopython.
   - Handle construct tag offsets (TERT: auth = label - 35, matching UniProt O14746; TPP1 7TRE: auth = label + 86; TPP1 2I46: auth = label + 89).
   - Implement explicit selection of 2I46 Chain A and handle Chain B DE-loop disorder.
   - Sequence-aware pairwise alignment to map equivalent residue correspondences across 7TRE, 7TRD, 7TRF, and 2I46 before geometric alignment.
2. Expose core mapping API in src/telotrade/__init__.py.
3. Write unit tests in tests/test_sequence_mapping.py, run pytest and ruff, verify 100% pass.
4. Write 5-component handoff.md and send completion message to parent.
