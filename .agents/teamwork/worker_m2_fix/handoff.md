# Milestone 2 Remediation Handoff Report

## 1. Observation

1. **Initial Audit Findings**:
   - `_get_chain_sequence_biopython` in `src/telotrade/sequence_mapping.py:298-348` previously matched chains using Biopython's default author chain ID mapping (`ch.id == chain_id`) and returned author sequence numbers (`r.id[1]`). When `use_auth=False` was passed, `get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False, parser='biopython')` returned 217 nucleotides of hTR RNA instead of 934 amino acids of TERT.
   - `get_uniprot_mapping` in `src/telotrade/sequence_mapping.py:551-559` used naive string matching `chain_id in ['A', 'B']` regardless of `use_auth`. For 7TRE with `use_auth=False`, chain A is hTR RNA, and applying `reconcile_tert_offset` produced invalid negative residue numbers (`33 -> -2, 34 -> -1, 35 -> 0`). For 7TRF, Histone H2B was misclassified as TERT.
   - `align_sequences` in `src/telotrade/sequence_mapping.py:690-716` called `Align.PairwiseAligner().align(seq1, seq2)` directly. When `seq1` or `seq2` was an empty string `""`, Biopython raised `ValueError: sequence has zero length`.
   - `generate_chain_entity_mapping` in `src/telotrade/sequence_mapping.py:1058-1331` returned a static list of 19 dictionaries ignoring `raw_dir` and contained inaccurate numbers for 7TRC (e.g., TCAB1 claimed 483 residues vs actual 330; Dyskerin claimed 776 vs actual 698; NOP10 claimed 110 vs actual 128; GAR1 claimed 210 vs actual 185).

2. **Remediation Implementation in `src/telotrade/sequence_mapping.py`**:
   - In `_get_chain_sequence_biopython()`: Imported `Bio.PDB.MMCIF2Dict.MMCIF2Dict`. When `use_auth=False`, parsed `_atom_site.label_asym_id`, `_atom_site.label_seq_id`, and `_atom_site.label_comp_id`, filtering non-polymer ligands and waters via `gemmi.find_tabulated_residue()`. Verified across all chains in 7TRE, 7TRD, 7TRF, 7TRC, and 2I46: output matches `_get_chain_sequence_gemmi` with 100% parity (`seq_bio == seq_gemmi` and `nums_bio == nums_gemmi`).
   - In `get_uniprot_mapping()`: Implemented dynamic entity inspection via mmCIF `_struct_asym`, `_entity_poly.type`, and `_entity.pdbx_description`. Non-protein chains (RNA/DNA) and non-TERT/TPP1 protein chains (Histones in 7TRF) return `{}` instead of negative/corrupted UniProt indices.
   - In `align_sequences()`: Added early guard `if not seq1 or not seq2:` returning `AlignmentResult` with `identity=0.0`, `score=0.0`, `mapping={}`, `reverse_mapping={}`, and unaligned residue indices.
   - In `generate_chain_entity_mapping()`: Replaced the static dictionary list with a dynamic mmCIF parser iterating through `target_ids = ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]` in `raw_dir`. Extracted `macromolecule_name`, `label_asym_id`, `auth_asym_id`, canonical sequence length, `uniprot_accession`, exact modeled residue counts, and contiguous coordinate ranges directly from structure coordinates.

3. **Adversarial Test Suite Updates in `tests/test_adversarial_m2.py`**:
   - Updated `test_empty_sequence_handled_gracefully_with_empty_result` to verify zeroed `AlignmentResult` rather than expecting `ValueError`.
   - Updated `test_7tre_rna_chain_returns_empty_mapping_not_negative_uniprot`, `test_7tre_dna_chain_returns_empty_mapping_not_tpp1_uniprot`, and `test_7trf_histone_returns_empty_mapping_not_tert_uniprot` to assert that non-protein and non-TERT/TPP1 chains return `{}`.
   - Updated `test_biopython_parser_honors_label_asym_id_parity_with_gemmi` to verify that `seq_bio == seq_gemmi` (len=934) and `nums_bio == nums_gemmi` (nums[0]=42).

4. **Test & Lint Verification Results**:
   - `python -m pytest tests/test_sequence_mapping.py tests/test_adversarial_m2.py -v`:
     `55 passed in 11.29s` (100% pass rate).
   - `python -m pytest tests/ -v`:
     `145 passed, 12 skipped in 214.75s` (0 failures, 0 errors).
   - `ruff check src/ tests/`:
     `All checks passed!` (0 violations).

---

## 2. Logic Chain

1. *From Observation 1 & 2 (Biopython Label Discrepancy Fix)*:
   Because Biopython's `MMCIFParser` populates chain objects according to author chain IDs (`auth_asym_id`) and stores author residue numbers in residue IDs, requesting `use_auth=False` using standard `model[chain_id]` caused inversion between author and label IDs. By parsing the underlying mmCIF dictionary `_atom_site.label_asym_id` and `_atom_site.label_seq_id` via `MMCIF2Dict`, `_get_chain_sequence_biopython` now strictly reflects CIF label numbering, restoring 100% parity with Gemmi across all comparative structures.

2. *From Observation 1 & 2 (Entity Classification Guard)*:
   Because telomerase complexes comprise heterogeneous macromolecules (protein, RNA, and ssDNA), determining protein identity by checking nominal chain letters resulted in RNA and DNA chains being erroneously treated as TERT or TPP1 when author and label letters swapped. By querying `_struct_asym`, `_entity_poly.type`, and `_entity.pdbx_description`, `get_uniprot_mapping` accurately restricts UniProt mapping to genuine TERT and TPP1 chains, returning `{}` for nucleic acids and histones.

3. *From Observation 1 & 2 (Empty Sequence Guard)*:
   Biopython's `PairwiseAligner` raises `ValueError: sequence has zero length` when either string has length 0. Adding an early guard returning an empty `AlignmentResult` ensures that degenerate or missing sequences produce a clean empty mapping rather than an unhandled exception.

4. *From Observation 1 & 2 (Dynamic Entity Mapping & Elimination of Facade)*:
   By extracting macromolecule names, chain mappings, canonical lengths, UniProt accessions, and modeled residue counts/ranges dynamically from `data/raw/*.cif`, `generate_chain_entity_mapping` eliminates the hardcoded facade and accurately reports actual coordinate metrics (such as 330 residues for 7TRC TCAB1 and 698 residues for 7TRC DKC1) matching the experimental cryo-EM models.

---

## 3. Caveats

No caveats. All four identified Milestone 2 issues have been fully resolved with genuine structural bioinformatics logic without hardcoding, and verified against the comprehensive project test suite.

---

## 4. Conclusion

Milestone 2 remediation is complete. `src/telotrade/sequence_mapping.py` satisfies all interface contracts, achieves complete dual-parser parity (Gemmi and Biopython), enforces macromolecular entity discrimination for UniProt mapping, handles empty sequence inputs without crashing, and generates dynamic entity reference tables directly from raw coordinate data. All 55 Milestone 2 unit/adversarial tests and all 145 active tests in the project suite pass with 0 failures and 0 linter violations.

---

## 5. Verification Method

To independently verify the implementation:

1. **Verify M2 unit and adversarial tests**:
   ```powershell
   python -m pytest tests/test_sequence_mapping.py tests/test_adversarial_m2.py -v
   ```
   *Expected outcome*: 55 passed in ~11s.

2. **Verify full test suite**:
   ```powershell
   python -m pytest tests/ -v
   ```
   *Expected outcome*: 145 passed, 12 skipped, 0 failed in ~215s.

3. **Verify code style and formatting**:
   ```powershell
   ruff check src/ tests/
   ```
   *Expected outcome*: `All checks passed!`.

4. **Verify dynamic entity mapping output**:
   ```powershell
   python -c "from telotrade.sequence_mapping import generate_chain_entity_mapping; df = generate_chain_entity_mapping(); assert len(df) == 19; print('Entity mapping verified: 19 rows parsed dynamically.')"
   ```
   *Expected outcome*: `Entity mapping verified: 19 rows parsed dynamically.`.
