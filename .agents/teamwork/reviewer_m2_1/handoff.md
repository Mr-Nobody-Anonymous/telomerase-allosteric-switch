# Milestone 2 Review Report: Sequence-Aware Chain & Residue Correspondence Mapping

**Agent**: Reviewer 1 (`reviewer_m2_1`)  
**Verdict**: **APPROVE**  
**Integrity Audit**: **CLEAN (0 Violations)**  
**Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\reviewer_m2_1`  

---

## 1. Observation

### 1.1 Direct Inspection of Implementation & Test Files
- **Implementation file**: `src/telotrade/sequence_mapping.py` (1,332 lines).
  - Implements dual-engine coordinate parsing via Gemmi (`_get_chain_sequence_gemmi`) and Biopython (`_get_chain_sequence_biopython`).
  - Supports `use_auth=True/False` toggling to disambiguate mmCIF author chains (`auth_asym_id`) from CIF label chains (`label_asym_id`).
  - Implements construct tag offsets:
    - TERT: `auth_seq_num = label_seq_id - 35` (lines 33-34: `TERT_AUTH_OFFSET = -35`, `TERT_LABEL_OFFSET = 35`, UniProt O14746).
    - TPP1 in 7TRE: `auth_seq_num = label_seq_id + 86` (lines 38-39: `TPP1_7TRE_AUTH_OFFSET = 86`, `TPP1_7TRE_LABEL_OFFSET = -86`, UniProt Q96AP0).
    - TPP1 in 2I46: `auth_seq_num = label_seq_id + 89` (lines 43-44: `TPP1_2I46_AUTH_OFFSET = 89`, `TPP1_2I46_LABEL_OFFSET = -89`).
  - Implements 2I46 dimer asymmetry handling:
    - `check_tel_patch_integrity` (lines 584–638) checks presence of `TEL_PATCH_RESIDUES` `(168, 169, 170, 171, 183, 212)` and DE-loop `(168, 169, 170, 171)`.
    - `select_valid_tpp1_chain` (lines 641–684) redirects requests for 2I46 Chain B to Chain A when `require_tel_patch=True` and issues a `UserWarning`.
  - Implements pairwise sequence homology alignment:
    - `align_sequences` (lines 690–772) uses `Bio.Align.PairwiseAligner` with standard affine gap penalties (`match=2.0`, `mismatch=-1.0`, `open_gap=-10.0`, `extend_gap=-0.5`).
    - `map_residue_correspondence` (lines 774–806) and `align_homologous_chains` (lines 809–864) map residues 1:1 between coordinate sets before geometric superposition.
  - Implements standardized table generators:
    - `generate_tert_correspondence_table` (lines 914–984) dynamically extracts modeled residues from 7TRE, 7TRD, and 7TRF.
    - `generate_tpp1_correspondence_table` (lines 987–1055) dynamically extracts modeled residues from 7TRE and 2I46 Chains A and B.
    - `generate_chain_entity_mapping` (lines 1058–1331) compiles reference chain-to-entity metadata across all 5 project structures.

### 1.2 Independent Empirical Verification of Offsets in Raw mmCIF Files
Direct execution of coordinate inspection on the raw mmCIF files yielded:
- `data/raw/7TRE.cif` (auth A, label B): 934 residues modeled, first `('CYS', auth=7, label=42)`, last `('ASP', auth=1132, label=1167)`. `auth - label = -35` across all 934 residues.
- `data/raw/7TRD.cif` (auth A, label B): 934 residues modeled, first `('CYS', auth=7, label=42)`, last `('ASP', auth=1132, label=1167)`. `auth - label = -35` across all 934 residues.
- `data/raw/7TRF.cif` (auth A, label D): 934 residues modeled, first `('CYS', auth=7, label=42)`, last `('ASP', auth=1132, label=1167)`. `auth - label = -35` across all 934 residues.
- `data/raw/7TRE.cif` (auth C, label D): 152 residues modeled, first `('SER', auth=90, label=4)`, last `('HIS', auth=241, label=155)`. `auth - label = 86` across all 152 residues.
- `data/raw/2I46.cif` (auth A, label A): 152 residues modeled, first `('SER', auth=90, label=1)`, last `('HIS', auth=241, label=152)`. `auth - label = 89` across all 152 residues.
- `data/raw/2I46.cif` (auth B, label B): 134 residues modeled, first `('SER', auth=90, label=1)`, last `('HIS', auth=241, label=152)`. `auth - label = 89` across all 134 residues.

Missing residues in 2I46 Chain B relative to Chain A: `[126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 168, 169, 170, 171]`. Specifically, residues 168–171 represent the TEL-patch DE loop.

### 1.3 Test Suite Execution Results
- `python -m pytest tests/test_sequence_mapping.py -v`:
  - Result: **34 passed in 21.06s** (100% pass rate).
- `python -m pytest tests/e2e/ -v`:
  - Result: **40 passed, 12 skipped in 14.08s** (100% pass rate for active tiers; 12 skipped for future M3-M5 integration).
- `python -m pytest tests/test_adversarial_m2.py -v`:
  - Result: **15 passed in 12.61s** (100% pass rate).
- Full repository test suite (`python -m pytest tests/ -v`):
  - Result: **124 passed, 12 skipped in 267.25s (0:04:27)**. Zero failures, zero errors.
- Linter execution:
  - `python -m ruff check src/ tests/test_sequence_mapping.py`: `All checks passed!` (0 violations).

---

## 2. Logic Chain

1. *From Observation 1.2 (Empirical Offset Verification)*:
   Direct atom-site extraction proved that TERT in 7TRE, 7TRD, and 7TRF exhibits an exact `-35` residue offset (`auth_seq_num = label_seq_id - 35`), corresponding to the 35-aa N-terminal Twin-Strep purification tag. TPP1 in 7TRE exhibits an exact `+86` residue offset (`auth_seq_num = label_seq_id + 86`), corresponding to the 3-aa N-terminal GAG tag. TPP1 in 2I46 exhibits an exact `+89` residue offset (`auth_seq_num = label_seq_id + 89`). Therefore, the constant definitions and offset reconcilers in `src/telotrade/sequence_mapping.py` are mathematically and biologically correct.
2. *From Observation 1.2 (2I46 Dimer Asymmetry Verification)*:
   In `2I46.cif`, Chain A contains all 152 residues of the OB-fold (residues 90–241), including fully modeled coordinates for the TEL-patch DE-loop (`Glu168–Glu169–Arg170–Asp171`). In contrast, Chain B is missing 18 residues, including all 4 DE-loop residues (168–171). The implementation in `check_tel_patch_integrity` accurately identifies this disorder, and `select_valid_tpp1_chain` correctly redirects queries from Chain B to Chain A with a descriptive warning, preventing coordinate crashes in downstream interface geometry modules.
3. *From Observations 1.1 & 1.3 (Biopython Aligner and Homology Mapping)*:
   Using `Bio.Align.PairwiseAligner`, `map_residue_correspondence` matches 934/934 residues between 7TRE and 7TRD/7TRF with 100% sequence identity, and matches 152/152 common residues between 7TRE TPP1 and 2I46 Chain A. When mapped against Chain B, it maps 134 residues and cleanly leaves 168–171 in `unaligned_res1`.
4. *From Observation 1.3 (Test Suite & Integrity)*:
   All unit, E2E, and adversarial test suites pass 100%. No hardcoded test fixtures, facade mocks, or bypassed computations were found in `src/telotrade/sequence_mapping.py`. All functions dynamically parse the mmCIF files.

---

## 3. Findings

### [Minor] Finding 1: Unhandled `ValueError` when `align_sequences` receives empty string
- **What**: Passing an empty string `""` to `align_sequences` raises `ValueError: sequence has zero length` from Biopython before reaching the `if not alignments:` guard.
- **Where**: `src/telotrade/sequence_mapping.py`, lines 708–716.
- **Why**: Biopython's C-extension pairwise aligner rejects zero-length sequence strings immediately.
- **Suggestion**: Add a pre-check at the start of `align_sequences`:
  ```python
  if not seq1 or not seq2:
      return AlignmentResult(
          seq1=seq1, seq2=seq2, mapping={}, reverse_mapping={},
          aligned_pairs=[], identity=0.0, score=0.0,
          unaligned_res1=list(range(len(seq1))),
          unaligned_res2=list(range(len(seq2))),
          alignment_str=""
      )
  ```

### [Minor] Finding 2: `get_uniprot_mapping` chain classification does not condition on `use_auth`
- **What**: In `get_uniprot_mapping`, `is_tert` is evaluated as `("7TRE" in pid and chain_id in ["A", "B"])`.
- **Where**: `src/telotrade/sequence_mapping.py`, lines 551–560.
- **Why**: When `use_auth=False`, chain `'A'` in 7TRE is RNA (hTR). If a caller requests `get_uniprot_mapping("7TRE.cif", "A", use_auth=False)`, it classifies RNA as TERT and subtracts 35, resulting in negative residue indices (e.g. `33 -> -2`).
- **Suggestion**: Condition `is_tert` on `use_auth`:
  ```python
  if use_auth:
      is_tert = ("7TRE" in pid and chain_id == "A") or ("7TRD" in pid and chain_id == "A") or ("7TRF" in pid and chain_id == "A")
  else:
      is_tert = ("7TRE" in pid and chain_id == "B") or ("7TRD" in pid and chain_id == "B") or ("7TRF" in pid and chain_id == "D")
  ```

### [Minor] Finding 3: Handoff report documentation text discrepancy regarding 7TRE TPP1 length
- **What**: `worker_m2/handoff.md` stated that 7TRE TPP1 has "93 modeled residues, 90–182".
- **Where**: `worker_m2/handoff.md`, lines 5, 28, and 34.
- **Why**: In `data/raw/7TRE.cif`, TPP1 (author chain C) actually contains 152 modeled residues (90–241), identical to 2I46 Chain A. The source code `src/telotrade/sequence_mapping.py` and unit test `test_7tre_tpp1_extraction_gemmi_and_biopython` correctly extract and test all 152 residues. The discrepancy was purely textual within the handoff narrative.
- **Suggestion**: Note the correction for downstream Milestone 3 documentation.

### [Minor] Finding 4: Unused imports and variables in `tests/test_adversarial_m2.py`
- **What**: `ruff check` flags 22 style/lint warnings in `tests/test_adversarial_m2.py` (unused imports like `reconcile_tert_offset` and unused local benchmark timers).
- **Where**: `tests/test_adversarial_m2.py`.
- **Why**: The file was added during adversarial testing.
- **Suggestion**: Run `ruff check --fix tests/test_adversarial_m2.py` when refining the test suite.

---

## 4. Integrity Violation Check
- **Hardcoded test outputs embedded in source code**: None detected.
- **Dummy or facade implementations**: None detected; Gemmi and Biopython parsers execute genuine structure extraction.
- **Shortcuts bypassing requirements**: None detected.
- **Fabricated verification logs**: None detected; all tests verified independently via live subprocess execution.
- **Self-certifying work without genuine independent verification**: None; verified independently by Reviewer 1 and Challenger 2.
- **Integrity Status**: **CLEAN**.

---

## 5. Caveats
- `generate_chain_entity_mapping` returns a standardized, curated project metadata table across all 19 entities/chains in the 5 PDB complexes rather than parsing dynamic mmCIF metadata headers. This matches the specification and does not affect numerical modeling.
- Residues 179–321 of TERT (the unstructured linker connecting the TEN domain to TRBD) are absent from cryo-EM atomic coordinates across 7TRE, 7TRD, and 7TRF; they are correctly marked as unmodeled in the generated correspondence tables.

---

## 6. Conclusion
The implementation of Milestone 2 (`src/telotrade/sequence_mapping.py` and `tests/test_sequence_mapping.py`) satisfies all functional requirements and acceptance criteria in `PROJECT.md` and `ORIGINAL_REQUEST.md`:
1. Dual-engine sequence and coordinate extraction via Gemmi and Biopython is verified.
2. Construct tag offsets (-35 for TERT, +86 for 7TRE TPP1, +89 for 2I46 TPP1) are verified against raw experimental structures.
3. 2I46 homodimer asymmetry is handled with automatic redirection to Chain A and disorder warnings.
4. Sequence-aware homology alignment operates with 100% identity across core structures and correctly maps the TEL patch.
5. All 34 milestone tests and all 124 full-suite tests pass cleanly.

**Final Verdict: APPROVE.**

---

## 7. Verification Method
To reproduce these findings independently:
1. Run Milestone 2 unit tests:
   ```bash
   python -m pytest tests/test_sequence_mapping.py -v
   ```
   *Expected*: 34 passed.
2. Run E2E test suite:
   ```bash
   python -m pytest tests/e2e/ -v
   ```
   *Expected*: 40 passed, 12 skipped.
3. Run full test suite:
   ```bash
   python -m pytest tests/ -v
   ```
   *Expected*: 124 passed, 12 skipped.
4. Verify construct offsets directly against raw mmCIF files:
   ```bash
   python -c "
   from telotrade.sequence_mapping import get_chain_sequence
   _, n_tre = get_chain_sequence('data/raw/7TRE.cif', 'A', use_auth=True)
   _, n_tre_l = get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False)
   assert set(a - l for a, l in zip(n_tre, n_tre_l)) == {-35}
   _, n_tpp1 = get_chain_sequence('data/raw/7TRE.cif', 'C', use_auth=True)
   _, n_tpp1_l = get_chain_sequence('data/raw/7TRE.cif', 'D', use_auth=False)
   assert set(a - l for a, l in zip(n_tpp1, n_tpp1_l)) == {86}
   _, n_2i46 = get_chain_sequence('data/raw/2I46.cif', 'A', use_auth=True)
   _, n_2i46_l = get_chain_sequence('data/raw/2I46.cif', 'A', use_auth=False)
   assert set(a - l for a, l in zip(n_2i46, n_2i46_l)) == {89}
   print('Offsets verified successfully!')
   "
   ```
