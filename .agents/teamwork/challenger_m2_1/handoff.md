# Handoff Report — Challenger 1: Milestone 2 Adversarial Stress Testing (`challenger_m2_1`)

**Target Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\challenger_m2_1`  
**Date**: 2026-10-16  
**Type**: Hard Handoff (Adversarial Empirical Evaluation Complete)  
**Verdict**: **REQUEST_CHANGES**

---

## 1. Observation

### 1.1 Test Suite & Verification Results
- **Milestone 2 Worker Suite**:
  `python -m pytest tests/test_sequence_mapping.py -v`:
  `34 passed in 20.96s` (100% pass on baseline developer tests).
- **Adversarial Stress Test Suite** (`tests/test_adversarial_m2.py`):
  `python -m pytest tests/test_adversarial_m2.py -v`:
  `21 passed in 9.70s` (100% pass across all 6 challenge dimensions).
- **Full Repository Test Suite**:
  `python -m pytest tests/ -v`:
  `145 passed, 12 skipped in 239.18s (0:03:59)`. Zero failures, zero errors, zero regressions.
- **Code Linter**:
  `python -m ruff check src tests/`:
  `All checks passed!` (0 lint errors).

---

### 1.2 Empirical Stress-Testing Findings

#### Finding 1 (High): Unhandled `ValueError: sequence has zero length` in `align_sequences()` on Empty Sequence
- **Location**: `src/telotrade/sequence_mapping.py:715`:
  ```python
  def align_sequences(seq1: str, seq2: str, mode: str = "global") -> AlignmentResult:
      ...
      alignments = aligner.align(seq1, seq2)
      if not alignments:
          return AlignmentResult(...)
  ```
- **Observed Exception**:
  When `seq1` or `seq2` has length 0 (e.g. unmodeled chain or missing coordinates in a mmCIF file):
  ```
  Traceback (most recent call last):
    File "<string>", line 1, in <module>
    File "src/telotrade/sequence_mapping.py", line 715, in align_sequences
      alignments = aligner.align(seq1, seq2)
    File "Bio/Align/__init__.py", line 4501, in align
      score, paths = super().align(sA, sB, strand)
  ValueError: sequence has zero length
  ```
- **Analysis**:
  The developer attempted to handle non-alignments on line 716 (`if not alignments:`), but calling Biopython's `aligner.align()` with an empty string raises `ValueError` before line 716 is reached.

---

#### Finding 2 (Critical): Cross-Entity Chain Conflation in `get_uniprot_mapping()` Producing Negative / Corrupted UniProt Numbers
- **Location**: `src/telotrade/sequence_mapping.py:551-576`:
  ```python
  is_tert = (
      ("7TRE" in pid and chain_id in ["A", "B"])
      or ("7TRD" in pid and chain_id in ["A", "B"])
      or ("7TRF" in pid and chain_id in ["A", "D"])
  )
  is_tpp1 = (
      ("7TRE" in pid and chain_id in ["C", "D"])
      or ("2I46" in pid and chain_id in ["A", "B"])
  )
  ```
- **Observed Behavior**:
  The function attempts to support both author and label chain IDs by checking set inclusion (`["A", "B"]`, `["C", "D"]`, `["A", "D"]`) without checking `use_auth`:
  1. In `7TRE`:
     - Author Chain A is TERT, Label Chain B is TERT.
     - Author Chain B is hTR RNA, Label Chain A is hTR RNA.
     - When `use_auth=False`, requesting Chain `A` (which is RNA) matches `chain_id in ["A", "B"]` -> flagged as `is_tert=True`.
     - TERT offset `num - 35` is applied to RNA residues:
       ```python
       get_uniprot_mapping("data/raw/7TRE.cif", "A", use_auth=False)
       # Returns: [(33, -2), (34, -1), (35, 0), (36, 1), ...]
       ```
       Generating illegal **negative** UniProt residue numbers for RNA!
  2. In `7TRE`:
     - Author Chain C is TPP1, Label Chain D is TPP1.
     - Author Chain D is DNA, Label Chain C is DNA.
     - When `use_auth=False`, requesting Chain `C` (DNA) matches `chain_id in ["C", "D"]` -> flagged as `is_tpp1=True`.
     - TPP1 offset `num + 86` is applied to DNA residues:
       ```python
       get_uniprot_mapping("data/raw/7TRE.cif", "C", use_auth=False)
       # Returns: [(12, 98), (13, 99), (14, 100), ..., (18, 104)]
       ```
       Corrupting telomeric DNA substrate residues 12-18 into TPP1 UniProt numbers 98-104!
  3. In `7TRF`:
     - Author Chain A is TERT, Label Chain D is TERT.
     - Label Chain A is Histone H2B.
     - When `use_auth=False`, requesting Chain `A` matches `chain_id in ["A", "D"]` -> flagged as `is_tert=True`, mapping Histone H2B residues (label 36..38) to TERT UniProt positions (1..3)!

---

#### Finding 3 (Critical): Biopython Parser Branch in `_get_chain_sequence_biopython()` Inverts Extracted Macromolecules Under `use_auth=False`
- **Location**: `src/telotrade/sequence_mapping.py:317-321`:
  ```python
  def _get_chain_sequence_biopython(cif_path: Path, chain_id: str, use_auth: bool = True):
      ...
      if use_auth:
          if chain_id in model:
              target_chain = model[chain_id]
      else:
          # Look up by label_asym_id in Biopython
          for ch in model:
              if ch.id == chain_id:
                  target_chain = ch
                  break
  ```
- **Observed Behavior**:
  In Biopython's `MMCIFParser`, `ch.id` is `auth_asym_id`, and `r.id[1]` is `auth_seq_num`.
  Comparing `parser='gemmi'` vs `parser='biopython'` on `7TRE.cif` with `use_auth=False` and `chain_id='B'`:
  - `gemmi` returns TERT protein: 934 residues, `seq='CRAVRSLLRS...'`, `res_nums=[42, 43, 44...]`.
  - `biopython` returns hTR RNA: 217 nucleotides, `seq='GGCCAUUUUU...'`, `res_nums=[33, 34, 35...]`!
  And for `chain_id='D'`, `gemmi` returns TPP1 protein (152 aa), while `biopython` returns DNA (7 nt)!
- **Analysis**:
  Biopython's `MMCIFParser` does not index chains or residues by `label_asym_id` / `label_seq_id`. Line 319 looks up `ch.id == chain_id`, which compares against `auth_asym_id`. Biopython completely inverts the extracted macromolecule and returns author residue numbers.

---

#### Finding 4 (Medium): Offset Reconcilers Accept Out-of-Range/Negative Values and Generate Illegal Negative CIF Labels
- **Location**: `src/telotrade/sequence_mapping.py:455-516`:
  `reconcile_tert_offset`, `reconcile_tpp1_offset`, `reconcile_label_from_auth`
- **Observed Behavior**:
  - `reconcile_label_from_auth(auth_seq_num=10, protein="TPP1", structure_id="7TRE")` returns `-76`.
  - `reconcile_label_from_auth(auth_seq_num=10, protein="TPP1", structure_id="2I46")` returns `-79`.
  - `reconcile_tert_offset(label_seq_id=10)` (in Twin-Strep tag 1..35) returns `-25`.
- **Analysis**:
  In mmCIF, `label_seq_id` is an index strictly $\ge 1$. The functions perform naive unchecked linear arithmetic without bounds checking or domain clamping.

---

#### Finding 5 (Medium): Global Alignment Mode Stretches Internal Gaps on Truncated Sequences with Matching Terminal Residues
- **Location**: `src/telotrade/sequence_mapping.py:708-715`:
- **Observed Behavior**:
  When aligning a truncated query whose C-terminal residue matches the C-terminal residue of the reference (e.g. `s_full = "MSTNPKPQRKWTHRLVRPWIEEVRAL"`, `s_trunc = s_full[:15]` which ends in 'L'), `align_sequences(mode="global")` aligns query index 14 ('L') to reference index 25 ('L'), opening an internal gap of 11 residues across `VRPWIEEVRA` and leaving reference index 14 unaligned.
- **Analysis**:
  In global alignment mode, end gaps are penalized equally to internal gaps. Matching the terminal character yields $+2.0$, which exceeds the end gap extension penalty, causing the aligner to stretch an internal gap.

---

#### Finding 6 (Low): `map_residue_correspondence()` Silently Returns Spurious Mappings on Non-Homologous Sequences
- **Location**: `src/telotrade/sequence_mapping.py:774-864`:
- **Observed Behavior**:
  Aligning 7TRE TERT (Chain A) against 7TRF Histone H2B (Chain F) yields 2.2% sequence identity and score $-477.5$. `map_residue_correspondence` returns a 90-residue dictionary of spurious pairings without warning or error.

---

### 1.3 Verified Robust Implementations
Despite the findings above, several critical components demonstrated exceptional robustness under adversarial stress:
1. **2I46 Chain A vs Chain B Dimer Asymmetry**:
   - `check_tel_patch_integrity` correctly identifies that Chain A has all TEL-patch residues modeled (152 residues, DE-loop 168-171 intact), while Chain B has DE-loop 168-171 missing coordinates (134 residues total, missing 126-139 and 168-171).
   - Pairwise alignment between 2I46 Chain A and Chain B correctly maps all 134 structured residues with **ZERO index drift** (`r1 == r2` for all 134 positions).
   - `select_valid_tpp1_chain` strictly redirects Chain B to Chain A with a `UserWarning` when `require_tel_patch=True`.
2. **Homologous Cross-Structure Correspondence**:
   - 7TRE TERT vs 7TRD TERT and 7TRF TERT align with 100% identity across all 934 modeled residues with zero drift.
   - 7TRE TPP1 (152 modeled residues) vs 2I46 Chain A (152 modeled residues) aligns with 100% identity across all positions 90-241.
3. **Resilience to Mutations and Gaps**:
   - Alanine scanning mutations (e.g. hexamutant in DE-loop) maintain 1:1 residue correspondence without spurious gaps.
   - Internal single and double insertions/deletions preserve flanking and downstream coordinate indices.

---

## 2. Logic Chain

1. **Premise 1 (Claimed Zero-Defect Sequence Mapping & Disambiguation)**:
   In `worker_m2/handoff.md`, the implementation is claimed to *"deliver robust, high-performance sequence extraction, chain disambiguation, construct tag reconcilers, 2I46 asymmetry handling, Biopython pairwise homology alignment, and standardized table generators"* with 100% pass rate.
2. **Premise 2 (Empirical Reproduction of Cross-Entity Conflation in `get_uniprot_mapping`)**:
   By querying `get_uniprot_mapping` with `use_auth=False` on RNA (Chain A in 7TRE) and DNA (Chain C in 7TRE), we empirically observed that RNA was mapped to negative TERT UniProt indices `(33, -2), (34, -1), (35, 0)` and DNA was mapped to TPP1 UniProt indices `(12, 98), ..., (18, 104)`.
3. **Premise 3 (Empirical Reproduction of Biopython Parser Inversion)**:
   By querying `get_chain_sequence('7TRE.cif', 'B', use_auth=False, parser='biopython')`, Biopython returned RNA nucleotides (`GGCCAUUUUU`, 217 nt), whereas `parser='gemmi'` returned TERT protein (`CRAVRSLLRS`, 934 aa), completely breaking cross-parser parity.
4. **Premise 4 (Empirical Reproduction of `ValueError` on Empty Sequence)**:
   Calling `align_sequences("", "ACDEF")` or `align_sequences("", "")` raised an uncaught `ValueError: sequence has zero length`.
5. **Conclusion from Logic Chain**:
   While the primary target protein chains (TERT Chain A and TPP1 Chain C in 7TRE, Chain A in 2I46) align cleanly when `use_auth=True` with Gemmi, `src/telotrade/sequence_mapping.py` contains severe logic bugs in `get_uniprot_mapping()` and `_get_chain_sequence_biopython()` under `use_auth=False`, as well as crash vulnerabilities on empty sequence inputs. Therefore, Milestone 2 cannot be approved without remediation.

---

## 3. Caveats

1. **Primary Target Analysis (Milestone 3 Unblocked for Baseline)**:
   Because Milestone 3 (Interface Geometry) primarily operates on `7TRE.cif` TERT (author chain A) and TPP1 (author chain C) using Gemmi with `use_auth=True`, the core calculations for contact networks and buried surface area are not directly broken by Finding 2 or Finding 3 as long as default author numbering and Gemmi are used.
2. **2I46 Chain A TEL-Patch**:
   2I46 Chain A is thoroughly validated and confirmed to have the complete TEL patch (168-171, 183, 212) resolved. Downstream comparison between bound TPP1 (7TRE) and unbound TPP1 (2I46 Chain A) is structurally sound.

---

## 4. Conclusion & Actionable Recommendations

**Verdict**: **REQUEST_CHANGES**

To resolve this verdict, the worker must apply the following targeted fixes:

1. **Fix `get_uniprot_mapping()` to Condition on `use_auth`**:
   In `src/telotrade/sequence_mapping.py:551-576`, replace naive letter lists with strict entity/chain mappings:
   - For 7TRE:
     - TERT: `chain_id == "A"` if `use_auth` else `chain_id == "B"`.
     - TPP1: `chain_id == "C"` if `use_auth` else `chain_id == "D"`.
   - For 7TRD:
     - TERT: `chain_id == "A"` if `use_auth` else `chain_id == "B"`.
   - For 7TRF:
     - TERT: `chain_id == "A"` if `use_auth` else `chain_id == "D"`.
   - For other chains (RNA, DNA, histones), return identity or raise `ValueError(f"Chain {chain_id} is not a supported TERT/TPP1 polymer")`.
2. **Fix `_get_chain_sequence_biopython()` for `use_auth=False`**:
   Biopython `MMCIFParser` structures do not store `label_asym_id` in `ch.id` or `label_seq_id` in `r.id[1]`.
   Either:
   - Read `_atom_site.label_asym_id` and `_atom_site.label_seq_id` from the parser dictionary, or
   - Delegate `use_auth=False` extraction to Gemmi, or emit a clear warning / error that Biopython MMCIFParser only supports author coordinates.
3. **Fix `align_sequences()` Empty Sequence Check**:
   Add an upfront check in `align_sequences`:
   ```python
   if len(seq1) == 0 or len(seq2) == 0:
       return AlignmentResult(
           seq1=seq1, seq2=seq2, mapping={}, reverse_mapping={},
           aligned_pairs=[], identity=0.0, score=0.0,
           unaligned_res1=list(range(len(seq1))),
           unaligned_res2=list(range(len(seq2))),
           alignment_str="",
       )
   ```
4. **Add Input Clamping / Validation to Offset Reconcilers**:
   Validate that `auth_seq_num` and `label_seq_id` are positive integers within construct bounds, or clamp to prevent negative CIF label generation.

---

## 5. Verification Method

To independently reproduce all observations and verify the challenge findings:

1. **Run Milestone 2 Adversarial Test Suite**:
   ```bash
   python -m pytest tests/test_adversarial_m2.py -v
   ```
   *Expected result*: 21 passed in ~10 seconds.
2. **Reproduce Finding 1 (Empty Sequence Crash)**:
   ```bash
   python -c "from telotrade.sequence_mapping import align_sequences; align_sequences('', 'ACDEF')"
   ```
   *Expected result*: `ValueError: sequence has zero length`.
3. **Reproduce Finding 2 (Cross-Entity Conflation in `get_uniprot_mapping`)**:
   ```bash
   python -c "from telotrade.sequence_mapping import get_uniprot_mapping; print(list(get_uniprot_mapping('data/raw/7TRE.cif', 'A', use_auth=False).items())[:3])"
   ```
   *Expected result*: `[(33, -2), (34, -1), (35, 0)]` (negative UniProt numbers on RNA).
4. **Reproduce Finding 3 (Biopython Parser Entity Inversion)**:
   ```bash
   python -c "from telotrade.sequence_mapping import get_chain_sequence; print('Gemmi:', get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False, parser='gemmi')[0][:10]); print('Biopython:', get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False, parser='biopython')[0][:10])"
   ```
   *Expected result*:
   `Gemmi: CRAVRSLLRS` (TERT protein) vs `Biopython: GGCCAUUUUU` (RNA nucleotide).
5. **Run Full Test Suite**:
   ```bash
   python -m pytest tests/ -v
   ```
   *Expected result*: `145 passed, 12 skipped in ~4 minutes`.
