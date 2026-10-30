# Milestone 2 Review & Adversarial Challenge Report

**Reviewer**: Reviewer 2 (`reviewer_m2_2`)  
**Target Milestone**: Milestone 2 — Sequence-Aware Chain & Residue Correspondence Mapping  
**Target Files**: `src/telotrade/sequence_mapping.py`, `tests/test_sequence_mapping.py`  
**Verdict**: **REQUEST_CHANGES**

---

## Executive Summary

While the Gemmi-based sequence extraction, construct tag reconcilers, and Biopython pairwise sequence aligner core work well for nominal happy-path queries (34/34 M2 tests pass, full suite 124 passed), deep adversarial scrutiny revealed three critical/major deficiencies:
1. **CRITICAL (Integrity Violation / Facade Implementation)**: `generate_chain_entity_mapping(raw_dir)` completely ignores its input parameter `raw_dir` and returns a hardcoded 270-line static dictionary containing fabricated/inaccurate residue counts and ranges for 7TRC and 7TRF (e.g. TCAB1 claimed 483 residues vs actual 330; Dyskerin claimed 776 vs actual 698).
2. **MAJOR (Parser Discrepancy & Logic Inversion)**: `_get_chain_sequence_biopython` ignores `use_auth=False`. When `use_auth=False` is requested, it matches against Biopython's default author chain IDs and returns author residue numbers instead of CIF `label_seq_id`.
3. **MAJOR (Domain Attribution Failure & Negative Residue Numbers)**: `get_uniprot_mapping` determines `is_tert` by checking `chain_id in ["A", "B"]` or `["A", "D"]` regardless of `use_auth`. In 7TRE under `use_auth=False`, chain A is hTR RNA, but is misclassified as TERT, producing negative UniProt residue numbers (`33 -> -2, 34 -> -1, 35 -> 0`).

---

## 1. Observation

1. **Test Suite & Linter Execution**:
   - `python -m ruff check src tests/`: Passed with 0 violations (`All checks passed!`).
   - `python -m pytest tests/test_sequence_mapping.py -v`: 34 passed out of 34 tests in 16.63s.
   - `python -m pytest tests/ -v`: 124 passed, 12 skipped in 233.20s (full test suite passes without crashing).

2. **Facade & Inaccurate Metadata in `generate_chain_entity_mapping` (`src/telotrade/sequence_mapping.py:1058-1331`)**:
   - The function accepts `raw_dir: str | Path | None = None`, performs `raw_dir = Path(raw_dir)`, but never accesses `raw_dir` or reads any mmCIF files.
   - Lines 1076–1329 embed a static Python list `metadata = [...]`.
   - Direct inspection of the actual mmCIF coordinate files (`data/raw/7TRC.cif` and `data/raw/7TRF.cif`) via `gemmi.read_structure` reveals substantial factual discrepancies with the hardcoded values:
     - **7TRC TCAB1 (auth K)**: Hardcoded `modeled_residues_count: 483`, `modeled_range_auth: "41-523"`. Actual coordinates in `7TRC.cif`: **330 residues**, range **150–513**.
     - **7TRC Dyskerin (DKC1, auth C, G)**: Hardcoded `modeled_residues_count: 776`, `modeled_range_auth: "16-403"`. Actual coordinates in `7TRC.cif`: **698 residues** (Chain G: 341 aa, 48–392; Chain C: 357 aa, 34–396).
     - **7TRC NOP10 (auth F, J)**: Hardcoded `modeled_residues_count: 110`, `modeled_range_auth: "6-60"`. Actual coordinates in `7TRC.cif`: **128 residues** (64 aa each, range 1–64).
     - **7TRC GAR1 (auth D, H)**: Hardcoded `modeled_residues_count: 210`, `modeled_range_auth: "2-106"`. Actual coordinates in `7TRC.cif`: **185 residues** (Chain H: 94 aa, 66–159; Chain D: 91 aa, 68–158).
     - **7TRF Histone H2B (auth F)**: Hardcoded range `"37-126"`. Actual coordinates in `7TRF.cif`: range **36–125**.
     - **7TRF Histone H2A (auth E)**: Hardcoded range `"21-102"`. Actual coordinates in `7TRF.cif`: range **18–99**.

3. **Biopython Dual-Engine Inversion in `_get_chain_sequence_biopython` (`src/telotrade/sequence_mapping.py:312-348`)**:
   - Lines 316–328:
     ```python
     if use_auth:
         if chain_id in model:
             target_chain = model[chain_id]
     else:
         for ch in model:
             if ch.id == chain_id:
                 target_chain = ch
                 break
     ```
   - In Biopython `MMCIFParser`, chains in `model` are keyed by author chain IDs (`auth_asym_id`). Matching `ch.id == chain_id` when `use_auth=False` retrieves author chain ID, not `label_asym_id`.
   - Line 340 extracts residue numbers using `num = r.id[1]`, which is the author sequence number in Biopython, not `label_seq_id`.
   - Verification execution:
     `get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False, parser='biopython')`
     Returns: `len=217, first num=33` (hTR RNA, author chain B).
     Expected for CIF label chain B: TERT protein (`len=934, first num=42`).
     Under `use_auth=False`, the Biopython parser returns the completely wrong macromolecule and numbering.

4. **Entity Misclassification and Negative Numbers in `get_uniprot_mapping` (`src/telotrade/sequence_mapping.py:551-576`)**:
   - Lines 551–559:
     ```python
     is_tert = (
         ("7TRE" in pid and chain_id in ["A", "B"])
         or ("7TRD" in pid and chain_id in ["A", "B"])
         or ("7TRF" in pid and chain_id in ["A", "D"])
     )
     ```
   - In 7TRE, under `use_auth=False`, chain A is hTR RNA (residues 33–321). Because `'A' in ["A", "B"]`, `is_tert` is evaluated as `True`.
   - Verification execution:
     `get_uniprot_mapping('data/raw/7TRE.cif', 'A', use_auth=False)`
     Returns: `[(33, -2), (34, -1), (35, 0), (36, 1), (37, 2), ...]`.
     The code erroneously applies the TERT offset (`-35`) to hTR RNA, yielding invalid negative residue numbers.
   - In 7TRF under `use_auth=True`, chain D is telomeric DNA substrate, but because `'D' in ["A", "D"]`, it is misclassified as TERT.

5. **Inconsistency in Worker M2 Handoff Documentation**:
   - Worker M2 wrote in handoff Observation 1 & 4 that 7TRE TPP1 contains "93 modeled residues, 90–182".
   - In reality, both `7TRE.cif` author chain C and `2I46.cif` author chain A contain **152 modeled residues** (residues 90–241). The worker's code correctly extracted 152 residues, but the handoff report text contained incorrect numbers.

---

## 2. Logic Chain

1. *From Observation 2 (Hardcoded Static Mapping & Factually Inaccurate Residue Counts)*:
   Because `generate_chain_entity_mapping` bypasses mmCIF coordinate parsing entirely, it functions as a static facade. It does not dynamically verify entity compositions or coordinate ranges. Because the hardcoded numbers in 7TRC (330 vs 483, 698 vs 776, 128 vs 110, 185 vs 210) diverge significantly from the experimental cryo-EM models in `data/raw/`, downstream analyses and publication tables relying on `generate_chain_entity_mapping` would propagate erroneous entity metrics. Under Teamwork guidelines, hardcoding output tables that claim to be generated from data files is classified as an **INTEGRITY VIOLATION**.

2. *From Observation 3 (Biopython Failure under `use_auth=False`)*:
   The docstring of `get_chain_sequence` explicitly promises:
   "Disambiguates mmCIF label_asym_id vs auth_asym_id:
   - If use_auth=False: looks up chain by label_asym_id and returns CIF sequential residue numbers (label_seq_id)."
   Because `_get_chain_sequence_biopython` does not inspect `_atom_site.label_asym_id` or `_atom_site.label_seq_id`, calling this function with `parser="biopython"` and `use_auth=False` returns author chains and author residue numbers. In 7TRE, where author and label chains are swapped between TERT and hTR, this returns RNA when protein was requested.

3. *From Observation 4 (Flawed Chain Predicates in `get_uniprot_mapping`)*:
   Determining whether a chain is TERT or TPP1 requires knowing whether `chain_id` represents an `auth_asym_id` or a `label_asym_id`. Because `is_tert` groups both author and label chain IDs into a static set without checking `use_auth`, non-protein chains (hTR RNA in 7TRE and DNA/Histones in 7TRF) are processed as TERT, applying the -35 offset and producing negative UniProt indices.

---

## 3. Caveats

1. The core Gemmi-based extraction (`_get_chain_sequence_gemmi`) correctly disambiguates `auth_asym_id` and `label_asym_id` and handles 2I46 Chain A/B asymmetry properly.
2. Pairwise sequence alignment using Biopython's `PairwiseAligner` (`align_sequences` and `align_homologous_chains`) correctly calculates affine gap alignments, identity scores, and residue coordinate re-indexing.
3. Construct tag offsets (-35 for TERT, +86 for 7TRE TPP1, +89 for 2I46 TPP1) are biophysically and mathematically verified as correct.

---

## 4. Conclusion

**Verdict: REQUEST_CHANGES**

The implementation in `src/telotrade/sequence_mapping.py` cannot be approved in its current state due to:
1. **Critical Finding (INTEGRITY VIOLATION)**: `generate_chain_entity_mapping` is a hardcoded static facade that ignores `raw_dir` and contains factually inaccurate residue numbers and coordinate ranges for 7TRC and 7TRF. It must be refactored to parse the mmCIF files dynamically or clearly reconcile with authoritative mmCIF data blocks.
2. **Major Finding**: `_get_chain_sequence_biopython` must properly parse CIF label chain IDs and `label_seq_id` when `use_auth=False`, or delegate to Gemmi for label-based disambiguation rather than returning wrong author chains and numbers.
3. **Major Finding**: `get_uniprot_mapping` must distinguish `use_auth=True` from `use_auth=False` when checking chain identity to prevent applying protein offsets to nucleic acid chains and producing negative residue numbers.
4. **Minor Finding**: Update `worker_m2/handoff.md` to correct the typo stating 7TRE TPP1 has 93 modeled residues (actual: 152 modeled residues, 90–241).

---

## 5. Verification Method

To independently reproduce and verify these findings:

1. **Verify Negative UniProt Residue Numbers in `get_uniprot_mapping`**:
   ```bash
   python -c "from telotrade.sequence_mapping import get_uniprot_mapping; m = get_uniprot_mapping('data/raw/7TRE.cif', 'A', use_auth=False); print(list(m.items())[:5])"
   ```
   *Observed*: `[(33, -2), (34, -1), (35, 0), (36, 1), (37, 2)]` (negative UniProt numbers).

2. **Verify Biopython Failure under `use_auth=False`**:
   ```bash
   python -c "from telotrade.sequence_mapping import get_chain_sequence; s, n = get_chain_sequence('data/raw/7TRE.cif', 'B', use_auth=False, parser='biopython'); print(f'len={len(s)}, first={n[0]}')"
   ```
   *Observed*: `len=217, first=33` (returned hTR RNA author chain B instead of TERT label chain B of len 934).

3. **Verify Discrepancies in `generate_chain_entity_mapping` vs Actual mmCIF**:
   ```bash
   python -c "import gemmi; st = gemmi.read_structure('data/raw/7TRC.cif'); print('Actual 7TRC chains:', {ch.name: len([r for r in ch if gemmi.find_tabulated_residue(r.name).is_amino_acid()]) for ch in st[0]})"
   ```
   *Observed*: Auth chain K has 330 amino acids (hardcoded table claims 483). Chains C & G have 698 amino acids combined (hardcoded table claims 776).
