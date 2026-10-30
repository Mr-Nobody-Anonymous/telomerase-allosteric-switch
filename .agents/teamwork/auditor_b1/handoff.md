# Forensic Audit Report: Milestone B1 Work Products

**Work Product**: Milestone B1 Deliverables (`src/telotrade/sequence_mapping.py`, `results/tables/residue_correspondence_tpp1.csv`, `results/audit/evidence_ledger.csv`, `results/audit/repository_inventory.md`, `results/literature/prior_art_matrix.csv`, `results/literature/claim_verification.md`, `results/literature/novelty_gap_report.md`)  
**Auditor**: `auditor_b1` (Milestone B1 Forensic Auditor)  
**Profile**: General Project (Integrity Mode: **Benchmark**)  
**Verdict**: **CLEAN**  

---

### Phase Results
- **Hardcoded Output Detection**: PASS — 0 instances of hardcoded test results, expected output literals, mock values, or dummy strings in `src/telotrade/`.
- **Facade Implementation Detection**: PASS — 0 facade functions, 0 empty stubs, 0 `NotImplementedError` placeholders. All calculations derive dynamically from mmCIF parsing and linear algebra.
- **Pre-populated Artifact Detection**: PASS — 0 stale `.log` files or unverified pre-generated outputs; CSV serialization matches dynamic in-memory generation 100%.
- **Cryptographic Raw Data Integrity**: PASS — 5/5 raw mmCIF files match authoritative SHA-256 hashes recorded in `data_manifest.csv` with zero divergence.
- **Empirical 3D Coordinate Concordance**: PASS — `residue_correspondence_tpp1.csv` matches 100% of modeled residues in `7TRE.cif` Chain C (152/152), `2I46.cif` Chain A (152/152), and `2I46.cif` Chain B (134/134, correctly identifying DE-loop disorder at residues 168–171).
- **Literature & Prior-Art Rigor**: PASS — `prior_art_matrix.csv`, `claim_verification.md`, and `novelty_gap_report.md` strictly designate canonical contacts (`Arg774–Glu168` and `Lys78–Glu215`) as Positive Controls, separate static RMSD (0.338 Å) from dynamic allostery, deconstruct sequential binding claims, and quarantine TELIN-1 as an unvalidated hypothetical proposal.
- **Dependency Audit (Benchmark Mode)**: PASS — Core algorithms built from scratch using permitted foundational numerical libraries (`numpy`, `pandas`, `gemmi`, `biopython`, `freesasa`); no prohibited third-party black-box tools or core execution delegation.
- **Build & Test Suite Execution**: PASS — `ruff check src/ tests/` (0 errors); `pytest tests/ -v` (309 passed, 0 failed, 2 warnings in 212.42s).

---

## 1. Observation

1. **Inspection of Bugfix in `src/telotrade/sequence_mapping.py`**:
   - Location: `src/telotrade/sequence_mapping.py:1145-1165`
   - Code:
     ```python
     # Canonical TPP1 OB-fold sequence from 7TRE Entity 4
     # Recombinant construct has a 3-aa N-terminal tag (Gly87-Ala88-Gly89 in entity numbering).
     # Stripping the 3-aa tag yields ob_seq starting at residue 90 (Ser90 = ob_seq[0]).
     # UniProt Q96AP0 residue 89 corresponds to full_poly[2] (Gly89).
     full_poly = get_polymer_entity_sequence(cif_7tre, entity_id_or_chain=4)
     ob_seq = full_poly[3:] if len(full_poly) >= 163 else full_poly

     rows: list[dict[str, Any]] = []

     # UniProt Q96AP0 residues 89 to 249 (161 positions)
     for resnum in range(89, 250):
         if resnum == 89:
             resname_1 = full_poly[2] if len(full_poly) >= 3 else "G"
         else:
             idx = resnum - 90
             resname_1 = ob_seq[idx] if 0 <= idx < len(ob_seq) else "X"
         resname_3 = one_to_three_aa(resname_1)
     ```
   - Grep searches across `src/telotrade/`:
     - Pattern `mock`: 0 results.
     - Pattern `dummy`: 0 results.
     - Pattern `hardcode`: 0 results.
     - Pattern `NotImplementedError`: 0 results.

2. **Empirical Concordance of `results/tables/residue_correspondence_tpp1.csv`**:
   - Tested directly against 3D atomic coordinates from raw mmCIF files:
     ```powershell
     # 7TRE Chain C (152 modeled residues)
     Checked 152 residues against 7TRE. Mismatches: 0
     All 7TRE modeled residues match perfectly!

     # 2I46 Chain A and Chain B
     2I46 Chain A: Checked 152 residues. Mismatches: 0
     2I46 Chain B: Checked 134 residues. Mismatches: 0
     ```
   - Inspection of canonical residues 165–174:
     - Residue 165: SER (auth 165, modeled=True)
     - Residue 166: ASP (auth 166, modeled=True)
     - Residue 167: TRP (auth 167, modeled=True)
     - Residue 168: GLU (auth 168, modeled=True, TEL-patch=True, 2I46 B=False)
     - Residue 169: GLU (auth 169, modeled=True, TEL-patch=True, 2I46 B=False)
     - Residue 170: LYS (auth 170, modeled=True, TEL-patch=True, 2I46 B=False)
     - Residue 171: GLU (auth 171, modeled=True, TEL-patch=True, 2I46 B=False)
     - Residue 172: PHE (auth 172, modeled=True, 2I46 B=True)
     - Residue 173: GLY (auth 173, modeled=True, 2I46 B=True)
     - Residue 174: PHE (auth 174, modeled=True, 2I46 B=True)
   - Dynamic re-execution vs. Disk CSV:
     - `generate_tpp1_correspondence_table()` was executed in memory and compared with `results/tables/residue_correspondence_tpp1.csv`.
     - Result: `CSVs are 100% IDENTICAL!`.

3. **Cryptographic Validation of Raw mmCIF Files**:
   - Recalculated SHA-256 for all 5 files in `data/raw/` vs `results/tables/data_manifest.csv`:
     - `7TRE.cif`: `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` -> MATCH (OK)
     - `7TRD.cif`: `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` -> MATCH (OK)
     - `7TRF.cif`: `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` -> MATCH (OK)
     - `7TRC.cif`: `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` -> MATCH (OK)
     - `2I46.cif`: `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` -> MATCH (OK)

4. **Audited Deliverables Content Verification**:
   - `results/audit/evidence_ledger.csv`:
     - 24 discrete claims cataloged.
     - Claims CLM-005, CLM-010, CLM-014, and CLM-017 classified as `uncertain` due to 3.5 Å coordinate precision / map compaction.
     - Claim CLM-022 (static 0.338 Å RMSD proves absence of allostery) classified as `contradicted`.
     - Claims CLM-021 (sequential induced-fit binding) and CLM-023 (TELIN-1 inhibitor) classified as `untested`.
     - Claim CLM-024 documents resolution of the TPP1 sequence naming defect as `established`.
   - `results/audit/repository_inventory.md`:
     - Accurately records Git repository status as uninitialized.
     - Comprehensive inventory of 8 source modules, 12 test files, 21 tables, and 4 figures.
   - `results/literature/prior_art_matrix.csv`:
     - 22 structured entries cross-referencing Liu et al. (Nature 2022), Sekne et al. (Science 2022), Wang et al. (Nature 2007), and Nandakumar et al. (Nature 2012).
     - PA-001 (Arg774–Glu168) and PA-002 (Lys78–Glu215) explicitly classified as `Positive Control`.
     - PA-008, PA-010, PA-019 classified as `Unvalidated Proposal`.
     - Genuine novelties (multi-threshold contact matrix PA-014, domain decomposition PA-015, Lys170 remodeling peak PA-012) accurately isolated.
   - `results/literature/claim_verification.md` & `results/literature/novelty_gap_report.md`:
     - Fully articulated academic critique establishing evidence boundaries, experimental error envelopes, and de-hyping static model extrapolations.

5. **Linting and Test Suite Execution**:
   - `ruff check src/ tests/`: `All checks passed!`
   - `python -m pytest tests/ -v`: `309 passed, 2 warnings in 212.42s (0:03:32)`.

---

## 2. Logic Chain

1. **Premise**: In Benchmark mode, all code must genuinely implement required functions without shortcuts, mocks, hardcoded test strings, or unverified outputs.
   - *Observation*: Source analysis revealed dynamic parsing via Gemmi mmCIF readers, explicit indexing math, and zero hardcoded test literals.
2. **Premise**: The Gen1 off-by-one residue correspondence error shifted `resname_canonical` relative to UniProt Q96AP0 and 3D coordinate numbers.
   - *Observation*: Adjusting the loop index to `idx = resnum - 90` for `resnum >= 90` and mapping residue 89 to `full_poly[2]` yields 0 mismatches across all 152 residues in 7TRE, 152 in 2I46 A, and 134 in 2I46 B.
3. **Premise**: Regenerated artifacts on disk must be genuinely generated by the verified codebase, not manually crafted or pre-populated.
   - *Observation*: Direct programmatic regeneration in memory and comparison against `results/tables/residue_correspondence_tpp1.csv` proved 100% identity.
4. **Premise**: Ground-truth user constraints in `ORIGINAL_REQUEST.md` (R1 and R2) require separating established positive controls from novel findings and tracing every numerical claim to source code and limitations.
   - *Observation*: `evidence_ledger.csv`, `repository_inventory.md`, `prior_art_matrix.csv`, `claim_verification.md`, and `novelty_gap_report.md` fulfill every acceptance criterion with complete factual rigor.
5. **Conclusion**: All Milestone B1 deliverables satisfy all integrity, structural, and behavioral requirements. The work product is authentic and uncompromised.

---

## 3. Caveats

1. **Unmodeled Disordered Loops**: Residues 89 and 242–249 of TPP1 are absent from electron density in 7TRE; their sequence identities derive from the recombinant entity sequence in mmCIF block `_entity_poly.pdbx_seq_one_letter_code`.
2. **Phenix Coordinate Compaction**: As noted in CLM-005, the 2.32 Å distance between TERT Arg774 and TPP1 Glu168 in 7TRE reflects cryo-EM real-space refinement compaction rather than a physically relaxed solution bond length.
3. **Uninitialized Git Workspace**: As documented in `repository_inventory.md`, the workspace lacks a `.git` tracking directory. Checksums and automated manifest checks serve as the primary cryptographic provenance baseline.

---

## 4. Conclusion

The work products submitted for Milestone B1 pass all forensic integrity, code quality, data consistency, and behavioral tests. 
- **Verdict**: **CLEAN**.
- **Action**: Milestone B1 is officially verified and approved. The project is cleared to proceed to Milestone B2 (Falsifiable Hypotheses & Independent Structural Recalculation).

---

## 5. Verification Method

To independently reproduce and verify this audit:

1. **Verify Raw Data Cryptographic Integrity**:
   ```powershell
   python -c "import hashlib, pandas as pd; df = pd.read_csv('results/tables/data_manifest.csv'); [(print(r['pdb_id'], 'OK') if hashlib.sha256(open(f'data/raw/{r.pdb_id}.cif', 'rb').read()).hexdigest() == r['sha256'] else print(r['pdb_id'], 'FAIL')) for _, r in df.iterrows()]"
   ```
   *Expected result*: All 5 entries print `OK`.

2. **Verify 3D Coordinate Concordance of TPP1 Table**:
   ```powershell
   python -c "import gemmi, pandas as pd; st = gemmi.read_structure('data/raw/7TRE.cif'); chain_c = st[0]['C']; res_dict = {r.seqid.num: r.name for r in chain_c if r.name != 'HOH'}; df = pd.read_csv('results/tables/residue_correspondence_tpp1.csv'); mismatches = [int(r['uniprot_resnum']) for _, r in df.iterrows() if int(r['uniprot_resnum']) in res_dict and r['resname_canonical'] != res_dict[int(r['uniprot_resnum'])]]; print('Mismatches:', len(mismatches))"
   ```
   *Expected result*: `Mismatches: 0`.

3. **Verify Linter and Pytest Suite**:
   ```powershell
   ruff check src/ tests/
   python -m pytest tests/ -v
   ```
   *Expected result*: 0 linter errors; 309 tests passing.

4. **Invalidation Conditions**:
   - Any mismatch between `residue_correspondence_tpp1.csv` and mmCIF atom coordinates.
   - Any test failure in `pytest tests/`.
   - Any SHA-256 hash discrepancy in `data/raw/`.
