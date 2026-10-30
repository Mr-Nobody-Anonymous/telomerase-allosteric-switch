# Handoff Report: Repository Inventory & Evidence Ledger Audit (Requirement R1 / Milestone B1)

**From**: `explorer_b1_repo` (Codebase & Artifact Auditor)  
**To**: `orchestrator_gen2` / Downstream Milestone Agents (`spec_miner_b1_lit`, `worker_b2`, etc.)  
**Date**: 2026-10-17T10:31:00Z  
**Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b1_repo`  
**Deliverables Published**:
- `results/audit/evidence_ledger.csv`
- `results/audit/repository_inventory.md`
- Local copies: `report.md`, `handoff.md`, `evidence_ledger.csv`, `repository_inventory.md` in working directory

---

## 1. Observation

1. **Git Repository Status**:
   - Command executed: `git status`
   - Direct output: `fatal: not a git repository (or any of the parent directories): .git`
   - Observation: Workspace is uninitialized under Git; no commit history or branch tracking exists.

2. **Raw Structural Data Integrity**:
   - Path: `data/raw/` containing `7TRE.cif`, `7TRD.cif`, `7TRF.cif`, `7TRC.cif`, and `2I46.cif`.
   - Manifest: `data_manifest.csv` and `data/manifest/data_manifest.csv`.
   - Hashes computed via `hashlib.sha256`:
     - `7TRE.cif` (1,518,438 bytes): `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` (Match)
     - `7TRD.cif` (1,394,494 bytes): `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` (Match)
     - `7TRF.cif` (1,542,948 bytes): `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` (Match)
     - `7TRC.cif` (1,647,704 bytes): `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` (Match)
     - `2I46.cif` (273,104 bytes): `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` (Match)
   - All 5 raw files match expected SHA-256 hashes 100%.

3. **Automated Test Execution**:
   - Command: `pytest`
   - Execution log output: `300 passed, 2 warnings in 598.86s (0:09:58)`
   - 100% of the 300 tests passed across E2E and adversarial test suites.

4. **Code Defect in Sequence Correspondence**:
   - File: `src/telotrade/sequence_mapping.py`, lines 1146–1156.
   - Code snippet:
     ```python
     full_poly = get_polymer_entity_sequence(cif_7tre, entity_id_or_chain=4)
     ob_seq = full_poly[3:] if len(full_poly) >= 163 else full_poly
     for resnum in range(89, 250):
         idx = resnum - 89
         resname_1 = ob_seq[idx] if 0 <= idx < len(ob_seq) else "X"
     ```
   - Inspection of `7TRE.cif` Entity 4: Slicing `[3:]` removes the 3-aa GAG tag, leaving `ob_seq[0]` as residue Gly90. Residue 89 was not part of the recombinant construct. For `resnum = 89`, `idx = 0` (Gly90), shifting the entire canonical sequence column `resname_canonical` in `results/tables/residue_correspondence_tpp1.csv` by +1 residue.
   - Example in `results/tables/residue_correspondence_tpp1.csv`:
     - Line 78: `165,ASP` (correct: SER165)
     - Line 79: `166,TRP` (correct: ASP166)
     - Line 80: `167,GLU` (correct: TRP167)
     - Line 81: `168,GLU` (correct: GLU168)
     - Line 82: `169,LYS` (correct: GLU169)
     - Line 83: `170,GLU` (correct: LYS170)
     - Line 84: `171,PHE` (correct: GLU171)

5. **Coordinate Anomaly in Primary Salt Bridge**:
   - File: `results/tables/polar_interactions_7tre.csv`, line 2:
     `Salt Bridge,774,ARG,NH1,TRAP_IFD,168,GLU,OE1,2.32,True`
   - Observation: Distance is 2.32 Å. Sum of Bondi vdW radii for N (1.55 Å) and O (1.52 Å) is 3.07 Å. Normal non-covalent salt bridge distances range between 2.7 and 3.2 Å. 2.32 Å indicates severe coordinate compaction in the deposited 3.5 Å cryo-EM model.

6. **Allostery and Kinetic Mechanism Claims in `report.md`**:
   - Section 1 & 4.4 claims: "catalytic core (residues 326–1132) is exceptionally rigid, exhibiting a Cα root-mean-square displacement (RMSD) of only 0.338 Å... without any propagating allosteric distortion into the RNA-template duplex".
   - Section 4.5 claims: "induced-fit recruitment mechanism: TPP1 presents a pre-organized, rigid hydrophobic docking knob (Leu183), which guides initial encounter docking... followed by conformational collapse of the flexible, acidic DE-loop".
   - Section 5 & Table 4 references: TELIN-1 candidate inhibitor.

---

## 2. Logic Chain

1. **Premise**: If raw mmCIF hashes match `data_manifest.csv`, raw coordinates are untampered and reproducible.
   - *Supported by Observation 2*: 100% hash match verified.
2. **Premise**: If `sequence_mapping.py` uses `idx = resnum - 89` on an array starting at residue 90, column `resname_canonical` has an off-by-one error.
   - *Supported by Observation 4*: Direct inspection of `7TRE.cif` entity 4 confirms `ob_seq[0]` is Gly90. Thus residue 166 was assigned residue 167's name (TRP).
3. **Premise**: Coordinate extraction in `interface_geometry.py` and `cross_state.py` accesses `r.name` and atom coordinates directly from the coordinate block, bypassing `residue_correspondence_tpp1.csv`.
   - *Supported by Inspection of `interface_geometry.py:209-236` and `cross_state.py:218-235`*: Both modules read directly from mmCIF atom records; therefore, interface contacts and RMSD calculations are unaffected by the sequence table naming defect.
4. **Premise**: A 2.32 Å N...O distance between heavy atoms is shorter than normal physical van der Waals boundaries.
   - *Supported by Observation 5*: At 3.5 Å resolution, cryo-EM map coordinate uncertainty is ~0.35–0.5 Å. The 2.32 Å distance is computationally present in the deposited CIF but represents a local coordinate compaction artifact during Phenix real-space refinement.
5. **Premise**: Static structural models with coordinate uncertainty ~0.35–0.5 Å cannot establish lack of allostery or prove a sequential kinetic mechanism.
   - *Supported by Observation 6 and physical principles*: Dynamic allostery (entropy, vibrational modes) can occur without static backbone displacement. Static crystal/cryo-EM coordinates represent frozen endpoints and cannot determine temporal kinetics without biophysical assays (e.g. stopped-flow, single-molecule FRET).
6. **Premise**: TELIN-1 has zero wet-lab data and zero citations in pharmacological databases.
   - *Supported by literature audits*: Searches in PubMed, ChEMBL, and PubChem yield 0 hits. TELIN-1 is purely an unvalidated computational candidate.

---

## 3. Caveats

1. **Coordinate Refinement Independence**: Models 7TRE, 7TRD, and 7TRF were refined by the same research group using identical starting models and Phenix real-space refinement protocols. Some coordinate identity in rigid cores (RMSD = 0.338 Å) is influenced by refinement starting-model inertia rather than strictly independent experimental observations.
2. **X-Ray Lattice Effects in 2I46**: 2I46 Chain A DE-loop coordinates are stabilized by crystal packing contacts in space group $I 4_1 2 2$, whereas Chain B DE-loop is disordered. Thus, the 4.293 Å displacement relative to 7TRE reflects a combination of binding adaptation, crystal packing, and high thermal mobility ($B > 80\ \text{Å}^2$).
3. **No Code Modification Undertaken**: In accordance with the explorer archetype's read-only mandate, no source code in `src/telotrade/` was modified. The off-by-one indexing defect was documented for downstream implementers.

---

## 4. Conclusion

1. **Requirement R1 is 100% complete**:
   - `results/audit/evidence_ledger.csv` has been generated and populated with 24 audited claims tracing metric names, reported values, source files, source scripts, execution parameters, limitations, and classifications (17 established, 4 uncertain, 2 contradicted, 2 untested).
   - `results/audit/repository_inventory.md` has been compiled documenting Git status, system environment, raw data hashes, source code inventory, 21 tables, 4 figures, and test suite baselines.
2. **Key Corrective Guidance for Downstream Agents**:
   - **Worker B2 / Independent Recalculation**: Recalculate BSA, contacts, and RMSDs using alternative libraries (e.g. BioPython vs Gemmi; alternative probe radii 1.2–1.6 Å; multi-cutoff thresholds 3.5/4.0/4.5 Å) to stress-test the 2.32 Å salt bridge and 0.338 Å RMSD.
   - **Worker B3 / Mechanistic Reports**: Explicitly reframe the DE-loop report to decouple static coordinate displacement from kinetic sequence, and reframe the allostery report to decouple static RMSD from dynamic/entropic allosteric coupling. Strictly quarantine TELIN-1 as an unvalidated hypothetical candidate.
   - **Worker B5 / Final Report Rewrite**: Correct the off-by-one sequence correspondence in `results/tables/residue_correspondence_tpp1.csv` and rewrite `report.md` adhering to the required 13-section structure with an explicit catalog of rejected/overstated claims.

---

## 5. Verification Method

1. **Verify Deliverable Files Exist**:
   ```powershell
   Get-Item results\audit\evidence_ledger.csv
   Get-Item results\audit\repository_inventory.md
   ```
2. **Verify SHA-256 Checksums and Baseline Tests**:
   ```powershell
   python -m telotrade.validate
   python -m pytest tests/e2e/test_tier1_features.py -v
   ```
3. **Verify Off-by-One Code Bug**:
   Inspect `src/telotrade/sequence_mapping.py` at line 1154 and compare row 79 of `results/tables/residue_correspondence_tpp1.csv` against `7TRE.cif` Chain C residue 166.
4. **Invalidation Conditions**:
   - Any raw CIF hash in `data/raw/` deviating from `data_manifest.csv`.
   - Any of the 24 claims in `evidence_ledger.csv` failing mathematical reproducibility against deposited CIFs.
