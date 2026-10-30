# Dispatch for Milestone 1 Worker (`worker_m1`)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Read PROJECT.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md`

Also read the survey reports:
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md`
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\handoff.md`
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_env\handoff.md`
- `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\handoff.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1`

## Mandatory Integrity Warning
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Milestone 1 Scope & Deliverables
1. **Download & Archive Raw Data** (`data/raw/`):
   - Download official mmCIF files for `7TRE`, `7TRD`, `7TRF`, `7TRC`, and `2I46` from `https://files.rcsb.org/download/{id}.cif` into `data/raw/{id}.cif`.
   - Download validation reports (PDBe full validation PDF from `https://www.ebi.ac.uk/pdbe/entry-files/download/{pid}_full_validation.pdf` and/or RCSB validation XML) into `data/raw/`.
   - Compute SHA-256 checksums and verify they match authoritative values:
     - 7TRE: `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a`
     - 7TRD: `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4`
     - 7TRF: `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd`
     - 7TRC: `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea`
     - 2I46: `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779`
   - Store raw files in `data/raw/` and ensure they remain immutable.

2. **Immutable Data Manifest** (`data_manifest.csv`):
   - Compile `data_manifest.csv` (both in project root `data_manifest.csv` and `data/manifest/data_manifest.csv` for safety).
   - Columns: `pdb_id,experimental_method,resolution,sha256,source_url,chain_mappings,unmodeled_segments,uncertainty_notes`
   - Populate with exact details documented in `survey_pdb/report.md` (e.g. 7TRE cryo-EM 3.5 Å, Q-score 0.392, clashscore 7.45; 7TRD 3.3 Å; 7TRF 3.5 Å; 7TRC 3.4 Å; 2I46 X-ray 2.70 Å).

3. **Validation Engine** (`src/telotrade/validate.py` & `src/telotrade/__init__.py`):
   - Implement `src/telotrade/__init__.py`.
   - Implement `src/telotrade/validate.py` with robust CLI and Python API:
     - Check file existence and verify SHA-256 hashes against `data_manifest.csv`.
     - Parse mmCIF structures via Biopython (`MMCIFParser(QUIET=True)`) and Gemmi (`gemmi.cif.Document()`).
     - Execute cleanly without unhandled exceptions on missing residues, unmodeled loops, disordered atoms, or alternate locations (`altloc`).
     - Provide summary metrics extraction (number of models, chains, residues, atoms, missing residue count).
   - CLI syntax: `python -m telotrade.validate` or `python src/telotrade/validate.py`.

4. **Testing & Verification**:
   - Write unit tests in `tests/test_validate.py`.
   - Run tests: `python -m pytest tests/test_validate.py` and ensure 100% pass.
   - Run `python src/telotrade/validate.py` and document output.

Write your 5-component `handoff.md` in your directory and send a completion message to parent.


## 2026-10-16T14:39:25Z
[Message] timestamp=2026-10-16T14:39:25Z sender=eb46fecb-7562-4ce7-aa98-998216cee2c4 priority=MESSAGE_PRIORITY_HIGH content=You are the Milestone 1 Implementation Worker.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Read PROJECT.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md
Read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m1\DISPATCH.md
Read survey reports:
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md
- C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\handoff.md

MANDATORY INTEGRITY WARNING:
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

Scope:
1. Download mmCIF files for 7TRE, 7TRD, 7TRF, 7TRC, 2I46 into data/raw/ and verify SHA-256 hashes. Download validation PDFs/XMLs into data/raw/.
2. Create immutable data_manifest.csv (at project root and in data/manifest/).
3. Implement src/telotrade/__init__.py and src/telotrade/validate.py.
4. Write unit tests in tests/test_validate.py, run pytest and validate.py directly to verify 100% pass.
5. Write your 5-component handoff.md and send message to parent.
