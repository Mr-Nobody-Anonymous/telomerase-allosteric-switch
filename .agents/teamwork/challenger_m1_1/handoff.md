# Handoff Report — Challenger 1: Milestone 1 Adversarial Stress Testing (`challenger_m1_1`)

**Target Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\challenger_m1_1`  
**Date**: 2026-10-16  
**Type**: Hard Handoff (Adversarial Empirical Evaluation Complete)  
**Verdict**: **REQUEST_CHANGES**

---

## 1. Observation

### 1.1 Independent Data Integrity & Checksum Verification
An independent SHA-256 verification was conducted across all files in `data/raw/` against `data_manifest.csv`:

```
=== INDEPENDENT SHA-256 VERIFICATION ===
7TRE: size=1518438 bytes, expected=2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a, actual=2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a, match=True
7TRD: size=1394494 bytes, expected=faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4, actual=faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4, match=True
7TRF: size=1542948 bytes, expected=e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd, actual=e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd, match=True
7TRC: size=1647704 bytes, expected=99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea, actual=99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea, match=True
2I46: size=273104 bytes, expected=4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779, actual=4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779, match=True

All mmCIF checksums match manifest: True
Root manifest matches data/manifest/data_manifest.csv identically: True
All 5 full validation PDFs and 5 validation XML files exist in data/raw/: True
```

### 1.2 Baseline Test Suite Execution
Running the existing unit test suite:
- `python -m pytest tests/test_validate.py -v`:
  `27 passed in 48.58s` (100% pass on baseline files).
- `python -m ruff check src tests/`:
  `All checks passed!` (0 lint errors).

### 1.3 Adversarial Stress Testing Observations (`tests/test_adversarial_m1.py`)
Executing the 23-test adversarial stress harness revealed two reproducible unhandled exceptions:

#### Finding 1: Unhandled `IndexError: invalid vector subscript` on Zero-Byte / Corrupted CIF Files
When `parse_structure_summary()` or `parse_validation_metrics()` encounters an empty (0-byte) CIF file, Gemmi's C++ binding raises `IndexError: invalid vector subscript`.

- **Location**: `src/telotrade/validate.py:367`:
  ```python
  structure = gemmi.read_structure(str(path))
  ```
  and `src/telotrade/validate.py:297`:
  ```python
  doc = gemmi.cif.read_file(str(cif_path))
  ```
- **Caught Exception Tuples**:
  - Line 326: `except (ValueError, KeyError, RuntimeError, OSError) as e:` (omits `IndexError`)
  - Line 424: `except (RuntimeError, ValueError, OSError, AttributeError) as e:` (omits `IndexError`)
  - Line 553: `except (FileNotFoundError, ValueError, RuntimeError, OSError) as e:` (omits `IndexError`)
- **Impact & Traceback**:
  When `validate_all()` or CLI `python -m telotrade.validate` is invoked on a dataset where any file is 0-byte or corrupted, `validate_all()` crashes completely:
  ```
  Traceback (most recent call last):
    File "src/telotrade/validate.py", line 552, in main
      results = validate_all(raw_dir=raw_dir, manifest_path=manifest_path)
    File "src/telotrade/validate.py", line 485, in validate_all
      sum_data = parse_structure_summary(cif_path)
    File "src/telotrade/validate.py", line 367, in parse_structure_summary
      structure = gemmi.read_structure(str(path))
  IndexError: invalid vector subscript
  ```

#### Finding 2: Unhandled `KeyError: '_atom_site.id'` on Header-Only / Coordinate-Empty CIF Files
When `parse_structure_summary()` is given a valid CIF format file that contains no atomic coordinate block (e.g. an empty data block `data_EMPTY`), Biopython's `MMCIFParser` raises `KeyError: '_atom_site.id'`.

- **Location**: `src/telotrade/validate.py:430`:
  ```python
  bio_struct = parser.get_structure(entry_id, str(path))
  ```
  which calls `Bio/PDB/MMCIFParser.py:154`:
  ```python
  atom_serial_list = mmcif_dict["_atom_site.id"]
  KeyError: '_atom_site.id'
  ```
- **Caught Exception Tuple in `validate.py:442`**:
  ```python
  except (Bio.PDB.PDBExceptions.PDBConstructionException, ValueError, OSError, AttributeError) as e:
  ```
  Notice: `KeyError` is omitted! Consequently, `KeyError` bubbles uncaught through `parse_structure_summary`, crashing the calling code.

---

## 2. Logic Chain

1. **Premise 1 (Claimed Zero-Crash Guarantee)**:
   In `worker_m1/handoff.md` and `PROJECT.md` Feature 3, the implementation is specified as having *"zero-crash handling of missing residues, unmodeled loops, disordered atoms, and alternate locations (altlocs)"* and that `src/telotrade/validate.py` executes cleanly without unhandled exceptions.
2. **Premise 2 (Empirical Reproduction of Crashes)**:
   By feeding edge-case inputs (0-byte `.cif` file and coordinate-empty `.cif` block) to `parse_structure_summary`, `parse_validation_metrics`, `validate_all`, and the CLI `main()`, we empirically observed unhandled `IndexError` (Gemmi vector subscript) and `KeyError` (Biopython atom table lookup) crashes.
3. **Premise 3 (Root Cause Analysis)**:
   Python's exception hierarchy does not map `IndexError` or `KeyError` (subclasses of `LookupError`) to `RuntimeError`, `ValueError`, `OSError`, or `AttributeError`. Because the `try/except` clauses in lines 326, 424, 442, and 553 explicitly specify narrow tuples that exclude `IndexError` and `KeyError`, the zero-crash guarantee is broken when non-ideal or corrupted files enter the pipeline.
4. **Conclusion from Logic Chain**:
   While the dataset files currently on disk (`data/raw/*.cif`) parse cleanly, the validation engine's defensive error-handling contract is fragile. A single corrupted download or 0-byte file will cause uncaught traceback termination rather than clean failure reporting.

---

## 3. Caveats

- **Authentic Dataset Integrity**: No data integrity flaws were found in the 5 downloaded target files (`7TRE`, `7TRD`, `7TRF`, `7TRC`, `2I46`), validation PDFs, or XML reports. All bit-level hashes match authoritative specs.
- **Normal Flow Operability**: For valid coordinate files, the parser correctly handles unmodeled loops (e.g. 7TRE TERT 179–321), missing TEL-patch loops (2I46 Chain B 168–171), and alternate conformations.
- **Scope Scoping**: The failure modes only manifest under malformed, empty, or truncated CIF inputs.

---

## 4. Conclusion & Verdict

**Verdict**: **REQUEST_CHANGES**

The baseline data acquisition and manifest compilation are sound and pass all integrity checks. However, `src/telotrade/validate.py` requires a narrow, essential robustness fix to fulfill its specification of zero-crash handling on malformed inputs before advancing to Milestone 2.

### Actionable Remediation Required from `worker_m1`:
In `src/telotrade/validate.py`:
1. In `parse_validation_metrics()` (around line 326):
   Add `IndexError` to the caught exception tuple:
   ```python
   except (ValueError, KeyError, RuntimeError, OSError, IndexError) as e:
       metrics["cif_error"] = str(e)
   ```
2. In `parse_structure_summary()` (around line 424):
   Add `IndexError` and `KeyError` (or catch `Exception`) for Gemmi:
   ```python
   except (RuntimeError, ValueError, OSError, AttributeError, IndexError, KeyError) as e:
       summary["gemmi_status"] = f"ERROR: {e}"
   ```
3. In `parse_structure_summary()` (around line 442):
   Add `KeyError` and `IndexError` (or catch `Exception`) for Biopython:
   ```python
   except (Bio.PDB.PDBExceptions.PDBConstructionException, ValueError, OSError, AttributeError, KeyError, IndexError) as e:
       summary["biopython_status"] = f"ERROR: {e}"
   ```
4. In `main()` (around line 553):
   Expand the catch to `(FileNotFoundError, ValueError, RuntimeError, OSError, IndexError, KeyError)` or general `Exception` so the CLI never crashes with an unhandled traceback.

---

## 5. Verification Method

To independently reproduce the findings and verify subsequent fixes:

1. **Run Full Adversarial Stress Suite**:
   ```powershell
   python -m pytest tests/test_adversarial_m1.py -v
   ```
   *Expected Failure before fix*: `2 failed, 21 passed` (fails on `test_zero_byte_file_fails_cleanly` and `test_valid_syntax_empty_datablock_handled`).  
   *Expected Pass after fix*: `23 passed in ~85s`.

2. **Direct CLI 0-Byte Crash Reproduction**:
   ```powershell
   python -c "
   import tempfile, shutil, subprocess, sys
   from pathlib import Path
   from telotrade.validate import TARGET_PDB_IDS

   with tempfile.TemporaryDirectory() as tmpdir:
       raw_dir = Path(tmpdir) / 'raw'
       raw_dir.mkdir()
       shutil.copy('data_manifest.csv', Path(tmpdir) / 'data_manifest.csv')
       for pid in TARGET_PDB_IDS:
           shutil.copy(Path('data/raw') / f'{pid}.cif', raw_dir / f'{pid}.cif')
       (raw_dir / '7TRE.cif').write_bytes(b'')
       p = subprocess.run([sys.executable, '-m', 'telotrade.validate', '--raw-dir', str(raw_dir), '--manifest', str(Path(tmpdir) / 'data_manifest.csv')], capture_output=True, text=True)
       print('Return code:', p.returncode)
       print('Traceback in stderr:', 'IndexError' in p.stderr)
   "
   ```
