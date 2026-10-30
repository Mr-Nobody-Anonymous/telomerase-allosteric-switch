# Environment & Tooling Survey Report

**Date**: 2026-10-16  
**Investigator**: Survey Agent 2 (`survey_env` — Environment & Tooling Explorer)  
**Workspace**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research`  
**Target Pipeline**: Telomerase TERT–TPP1 Structural Bioinformatics Pipeline (7TRE, 7TRD, 7TRF, 7TRC, 2I46)

---

## Executive Summary

The local environment is fully capable of executing the computational structural bioinformatics study end-to-end. Python 3.10.11 (64-bit AMD64) is installed with core numerical and scientific libraries (`numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `networkx`, `pytest`, `requests`). All 5 target mmCIF structures (7TRE, 7TRD, 7TRF, 7TRC, 2I46) and their validation metadata were tested and verified to be immediately downloadable over HTTPS from RCSB PDB and PDBe. The three structural bioinformatics libraries previously absent (`biopython`, `gemmi`, `freesasa`) have pre-compiled Windows wheels and were successfully installed and verified. Disk space is 3.84 GB free, which is more than sufficient for the ~5 MB of raw coordinate data and intermediate analysis matrices, though large external downloads should be managed conscientiously.

---

## 1. Python Runtime & Core Tooling

### 1.1 Python Binary & System Architecture
- **Python Executable**: `C:\Users\hp\AppData\Local\Programs\Python\Python310\python.exe`
- **Python Version**: `3.10.11 (tags/v3.10.11:7d4cc5a, Apr  5 2023, 00:38:17) [MSC v.1929 64 bit (AMD64)]`
- **Pip Version**: `26.2.1`
- **Operating System**: Windows 11 / Windows Server AMD64 (64-bit)
- **Version Control**: `git version 2.49.0.windows.1`
- **Test Runner**: `pytest 9.1.1`, `pytest-cov 7.1.0`

### 1.2 Package Status Matrix

| Package | Status | Version | Installation Location / Source | Pipeline Role |
| :--- | :--- | :--- | :--- | :--- |
| **numpy** | Pre-installed | `1.26.4` | `...\Python310\lib\site-packages\numpy` | Vectorized coordinate geometry & distance matrices |
| **scipy** | Pre-installed | `1.15.3` | `...\Roaming\Python\Python310\site-packages\scipy` | KDTree spatial search, Kabsch SVD / rotation |
| **pandas** | Pre-installed | `2.3.3` | `...\Roaming\Python\Python310\site-packages\pandas` | Data manifest, residue contacts, CSV outputs |
| **matplotlib** | Pre-installed | `3.10.7` | `...\Roaming\Python\Python310\site-packages\matplotlib` | Publication-grade plots & RMSD profiles |
| **seaborn** | Pre-installed | `0.13.2` | `...\Roaming\Python\Python310\site-packages\seaborn` | Heatmaps & contact frequency matrices |
| **requests** | Pre-installed | `2.32.5` | `...\Roaming\Python\Python310\site-packages\requests` | HTTPS downloads from RCSB PDB & REST API |
| **networkx** | Pre-installed | `3.1` | `...\Roaming\Python\Python310\site-packages\networkx` | Residue-residue interaction network analysis |
| **pytest** | Pre-installed | `9.1.1` | `...\Roaming\Python\Python310\site-packages\pytest` | Automated unit testing & regression checks |
| **biopython** | Verified / Installed | `1.88` | `...\Python310\lib\site-packages\Bio` | mmCIF parsing, sequence alignment, Superimposer |
| **gemmi** | Verified / Installed | `0.7.5` | `...\Python310\lib\site-packages\gemmi` | Fast mmCIF parsing, residue correspondence, symmetry |
| **freesasa** | Verified / Installed | `2.2.1` | `...\Python310\lib\site-packages\freesasa` | Buried Surface Area (BSA / SASA) calculation |

### 1.3 Pip Warning Note
During pip operations, a non-fatal warning is reported:
`WARNING: Ignoring invalid distribution -astmcp (c:\users\hp\appdata\local\programs\python\python310\lib\site-packages)`
Inspection confirmed this is an orphaned directory `~astmcp-4.0.3.dist-info` from a prior package update. It does not impede pip operations or library functionality.

---

## 2. Network Connectivity & RCSB PDB Data Ingestion

### 2.1 mmCIF Download Verification
Network endpoints for all 5 required structures were verified via direct HTTPS requests to `https://files.rcsb.org/download/{id}.cif`:

| Entry ID | File Type | HTTP Status | File Size (Bytes) | Verification Status |
| :--- | :--- | :--- | :--- | :--- |
| **7TRE** | Primary TERT–TPP1 complex (3.5 Å cryo-EM) | `200 OK` | 328,668 (~321 KB) | Verified |
| **7TRD** | Human telomerase catalytic reference (3.3 Å) | `200 OK` | 301,812 (~295 KB) | Verified |
| **7TRF** | Catalytic core RNP with H2A–H2B | `200 OK` | 338,967 (~331 KB) | Verified |
| **7TRC** | Human telomerase H/ACA RNP lobe | `200 OK` | 362,679 (~354 KB) | Verified |
| **2I46** | Unbound TPP1 OB-fold reference (X-ray) | `200 OK` | 69,582 (~68 KB) | Verified |

Total download payload for all 5 mmCIF files is ~1.4 MB, requiring less than 2 seconds of network transfer.

### 2.2 RCSB PDB Data REST API & Validation Reports
- **REST API (`https://data.rcsb.org/rest/v1/core/entry/{id}`)**:
  - Responds `200 OK` with complete structured metadata.
  - Contains embedded validation metrics:
    - `pdbx_vrpt_summary`: RNA suiteness, validation step history, report date.
    - `pdbx_vrpt_summary_em`: Cryo-EM Q-scores (`0.392` for 7TRE), all-atom inclusion (`0.794`), backbone inclusion (`0.866`), contour levels.
    - `pdbx_vrpt_summary_geometry`: Clashscore (`7.45` for 7TRE), Ramachandran outliers (`0.0%`), rotamer outliers (`0.11%`), bond/angle RMSZ (`0.29` / `0.43`).
- **Validation Report PDFs (`https://www.ebi.ac.uk/pdbe/entry-files/download/{id}_validation.pdf`)**:
  - Responds `200 OK` for all 5 structures:
    - `7tre`: 1,119,252 bytes
    - `7trd`: 1,165,270 bytes
    - `7trf`: 1,035,543 bytes
    - `7trc`: 1,258,368 bytes
    - `2i46`: 516,343 bytes

---

## 3. Storage & Workspace Permissions

### 3.1 Disk Capacity (Drive C:)
- **Total Capacity**: 476.07 GB
- **Used Space**: 472.23 GB
- **Available Free Space**: 3.84 GB (~3,932 MB)

**Capacity Assessment**: The entire research pipeline (5 raw mmCIF files, validation reports, parsed contact tables, RMSD matrices, and figures) will consume less than 20–30 MB total (< 1% of available free disk space). However, because free space is under 4 GB, caching bulky EM density maps (`.mrc` files are typically 200–500 MB each) should be avoided unless strictly necessary.

### 3.2 Directory Permissions
- Confirmed full read, write, and deletion privileges in `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research`.
- Confirmed write permissions in Python user site-packages (`C:\Users\hp\AppData\Roaming\Python\Python310\site-packages`) and global site-packages (`C:\Users\hp\AppData\Local\Programs\Python\Python310\lib\site-packages`).

---

## 4. Pipeline Execution Readiness Checklist

| Requirement | Tooling / Asset | Status | Notes |
| :--- | :--- | :--- | :--- |
| **R1. Data Acquisition** | `requests`, `hashlib`, RCSB PDB HTTPS | **READY** | URLs verified, mmCIF and validation data accessible |
| **R2. Sequence Alignment** | `biopython.Align`, `gemmi` | **READY** | `biopython==1.88` & `gemmi==0.7.5` verified |
| **R3. Interface Geometry** | `scipy.spatial.cKDTree`, `freesasa`, `numpy` | **READY** | `freesasa==2.2.1` verified; contact cutoff math tested |
| **R4. Cross-State RMSD** | `Bio.PDB.Superimposer`, `scipy.spatial` | **READY** | SVD Kabsch & Cα alignment ready |
| **R5. Synthesis & Figures** | `matplotlib`, `seaborn`, `pandas` | **READY** | Vector/PNG plotting and CSV table exporters ready |
| **Automated Testing** | `pytest` | **READY** | `pytest 9.1.1` operational for unit tests |

---

## 5. Recommended Implementation Commands

For downstream agents running the pipeline:

1. **Verify environment**:
   ```powershell
   python -c "import Bio, gemmi, freesasa, numpy, scipy, matplotlib, pandas, requests; print('Environment Ready')"
   ```
2. **Execute tests**:
   ```powershell
   pytest tests/
   ```
3. **Data Acquisition**:
   Use standard `requests.get` with streaming write to `data/raw/{id}.cif` and compute SHA-256 with `hashlib.sha256()`.
