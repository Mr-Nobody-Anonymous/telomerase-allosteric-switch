# Authoritative Structural Specification & PDB Survey Report
**Target Entries**: PDB `7TRE`, `7TRD`, `7TRF`, `7TRC`, and `2I46`  
**Date**: 2026-10-16  
**Investigator**: PDB & Structural Specification Miner (`survey_pdb`)  
**Mission**: Authoritative specification mining of mmCIF schemas, chain entity mappings, sequence numbering conventions, TEL-patch residue correspondences, experimental parameters, and validation quality metrics to establish baseline data integrity for the human telomerase structural bioinformatics study.

---

## Executive Summary

An exhaustive empirical survey was conducted by directly downloading and probing the official mmCIF coordinate files, wwPDB validation reports, EMDB map records, and RCSB/PDBe REST API endpoints for all five target structures (`7TRE`, `7TRD`, `7TRF`, `7TRC`, `2I46`).

Key findings with critical pipeline impact:
1. **Discrepancy in Chain ID Conventions (`label_asym_id` vs `auth_asym_id`)**:
   In `7TRE`, the primary catalytic core complex, the canonical mmCIF chain ID (`label_asym_id`) for TERT is **`B`**, but the author chain ID (`auth_asym_id`) is **`A`**. Conversely, telomerase RNA (hTR) has `label_asym_id` = **`A`** and `auth_asym_id` = **`B`**. TPP1 has `label_asym_id` = **`D`** and `auth_asym_id` = **`C`**. Structural alignment algorithms and contact analysis pipelines **must strictly declare** whether they operate on `label_asym_id` or `auth_asym_id` to prevent swapping protein and RNA coordinates.
2. **TERT Sequence Numbering Offset (+35 aa Twin-Strep Tag)**:
   In `7TRE`, `7TRD`, and `7TRF`, TERT was expressed with an N-terminal 35-residue Twin-Strep purification tag (`pdb_seq_num` -34 to 0, `seq_id` 1 to 35). As a result, canonical UniProt human TERT (O14746) numbering corresponds exactly to author residue numbering (`auth_seq_num`), whereas CIF `label_seq_id` is offset by **+35**:
   $$\text{label\_seq\_id} = \text{auth\_seq\_num} + 35$$
   The first modeled residue in TERT is `Cys7` (`auth_seq_num = 7`, `label_seq_id = 42`).
3. **TPP1 TEL-Patch Residue Modeling & 2I46 Asymmetry**:
   - In **`7TRE`** (bound TPP1), the entire TEL-patch region (`Glu168`, `Glu169`, `Lys170`, `Glu171`, `Leu183`, `Leu212`) is fully resolved and modeled (`auth_seq_num = label_seq_id + 86`).
   - In **`2I46`** (unbound TPP1 crystal structure, dimer in ASU), **Chain A has the TEL-patch DE-loop (`Glu168–171`) fully modeled**, whereas **in Chain B, residues 168–171 are completely disordered/unobserved** (`auth_seq_num = '?'`). All comparative structural analysis between bound TPP1 (7TRE) and unbound TPP1 (2I46) **must strictly target Chain A** of 2I46.
4. **Resolution Clarification for `2I46`**:
   Authoritative crystallographic data confirm `2I46` was determined at **2.70 Å resolution** (R_work = 0.2135, R_free = 0.2509; space group $P2_12_12_1$), correcting the 1.7 Å resolution mention in early documentation drafts.
5. **Validation PDF Endpoints**:
   Direct programmatic access to wwPDB Full Validation Report PDFs is blocked (HTTP 403) on `files.rcsb.org` for automated user agents, but is **100% accessible with HTTP 200 via the official PDBe endpoint** (`https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_full_validation.pdf`). Official validation XMLs are accessible via `files.rcsb.org`.

---

## Features Discovered

| # | Category | Feature | Description | Inputs | Outputs | Error Behavior | Discovered Via |
|---|----------|---------|-------------|--------|---------|----------------|----------------|
| 1 | Download Endpoints | RCSB mmCIF Download | Direct HTTP retrieval of full coordinate file in mmCIF format | PDB ID (e.g. `7TRE`) | mmCIF plain text stream (`.cif`) | HTTP 404 if invalid PDB ID | RCSB REST API / HTTP HEAD/GET probes |
| 2 | Download Endpoints | RCSB Compressed mmCIF | Gzip-compressed mmCIF archive for rapid download | PDB ID | Gzipped mmCIF archive (`.cif.gz`) | HTTP 404 if invalid PDB ID | RCSB REST API probe |
| 3 | Download Endpoints | PDBe Full Validation Report PDF | Complete official wwPDB validation report PDF document | PDB ID (lowercase) | PDF binary stream | HTTP 404 if invalid PDB ID | PDBe REST API probe |
| 4 | Download Endpoints | RCSB Validation XML | Machine-readable validation metrics (XML format) | PDB ID, middle 2 chars | Gzipped XML (`_validation.xml.gz`) | HTTP 404 if invalid PDB ID | RCSB FTP/HTTP server probe |
| 5 | Download Endpoints | EMDB Map Volume | 3D Cryo-EM primary reconstructed electron density map | EMDB ID (e.g. `EMD-26087`) | Gzipped CCP4/MRC map (`.map.gz`) | HTTP 404 if invalid EMDB ID | EBI EMDB FTP/HTTP server probe |
| 6 | Entity Mapping | 7TRE Macromolecular Entities | 4 distinct entities: hTR (1), TERT (2), DNA substrate (3), TPP1 (4) | Entry 7TRE | Entities with lengths 451 nt, 1167 aa, 18 nt, 163 aa | Malformed entity if index > 4 | RCSB Core Polymer Entity API & `_entity` CIF loop |
| 7 | Entity Mapping | 7TRD Macromolecular Entities | 3 entities: hTR (1), TERT (2), DNA substrate (3); TPP1 is ABSENT | Entry 7TRD | Entities with lengths 451 nt, 1167 aa, 18 nt | KeyError if querying TPP1 in 7TRD | RCSB Core Polymer Entity API & `_entity` CIF loop |
| 8 | Entity Mapping | 7TRF Macromolecular Entities | 5 entities: Histone H2B (1), Histone H2A (2), hTR (3), TERT (4), DNA (5) | Entry 7TRF | H2B (126 aa), H2A (129 aa), hTR (451 nt), TERT (1167 aa), DNA (18 nt) | Chain mismatch if assuming TERT is entity 2 | RCSB Core Polymer Entity API & `_entity` CIF loop |
| 9 | Entity Mapping | 7TRC Macromolecular Entities | 6 entities: TCAB1 (1), DKC1 (2), NOP10 (3), hTR (4), GAR1 (5), NHP2 (6) | Entry 7TRC | H/ACA lobe ribonucleoprotein subunits | TERT/TPP1 absent (only H/ACA lobe) | RCSB Core Polymer Entity API & `_entity` CIF loop |
| 10 | Entity Mapping | 2I46 Macromolecular Entities | 1 entity: TPP1 / ACD OB-fold domain (homodimer in ASU: chains A & B) | Entry 2I46 | Polypeptide (161 aa, UniProt Q96AP0) | Single entity with 2 instances (chains A, B) | RCSB Core Polymer Entity API & `_entity` CIF loop |
| 11 | Chain Mapping | CIF vs Author Chain Discrepancy | Non-identical chain IDs between mmCIF canonical (`label_asym_id`) and legacy PDB (`auth_asym_id`) | Structure and chain | Exact mapping pair (e.g. 7TRE TERT: label `B` = auth `A`) | Silent coordinate mismatch if using wrong ID | `_struct_asym` and `_atom_site` CIF analysis |
| 12 | Numbering Schemes | TERT Twin-Strep Tag Offset | 35-aa N-terminal purification tag in recombinant construct | TERT residue query | `auth_seq_num = label_seq_id - 35`; `auth_seq_num` matches UniProt O14746 | Index shift of 35 if using `label_seq_id` | `_pdbx_poly_seq_scheme` in 7TRE, 7TRD, 7TRF |
| 13 | Numbering Schemes | 7TRE TPP1 GAG Cloning Offset | 3-aa N-terminal GAG cloning artifact (residues 87-89 in construct) | TPP1 residue query in 7TRE | `auth_seq_num = label_seq_id + 86`; `auth_seq_num` matches UniProt Q96AP0 | 3-residue shift relative to 2I46 label numbering | `_pdbx_poly_seq_scheme` in 7TRE |
| 14 | Numbering Schemes | 2I46 TPP1 Numbering Scheme | Unbound crystal construct starts directly at Ser90 | TPP1 residue query in 2I46 | `auth_seq_num = label_seq_id + 89`; `auth_seq_num` matches UniProt Q96AP0 | Label seq IDs differ between 7TRE and 2I46 | `_pdbx_poly_seq_scheme` in 2I46 |
| 15 | Structural Epitope | TEL-Patch Residue Coordinates | Critical recruitment residues in TPP1: Glu168, Glu169, Lys170, Glu171, Leu183, Leu212 | Residue identifier | 3D atomic coordinates in 7TRE chain D (auth C) and 2I46 chain A | Chain B in 2I46 has missing coordinates | `_atom_site` coordinate parsing |
| 16 | Structural Epitope | TERT Receptor Cradle Contacts | Direct atomic contact network between TPP1 TEL-patch and TERT | Distance threshold ≤ 4.5 Å | 32 residue-residue pairs (Glu168-Arg774, Glu171-Ala45, Leu183-Leu769, etc.) | Empty set if chains misidentified | 3D Euclidean distance calculation on 7TRE |
| 17 | Validation Metrics | Cryo-EM Map-Model Inclusion | Atom inclusion and real-space correlation scores (Q-score, all-atom / backbone inclusion) | Cryo-EM mmCIF & EMDB volume | 7TRE: Q=0.392, atom inclusion=0.794; 7TRD: Q=0.447; 7TRF: Q=0.369; 7TRC: Q=0.498 | N/A for X-ray (2I46) | `pdbx_vrpt_summary_em` validation records |
| 18 | Validation Metrics | Geometric Quality & Outliers | Clashscore, Ramachandran favored/outliers, rotamer outliers, RNA suiteness | PDB structure | 7TRE: Clash=7.45, Rama=0.0%, Rot=0.11%, RNA=0.36; 2I46: Clash=12.97, Rama=0.36% | N/A | `pdbx_vrpt_summary_geometry` validation records |
| 19 | RNA Architecture | hTR Transcript Domain Segregation | Catalytic core structures (7TRE/D/F) model pseudoknot & CR4/CR5; 7TRC models H/ACA lobe | Structure entry | 7TRE/D/F: nt 33-149, 163-191, 251-321 (217 nt); 7TRC: nt 211-218, 362-392, 397-449 (92 nt) | Zero overlap between catalytic core and H/ACA models | `_pdbx_poly_seq_scheme` RNA segment analysis |
| 20 | Substrate Complex | Telomeric DNA Substrate Alignment | 18-nt single-stranded DNA substrate with 7 resolved base pairs at 3' end | Structures 7TRE, 7TRD, 7TRF | 7 modeled nucleotides (`12:DT` to `18:DG`) paired with RNA template | Author chain ID is `D` in 7TRE/7TRF, but `N` in 7TRD | `_atom_site` and `_pdbx_poly_seq_scheme` |

---

## Edge Cases

| # | Feature | Input | Observed Behavior |
|---|---------|-------|-------------------|
| 1 | 2I46 TEL-patch Disordered Asymmetry | Chain B of `2I46` queried for residues 168–171 | **Residues are unobserved / missing** (`auth_seq_num = '?'`). Residues Glu168, Glu169, Lys170, Glu171 are only resolved in **Chain A**. Downstream comparative alignment must select Chain A or fall back to an exception. |
| 2 | TERT Recombinant Tag Numbering | Querying TERT catalytic residues using mmCIF `label_seq_id` | **35-residue index shift**: The Twin-Strep tag occupies `seq_id` 1–35 (`pdb_seq_num` -34 to 0). TERT Asp712 has `label_seq_id = 747`, Asp868 has `712 + 35 = ?` (Asp868 is `label_seq_id = 903`). Using `auth_seq_num` matches canonical UniProt O14746. |
| 3 | TPP1 N-Terminal Offset Between 7TRE and 2I46 | Comparing `label_seq_id` of TPP1 between 7TRE and 2I46 | **3-residue discrepancy**: Glu168 is `label_seq_id = 82` in 7TRE (due to Gly87-Ala88-Gly89 cloning artifact), but `label_seq_id = 79` in 2I46. In both files, `auth_seq_num = 168`. Sequence alignment must use canonical UniProt sequence or `auth_seq_num`. |
| 4 | Chain ID Transposition in 7TRE | Requesting TERT coordinates using CIF chain `A` | **Yields RNA instead of Protein**: In 7TRE, `label_asym_id` 'A' is Telomerase RNA (451 nt), while `label_asym_id` 'B' is TERT. In author naming, TERT is 'A' and RNA is 'B'. Pipelines mixing the two conventions will pass RNA coordinates into protein alignment algorithms. |
| 5 | DNA Substrate Author Chain ID Inconsistency | Querying DNA substrate in 7TRD vs 7TRE | In 7TRE and 7TRF, the DNA chain has `auth_asym_id = 'D'`. In 7TRD, the DNA chain has `auth_asym_id = 'N'` (nucleic acid), while CIF `label_asym_id` is consistently 'C' (7TRE, 7TRD) and 'E' (7TRF). |
| 6 | Absence of TPP1 in 7TRD, 7TRF, and 7TRC | Attempting to calculate TPP1-TERT contacts or RMSD on 7TRD/7TRF/7TRC | **Entity not present**: TPP1 is present **only in 7TRE** (catalytic core complex) and **2I46** (isolated OB-fold). Cross-state comparison (R4) must align homologous TERT domains across 7TRD/7TRE/7TRF to measure TPP1-induced conformational changes. |
| 7 | Disordered TEN Domain Linker in TERT | Tracing continuous Cα backbone across TERT | In all three cryo-EM TERT models (7TRE, 7TRD, 7TRF), the flexible linker connecting the TEN domain to the TRBD domain (residues 179–321 in UniProt numbering) is **completely unmodeled (143 residues missing)**. RMSD calculations must segment TERT by domain. |
| 8 | Validation Report PDF Download 403 on RCSB | Script attempting automated `requests.get()` to `files.rcsb.org/..._full_validation.pdf` | **Returns HTTP 403 Forbidden**: RCSB Cloudflare firewall blocks automated PDF requests without specific session cookies. The identical official wwPDB validation PDF is downloadable with **HTTP 200** via PDBe (`https://www.ebi.ac.uk/pdbe/entry-files/download/{pid}_full_validation.pdf`). |
| 9 | 2I46 Reported Resolution Discrepancy | Ingesting 2I46 assuming 1.7 Å resolution from informal documentation | **Actual resolution is 2.70 Å**: In mmCIF `2I46.cif`, `_refine.ls_d_res_high` and `_reflns.d_resolution_high` are both 2.70 Å. Wilson B-factor is 43.36 Å². Data completeness is 97.2%. |
| 10 | DNA-RNA Template Pairing Geometry | Selecting full 18-nt DNA substrate for contact analysis | Only **7 nucleotides** (`12:DT` through `18:DG`) are modeled in electron density; nucleotides 1–11 are disordered in bulk solvent. |

---

## Detailed Structural Dossier by PDB Entry

### 1. PDB ID: `7TRE`
- **Title**: Human telomerase catalytic core with shelterin protein TPP1
- **Citation**: He, Q., Ho, F.N., Radovancevic, R., et al. *Structure of active human telomerase with telomere shelterin protein TPP1*. Nature (2022). [DOI: 10.1038/s41586-022-04661-3]
- **Experimental Method**: Single particle cryo-electron microscopy (`ELECTRON MICROSCOPY`)
- **Resolution**: **3.50 Å** (determined by gold-standard Fourier shell correlation FSC = 0.143 cutoff)
- **Primary Accession**: PDB `7TRE`, EMDB `EMD-26087`
- **Microscope & Imaging**:
  - Microscope: FEI Titan Krios
  - Accelerating Voltage: 300 kV
  - Symmetry: Point group C1 (asymmetric)
- **Validation Quality Summary**:
  - Cryo-EM Q-score: **0.392**
  - Atom Inclusion (all atoms): **79.4%** (0.794)
  - Atom Inclusion (backbone): **86.6%** (0.866)
  - Primary Map Recommended Contour Level: **0.018**
  - MolProbity Clashscore: **7.45**
  - Ramachandran Outliers: **0.00%** (Favored: 97.2%)
  - Rotamer Outliers: **0.11%**
  - RNA Suiteness: **0.36**
  - Geometry RMSZ: Bonds = 0.29, Angles = 0.43
- **Macromolecular Entity & Chain Architecture**:
  1. **Entity 1: Telomerase RNA (hTR)**
     - Molecule Type: polyribonucleotide (RNA)
     - Construct Length: 451 nt | Modeled Count: 217 nt
     - Canonical mmCIF Chain (`label_asym_id`): **`A`**
     - Author Chain (`auth_asym_id`): **`B`**
     - Modeled Segments: nt 33–149 (Pseudoknot & Template `5'-CUAACCCUAAC-3'`), nt 163–191, nt 251–321 (CR4/CR5 domain).
     - Disordered Segments: nt 1–32, 150–162, 192–250, 322–451.
  2. **Entity 2: Telomerase Reverse Transcriptase (TERT)**
     - Molecule Type: polypeptide(L) (Protein)
     - UniProt ID: **O14746** (Human TERT)
     - Construct Length: 1167 aa | Modeled Count: 934 aa | Disordered Count: 233 aa
     - Canonical mmCIF Chain (`label_asym_id`): **`B`**
     - Author Chain (`auth_asym_id`): **`A`**
     - Cloning Tag: N-terminal Twin-Strep tag (`GHMSAWSHPQFEKGGGSGGGSGGSAWSHPQFEKSA`, 35 aa, `seq_id` 1–35, `pdb_seq_num` -34 to 0, completely unmodeled).
     - First Modeled Residue: `Cys7` (`seq_id = 42`, `auth_seq_num = 7`).
     - Last Modeled Residue: `Asp1132` (`seq_id = 1167`, `auth_seq_num = 1132`).
     - Numbering Rule: $\text{auth\_seq\_num} = \text{label\_seq\_id} - 35$.
     - Major Unmodeled Loops: aa 1–6 (N-term), aa 105–116 (TEN domain loop), aa 179–321 (TEN–TRBD linker), aa 417–443 (TRBD loop), aa 641–650 (RT loop).
  3. **Entity 3: Telomeric Repeat Substrate (DNA)**
     - Molecule Type: polydeoxyribonucleotide (DNA)
     - Sequence: `5'-d(TTAGGGTTAGGGTTAGGG)-3'` (18 nt)
     - Canonical mmCIF Chain (`label_asym_id`): **`C`**
     - Author Chain (`auth_asym_id`): **`D`**
     - Modeled Bases: 7 nt (`DT12`, `DT13`, `DT14`, `DA15`, `DG16`, `DG17`, `DG18`).
  4. **Entity 4: Adrenocortical Dysplasia Protein Homolog (TPP1 / ACD)**
     - Molecule Type: polypeptide(L) (Protein)
     - UniProt ID: **Q96AP0** (Human TPP1, OB-fold domain)
     - Construct Length: 163 aa | Modeled Count: 152 aa | Disordered Count: 11 aa
     - Canonical mmCIF Chain (`label_asym_id`): **`D`**
     - Author Chain (`auth_asym_id`): **`C`**
     - Construct Range: Gly87 to Asn249 (includes cloning tripeptide `Gly87-Ala88-Gly89`).
     - Modeled Range: `Ser90` to `His241` (`auth_seq_num` 90 to 241; `label_seq_id` 4 to 155).
     - Unmodeled Residues: `Gly87`, `Ala88`, `Gly89` (N-terminal cloning tag, `seq_id` 1–3); `Glu242` to `Asn249` (C-terminal tail, `seq_id` 156–163).
     - Numbering Rule: $\text{auth\_seq\_num} = \text{label\_seq\_id} + 86$.
- **Authoritative File URLs & Checksums**:
  - mmCIF URL: `https://files.rcsb.org/download/7TRE.cif`
  - mmCIF Size: 1,518,438 bytes | **SHA-256**: `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a`
  - Full Validation PDF URL: `https://www.ebi.ac.uk/pdbe/entry-files/download/7tre_full_validation.pdf` (1,139,269 bytes, HTTP 200)
  - Validation XML URL: `https://files.rcsb.org/pub/pdb/validation_reports/tr/7tre/7tre_validation.xml.gz` (37,751 bytes, HTTP 200)
  - Primary EMDB Map URL: `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26087/map/emd_26087.map.gz` (62,724,891 bytes, HTTP 200)

---

### 2. PDB ID: `7TRD`
- **Title**: Human telomerase catalytic core structure at 3.3 Angstrom
- **Citation**: He, Q., Ho, F.N., Radovancevic, R., et al. Nature (2022).
- **Experimental Method**: Single particle cryo-electron microscopy (`ELECTRON MICROSCOPY`)
- **Resolution**: **3.30 Å** (FSC = 0.143 cutoff)
- **Primary Accession**: PDB `7TRD`, EMDB `EMD-26086`
- **Validation Quality Summary**:
  - Cryo-EM Q-score: **0.447**
  - Atom Inclusion (all atoms): **58.7%** (0.587)
  - Atom Inclusion (backbone): **66.8%** (0.668)
  - Primary Map Recommended Contour Level: **0.030**
  - MolProbity Clashscore: **7.59**
  - Ramachandran Outliers: **0.00%**
  - Rotamer Outliers: **0.12%**
  - RNA Suiteness: **0.38**
  - Geometry RMSZ: Bonds = 0.21, Angles = 0.39
- **Macromolecular Entities & Chains**:
  1. **Entity 1: Telomerase RNA (hTR)**
     - CIF Chain (`label_asym_id`): **`A`** | Author Chain (`auth_asym_id`): **`B`**
     - Construct: 451 nt | Modeled: 217 nt (segments 33–149, 163–191, 251–321)
  2. **Entity 2: Telomerase Reverse Transcriptase (TERT)**
     - UniProt: **O14746**
     - CIF Chain (`label_asym_id`): **`B`** | Author Chain (`auth_asym_id`): **`A`**
     - Construct: 1167 aa | Modeled: 934 aa | Disordered: 233 aa
     - Twin-Strep Tag: aa 1–35 (`pdb_seq_num` -34 to 0)
     - Numbering Rule: $\text{auth\_seq\_num} = \text{label\_seq\_id} - 35$.
  3. **Entity 3: Telomeric DNA Substrate**
     - Sequence: 18 nt
     - CIF Chain (`label_asym_id`): **`C`** | Author Chain (`auth_asym_id`): **`N`** (Note author chain ID is `N`!)
     - Modeled Bases: 7 nt (`DT12` to `DG18`)
- **Note on TPP1**: **Absent** (unbound catalytic core reference state).
- **Authoritative File URLs & Checksums**:
  - mmCIF URL: `https://files.rcsb.org/download/7TRD.cif`
  - mmCIF Size: 1,394,494 bytes | **SHA-256**: `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4`
  - Full Validation PDF URL: `https://www.ebi.ac.uk/pdbe/entry-files/download/7trd_full_validation.pdf` (1,182,982 bytes, HTTP 200)
  - Validation XML URL: `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trd/7trd_validation.xml.gz` (34,865 bytes, HTTP 200)
  - Primary EMDB Map URL: `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26086/map/emd_26086.map.gz` (62,791,352 bytes, HTTP 200)

---

### 3. PDB ID: `7TRF`
- **Title**: Human telomerase catalytic core RNP with H2A/H2B
- **Citation**: He, Q., Ho, F.N., Radovancevic, R., et al. Nature (2022).
- **Experimental Method**: Single particle cryo-electron microscopy (`ELECTRON MICROSCOPY`)
- **Resolution**: **3.70 Å** (FSC = 0.143 cutoff)
- **Primary Accession**: PDB `7TRF`, EMDB `EMD-26088`
- **Validation Quality Summary**:
  - Cryo-EM Q-score: **0.369**
  - Atom Inclusion (all atoms): **75.4%** (0.754)
  - Atom Inclusion (backbone): **85.0%** (0.850)
  - Primary Map Recommended Contour Level: **0.017**
  - MolProbity Clashscore: **7.52**
  - Ramachandran Outliers: **0.09%**
  - Rotamer Outliers: **0.10%**
  - RNA Suiteness: **0.38**
  - Geometry RMSZ: Bonds = 0.32, Angles = 0.48
- **Macromolecular Entities & Chains**:
  1. **Entity 1: Histone H2B (type 1-C/E/F/G/I)**
     - UniProt: **P62807** | Length: 126 aa
     - CIF Chain (`label_asym_id`): **`A`** | Author Chain (`auth_asym_id`): **`F`**
  2. **Entity 2: Histone H2A**
     - UniProt: **Q9BTM1** | Length: 129 aa
     - CIF Chain (`label_asym_id`): **`B`** | Author Chain (`auth_asym_id`): **`E`**
  3. **Entity 3: Telomerase RNA (hTR)**
     - CIF Chain (`label_asym_id`): **`C`** | Author Chain (`auth_asym_id`): **`B`**
     - Construct: 451 nt | Modeled: 217 nt (segments 33–149, 163–191, 251–321)
  4. **Entity 4: Telomerase Reverse Transcriptase (TERT)**
     - UniProt: **O14746**
     - CIF Chain (`label_asym_id`): **`D`** | Author Chain (`auth_asym_id`): **`A`**
     - Construct: 1167 aa | Modeled: 934 aa | Disordered: 233 aa
     - Numbering Rule: $\text{auth\_seq\_num} = \text{label\_seq\_id} - 35$.
  5. **Entity 5: Telomeric DNA Substrate**
     - CIF Chain (`label_asym_id`): **`E`** | Author Chain (`auth_asym_id`): **`D`**
     - Modeled Bases: 7 nt (`DT12` to `DG18`)
- **Note on TPP1**: **Absent** (complexed with nucleosomal histones H2A–H2B instead).
- **Authoritative File URLs & Checksums**:
  - mmCIF URL: `https://files.rcsb.org/download/7TRF.cif`
  - mmCIF Size: 1,542,948 bytes | **SHA-256**: `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd`
  - Full Validation PDF URL: `https://www.ebi.ac.uk/pdbe/entry-files/download/7trf_full_validation.pdf` (1,055,729 bytes, HTTP 200)
  - Validation XML URL: `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trf/7trf_validation.xml.gz` (39,193 bytes, HTTP 200)
  - Primary EMDB Map URL: `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26088/map/emd_26088.map.gz` (62,350,199 bytes, HTTP 200)

---

### 4. PDB ID: `7TRC`
- **Title**: Human telomerase H/ACA RNP at 3.3 Angstrom
- **Citation**: He, Q., Ho, F.N., Radovancevic, R., et al. Nature (2022).
- **Experimental Method**: Single particle cryo-electron microscopy (`ELECTRON MICROSCOPY`)
- **Resolution**: **3.30 Å** (FSC = 0.143 cutoff)
- **Primary Accession**: PDB `7TRC`, EMDB `EMD-26085`
- **Validation Quality Summary**:
  - Cryo-EM Q-score: **0.498**
  - Atom Inclusion (all atoms): **79.4%** (0.794)
  - Atom Inclusion (backbone): **83.9%** (0.839)
  - Primary Map Recommended Contour Level: **0.025**
  - MolProbity Clashscore: **8.18**
  - Ramachandran Outliers: **0.00%**
  - Rotamer Outliers: **0.21%**
  - RNA Suiteness: **0.51**
  - Geometry RMSZ: Bonds = 0.14, Angles = 0.37
- **Macromolecular Entities & Chains (H/ACA Lobe)**:
  1. **Entity 1: TCAB1 (Telomerase Cajal body protein 1 / WRAP53)**
     - UniProt: **Q9BUR4** | Length: 548 aa
     - CIF Chain: **`A`** | Author Chain: **`K`**
  2. **Entity 2: Dyskerin (DKC1)**
     - UniProt: **O60832** | Length: 514 aa
     - CIF Chains: **`B`**, **`J`** | Author Chains: **`C`**, **`G`**
  3. **Entity 3: NOP10 (H/ACA complex subunit 3)**
     - UniProt: **Q9NPE3** | Length: 64 aa
     - CIF Chains: **`C`**, **`H`** | Author Chains: **`F`**, **`J`**
  4. **Entity 4: Telomerase RNA (hTR H/ACA domain)**
     - CIF Chain: **`D`** | Author Chain: **`B`**
     - Construct: 451 nt | Modeled: 92 nt (segments 211–218, 362–392, 397–449)
  5. **Entity 5: GAR1 (H/ACA complex subunit 1)**
     - UniProt: **Q9NY12** | Length: 217 aa
     - CIF Chains: **`E`**, **`G`** | Author Chains: **`D`**, **`H`**
  6. **Entity 6: NHP2 (H/ACA complex subunit 2)**
     - UniProt: **Q9NX24** | Length: 153 aa
     - CIF Chains: **`F`**, **`I`** | Author Chains: **`E`**, **`I`**
- **Note on TERT and TPP1**: **Absent** (represents the distinct non-catalytic H/ACA ribonucleoprotein lobe).
- **Authoritative File URLs & Checksums**:
  - mmCIF URL: `https://files.rcsb.org/download/7TRC.cif`
  - mmCIF Size: 1,647,704 bytes | **SHA-256**: `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea`
  - Full Validation PDF URL: `https://www.ebi.ac.uk/pdbe/entry-files/download/7trc_full_validation.pdf` (1,279,593 bytes, HTTP 200)
  - Validation XML URL: `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trc/7trc_validation.xml.gz` (49,084 bytes, HTTP 200)
  - Primary EMDB Map URL: `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26085/map/emd_26085.map.gz` (122,175,885 bytes, HTTP 200)

---

### 5. PDB ID: `2I46`
- **Title**: Crystal structure of human TPP1
- **Citation**: Wang, F., Podell, E.R., Zaug, A.J., Batey, R.T., Cech, T.R., Lei, M. *The POT1-TPP1 telomere complex is a telomerase processivity factor*. Nature 445:506-510 (2007). [DOI: 10.1038/nature05454]
- **Experimental Method**: X-ray diffraction (`X-RAY DIFFRACTION`)
- **Resolution**: **2.70 Å** (Correct authoritative value; Space Group: $P2_12_12_1$)
- **Unit Cell**: $a = 48.24\text{ \AA},\; b = 86.30\text{ \AA},\; c = 91.43\text{ \AA},\; \alpha = \beta = \gamma = 90.00^\circ$
- **Crystallographic Refinement Statistics**:
  - High Resolution Limit: **2.70 Å** (Low: 50.0 Å)
  - $R_{\text{work}}$: **0.2135** (DCC = 0.2101)
  - $R_{\text{free}}$: **0.2509** (DCC = 0.2474)
  - Completeness: **97.2%**
  - Wilson B-factor Estimate: **43.36 Å²**
  - Data Anisotropy: **0.569**
  - MolProbity Clashscore: **12.97**
  - Ramachandran Outliers: **0.36%** (Favored: 92.5%)
  - Rotamer Outliers: **4.12%**
  - RSRZ Outliers: **13.29%**
- **Macromolecular Entity & Chains**:
  - **Entity 1: Adrenocortical Dysplasia Protein Homolog (TPP1 / ACD)**
    - UniProt: **Q96AP0** | Construct Length: 161 aa (`Ser90` to `Ala250`)
    - Asymmetric Unit Contents: Homodimer composed of **Chain A** and **Chain B**
    - Canonical CIF Chains: **`A`**, **`B`** | Author Chains: **`A`**, **`B`**
    - Numbering Rule: $\text{auth\_seq\_num} = \text{label\_seq\_id} + 89$.
  - **Chain A Conformational Status**:
    - Modeled Residues: `Ser90` to `Glu243` (154 aa)
    - **TEL-Patch DE-Loop (`Glu168–171`)**: **FULLY MODELED**
      - `Glu168`: `label_seq_id = 79`, `auth_seq_num = 168`
      - `Glu169`: `label_seq_id = 80`, `auth_seq_num = 169`
      - `Lys170`: `label_seq_id = 81`, `auth_seq_num = 170`
      - `Glu171`: `label_seq_id = 82`, `auth_seq_num = 171`
      - `Leu183`: `label_seq_id = 94`, `auth_seq_num = 183`
      - `Leu212`: `label_seq_id = 123`, `auth_seq_num = 212`
  - **Chain B Conformational Status**:
    - Modeled Residues: `Ser90` to `Glu243` EXCEPT loop 168–171
    - **TEL-Patch DE-Loop (`Glu168–171`)**: **DISORDERED / UNOBSERVED** (`auth_seq_num = '?'`)
    - `Leu183`: `label_seq_id = 94`, `auth_seq_num = 183` (Modeled)
    - `Leu212`: `label_seq_id = 123`, `auth_seq_num = 212` (Modeled)
- **Authoritative File URLs & Checksums**:
  - mmCIF URL: `https://files.rcsb.org/download/2I46.cif`
  - mmCIF Size: 273,104 bytes | **SHA-256**: `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779`
  - Full Validation PDF URL: `https://www.ebi.ac.uk/pdbe/entry-files/download/2i46_full_validation.pdf` (524,158 bytes, HTTP 200)
  - Validation XML URL: `https://files.rcsb.org/pub/pdb/validation_reports/i4/2i46/2i46_validation.xml.gz` (16,330 bytes, HTTP 200)

---

## TEL-Patch & TERT Receptor Cradle Contact Architecture

Analysis of the 7TRE bound interface between TPP1 (CIF chain D / Author chain C) and TERT (CIF chain B / Author chain A) reveals a bipartite contact cradle. The acidic DE-loop and hydrophobic anchor of TPP1 nestle into a cleft formed between the **TEN domain** and the **TRAP/IFD domain** of TERT:

### Direct Contacts of the TEL-Patch Residues:
1. **Glu168 (TPP1)**:
   - Primary Salt Bridge / H-Bond: **TERT Arg774** (min heavy-atom distance **2.32 Å**)
   - Hydrophobic Capping Contact: **TERT Ala46** (min distance **3.96 Å**)
2. **Glu169 (TPP1)**:
   - Polar / Packing Contacts: **TERT Ala46** (**3.35 Å**), **TERT Ala136** (**4.21 Å**), **TERT Gly135** (**4.25 Å**)
3. **Lys170 (TPP1)**:
   - Solvent-exposed basic sidechain projecting outward from the cradle, allowing steric tolerance for electrostatic steering.
4. **Glu171 (TPP1)**:
   - Specific Polar / Backbone Contacts: **TERT Ala45** (**3.45 Å**), **TERT Pro44** (**3.51 Å**), **TERT Asp43** (**3.97 Å**)
5. **Leu183 (TPP1)**:
   - Core Hydrophobic Pocket: **TERT Leu769** (**3.67 Å**)
6. **Leu212 (TPP1)**:
   - Hydrophobic Anchor: **TERT Pro771** (**3.52 Å**), **TERT Leu769** (**3.75 Å**)

### Peripheral Reinforcing Contacts:
- **Arg92 (TPP1)** forms polar contacts with **TERT Asn799** (3.24 Å) and **Ser121** (3.85 Å).
- **Leu93 (TPP1)** packs against **TERT Leu766** (3.64 Å) and **Leu798** (3.98 Å).
- **Leu211 (TPP1)** packs against **TERT Tyr772** (3.36 Å).
- **Pro213 (TPP1)** packs against **TERT Asp129** (3.55 Å), **Tyr772** (3.56 Å), **Ala130** (3.62 Å), and **Pro771** (3.63 Å).
- **Glu215 (TPP1)** forms a salt bridge with **TERT Lys78** (3.58 Å) and contacts **Arg132** (4.07 Å).

This demonstrates that the **TEN domain** (residues 43–46, 78, 121–136) coordinates the acidic residues `Glu169`, `Glu171`, and `Glu215`, while the **TRAP/IFD domain** (residues 766–775, 798–799) accommodates the salt bridge with `Glu168` and the hydrophobic anchors `Leu183`, `Leu211`, and `Leu212`.

---

## Authoritative Master Data Manifest Specification

The following tabular specification defines the schema and exact values required for generating `data_manifest.csv` to satisfy Requirement R1:

| PDB_ID | Method | Resolution_A | EMDB_ID | CIF_Size_Bytes | CIF_SHA256 | Val_PDF_URL | Val_XML_URL | Map_URL | TERT_CIF_Chain | TERT_Auth_Chain | TPP1_CIF_Chain | TPP1_Auth_Chain | RNA_CIF_Chain | RNA_Auth_Chain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `7TRE` | ELECTRON MICROSCOPY | 3.50 | EMD-26087 | 1518438 | `2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a` | `https://www.ebi.ac.uk/pdbe/entry-files/download/7tre_full_validation.pdf` | `https://files.rcsb.org/pub/pdb/validation_reports/tr/7tre/7tre_validation.xml.gz` | `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26087/map/emd_26087.map.gz` | B | A | D | C | A | B |
| `7TRD` | ELECTRON MICROSCOPY | 3.30 | EMD-26086 | 1394494 | `faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4` | `https://www.ebi.ac.uk/pdbe/entry-files/download/7trd_full_validation.pdf` | `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trd/7trd_validation.xml.gz` | `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26086/map/emd_26086.map.gz` | B | A | N/A | N/A | A | B |
| `7TRF` | ELECTRON MICROSCOPY | 3.70 | EMD-26088 | 1542948 | `e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd` | `https://www.ebi.ac.uk/pdbe/entry-files/download/7trf_full_validation.pdf` | `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trf/7trf_validation.xml.gz` | `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26088/map/emd_26088.map.gz` | D | A | N/A | N/A | C | B |
| `7TRC` | ELECTRON MICROSCOPY | 3.30 | EMD-26085 | 1647704 | `99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea` | `https://www.ebi.ac.uk/pdbe/entry-files/download/7trc_full_validation.pdf` | `https://files.rcsb.org/pub/pdb/validation_reports/tr/7trc/7trc_validation.xml.gz` | `https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-26085/map/emd_26085.map.gz` | N/A | N/A | N/A | N/A | D | B |
| `2I46` | X-RAY DIFFRACTION | 2.70 | N/A | 273104 | `4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779` | `https://www.ebi.ac.uk/pdbe/entry-files/download/2i46_full_validation.pdf` | `https://files.rcsb.org/pub/pdb/validation_reports/i4/2i46/2i46_validation.xml.gz` | N/A | N/A | N/A | A, B | A, B | N/A | N/A |

---

## Technical Guidance for Downstream Engineering Agents

1. **Biopython / Gemmi Parser Conventions**:
   - When loading mmCIF files in Gemmi or Biopython (`Bio.PDB.MMCIFParser`), Gemmi defaults to `subchain_id` or `label_asym_id`, whereas Biopython often indexes chains by `auth_asym_id`.
   - **Recommendation**: Always explicitly extract both `_atom_site.label_asym_id` and `_atom_site.auth_asym_id`. To align TERT between 7TRE, 7TRD, and 7TRF, select:
     - In 7TRE: `label_asym_id == 'B'` (or `auth_asym_id == 'A'`)
     - In 7TRD: `label_asym_id == 'B'` (or `auth_asym_id == 'A'`)
     - In 7TRF: `label_asym_id == 'D'` (or `auth_asym_id == 'A'`)
2. **Residue Alignment Logic (Requirement R2)**:
   - Do **NOT** use `label_seq_id` for direct residue matching between TERT and TPP1 constructs.
   - Use `auth_seq_num` or perform an automated pairwise sequence alignment against the reference UniProt FASTA sequence (`O14746` for TERT, `Q96AP0` for TPP1).
3. **TPP1 OB-Fold Comparison (7TRE vs 2I46)**:
   - Compare 7TRE TPP1 (`label_asym_id == 'D'`) specifically with **`2I46` Chain A** (`label_asym_id == 'A'`).
   - If aligning to 2I46 Chain B, mask out residues 168–171 to prevent alignment errors caused by missing coordinates.
