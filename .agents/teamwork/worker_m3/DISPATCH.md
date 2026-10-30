# DISPATCH — worker_m3

## Mission
Implement Milestone 3: Quantitative Interface Geometry & Interaction Profiling (R3) in `src/telotrade/interface_geometry.py` and unit tests in `tests/test_interface_geometry.py`.

## Working Directory
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m3`

## Mandatory Reference Documents
Read these files thoroughly before starting:
1. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md` (authoritative user requirements)
2. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\orchestrator\PROJECT.md` (architecture, milestones, interface contracts)
3. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md` (Section 2: Contact networks, BSA formulas, polar networks, TEL-patch mapping, exact benchmarks; Section 4: CSV table schemas)
4. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\report.md` (Chain IDs, offsets)
5. `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\src\telotrade\sequence_mapping.py` (Residue mapping, offsets, chain conversions)

## MANDATORY INTEGRITY WARNING
DO NOT CHEAT. All implementations must be genuine. DO NOT hardcode test results, create dummy/facade implementations, or circumvent the intended task. A teamwork_preview_auditor will independently verify your work. Integrity violations WILL be detected and your work WILL be rejected.

## Detailed Scope & Technical Requirements
1. **Module Implementation in `src/telotrade/interface_geometry.py`**:
   - Contact network analysis:
     - `compute_interface_contacts(cif_path, chain_tert='A', chain_tpp1='C', cutoffs=[3.5, 4.0, 4.5], use_auth=True)`:
       - Uses Gemmi or NumPy grid hashing to find all heavy-atom contacts between TERT and TPP1 in 7TRE at 3.5 Å, 4.0 Å, and 4.5 Å.
       - Supports both author chain IDs (`A` for TERT, `C` for TPP1) and label chain IDs (`B` for TERT, `D` for TPP1).
       - Classifies interaction types:
         - `Salt Bridge`: Anionic oxygen (Asp OD1/OD2, Glu OE1/OE2) to basic nitrogen (Arg NH1/NH2/NE, Lys NZ, His ND1/NE2) <= 4.0 Å.
         - `Polar/H-Bond`: Non-salt-bridge polar atom contacts (N, O, S) <= 3.5 Å (or <= 4.0 Å if polar).
         - `Hydrophobic`: Aliphatic/aromatic carbon-carbon contacts (Ala, Val, Leu, Ile, Pro, Phe, Trp, Tyr, Met) <= 4.5 Å.
         - `Van der Waals`: Other heavy atom pairs <= cutoff.
       - Annotates TERT structural domain: `TEN` (residues 7-178), `TRAP_IFD` (residues 766-800), `RT` (residues 540-765, 801-930), etc.
       - Flags `is_tel_patch` on TPP1: `Glu168, Glu169, Lys170, Glu171, Leu183, Leu212` (and optionally `Arg180, Glu215`).
   - Buried Surface Area (BSA) calculation:
     - `compute_bsa(cif_path, chain_tert='A', chain_tpp1='C', probe_radius=1.4, n_sphere_points=960, use_auth=True)`:
       - Computes SASA of isolated partner A ($\text{SASA}_A$), isolated partner B ($\text{SASA}_B$), and binary complex ($\text{SASA}_{AB}$).
       - Uses FreeSASA or Biopython `ShrakeRupley` ($r_p = 1.4\ \text{Å}$, $N=960$).
       - Formula: $\text{BSA} = (\text{SASA}_A + \text{SASA}_B - \text{SASA}_{AB}) / 2$.
       - Total $\Delta\text{SASA} = \text{SASA}_A + \text{SASA}_B - \text{SASA}_{AB}$.
       - Expected 7TRE empirical benchmark: $\text{BSA} \approx 856.82\ \text{Å}^2$ ($\Delta\text{SASA} \approx 1,713.65\ \text{Å}^2$).
   - Polar interaction profiling:
     - `compute_polar_interactions(cif_path, chain_tert='A', chain_tpp1='C', use_auth=True)`:
       - Identifies all salt bridges and hydrogen bonds.
       - Specifically validates primary salt bridge: TERT Arg774 (NH1/NE) to TPP1 Glu168 (OE1/OE2) at ~2.32 Å.
       - Specifically validates secondary clamp: TERT Lys78 (NZ) to TPP1 Glu215 (OE2) at ~3.58 Å.
   - TEL-patch profiling:
     - `map_tel_patch_contacts(contacts_df)`:
       - Generates dedicated table for TEL-patch residues (`Glu168–171`, `Leu183`, `Leu212`, `Arg180`, `Glu215`), documenting partners, interaction types, and distances.
       - Specifically validates hydrophobic anchors: TERT Leu769 to TPP1 Leu183 (3.67 Å) and TERT Pro771 to TPP1 Leu212 (3.52 Å).
   - Cutoff robustness analysis:
     - `compute_contact_robustness(cif_path, chain_tert='A', chain_tpp1='C', cutoffs=[3.5, 4.0, 4.5])`:
       - Tracks persistence across cutoffs: `Core Robust` (present at 3.5, 4.0, 4.5 Å), `Intermediate` (present at 4.0, 4.5 Å), `Threshold Sensitive` (present only at 4.5 Å).
   - Table Generation:
     - Functions to generate and save CSV files to `results/tables/`:
       1. `interface_contacts_7tre_3.5A.csv` (and `interface_contacts_3.5A.csv`)
       2. `interface_contacts_7tre_4.0A.csv` (and `interface_contacts_4.0A.csv`)
       3. `interface_contacts_7tre_4.5A.csv` (and `interface_contacts_4.5A.csv`)
       4. `interface_contacts_robustness.csv`
       5. `buried_surface_area.csv`
       6. `polar_interactions_7tre.csv` (and `polar_interactions.csv`)
       7. `tel_patch_contacts.csv`
     - Column schemas must match `survey_methodology/report.md` Section 4.
2. **Expose Core APIs in `src/telotrade/__init__.py`**:
   - Export `compute_interface_contacts`, `compute_bsa`, `compute_polar_interactions`, `map_tel_patch_contacts`, `compute_contact_robustness`.
3. **Unit Tests in `tests/test_interface_geometry.py`**:
   - Comprehensive tests validating contact counts across cutoffs, benchmark distances (Arg774-Glu168, Lys78-Glu215, Leu769-Leu183, Pro771-Leu212), BSA calculation ($\approx 856.82\ \text{Å}^2 \pm 50\ \text{Å}^2$), robustness scoring, CSV export, and handling both author and label chain IDs.
4. **Verification**:
   - Run `pytest tests/test_interface_geometry.py -v`.
   - Run full test suite `pytest tests/ -v`.
   - Run `ruff check src/ tests/`.
5. **Handoff**:
   - Write standard 5-component `handoff.md` (Observation, Logic Chain, Caveats, Conclusion, Verification Method) in your working directory and notify orchestrator.


## 2026-10-16T16:46:21Z
Received dispatch from parent eb46fecb-7562-4ce7-aa98-998216cee2c4:
You are Milestone 3 Implementation Worker.
Working directory: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\worker_m3
Scope:
1. Implement src/telotrade/interface_geometry.py with contact networks (3.5 Å, 4.0 Å, 4.5 Å), Buried Surface Area (BSA) via Shrake-Rupley / FreeSASA, polar interaction profiling (salt bridges including Arg774-Glu168 and Lys78-Glu215, H-bonds), TEL-patch mapping (hydrophobic anchors Leu769-Leu183, Pro771-Leu212), cutoff robustness scoring, and table generators for results/tables/.
2. Expose core APIs in src/telotrade/__init__.py.
3. Write unit tests in tests/test_interface_geometry.py.
4. Run tests: pytest tests/test_interface_geometry.py tests/ -v and ruff check. Verify 100% pass.
5. Write your 5-component handoff.md and send completion message to parent.
