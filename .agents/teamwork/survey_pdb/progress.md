# Progress — survey_pdb

Last visited: 2026-10-16T14:32:00Z
Status: Completed

## Completed
- Received dispatch instructions and appended to DISPATCH.md.
- Initialized BRIEFING.md, progress.md, and working directories.
- Queried RCSB REST, Entry, and Polymer Instance APIs for 7TRE, 7TRD, 7TRF, 7TRC, and 2I46.
- Downloaded and locally verified all 5 mmCIF coordinate files into `cifs/`.
- Computed exact SHA-256 hashes and byte lengths for all mmCIF files.
- Fully mapped macromolecular entities, CIF `label_asym_id`, and author `auth_asym_id` across all entries.
- Uncovered author vs label sequence numbering offsets and cloning tags:
  - TERT: N-terminal 35-residue Twin-Strep tag (`pdb_seq_num` -34 to 0; `auth_seq_num = label_seq_id - 35` matches UniProt O14746).
  - TPP1 in 7TRE: 3-residue cloning tag (`Gly87-Ala88-Gly89`; `auth_seq_num = label_seq_id + 86` matches UniProt Q96AP0).
  - TPP1 in 2I46: `auth_seq_num = label_seq_id + 89`. Discovered that Chain A has the TEL-patch DE-loop (168-171) modeled, while in Chain B the DE-loop is disordered/unmodeled.
- Mapped TEL patch residues (Glu168-171, Leu183, Leu212) and their atomic contact cradle in TERT (TEN domain: Asp43, Pro44, Ala45, Ala46; TRAP/IFD domain: Leu766, Leu769, Pro771, Tyr772, Arg774, Asn799).
- Tested and verified HTTP 200 download endpoints for mmCIFs, PDBe validation PDFs, RCSB validation XMLs, and EMDB map volumes.
- Extracted official wwPDB validation report metrics (Q-score, Clashscore, Ramachandran outliers, Rotamer outliers, RNA suiteness, R_work/R_free).
- Generated publication-grade specification report `report.md` with complete Features Discovered table (20 features), Edge Cases table (10 edge cases), and master data manifest specification.
- Generated self-contained 5-Component handoff report `handoff.md`.
- Updated BRIEFING.md.
