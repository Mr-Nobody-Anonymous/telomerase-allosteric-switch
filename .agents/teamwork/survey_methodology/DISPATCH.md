# Dispatch for Survey Agent 3 (Methodology & Pipeline Explorer)

## Task Description
Read ORIGINAL_REQUEST.md at:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md`

Your working directory is:
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology`

## Mission
Investigate and define the structural bioinformatics algorithms, mathematical formulas, and pipeline requirements for R1–R5:
1. R2 Sequence-aware residue correspondence mapping: Gemmi vs Biopython mmCIF parsing, sequence alignment (pairwise / Smith-Waterman or Needleman-Wunsch or coordinate sequence extraction) before geometric alignment. Handling author vs label sequence numbering.
2. R3 Interface geometry:
   - Contact networks at 3.5 Å, 4.0 Å, 4.5 Å (atom selection: all heavy atoms vs Cα).
   - BSA: formula $(\mathrm{SASA}_A + \mathrm{SASA}_B - \mathrm{SASA}_{AB}) / 2$. Implementation via Shrake-Rupley (Biopython `ShrakeRupley` or custom/gemmi/freesasa).
   - Polar networks (H-bonds: donor-acceptor distances and angles; salt bridges: acidic vs basic side chains <= 4.0 Å).
   - Hydrophobic packing and TEL-patch residues mapping.
3. R4 Cross-state comparative structural analysis & robustness:
   - Superposition / Cα RMSD between 7TRE, 7TRD, 7TRF (global, catalytic core, TEN domain, TRAP domain).
   - Contact persistence / sensitivity across cutoffs (robustness table).
   - TPP1 conformational changes between 7TRE bound vs 2I46 unbound.
4. R5 Synthesis report and table/figure schemas:
   - Structure of `report.md`, `results/tables/`, `results/figures/`.
   - Clear distinction between Axis A (molecular mechanism) and Axis B (therapeutic selectivity: cancer vs healthy aging / stem cells).

Write your findings to `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md` and complete with `handoff.md`.

## 2026-10-16T14:16:27Z
Sender: eb46fecb-7562-4ce7-aa98-998216cee2c4
Message:
You are the Methodology & Pipeline Explorer.
Your working directory is: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology
Read ORIGINAL_REQUEST.md at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\ORIGINAL_REQUEST.md
Also read your dispatch file at: C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\DISPATCH.md

Your mission is to investigate:
1. Algorithmic details for sequence-aware residue correspondence mapping between homologous chains (TERT chains across 7TRE, 7TRD, 7TRF; TPP1 between 7TRE and 2I46) using Gemmi/Biopython.
2. Algorithmic details for contact network computation (3.5 Å, 4.0 Å, 4.5 Å cutoffs; heavy atom vs Cα), Buried Surface Area (BSA formula), polar interactions (H-bonds, salt bridges <= 4.0 Å), and TEL-patch mapping.
3. Cross-state comparative structural analysis: superposition algorithms, domain definitions (TEN: 1-199/330, TRAP, Catalytic core), RMSD metrics, and cutoff robustness testing.
4. Schema for all CSV tables in results/tables/ and outline of report.md (Axis A vs Axis B).

Save your detailed report to C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_methodology\report.md and write handoff.md in your directory.
Maintain progress in progress.md in your directory. Send a message to parent when completed.
