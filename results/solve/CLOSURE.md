# Closure ledger — what of "all of it" is actually solved

    python scripts/solve/wp46_closure.py
      -> results/solve/CLOSURE.md  (this file)
      -> results/solve/wp46_closure.json / .csv

Generated from the artefacts, not from memory. Every status below is checked
against the JSON of the work package that produced the number.

## The answer

**No. 7 of 15 sub-problems are solved, and all 7 are solved as CALCULATIONS. Zero are solved as evidence, because nothing has been made or measured. 4 are not solved, 2 are partial, and 2 cannot be answered by any computation here. The single largest open item is P6 - the presenting scaffold - which blocks the only route that reached the full triad.**

| status | count |
|---|---|
| SOLVED | **7** |
| PARTIAL | **2** |
| NOT SOLVED | **4** |
| UNMEASURABLE | **2** |

## Every sub-problem

| id | area | problem | status | kind |
|---|---|---|---|---|
| P1 | Interface | Is the TPP1 TEL patch one binding surface? | **SOLVED** | geometric |
| P2 | Interface | Can a drug-like small molecule span both halves? | **SOLVED** | geometric |
| P3 | Interface | What is the shortest path across the gulf? | **SOLVED** | geometric |
| P4 | Triad | Can an alpha-helix present the E168/E169/E171 triad? | **SOLVED** | geometric |
| P5 | Triad | Can the triad be presented at all, in the native geometry? | **SOLVED** | geometric |
| P6 | Triad | Does a 19-mer peptide HOLD that geometry? | **NOT SOLVED** | design-level |
| P7 | Modularity | Can each module solve only its own site, fully? | **PARTIAL** | geometric |
| P8 | Modularity | Do merged modules stay stable together? | **SOLVED** | design-level |
| P9 | Modularity | Can one groove be blocked and repaired at once? | **SOLVED** | logical |
| P10 | Screening | Do the designs pass computed filters? | **PARTIAL** | computed |
| P11 | Screening | What is the experimental success rate? | **UNMEASURABLE** | wet-lab |
| P12 | Screening | Do any of these designs bind telomerase? | **UNMEASURABLE** | wet-lab |
| P13 | Mechanism | What does TPP1 change mechanically (M5)? | **NOT SOLVED** | wet-lab |
| P14 | Mechanism | Is n measured? | **NOT SOLVED** | wet-lab |
| P15 | Drug | Is there a drug? | **NOT SOLVED** | design-level |

## What each status rests on

**P1 — Is the TPP1 TEL patch one binding surface?** — *SOLVED*

No - bipartite. E215 grips K78/R132, E168 grips R774, 20.26 A apart.

**P2 — Can a drug-like small molecule span both halves?** — *SOLVED*

No - and it is a class exclusion. 7.8 phenyls needed: MW ~683, logP ~11.

**P3 — What is the shortest path across the gulf?** — *SOLVED*

25.54 A with a 119 deg bend and tortuosity 1.26; a straight line is blocked.

**P4 — Can an alpha-helix present the E168/E169/E171 triad?** — *SOLVED*

No - 2 of 3 at best. The segment is a kinked turn (psi +24.4 deg, CA168-CA171 7.68 A), not a helix.

**P5 — Can the triad be presented at all, in the native geometry?** — *SOLVED*

Yes - graft the two native segments. 3 of 3 to 0.02 A: DWEEKEF-GGGGG-FSLLPTE, E215 at 4.10 A from K78.

**P6 — Does a 19-mer peptide HOLD that geometry?** — *NOT SOLVED*

No covalent route exists: closest disulfide Cb-Cb 6.21 A (needs 3.5-5.0), head-to-tail 14.28 A CA-CA (needs ~6), and freezing the linker costs 12-20 kcal/mol against 5.4 kcal/mol of contact energy.

Blocker: a rigid presenting scaffold, which needs PDB scaffold screening + loop modelling (no Rosetta, no loop engine, no offline scaffold set here).

**P7 — Can each module solve only its own site, fully?** — *PARTIAL*

4 of 5 do. One site (S3) was declared unsolvable on a helix model that was wrong - corrected by P4/P5 - so the module result is retained but its premise is superseded.

Blocker: the MD3 module needs rebuilding on the graft (P5), not on a helix.

**P8 — Do merged modules stay stable together?** — *SOLVED*

4 of 5 merges STABLE, 1 COLLAPSE-RISK, 1 REFUSED. The linker alone decides it: same two modules, -6.04 kcal/mol attraction when flexible.

**P9 — Can one groove be blocked and repaired at once?** — *SOLVED*

No - the merge engine refuses it by construction, with a stated reason.

**P10 — Do the designs pass computed filters?** — *PARTIAL*

6/10 clear every gate (60 %), or 6/15 = 40 % with the controls counted. The battery does NOT fully validate on its controls: no single alert configuration accepts both approved drugs and rejects all three bad ones.

Blocker: no alert configuration is both sensitive and specific on the control set.

**P11 — What is the experimental success rate?** — *UNMEASURABLE*

Unknown. N = 0 designs have been made, so the rate has no denominator.

Blocker: synthesis and assay; no amount of computation supplies this.

**P12 — Do any of these designs bind telomerase?** — *UNMEASURABLE*

Unknown - and it is reported as unknown rather than scored.

Blocker: no docking engine, no OpenMM, no FEP; WP39 already refuted the in-repo rigid-body scorer, so a score would be a fabricated number.

**P13 — What does TPP1 change mechanically (M5)?** — *NOT SOLVED*

M5 carries P = 0.35. Two branches separated by one DNA-trap experiment.

Blocker: the P2 two-chase experiment has not been run.

**P14 — Is n measured?** — *NOT SOLVED*

No. The two-chase inversion is exact and unrun.

Blocker: wet lab.

**P15 — Is there a drug?** — *NOT SOLVED*

No, and the small-molecule route for the bipartite patch is excluded (P2). The protein route is unmeasured (P11, P12).

Blocker: P6 for the graft; P11/P12 for everything else.

## The distinction this ledger exists to protect

`SOLVED` in every row above means solved **as a calculation**. Six of fifteen
sub-problems are closed, and all six are closed against the deposit or against
arithmetic. **Zero are closed as evidence.** Nothing in this repository has been
synthesised, expressed or assayed, so no design here has been shown to bind
anything. A chain of improvements is not a chain of closures, and the six open
items are the honest measure of the distance left.

## What would move it

- P6: a scaffold search (needs PDB screening + loop modelling, not available here)
- P10: an alert configuration that is sensitive AND specific on the controls
- P11/P12: synthesis, then one binding and one activity measurement
- P13/P14: the P2 two-chase experiment
