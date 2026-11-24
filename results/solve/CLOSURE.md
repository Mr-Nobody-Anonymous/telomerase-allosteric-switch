# Closure ledger — what of "all of it" is actually solved

    python scripts/solve/wp46_closure.py
      -> results/solve/CLOSURE.md  (this file)
      -> results/solve/wp46_closure.json / .csv

Generated from the artefacts, not from memory. Every status below is checked
against the JSON of the work package that produced the number.

## The answer

**No. 8 of 16 sub-problems are solved, and all 8 are solved as CALCULATIONS. Zero are solved as evidence, because nothing has been made or measured. 3 are not solved, 3 are partial, and 2 cannot be answered by any computation here. The largest open item is still P6: WP47 found a real covalent constraint (W167-L211, 6.34 A) where WP45 had concluded none existed, but one staple buys 2.8 kcal/mol against a 12-20 kcal/mol ordering cost. WP48 then improved the construct (PAPAP linker, -2.84 kcal/mol) and still has 6.36 kcal/mol unpaid. A scaffold is still required.**

| status | count |
|---|---|
| SOLVED | **8** |
| PARTIAL | **3** |
| NOT SOLVED | **3** |
| UNMEASURABLE | **2** |

## Every sub-problem

| id | area | problem | status | kind |
|---|---|---|---|---|
| P1 | Interface | Is the TPP1 TEL patch one binding surface? | **SOLVED** | geometric |
| P2 | Interface | Can a drug-like small molecule span both halves? | **SOLVED** | geometric |
| P3 | Interface | What is the shortest path across the gulf? | **SOLVED** | geometric |
| P4 | Triad | Can an alpha-helix present the E168/E169/E171 triad? | **SOLVED** | geometric |
| P5 | Triad | Can the triad be presented at all, in the native geometry? | **SOLVED** | geometric |
| P6 | Triad | Does a peptide HOLD that geometry? | **PARTIAL** | design-level |
| P7 | Modularity | Can each module solve only its own site, fully? | **PARTIAL** | geometric |
| P8 | Modularity | Do merged modules stay stable together? | **SOLVED** | design-level |
| P9 | Modularity | Can one groove be blocked and repaired at once? | **SOLVED** | logical |
| P10 | Screening | Is a structural alert a usable filter at all? | **SOLVED** | computed |
| P10b | Screening | Does the full gate stack accept approved drugs? | **PARTIAL** | computed |
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

**P6 — Does a peptide HOLD that geometry?** — *PARTIAL*

A covalent constraint EXISTS, and WP45 missed it by screening disulfides only. Over six chemistries spanning 3.5-12.5 A, W167-L211 passes reach, direction, clash and payload tests at 6.34 A with cos 0.98/0.09 - and the native indole/Leu contact there is already 3.50 A. But the staple buys only 2.8 kcal/mol of ordering against a 12-20 kcal/mol cost: residual 9.2 kcal/mol unpaid. Constraint found, energy still not balanced.

Blocker: the remaining 6.36 kcal/mol even after a structured PAPAP linker and a second staple attempt - neither closes it, so a presenting scaffold is still required (PDB screening + loop modelling, unavailable here).

**P7 — Can each module solve only its own site, fully?** — *PARTIAL*

4 of 5 do. One site (S3) was declared unsolvable on a helix model that was wrong - corrected by P4/P5 - so the module result is retained but its premise is superseded.

Blocker: the MD3 module needs rebuilding on the graft (P5), not on a helix.

**P8 — Do merged modules stay stable together?** — *SOLVED*

4 of 5 merges STABLE, 1 COLLAPSE-RISK, 1 REFUSED. The linker alone decides it: same two modules, -6.04 kcal/mol attraction when flexible.

**P9 — Can one groove be blocked and repaired at once?** — *SOLVED*

No - the merge engine refuses it by construction, with a stated reason.

**P10 — Is a structural alert a usable filter at all?** — *SOLVED*

No as a veto - proven unfixable, alert count caps at 0.70 balanced accuracy. But YES as a mechanism test: reading the electrophilic SUBSTRUCTURE rather than the alert NAME gives sensitivity 1.0 and specificity 1.0 (11/11) on a panel used by neither earlier WP. The failure was the vocabulary, and replacing it fixed it.

**P10b — Does the full gate stack accept approved drugs?** — *PARTIAL*

No. 10/11 end-to-end on the third panel; the one miss is Lipinski rejecting atorvastatin (MW 559, logP > 5). That is a pre-existing gate defect inherited from WP43, not a reactivity failure - it is attributed, not blamed on the layer that was fixed.

Blocker: Lipinski false rejection of large approved drugs; needs a rule that does not reject a marketed molecule at MW 559.

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
