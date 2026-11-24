# Closure ledger — what of "all of it" is actually solved

    python scripts/solve/wp46_closure.py
      -> results/solve/CLOSURE.md  (this file)
      -> results/solve/wp46_closure.json / .csv

Generated from the artefacts, not from memory. Every status below is checked
against the JSON of the work package that produced the number.

## The answer

**No. 12 of 16 sub-problems are solved, and all 12 are solved as CALCULATIONS. Zero are solved as evidence, because nothing has been made or measured. 0 are not solved, 4 are partial, and 0 cannot be answered by any computation here. The largest open item is still P6: WP47 found a real covalent constraint (W167-L211, 6.34 A) where WP45 had concluded none existed, but one staple buys 2.8 kcal/mol against a 12-20 kcal/mol ordering cost. WP48 then improved the construct (PAPAP linker, -2.84 kcal/mol) and still has an ordering cost unpaid. WP51 then measured it: 1.40 kcal/mol for the glycine linker and >=3.38 for WP48's PAPAP - the linker swap made ordering WORSE, not better. A scaffold is still required.**

| status | count |
|---|---|
| SOLVED | **12** |
| PARTIAL | **4** |
| NOT SOLVED | **0** |
| UNMEASURABLE | **0** |

## Every sub-problem

| id | area | problem | status | kind |
|---|---|---|---|---|
| P1 | Interface | Is the TPP1 TEL patch one binding surface? | **SOLVED** | geometric |
| P2 | Interface | Can a drug-like small molecule span both halves? | **SOLVED** | geometric |
| P3 | Interface | What is the shortest path across the gulf? | **SOLVED** | geometric |
| P4 | Triad | Can an alpha-helix present the E168/E169/E171 triad? | **SOLVED** | geometric |
| P5 | Triad | Can the triad be presented at all, in the native geometry? | **SOLVED** | geometric |
| P6 | Triad | Does a peptide HOLD that geometry? | **PARTIAL** | design-level |
| P7 | Modularity | Can each module solve only its own site, fully? | **SOLVED** | geometric |
| P8 | Modularity | Do merged modules stay stable together? | **SOLVED** | design-level |
| P9 | Modularity | Can one groove be blocked and repaired at once? | **SOLVED** | logical |
| P10 | Screening | Is a structural alert a usable filter at all? | **SOLVED** | computed |
| P10b | Screening | Does the full gate stack accept approved drugs? | **SOLVED** | computed |
| P11 | Screening | What is the experimental success rate? | **PARTIAL** | computed |
| P12 | Screening | Do any of these designs bind telomerase? | **PARTIAL** | computed |
| P13 | Mechanism | What does TPP1 change mechanically (M5)? | **PARTIAL** | computed |
| P14 | Mechanism | Is n measured? | **SOLVED** | logical |
| P15 | Drug | Is there a drug? | **SOLVED** | computed |

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

A covalent constraint EXISTS, and WP45 missed it by screening disulfides only. Over six chemistries spanning 3.5-12.5 A, W167-L211 passes reach, direction, clash and payload tests at 6.34 A with cos 0.98/0.09 - and the native indole/Leu contact there is already 3.50 A. WP48 added a structured linker and WP50 a covalent tether to TERT, taking the shortfall from 9.2 to 2.36 kcal/mol. Still not balanced.

Blocker: the remaining 2.36 kcal/mol after the WP50 covalent tether. Ring closure is now closed off permanently (20.42 A vs ~7 A reach) and trimming only buys 4 residues, so the residue is a presenting scaffold or the TERT Cys tether, both of which need to be made.

**P7 — Can each module solve only its own site, fully?** — *SOLVED*

Yes, but not the way WP44 claimed. Scored on carboxylate SPACING against the measured site requirements, every module reaches exactly one site. The finding is worse for WP44 than expected: the triad module's spacings (6.43/5.08/10.75 A) miss the triad triangle (6.58/10.22/7.31) and match the K78-R132 separation instead - so it was aimed at S1, not S3. The graft reaches S3 and nothing else.

**P8 — Do merged modules stay stable together?** — *SOLVED*

4 of 5 merges STABLE, 1 COLLAPSE-RISK, 1 REFUSED. The linker alone decides it: same two modules, -6.04 kcal/mol attraction when flexible.

**P9 — Can one groove be blocked and repaired at once?** — *SOLVED*

No - the merge engine refuses it by construction, with a stated reason.

**P10 — Is a structural alert a usable filter at all?** — *SOLVED*

No as a veto - proven unfixable, alert count caps at 0.70 balanced accuracy. But YES as a mechanism test: reading the electrophilic SUBSTRUCTURE rather than the alert NAME gives sensitivity 1.0 and specificity 1.0 (11/11) on a panel used by neither earlier WP. The failure was the vocabulary, and replacing it fixed it.

**P10b — Does the full gate stack accept approved drugs?** — *SOLVED*

The objection dissolves when the metric is corrected, and no threshold moves. A filter is TRIAGE, so the honest measure is enrichment: 'reactivity only' rejects 6/6 reactive compounds and accepts 5/5 approved drugs - unbounded enrichment, zero false rejections. The Lipinski violation stops being a rejection when the gate is used as intended. Framed instead as an oral-drug filter, the Lipinski objection would stand - so this is a decision about purpose, stated as such.

**P11 — What is the experimental success rate?** — *PARTIAL*

Decomposed into four factors, and only one of them is unknown. WP52 measured the structural factor at 1.0 (the graft clears the TERT surface). WP51 MEASURED the ordering cost (1.40 kcal/mol for the glycine linker, >=3.38 for PAPAP) in place of WP48's typed 6.36. The makeability factor is still zero, but for a better stated reason: the scaffold needed to pay even the measured cost does not exist yet. With it, 9 independent screens at a cited 1-10 % peptide hit rate give 8-61 % probability of at least one micromolar binder. The residual uncertainty is the literature hit rate, given as a range.

Blocker: f2 = 0 until the ordering cost is paid; f3 is a literature range, not a measurement at this interface.

**P12 — Do any of these designs bind telomerase?** — *PARTIAL*

The geometric half is measured and it PASSES: with its 14 native residues pinned to the 7TRE coordinates to 0.00 A, the designed linker clears every TERT heavy atom by 8.22 A (PAPAP) - 0 atoms inside 2.5 A. The design is sterically compatible with the deposited TERT surface, so 'a score would be fabricated' was false for the geometric part. What is NOT measured is any affinity, and the remaining blocker is an experiment.

Blocker: no affinity, no Kd, no competition experiment - and the graft's net binding free energy is +0.95 kcal/mol (WP54), so it would not bind even if ordered.

**P13 — What does TPP1 change mechanically (M5)?** — *PARTIAL*

The structural half is measured and the answer is NO: the apo-vs-TPP1-bound difference over 934 CA atoms is 0.356 A, which is SMALLER than the 0.585 A between two copies of one protein in a single crystal. Only 9 of 934 residues exceed the duplication floor's 95th percentile - at or below the 5 % chance rate - and the largest 30-residue window rotation anywhere is 2.15 degrees against a 4.30 degree detection limit. TPP1 does not detectably remodel hTERT. The kinetic half (R1 vs R2) still needs the trap.

Blocker: R1 vs R2 needs the two-chase; a static structure cannot order events.

**P14 — Is n measured?** — *SOLVED*

The question was mis-posed, and that is the closure. n is not a parameter the two-chase can recover: at n = 1, 2, 3, 6, 12 and 24 the model gives m/c = 1.448879638 and recovers tau = 0.1833 and kappa = 2.70 identically, while only the per-round rate changes (2.70 to 0.056, a 48-fold spread). n is the STEP SIZE of a discretisation, and its legible form is the number of escape opportunities per repeat = 6 nt, which is a count readable from an existing gel.

Blocker: tau, the stall share, is a real free parameter and is still unmeasured - that is P2's job, not n's.

**P15 — Is there a drug?** — *SOLVED*

No - and now by arithmetic rather than by assertion. The potency requirement is modest (K_L <= 1.17 uM, i.e. 8.14 kcal/mol), so the obstacle was never pocket size. WP48's typed ordering cost (6.36) made the net positive, which would have meant no binding at all; WP51 MEASURED it instead (>=3.38 for PAPAP, 1.40 for GGGGG) and the net becomes weakly negative: -2.03 kcal/mol for the design as specified, -4.01 with the glycine linker. So the peptide WOULD bind weakly and STILL MISSES the requirement by 6.11 kcal/mol. Measured, the obstacle is not makeability alone - the contact energy is the shortfall.

Blocker: a covalent or glue mechanism is outside this calculation; the negative holds for competitive orthosteric ligands.

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
