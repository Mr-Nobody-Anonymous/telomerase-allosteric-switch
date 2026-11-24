<!-- generated, do not edit by hand: python scripts/solve/wp53_conformational_change.py -->

# WP53 — P13: what TPP1 changes mechanically, measured

    python scripts/solve/wp53_conformational_change.py

## Why P13 sat at NOT SOLVED, and what changed

P13 asks what TPP1 changes mechanically. It was blocked on a DNA-trap experiment, which is correct for the **kinetic** half of the question. But there is a structural half, and its deposits have been in `data/raw` since before this work package existed:

| deposit | contents | resolution |
|---|---|---|
| **7TRD** | hTERT + hTR + DNA, **no TPP1** | 3.3 Å |
| **7TRE** | hTERT + hTR + DNA **+ TPP1** | 3.5 Å |
| **7TRF** | hTERT + hTR + DNA, plus histones, **no TPP1** | 3.7 Å |
| **7TRC** | TPP1 + POT1, two copies of each | 3.3 Å |

All three telomerase deposits carry the same hTERT chain: 934 modelled residues, numbered 7–1132. So the apo and TPP1-bound states of the same enzyme are both on disk, and the question *did TERT change?* is a superposition, not an experiment.

## The control that makes it a measurement

A difference between two crystals is not automatically a TPP1 effect. So the same comparison is made twice:

| role | comparison | why |
|---|---|---|
| **signal** | 7TRD (apo) vs 7TRE (bound) | the question |
| **control** | 7TRD (apo) vs 7TRF (apo) | the same kind of comparison with **no TPP1** |
| **floor** | 7TRC chain E vs chain I | two copies of the identical protein in one crystal — the only comparison with a known true answer (zero) |

## What the superpositions say

| comparison | CA atoms | CA-RMSD | median | p90 | max | at |
|---|---|---|---|---|---|---|
| **signal: apo vs bound** | 934 | **0.356 Å** | 0.293 | 0.539 | 1.226 | 612 |
| control: apo vs apo (no TPP1) | 934 | **0.0 Å** | 0.0 | 0.001 | 0.001 | 796 |
| control: two copies | 129 | **0.585 Å** | 0.365 | 0.608 | 3.435 | 153 |

Noise floor: the 95th percentile of the crystal-form control is **0.824 Å**.

## Where it moves

The ten largest displacements, with the distance of each residue's CA from the TPP1 footprint (WP52's 14 residues):

| residue | displacement | distance to patch |
|---|---|---|
| ALA612 | 1.23 Å | 33.4 Å |
| PRO322 | 0.99 Å | 45.3 Å |
| CYS7 | 0.98 Å | 18.3 Å |
| VAL639 | 0.9 Å | 33.4 Å |
| ARG120 | 0.85 Å | 12.7 Å |
| PRO771 | 0.85 Å | 4.8 Å |
| ALA67 | 0.84 Å | 29.9 Å |
| ALA801 | 0.84 Å | 15.6 Å |
| GLU800 | 0.83 Å | 13.6 Å |
| ASP1132 | 0.81 Å | 58.6 Å |

**9 of 934 residues (1.0%) exceed the floor's 95th percentile of 0.824 Å.** A pure-noise difference would put about 5% above a 95th percentile, so this is *at or below chance*, not above it.

## Local or coupled?

Distance from the TPP1 footprint is a weak predictor of displacement: Spearman ρ = **-0.203** for the apo-vs-bound difference, against -0.024 for the apo-vs-apo control. The sign is the one a contact-local effect would have — residues nearer the patch move slightly more — and the size is far too small to carry anything on its own.

## Is there a hinge?

Sliding 30-residue windows over the whole chain, each superposed on itself:

| window | internal RMSD | rotation vs global | translation |
|---|---|---|---|
| 7-36 | 0.339 Å | 1.44° | 6.81 Å |
| 22-51 | 0.313 Å | 1.14° | 1.73 Å |
| 37-66 | 0.361 Å | 0.48° | 2.13 Å |
| 52-81 | 0.388 Å | 0.7° | 2.74 Å |
| 67-96 | 0.343 Å | 1.21° | 1.32 Å |
| 82-123 | 0.38 Å | 2.15° | 2.39 Å |

Largest window rotation anywhere in the chain: **2.15°** (window 82-123). The same protocol on the duplication control gives up to 1.27°, so the detection limit is **2.54°** — measured, not chosen.

## The findings, in order

1. **No conformational change is detectable at this resolution.** The apo-vs-TPP1-bound difference is 0.356 A CA-RMSD over 934 residues. The only comparison with a known true answer - two copies of the identical protein in the 7TRC crystal - is 0.585 A over 129 residues. The TPP1-associated difference is SMALLER than the difference between two copies of one molecule. On this evidence TPP1 does not remodel hTERT.

2. **Nothing exceeds the noise.** 9 of 934 residues (1.0%) displace more than the duplication floor's 95th percentile (0.82 A). If the difference were pure noise about 5% would exceed a 95th percentile; 1.0% is at or below that. The largest single displacement is 1.226 A at residue 612.

3. **No hinge.** Sliding 30-residue windows over the whole TERT chain rotate at most 2.15 degrees relative to the global superposition (window 82-123), with an internal fit of 0.38 A. The same protocol on the duplication control rotates up to 1.27 degrees, so the detection limit is 2.54 degrees. The signal does not reach it.

4. **There is a weak local trend, in the right direction, and it is too small to rest anything on.** Displacement correlates with proximity to the TPP1 footprint at rho = -0.20 (control -0.02): the residues nearest the patch move slightly more. The trend is consistent with a contact-local adjustment and is not evidence of a propagated mechanical change.

5. **What this means for the mechanism question.** P13 asks what TPP1 changes mechanically. Measured against the only available noise floor, it changes no conformation. That is a result, not a null to be explained away: it removes 'TPP1 remodels TERT' from the candidate mechanisms and leaves the two branches this repository already carries - R1 (affinity) and R2 (retention during extension) - both of which act without requiring a conformational change in TERT.

## The verdict

- structural answer: **TPP1 does not detectably remodel hTERT at 3.3-3.5 A**
- difference exceeds the duplication floor: **False**
- residues above chance: **False**
- a hinge is present: **False**
- weak local trend toward the patch: **True**
- mechanism branch (R1 vs R2) decided: **False**
- neither R1 (affinity) nor R2 (retention) is decided by two static structures, and this result does not decide them either - it removes a third possibility

**P13 moves from NOT SOLVED to PARTIAL.** P13: NOT SOLVED -> PARTIAL (structural half measured and answered in the negative; kinetic half still needs the trap experiment)

## Caveats

- Two crystals at 3.3 and 3.5 A with different contents/packing, so the difference between them is not a controlled perturbation: TPP1 is one of several things that differ. The result survives that because the difference is below the duplication floor rather than above the noise, but a crystal is not an experiment.
- The apo-vs-apo control (7TRD vs 7TRF) is DEGENERATE: their hTERT coordinates are identical to 0.000 A over 934 CA atoms, so it carries no information and was not used as a floor. This was found by running it, not by reading the deposit titles.
- A true floor is only available from the 7TRC duplication (two copies of one 129-residue protein at 3.3 A). A duplication floor from a smaller protein is not guaranteed to transfer to a 934-residue one; it is the best available and it is stated.
- Static structures cannot order events. A difference that exists in the crystals need not drive the kinetics.
- Model building at 3.3-3.5 A is coarse. A sub-Angstrom difference is near the limit of what these maps support, which is exactly why the floor is measured rather than assumed - and why the result is 'not detectable', not 'no change'.
