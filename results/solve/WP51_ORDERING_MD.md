<!-- generated, do not edit by hand: python scripts/solve/wp51_ordering_md.py -->

# WP51 — the ordering cost, measured

    python scripts/solve/wp51_ordering_md.py

## Why this work package exists

The ordering cost has been a **model** through five packages:

| package | the number | how it was obtained |
|---|---|---|
| WP45 | 12–20 kcal/mol | typed |
| WP47 | 0.6 kcal/mol per torsion | typed, then multiplied |
| WP48 | same price, new linker | subtracted from itself |
| WP50 | 4 kcal/mol tether gain | a conventional constant |

**Every comparison between designs subtracted the same unmeasured quantity
from itself.** OpenMM 8.6.1 and RDKit are installed here, so the linker can
be built as atoms and its conformational distribution sampled.

## The first attempt at this measurement was wrong

The first version of this script replaced the linker with a harmonic
restraint whose stiffness it **set by hand** (5 / 40 / 200 kcal/mol/Å² for
GGGGG / PAPAP / AAAAA) and then reported the standard deviation of the
pharmacophore distances. It announced that AAAAA was stiffest — which was
the value that had been typed in. The ranking was installed, not measured.

That is exactly the failure mode this package exists to end — reading the input back out as if it were output — so it is recorded here rather than quietly deleted. What follows does not assign a stiffness. It measures a distribution.

## The restraint was INVERTED, and that produced the first 'result'

The sampler uses a flat-bottom restraint so the free peptide cannot unfold
out of the box. The first version read:

    0.5*k*max(r - r_on, 0)^2  +  0.5*k*max(r_off - r, 0)^2

The second term is non-zero for **every** r below `r_off`. That is not a
flat bottom: it is a harmonic well centred at (r_on + r_off)/2 = 14 A.

The consequence was unmistakable in hindsight: **nine independent builds of
three chemically different linkers all returned r = 14.00 ± 0.12 A**, with a
spread of 0.11 A for a flexible pentapeptide — which is not a
conformational distribution at all. The run was measuring its own restraint.

The correct form is `max(r - r_off, 0)^2 + max(r_on - r, 0)^2`, which is
exactly zero inside the window and walls only outside it. It is now
unit-checked at 5, 14 and 20 A in `--self-test` rather than read for
correctness.

## The quantity that needs pricing

The linker occupies the junction between TPP1 residues 166–172 and
209–215. In 7TRE that junction spans:

| measurement | value |
|---|---|
| C(172)–N(209), the gap the linker fills | **16.31 Å** |
| CA(172)–CA(209) | 16.66 Å |
| N(166)–N(209) | 7.35 Å |
| native residues deleted | 36 |

This is **not** WP50's ring-closure distance C(215)–N(166) = 20.42 Å,
which asks a different question: closing the whole graft into a macrocycle.
The linker only has to span the junction, and the junction is narrower.

So the quantity to measure is: **how much free energy does it cost to hold
this chain's two ends at 16.31 Å?** That is a conformational entropy, and it is
read straight off the free chain's own distribution:

    ΔG_order = −kT · ln( p_window )

where `p_window` is the fraction of free conformations with an end-to-end
distance inside a window around 16.31 Å.

## Stage 1 — the linkers, all atoms, free in solvent

Each linker built by RDKit from sequence, hydrogenated, GB-Neck2, 300 K,
Langevin, 2 fs, 0.4 ns × 3 seeds. A flat-bottom restraint (zero
force between 3 and 25 Å) keeps the peptide in the box; the fraction of
frames touching it is reported, because a non-zero value would mean the
restraint is perturbing the measurement.

| linker | mean r | **σ** | 5–95 % | frames at boundary |
|---|---|---|---|---|
| GGGGG | 10.91 Å | **3.61 Å** | 4.16–16.02 Å | 0.000 |
| PAPAP | 10.11 Å | **2.26 Å** | 6.17–13.26 Å | 0.000 |
| AAAAA | 9.3 Å | **2.16 Å** | 5.59–12.36 Å | 0.000 |

### The ordering cost

| linker | p(±0.5 Å) | p(±1.0 Å) | p(±2.0 Å) | **ΔG_order (±1.0 Å)** |
|---|---|---|---|---|
| GGGGG | 0.060 | 0.093 | 0.190 | **1.404 kcal/mol** |
| PAPAP | 0.000 | 0.000 | 0.003 | **3.378 kcal/mol** (lower bound) |
| AAAAA | 0.000 | 0.000 | 0.000 | **3.378 kcal/mol** (lower bound) |

Cheapest: **GGGGG**.

## Stage 2 — does the pharmacophore survive a soft linker?

The two native segments, actual crystal coordinates, with the linker
reduced to one harmonic term whose stiffness is **derived** from stage 1
(k = kT/σ²) rather than chosen. This is the question the earlier packages
never asked: not *is the linker cheap*, but *does the pharmacophore stay
together at all* once a realistic tether is in place.

| linker | k (kcal/mol/Å²) | mean d(168–215) | **σ** | mean d(171–210) | **σ** | RMSD |
|---|---|---|---|---|---|---|
| GGGGG | 0.0454 | 19.14 Å | **0.57** | 13.71 Å | **1.01** | 3.75 Å |
| PAPAP | 0.1159 | 17.96 Å | **0.98** | 14.33 Å | **0.98** | 3.56 Å |
| AAAAA | 0.1269 | 18.27 Å | **0.89** | 14.28 Å | **0.9** | 3.83 Å |

**H1 (pharmacophore retained).** Worst σ over all constructs: 0.98 Å and 1.01 Å against a 3.0 Å threshold → **retained**.

**H2 (the native fold does not melt).** Worst mean all-atom RMSD from
the crystal start: 3.83 Å → **drifts**.

## What this can and cannot say

**It can say** what each linker's conformational distribution actually is,
what it costs to confine it to the native junction distance, and whether
the two native segments stay in register once it is in place. Those are
measured numbers with a measured spread, not conventions.

**It cannot say** that any construct folds, and the caveats are not
decoration:

- The stapled construct cannot be built: the non-natural residue template needs files.rcsb.org, which is unreachable. Staples are absent from every run; this compares peptide linkers only.
- No folding is sampled. The linker runs are 1 ns per replica on 2 CPU threads; the constructs sit near their crystal start.
- Stage 2 is 0.01 ns per replica. That was set by the measured cost of the 1547-atom construct (63.3 ms/step, 17.6 h per ns on 2 threads), not by what the question deserves. It is long enough to see whether the tether holds the two segments together and far too short to say anything about relaxation.
- dG_order = -kT ln p_window counts only the radial confinement. The angular part is omitted, so every value is a LOWER BOUND.
- Proline cis/trans isomerism is not sampled; PAPAP's repertoire is whatever RDKit's embedder chose.
- Stage 2 replaces the linker with ONE harmonic term. Its stiffness comes from the stage-1 measurement, but its functional form is an assumption.

## The honest verdict

The measured ordering costs span **1.404–3.378 kcal/mol** at the ±1 Å
window, and the spread between linkers is real rather than a tie.

**GGGGG is the cheapest construct measured** (p = 0.093 inside ±1 Å of the 16.31 Å junction, 1.404 kcal/mol), while PAPAP never came within ±1 Å in 300 frames, so its cost is a lower bound at 3.378 kcal/mol.

That ordering is **the opposite of what WP48 assumed.** WP48 replaced
WP45's flexible GGGGG linker with PAPAP and credited the change with
−2.84 kcal/mol on the grounds that proline is more rigid. Measured, the
rigid linker is *worse* at reaching the 16.31 Å junction: a stiffer chain
cannot span far, and PAPAP's population sits at 10.11 Å. The
WP48 saving was subtracted from itself rather than measured, and when
measured it changes sign.

The same applies to WP47's 0.6 kcal/mol per torsion, six torsions per
staple: 3.6 kcal/mol for a construct nobody measured. The measurement
does not confirm that price.

Because only the radial part of the confinement is counted, every value
here is a **lower bound**. It is the right quantity to compare linkers
with, and the wrong quantity to add to a binding free energy as if it were
the whole price.
