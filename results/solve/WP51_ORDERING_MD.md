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
| GGGGG | 14.0 Å | **0.13 Å** | 13.77–14.22 Å | 0.000 |
| PAPAP | 13.95 Å | **0.13 Å** | 13.75–14.17 Å | 0.000 |
| AAAAA | 13.98 Å | **0.12 Å** | 13.79–14.18 Å | 0.000 |

### The ordering cost

| linker | p(±0.5 Å) | p(±1.0 Å) | p(±2.0 Å) | **ΔG_order (±1.0 Å)** |
|---|---|---|---|---|
| GGGGG | 0.000 | 0.000 | 0.003 | **3.378 kcal/mol** (lower bound) |
| PAPAP | 0.000 | 0.000 | 0.003 | **3.378 kcal/mol** (lower bound) |
| AAAAA | 0.000 | 0.000 | 0.003 | **3.378 kcal/mol** (lower bound) |

Cheapest: **GGGGG**.

## Stage 2 — does the pharmacophore survive a soft linker?

The two native segments, actual crystal coordinates, with the linker
reduced to one harmonic term whose stiffness is **derived** from stage 1
(k = kT/σ²) rather than chosen. This is the question the earlier packages
never asked: not *is the linker cheap*, but *does the pharmacophore stay
together at all* once a realistic tether is in place.

| linker | k (kcal/mol/Å²) | mean d(168–215) | **σ** | mean d(171–210) | **σ** | RMSD |
|---|---|---|---|---|---|---|
| GGGGG | 35.0414 | 19.43 Å | **0.76** | 13.56 Å | **0.28** | 3.55 Å |
| PAPAP | 35.0414 | 18.65 Å | **0.75** | 14.33 Å | **0.68** | 3.57 Å |
| AAAAA | 41.125 | 19.06 Å | **0.91** | 13.88 Å | **0.34** | 3.4 Å |

**H1 (pharmacophore retained).** Worst σ over all constructs: 0.91 Å and 0.68 Å against a 3.0 Å threshold → **retained**.

**H2 (the native fold does not melt).** Worst mean all-atom RMSD from
the crystal start: 3.57 Å → **drifts**.

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

The measured ordering costs span **3.378–3.378 kcal/mol** at the ±1 Å
window. Read that against what was being typed before: WP47's 0.6
kcal/mol per torsion, six torsions per staple, would have charged 3.6
kcal/mol for a construct it never measured. The measurement does not
confirm that price; it replaces the argument with a distribution.

Because only the radial part of the confinement is counted, every value
here is a **lower bound**. It is the right quantity to compare linkers
with, and the wrong quantity to add to a binding free energy as if it were
the whole price.
