# WP26 — What is missing for a **simulation-only** resolution

The question: what must simulation supply to settle this before any wet-lab
transfer? Answer: **the structural half is closable; the kinetic half is not,
and no amount of compute changes that.**

    python scripts/solve/wp26_simulation_readiness.py
    -> results/solve/wp26_simulation_readiness.csv

## First, the uncomfortable part

The repository does contain molecular dynamics, and `docs/METHODOLOGY.md` states
it honestly — but the numbers are worth putting in one place:

| | on disk | needed for a converged answer | gap |
|---|---|---|---|
| simulated time / replica | **3.0 ps** (1000 steps × 2 fs) | ~100 µs | **3.3 × 10⁷×** |
| system size | **1,547 atoms** | ~150,000 | 97× |
| residue coverage | **97 of 1132 (8.6 %)** | 100 % | — |

The catalytic aspartates (**D868** and partners) are **not in the simulated
system**. The active site cannot be observed at all.

The timescale gap is **~8 orders of magnitude**. No re-analysis of those frames
can support a statement about TPP1's mechanism. The `p < 0.001` in
`md_replicate_summary.txt` is a real statistic about a 3 ps relaxation from a
starting structure, not about an equilibrium ensemble. With 3 replicas
(df = 4) the test also cannot resolve anything smaller than ~3× its own noise.

This does not undermine the manuscript — `paper/MANUSCRIPT.md` does not use MD
for any claim. It does mean MD cannot be recruited retroactively to close a gap.

## The four routes a simulation could take

| route | time | atoms | would settle |
|---|---|---|---|
| **S1** conformational equilibrium | ~100 µs | 150 k | M1 vs M4 |
| **S2** communication pathway | static | 40 k | allostery yes/no |
| **S3** alchemical free energy | ~1 µs/window | 150 k | is the interface real? |
| **S4** kinetics (`k_off`, `k_tr`) | ~10 ms effective | 150 k | M5 vs data |

## Ranked by (what it settles) ÷ (cost)

**1. S3 — alchemical free energy on the TPP1–TEN interface.** The only route
that yields a number a lab can *directly* falsify: ΔΔG for G100A, K757E and the
TEL-patch mutants Tesmer et al. already measured. ~1 µs per λ-window is routine
on a GPU cluster. This tests the foundation every other claim rests on.

**2. S2 — communication, with an all-atom distance-decayed graph.** Cheap, and
WP4 already demonstrated how to build the null correctly. But note WP23's
measured null: **no domain moves more than 0.53 Å** on TPP1 binding, so a
large-scale conduit is already disfavoured by experiment.

**3. S1 — equilibrium sampling.** Expensive *and* largely pre-empted: WP23
measured the endpoint difference at 0.435 Å, i.e. at the coordinate noise floor.
Do this only after S3 says the interface carries real energy.

**4. S4 — kinetics. Not feasible at the required accuracy.** Simulation delivers
`k_off` to perhaps a factor of 2–3; the tension to resolve is a factor of
**1.46**. *This is the key limit: no simulation will settle gap 2.*

## What simulation-only can and cannot close

| | closable by simulation? |
|---|---|
| Is the TPP1–TEN interface energetically real? | **yes** (S3) |
| Is a mechanical conduit plausible? | **yes** (S2, better graph) |
| The 2.73× vs ≥4× tension | **no** — needs precision nothing delivers |
| Gate 3 | **no** — it is an observation-model question about gel data |
| M2 vs M4 | **no** — kinetic and degenerate by construction |

## The honest summary

A simulation-only resolution is **not** achievable for the whole question. It
**is** achievable for the structural half: S3 plus a properly-controlled S2 would
move the interface claim from *geometrically plausible* to *energetically
quantified* — which is the strongest statement these methods support.

Beyond that, the binding constraint is not compute. It is that the two open
kinetic questions are **degenerate in the observable**, so no simulation of any
length can separate M2 from M4 without a direct `k_tr` measurement, and the
2.73-vs-4× tension needs 1.46× precision on an off-rate.

**The single highest-value simulation is S3**, and it answers a question the lab
can check independently in the same session — which is exactly what makes it the
right thing to do before transferring to a bench.
