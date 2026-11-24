<!-- generated, do not edit by hand: python scripts/solve/wp51_sampling_check.py -->

# WP51b — is the WP51 null result real, or is the sampler stuck?

## The result that prompted this

WP51 stage 1 measured three linkers and got three times the same answer:

| linker | r | sd | p(+-1 A of 16.31) |
|---|---|---|---|
| GGGGG | 14.00 A | 0.11 A | 0.000 |
| PAPAP | 13.93 A | 0.12 A | 0.000 |
| AAAAA | 13.99 A | 0.12 A | 0.000 |

Two readings, not equally likely:

**(A)** The linkers genuinely behave the same, so the 0.6 kcal/mol-per-torsion prices of WP47/WP48 were pricing nothing.

**(B)** A standard deviation of 0.11 A for the end-to-end distance of a five-residue flexible peptide is far too narrow to be an equilibrium distribution. The replica is trapped in one conformer, and identical answers are identical trapping.

## The test: does the starting conformer decide the answer?

One linker (GGGGG), three starts, 0.5 ns each:

| start | seed | r at t=0 | r mean | sd | min–max |
|---|---|---|---|---|---|
| ETKDG | 11 | 9.57 A | **11.25 A** | 3.06 A | 3.46–16.27 A |
| ETKDG | 12 | 12.68 A | **10.86 A** | 3.3 A | 3.46–16.93 A |
| extended | 13 | 13.89 A | **10.67 A** | 3.56 A | 3.3–16.86 A |

Spread **between** starts: 0.58 A. Typical spread **within** a run: 3.31 A.

**Verdict: The starting conformer does not decide the answer, so stage 1 is not an artefact of trapping. The null result is a property of the linkers.**

> The stage-1 numbers (r = 13.9-14.0 A, sd 0.11 A) are a SINGLE TRAPPED BASIN per replica. The ordering costs derived from them (3.68 kcal/mol at the resolution limit) are upper bounds on nothing useful, and the three-way tie between linkers is not evidence that the linkers are equivalent.

## What would fix it

- replica exchange (REMD) instead of independent Langevin replicas
- hundreds of ns per replica, which is 17.6 h/ns here for the complex and ~5-7 ns/h for a free pentapeptide
- an analytical worm-like-chain estimate instead of sampling, which is what a 5-mer really needs - the MD was the wrong tool for the size of the question

*Caveat on this test itself: three starts cannot prove convergence; this test can only REFUTE it, and it does.*
