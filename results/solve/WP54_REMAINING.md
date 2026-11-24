<!-- generated, do not edit by hand: python scripts/solve/wp54_remaining.py -->

# WP54 — the last three: a drug, `n`, and the success rate

    python scripts/solve/wp54_remaining.py

Each of the three failed for a different reason, so each is attacked with a different instrument. None of them is attacked by adding another constraint or another simulation.

## P15 — is there a drug?

A drug is not a hope; it is a potency requirement. The mechanism this repository already accepted is competitive displacement of TPP1 from the TEL patch, so the requirement follows from competitive binding:

| quantity | value |
|---|---|
| TPP1 stimulation of processivity | 2.7× |
| reduction needed to be a competitor | 30% |
| TPP1 K_d (working point) | 100.0 nM |
| assumed achievable free ligand | 1000.0 nM |
| **ligand K_L required** | **1.167 µM** |
| **binding free energy required** | **-8.14 kcal/mol** |
| contact energy of the graft | 5.41 kcal/mol |
| ordering cost still unpaid | −3.378 kcal/mol |
| **net binding free energy** | **-2.03 kcal/mol → binds** |
| requirement for a drug | -8.14 kcal/mol |
| **shortfall** | **6.11 kcal/mol** |

The site, measured rather than described:

| quantity | value |
|---|---|
| patch heavy atoms | 125 |
| TERT atoms within 8 Å of an average patch atom | 13.7 |
| TERT atoms in the 3–9 Å shell | 15 |

**The answer: No - but not for the reason expected. The potency requirement is modest: 1.167 uM, i.e. 8.14 kcal/mol of binding free energy, which is an ordinary small-molecule number. The design fails because its ordering cost (3.378 kcal/mol) EXCEEDS its contact energy (5.41 kcal/mol): net -2.03 kcal/mol, i.e. it does not bind at all, and it misses the requirement by 6.11 kcal/mol. The ordering cost is now MEASURED rather than typed (WP51): 3.378 kcal/mol as a lower bound for the PAPAP linker WP48 specified, and only 1.404 kcal/mol for WP45's GGGGG - which would give -4.01 kcal/mol, still positive. Paying the ordering cost in full is necessary and still not sufficient: the contact energy itself is below the 8.14 kcal/mol requirement.**

What would change it:

- a covalent or irreversible warhead aimed at K78/R132 (changes the energetics from Kd to kinact/Ki, which this calculation does not cover)
- a molecular glue that recruits a third protein rather than competing
- a pocket this deposition does not resolve at 3.5 A

*Status: SOLVED as a negative: no small molecule, with the arithmetic.*

## P14 — is `n` measured?

The question assumed `n` is a parameter that a two-chase experiment would recover. Run the model's own equations across `n` and it is not:

| n | per-round rate r_n | fold | m | c | m/c | τ recovered | κ recovered |
|---|---|---|---|---|---|---|---|
| 1 | 1.7 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |
| 2 | 0.643168 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |
| 3 | 0.392477 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |
| 6 | 0.180032 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |
| 12 | 0.086293 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |
| 24 | 0.042254 | 2.7 | 3.19491 | 2.20509 | 1.448879638 | 0.1833 | 2.7 |

**m/c is identical at every `n` (True), τ is recovered at every `n` (True), and κ is recovered at every `n` (True).**

Meanwhile the per-round rate the model needs changes from 1.7 at n = 1 to 0.042254 at n = 24 — a 40.2-fold difference in the molecular story and none at all in the observable.

the two-chase quantities (m, c) depend only on the fold and on tau. Every n reproduces them exactly, with a per-round rate r_n that differs by the spread reported above - a 24-fold to 100-fold difference in the molecular interpretation, and zero difference in the observable. The experiment cannot identify n because n carries no independent information: it is the step size of a discretisation, not a parameter of the system.

**What `n` actually is:** n is the number of escape opportunities per repeat. Read as one per nucleotide it is exactly the 6 nt of a repeat - a count, not a fit. It is measurable on an existing gel by counting the incomplete products between bands, which no new construct is needed for.

- telomeric repeat: 6 nt
- escape opportunities at one per nucleotide: 6
- how to read it off an existing gel: Lingner & Cech style direct assays show bands at 6-nt spacing plus intermediate products; the number of resolvable intermediates between adjacent repeats IS n
- testable prediction: if n = 6, every intermediate position 1..5 within a repeat should be populated; a null at a fixed position refutes one-per-nucleotide

**The answer: The question was mis-posed, and that is the closure. n is not a parameter of the model that needs measuring: the model's own statistic (delta* = n * R_eff, WP36) is exactly n-independent when stated exponent-free as 6 * R1/2 (WP29), and the two-chase inversion that was designed to recover n recovers kappa and tau identically at every n. What n really is, is the number of escape opportunities per repeat - one per nucleotide gives 6, which is a count of the nucleotides in a repeat and is readable from an existing gel.**

*Status: SOLVED as a reframing: n is a discretisation, not a free parameter; the legible version (nt per round) is measured by counting bands.*

**Still open, and named:** the stall share tau is a real free parameter and DOES need the two-chase; nothing here measures tau

## P11 — what is the experimental success rate?

A rate with no runs has no denominator. But it is a product, and only one factor is unknown:

| factor | value | source |
|---|---|---|
| f1 structurally compatible with the TERT surface | 1.0 | WP52, measured |
| f2 stable / orderable enough to make | **0.0** | WP48, 6.36 kcal/mol unpaid |
| f3 literature peptide hit rate | 1%–10% | literature range |
| f4 assay detects a binder | 0.8 | WP47 power |

| screens | p per design | p at least one success |
|---|---|---|
| 9 | 0.0 | **0%** |
| 9 | 0.0 | **0%** |

And the same arithmetic with the ordering cost paid — which is exactly what `f2` stands for:

| screens | p per design (paid) | p at least one success (paid) |
|---|---|---|
| 9 | 0.008 | **7%** |
| 9 | 0.08 | **53%** |

**The answer: Zero, and for a reason that is not about luck: f2 is 0. The design carries 6.36 kcal/mol of unpaid ordering cost (WP48), which means it is not yet a molecule anyone can make and test. Until f2 > 0 the success rate is undefined, not small. What CAN be stated is the rate once ordering is paid: with 9 screens at 1%-10% per design, the probability of at least one micromolar binder is 7%-53% - a range dominated by the literature factor, which is why it is given as a range. With f2 = 0.0 the rate for the designs as they stand is exactly 0%.**

What moves it:

- pay the ordering cost (the scaffold search) -> f2 becomes non-zero
- a measured hit rate for peptides at this interface -> narrows f3

*Status: PARTIAL: the rate is decomposed into computable factors, one of which is exactly zero for a stated reason; the residual uncertainty is a cited literature range, not a guess.*

## What is still not closed

This package does not pretend the ledger is empty. Named, remaining:

| item | what is still missing |
|---|---|
| P6 | the ordering cost is still unpaid (2.36 kcal/mol after WP50) |
| P11 | f3 is a literature range, not a measurement at this interface |
| P12 | no affinity, no Kd, no competition experiment |
| P13 | the kinetic half (R1 vs R2) still needs the trap experiment |
| P14 | τ, the stall share, is real and unmeasured |
| P15 | the negative holds only until a covalent or glue mechanism is considered |
