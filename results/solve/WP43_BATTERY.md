# WP43 — Testing all ten designs, and what "success rate" can and cannot mean

    python scripts/solve/wp43_test_battery.py
      -> results/solve/wp43_test_battery.json
      -> results/solve/wp43_test_battery.csv
    python scripts/solve/wp43_test_battery.py --self-test   -> 48 checks, 0 failed

The request was to test every design and report the success rate. Ten designs
were run through a computed battery. **The honest answer has no single number,**
and the reason is worth more than the number would have been.

## The three rates, kept apart

| rate | value | what it is |
|---|---|---|
| **In-silico pass, designs only** | **6/10 = 60 %** | computed filters cleared. Not evidence of binding |
| **In-silico pass, designs + controls** | **6/15 = 40 %** | the same, with the calibration set in the denominator |
| **Experimental success rate** | **not measurable** | N = 0 designs synthesised, expressed or assayed |

Reporting the first as "the success rate" would be the single most misleading
thing this repository could do, because it is computed, not observed, and it is
always higher than the number anyone actually wants. The gap between 60 % and
the real rate is the entire content of this work package.

## The battery first had to prove it can say no

A filter that passes everything has a 100 % pass rate and predicts nothing. So
five controls were run before any design was scored: two approved drugs that must
pass, three known-bad molecules that must fail.

Two configurations were tested, because the standard alert catalogs **do not
agree**:

| control | expected | PAINS only | PAINS + BRENK |
|---|---|---|---|
| imatinib (approved) | pass | pass | pass |
| aspirin (approved) | pass | **pass** | **REJECT** ← false positive |
| curcumin (frequent hitter) | fail | **pass** ← miss | REJECT |
| p-benzoquinone (redox cycler) | fail | REJECT | REJECT |
| nitrogen mustard (alkylator) | fail | **pass** ← miss | REJECT |
| **accuracy** | | **0.60** | **0.80** |

Neither is free. PAINS alone misses two of three genuinely bad molecules;
BRENK catches them and rejects aspirin. The battery therefore reports **both**
and uses the strict configuration as the headline, while keeping the permissive
result beside it — because choosing one silently is how a 0.80 accuracy gets
quoted as if it were 1.00.

## Two gates were wrong, and the controls found them

**QED was going to be a hard gate at ≥ 0.50.** Imatinib scores **0.389**. QED was
calibrated on lead-like oral drugs and penalises the basic solubilising side chain
that every kinase inhibitor carries, so a 0.50 cutoff rejects part of the
approved-drug world before it looks at a design. QED is now *reported*, not gated.

**A single fractional score was an artefact.** With one "passed/total" number,
adding a gate changes the **denominator** — and adding BRENK flipped KP-2 from
FAIL to CONDITIONAL without any property of the molecule changing. Chemistry and
geometry are now scored on two independent axes, and neither can be improved by
counting more gates.

## Every design, under both configurations

| id | MW | logP | QED | SA | span needed | chemistry | geometry | tier | strict |
|---|---|---|---|---|---|---|---|---|---|
| KP-1 | 318.3 | 4.42 | 0.747 | 1.62 | 14.26 / 25.5 | PASS | **cannot span** | GEOMETRY-FAIL | same |
| KP-2 | 394.4 | 6.08 | 0.430 | 1.64 | 17.95 / 25.5 | 1 violation | **cannot span** | GEOMETRY-FAIL | same |
| TEL-1 | 270.3 | 2.61 | 0.875 | 1.58 | 14.27 / 20.2 | PASS | **cannot span** | GEOMETRY-FAIL | same |
| TEL-2 | 270.3 | 2.61 | 0.875 | 1.68 | 12.81 / 20.2 | PASS | **cannot span** | GEOMETRY-FAIL | same |
| **TEL-3** | 179.2 | 0.41 | 0.705 | 1.52 | **8.68 / 6.0** | PASS | **reaches** | **MODEL-PASS** | **same** |

| id | aa | GRAVY | charge | helix | stability | gates | tier |
|---|---|---|---|---|---|---|---|
| **TELM-1** | 48 | −0.02 | −3.41 | 0.71 | stable | 8/8 | MODEL-PASS |
| **TELM-2** | 28 | +0.37 | −3.38 | 0.71 | predicted unstable | 8/8 | MODEL-PASS |
| **TELM-3** | 50 | −0.23 | −3.41 | 0.66 | stable | 8/8 | MODEL-PASS |
| **TENB-1** | 48 | −0.05 | +4.58 | 0.71 | stable | 8/8 | MODEL-PASS |
| **TENB-2** | 40 | +0.44 | −1.41 | 0.60 | predicted unstable | 8/8 | MODEL-PASS |

**The chemistry is not the problem and never was.** Every ligand is drug-like on
In-silico terms — QED 0.43–0.88, SA 1.52–1.68 (all "easy to make"), Lipinski
satisfied. The four failures fail *only* on geometry, which is exactly the WP42
result arriving from the other direction: you can make these molecules, you just
cannot make them long enough.

## The rankings, with the tiebreak stated

All five proteins clear all eight gates, so a "best protein" cannot be decided by
the battery. The tiebreak is applied in a fixed, declared order — predicted
stability, then how many distinct measured contacts the design addresses, then
|GRAVY| as the crude solubility proxy — rather than hidden in a weighted sum.

**Ligands:** 1. **TEL-3** (only pass, +2.68 Å margin) · 2. TEL-1 · 3. TEL-2 ·
4. KP-2 · 5. KP-1.

**Proteins:** 1. **TELM-1** (stable, 3 contacts, |GRAVY| 0.02) · 2. TENB-1 ·
3. TELM-3 · 4. TELM-2 · 5. TENB-2.

Two of those orderings are worth arguing with:

- **TELM-1 over TENB-1** is decided entirely by |GRAVY|, 0.02 versus 0.05. That is
  not a meaningful difference and should not be quoted as one.
- **TENB-2 ranks last** despite being the most novel proposal in the set. It
  addresses one measured contact and carries the highest |GRAVY|. If the
  groove-repair idea is right, **this ranking is wrong**, and it is wrong because
  the battery cannot see mechanism — only composition.

## What 60 % is comparable to, and why the comparison is weak

| method | published experimental hit rate |
|---|---|
| AlphaProteo, de novo binders | 9–88 % per target (Zambaldi 2024) |
| RFdiffusion + ProteinMPNN binders, unfiltered | ~1–5 % |
| de novo small-molecule design | < 1 % |

The designed proteins here are the class those methods report on, so **the
relevant expectation is single digits to tens of percent, not 60 %.** The 60 % is
a filter pass rate. Quoting it next to a hit rate would be a category error, and
the table exists so that nobody makes it by accident.

## What could not be tested, and what that costs

No OpenMM, no docking engine, no FEP. So there is **no binding score, no pose and
no affinity in this file** — WP39 already refuted this repository's only
rigid-body scorer, and running it again would produce a number without evidence.

Untested, in order of how much it matters:

1. **Binding, of any kind.** The first experiment is the only evidence that counts.
2. **Whether the peptides fold as designed.** A 28-residue amphipathic peptide can
   be a helix or a molten globule; the ideal-helix geometry says nothing about which.
3. **Selectivity.** The K78/R132 groove is basic and basic surfaces are common, so
   off-target binding is the first thing an assay will find.
4. **Permeability.** TEL-3 carries one acid; the others carry two. That is a
   liability the battery records but cannot quantify.
5. **Immunogenicity** of the designed proteins.

The single most informative next measurement remains **TELM-2 versus TELM-2-E→A**
(WP42 §7, test 1). Everything else in this document is downstream of whether that
acid face is required.

## Tests

48 checks, 0 failed. They pin the things that would otherwise drift:

- the battery rejects a quinone in both configurations and catches all three bad
  controls in strict;
- **permissive does not reject aspirin and strict does** — the two configurations
  are asserted to differ, so the comparison cannot silently become vacuous (it did,
  on the first implementation);
- imatinib's QED is asserted below 0.50, which is why QED is not a gate;
- the experimental rate is `None`, not a number, and the reason is recorded;
- the pass rate with controls is strictly lower than the designs-only rate;
- the four designs that cannot span are `GEOMETRY-FAIL`, and chemistry is *not*
  the axis they fail on;
- the best ligand is the only one that reaches, and no geometry failure outranks it;
- every design carries a stated rank basis.
