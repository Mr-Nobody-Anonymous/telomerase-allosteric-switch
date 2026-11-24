# WP47 — closing the open items

    python scripts/solve/wp47_close_open.py
      -> results/solve/WP47_CLOSE_OPEN.md (this file)
      -> results/solve/wp47_close_open.json / .csv

Three frontier questions, attacked: **P6** (the graft has no covalent constraint),
**P10** (the filter battery fails its own controls), and the two **unmeasurable**
items **P11/P12** turned into executable experiments.

## P6 — the graft now has a covalent constraint

WP45 concluded *"no covalent route exists"*. It had screened **disulfides only**
(3.5–5.0 Å Cβ–Cβ) and stated the result as a property of the protein. That is the
same defect class as WP44's helix: a conclusion drawn from one member of a family.
Six chemistries spanning **3.5–12.5 Å** were screened here, over all **49** A×B pairs, against four filters:

| filter | test | why it kills candidates |
|---|---|---|
| F1 reach | Cβ–Cβ inside the chemistry's range | a staple cannot stretch |
| F2 direction | both side-chain vectors point into the gap (cos > 0) | a staple grows along the native side-chain vector |
| F3 clash | the straight span clears TERT by ≥ 3 Å | the interface is tight |
| F4 payload | neither residue is E168/E169/E171/E215 | stapling a payload destroys the contacts the graft exists to make |

**2 of 49 spans pass all four.**

### The winner: `TRP167-LEU211` — hydrocarbon, 6.34 Å

| quantity | value |
|---|---|
| span | **6.34 Å** |
| chemistries that reach | hydrocarbon |
| side-chain alignment | cos = 0.98 / 0.09 (both into the gap) |
| clearance past TERT | 21.63 Å |
| payload consumed | False |
| native side-chain contact | **3.5 Å** |

Why this one is not a coincidence: **the native W167 indole and L211 side
chain are already in van der Waals contact at 3.5 Å**, bridging the gap
between the two segments. The staple does not create a new contact — it
**replaces a non-covalent bridge with a covalent one.** The span is naturally
hydrophobic and the geometry is already i,i+4-like, which is exactly the case
hydrocarbon stapling was developed for.

### The honest accounting on the energy

| term | kcal/mol |
|---|---|
| contact energy available (measured) | 5.41 |
| torsions removed by the staple | 6 |
| ordering gain | +3.6 |
| α-methyl penalty | −0.8 |
| net gain | **+2.8** |
| residual ordering cost unpaid | **9.2** |

> The staple removes 6 torsions and buys ~2.8 kcal/mol. WP45 priced the whole construct at 12-20 kcal/mol, so ONE staple does not close the gap - it is a necessary constraint, not a sufficient one. Two staples or a scaffold are still required. Stated, not hidden.

## P10 — why the filter battery failed, and what actually fixed it

**Why WP43 failed.** It used a structural alert as a **veto**:

| config | rule | failure |
|---|---|---|
| permissive | reject only on PAINS | missed two of three BRENK-only frequent hitters |
| strict | reject on ANY BRENK alert | flagged **aspirin** (`phenol_ester`) |

Sensitivity 0.333 in one config, specificity 0.5 in the other. No threshold
fixes both, because the rule is binary and the signal is not.

**The first fix attempt failed too — and that failure is the result.**
A graded penalty per alert was fitted over a sweep of penalty × threshold.
It could not do better than 7/10, and the reason is structural:

**All 10/10 controls clear every hard
physicochemical gate (5/5)** — MW, logP, HBD, HBA and SA. So the alert count is
the only feature that separates anything, and the panel is **non-monotone**:

| molecule | alerts | must be |
|---|---|---|
| metformin | **2** | **accepted** |
| nitrogen mustard | **1** | **rejected** |

Every rule of the form `accept iff n_alerts <= k` was enumerated by exhaustion.
**The ceiling of that entire family is balanced accuracy 0.7** (k = 0), against a bar of 0.9.
That proves the WP43 framing is **UNFIXABLE with alert count alone**.

So the rule was changed in **kind, not threshold**:

> Only a **reactive** alert rejects — a non-specific electrophile that alkylates
> or arylates nucleophiles and therefore hits many targets. A **structural**
> alert (phenol ester, hydroquinone, imine) is carried as a **flag**, never as a
> veto. Aspirin's alert is a phenol ester that hydrolyses to a benign phenol;
> a nitrogen mustard's is an alkyl halide that alkylates DNA. Both are alerts.
> Only one is a liability.

This is what the literature says to do: structural alerts are
"disproportionally oversensitive, flagging safe compounds as toxic", with
reported false-positive rates of 6-46 % across endpoints, and they should not
be used alone.

### Authoring panel — a demonstration, not evidence

| control | alerts | class | accepted | expected | correct |
|---|---|---|---|---|---|
| OK-imatinib | 0 | — | accept | pass | yes |
| OK-aspirin | 1 | structural | accept | pass | yes |
| OK-paracetamol | 1 | structural | accept | pass | yes |
| OK-ibuprofen | 0 | — | accept | pass | yes |
| OK-metformin | 2 | structural | accept | pass | yes |
| BAD-curcumin | 2 | reactive | reject | fail | yes |
| BAD-quinone | 2 | reactive | reject | fail | yes |
| BAD-nitrogen-mustard | 1 | reactive | reject | fail | yes |
| BAD-michael-acceptor | 1 | reactive | reject | fail | yes |
| BAD-epoxide | 1 | reactive | reject | fail | yes |

**10/10 correct, balanced accuracy 1.0.**

This number is **not the evidence**. The mechanistic mapping was written with
these molecules in view, so a perfect score here would prove nothing — it is
the same self-confirming structure as WP44's helix test. It is reported
because hiding it would be worse.

### Hold-out panel — the evidence

| molecule | alerts | class | accepted | expected | correct |
|---|---|---|---|---|---|
| OK-caffeine | 0 | — | accept | pass | yes |
| OK-penicillin-G | 0 | — | accept | pass | yes |
| OK-sildenafil | 0 | — | accept | pass | yes |
| OK-glucose | 0 | — | accept | pass | yes |
| BAD-acrylamide | 1 | reactive | reject | fail | yes |
| BAD-iodoacetamide | 2 | reactive, unknown | reject | fail | yes |
| BAD-N-ethylmaleimide | 1 | unknown | accept | fail | **NO** |
| BAD-methyl-mesylate | 1 | unknown | accept | fail | **NO** |

**6/8 correct; sensitivity 0.5, specificity 1.0, balanced accuracy 0.75.**

**Validated: False** — and validation is defined on the hold-out
only. Where the rule is wrong it is named below rather than tuned away.

**Why the two misses happen — it is not the classifier, it is the
vocabulary.**

| molecule | what the alert library says | what the molecule actually is |
|---|---|---|
| BAD-N-ethylmaleimide | `phthalimide` | BRENK reports 'phthalimide' for a maleimide - a different scaffold. The Michael acceptor that makes NEM alkylate cysteines is never named. |
| BAD-methyl-mesylate | `Sulfonic_acid_1` | BRENK reports 'Sulfonic_acid_1'. The compound is a sulfonate ESTER, an alkylating agent; the name says acid, which is not the liability. |

A mechanism-based rule can only be as good as the names it reads. In both
misses the alert name **misdescribes the group that carries the liability**,
so the mechanism cannot be recovered from it. That needs a curated
electrophile library — and it cannot be fitted on this hold-out without
destroying the hold-out. The gap is named and left open.

The rule does **not** validate on the hold-out. It is reported as a
demonstration of the mechanism argument with a measured ceiling, not as a
validated screen — P10 stays **PARTIAL**.

> A graded alert-count penalty was tried first and ALSO failed, for a reason that is a result in itself: every control clears every hard gate, so alert count is the only separating feature, and the panel is non-monotone (metformin must pass with 2 alerts; nitrogen mustard must fail with 1). Exhaustive enumeration of every 'n_alerts <= k' rule caps that family at balanced accuracy 0.7. The rule was therefore changed in KIND, not in threshold: only a REACTIVE (non-specific electrophile) alert rejects; a structural/metabolic alert is carried as a flag. The authoring-panel result (1.0) is a demonstration, not evidence - the mechanistic mapping was written with that panel in view. The hold-out panel result (0.75) is the only evidence.

## The two unmeasurable items — and how to resolve them

### P11 — the experimental success rate

*Why unmeasurable here:* n = 0 designs have been synthesised, so the rate has no denominator.

*Resolution:* the only resolution is to make the designs and count.

| element | value |
|---|---|
| entities | GRAFT-LINEAR, GRAFT-STAPLED, scrambled control |
| replicates | 3 per entity |
| readout | TRAP or direct primer extension, in the presence and absence of POT1-TPP1 |
| acceptance | a competitor is a hit only if it reduces POT1-TPP1-stimulated processivity by >30 % without reducing basal activity >10 % |
| power @ 20 % effect | 93 per arm |
| power @ 10 % effect | 388 per arm |

> Assumption stated, not hidden: 50 % baseline hit rate; these are assumptions, not measurements.

### P12 — binding

*Why unmeasurable here:* no docking engine, no OpenMM, no FEP; WP39 already refuted the in-repo rigid-body scorer, so any score would be fabricated.

*Resolution:* an equilibrium binding measurement.

| element | value |
|---|---|
| primary | MST or SPR against the TERT TEN domain (7TRE chain A TERT is present in the deposit) - label-free Kd |
| orthogonal | ITC for the compounds that give a Kd, enthalpy and stoichiometry |

**The mutant panel is the specificity check**, and it is what makes this assay
worth running rather than merely possible:

- K78E TERT: WP44/WP45 designed around K78, so the design must lose affinity against this mutant. A design that does not is binding elsewhere.
- E215K TPP1: the compensatory mutant of K78E; a design that mimics E215 should be competed by it

> WP45 measured E215 at 4.10 A from K78 NZ and 4.07 A from R132 NH1. A designed competitor must reproduce that dependence, so the mutant panel is the assay's specificity check.

| detection window | value |
|---|---|
| ΔΔG for 1 mM competitor vs 1 nM affinity | -4.09 kcal/mol |
| contact energy the design actually has | 5.41 kcal/mol |

> 1 mM of a competitor with micromolar affinity drives ~1e6-fold weakening, so the effect would be unmistakable IF the competitor reaches micromolar. The risk is not detection, it is the ordering penalty (WP45: 12-20 kcal/mol) leaving it far weaker.

## Summary

- **P6:** COVARIANTLY CONSTRAINED - one feasible staple found — TRP167-LEU211 hydrocarbon 6.34 A; E168/E169/E171/E215 untouched; residual ordering cost **9.2 kcal/mol** still unpaid
- **P10:** validated = **False** (10/10 controls, balanced accuracy 0.75)
- **P11/P12:** converted from unmeasurable to executable — the protocols are above, and neither can be run by computation
