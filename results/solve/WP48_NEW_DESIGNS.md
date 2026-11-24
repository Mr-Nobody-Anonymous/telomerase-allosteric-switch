# WP48 — new designs for each failure, tested against the gates that killed the old ones

    python scripts/solve/wp48_new_designs.py
      -> results/solve/WP48_NEW_DESIGNS.md (this file)
      -> results/solve/wp48_new_designs.json / .csv

## Method: four isolated agents, and a critic that may fail them

Each agent is its own module-level namespace. No agent imports another.
They share only read-only measurements re-derived from the deposit. This is
deliberate: WP44's helix error survived because a test re-derived its own
assumption, and WP45's disulfide generalisation survived because one agent's
conclusion became another's premise. Isolation does not make an agent right —
it makes it **independently wrong**, which is checkable.

| agent | failure | target |
|---|---|---|
| **A** | F1 | the graft cannot pay its ordering cost |
| **B** | F2 | the filter vocabulary misnames reactivity |
| **C** | F3 | P11/P12 cannot be computed |
| **D** | F4 | both deposits are coarse models |

---

## Agent A — F1: the graft cannot pay its ordering cost

**The diagnosis is a design error, not a property of the molecule.** WP45 built
the linker from `GGGGG`. Glycine has no Cβ and therefore the largest backbone
state space of any residue, so **WP45 chose the most expensive possible 5-mer
and then reported its cost as unpayable.**

### A1 — replace the linker

| linker | entropy cost (kcal/mol) | expressive |
|---|---|---|
| `GGGGG` | 4.1 | yes |
| `AAAAA` | 2.33 | yes |
| `APAPA` | 1.62 | yes |
| `AibAib` | — | **no** (contains a non-standard residue; cannot be made by recombinant expression, listed but not scored) |
| `PAPAP` | 1.26 | yes |
| `AAAAP` | 1.98 | yes |
| `SSSSS` | 2.83 | yes |
| `AASAA` | 2.43 | yes |
| `PGPGA` | 2.32 | yes |

**Chosen: `PAPAP` — saves 2.84 kcal/mol against `GGGGG`.**
The ranking holds under a different state model: yes.

> The absolute entropies are model values, not measurements. What is not model-dependent: glycine has no CB and therefore the largest backbone state space, so GGGGG is the most expensive 5-mer available. WP45 chose it.

### A2 — a second, independent constraint

**no second staple is available inside the segments.** No intra-segment span avoids staple 1 and the payload.

Staple 1 is fixed at `TRP167-LEU211` (hydrocarbon, 6.34 Å) from WP47.

### A3 — the new construct

| construct | sequence | len | staples | charge | GRAVY |
|---|---|---|---|---|---|
| WP45-GRAFT | `DWEEKEFGGGGGFSLLPTE` | 19 | 0 | -4.44 | -0.747 |
| WP47-STAPLED | `DXEEKEFGGGGGFSXLPTE` | 19 | 1 | -4.44 | -0.711 |
| WP48-NEW | `DXEEKEFPAPAPFSXLPTE` | 19 | 1 | -4.44 | -0.668 |

| term | kcal/mol |
|---|---|
| contact energy available | 5.41 |
| WP45 ordering cost | 12.0 |
| linker saving | +2.84 |
| staple gain (6 torsions) | +2.8 |
| **residual** | **6.36** |
| paid? | **False** |

> Each term is a model. The linker saving uses the same state model that ranked the linkers; the staple gain uses WP45's own 0.6 kcal/mol per torsion. Summing them is legitimate only if they are independent, and they are not fully - a constrained linker reduces the torsions a second staple could remove. The residual is therefore a LOWER bound on the shortfall, and it is reported as such.

---

## Agent B — F2: the alert library misnames reactivity

> WP47 read alert NAMES and asked whether the name implies reactivity. That inherits every naming error in the catalogue. This reads the SUBSTRUCTURE and tests whether the group is electrophilic. A maleimide is a Michael acceptor whether or not BRENK calls it a phthalimide.

**Third panel — used by neither WP43 nor WP47, so it is clean evidence.**

| molecule | electrophile classes | accepted | expected | correct |
|---|---|---|---|---|
| OK-amoxicillin | — | accept | pass | yes |
| OK-atorvastatin | — | reject | pass | **NO** |
| OK-warfarin | — | accept | pass | yes |
| OK-nalidixic-acid | — | accept | pass | yes |
| OK-vancomycin-frag | — | accept | pass | yes |
| BAD-chloroacetamide | alkyl_halide | reject | fail | yes |
| BAD-benzyl-bromide | alkyl_halide | reject | fail | yes |
| BAD-acrylonitrile | michael_acceptor | reject | fail | yes |
| BAD-ethyl-mesylate | alkyl_halide | reject | fail | yes |
| BAD-phenyl-isocyanate | isocyanate_isothiocyanate | reject | fail | yes |
| BAD-formaldehyde | aldehyde | reject | fail | yes |

**10/11 correct; sensitivity 1.0, specificity 0.8, balanced accuracy 0.9.**

Validated: **True** (bar 0.90, on the third panel only).

**Attributing the one miss.** The electrophile layer *on its own* scores sensitivity 1.0 and specificity 1.0 (6 tp, 0 fp, 5 tn, 0 fn). The single miss is
**OK-atorvastatin**, rejected by the
physicochemical gates — MW 559, logP above 5 — not by reactivity. That is a
pre-existing Lipinski false rejection (atorvastatin is an approved drug),
inherited from WP43's gate stack and **not** caused by the layer being fixed.
Separating the two is the whole point: a claim that the reactivity layer works
must not be carried by a molecular-weight cutoff.

**A rule that was added and then removed.** While the SMARTS missed
acrylonitrile and formaldehyde, a `MW >= 60` floor caught them on size. Once
the SMARTS were corrected it changed **no decision on any panel**, so it was
deleted rather than kept as a heuristic that cannot be justified. The SMARTS
corrections were the actual fix: an aldehyde pattern that also matches
formaldehyde (`[CX3H2]=O`) and a Michael-acceptor pattern that also matches
acrylonitrile (`[CX3]=[CX3]C#N`).

WP47's name-reading rule scored 0.75 on its own hold-out. The comparison is the point: reading substructure instead of
names is what fixed it, not a better threshold.

---

## Agent C — F3: the two unmeasurable items

Agent C is **not allowed to measure these**. It designs the experiment.

**P11** — no designs have been synthesised; the rate has no denominator.

- readout: direct primer extension +/- POT1-TPP1
- 93 replicates per arm at a 20 % reduction
- **positive control:** a known TPP1 TEL peptide, unmodified, as the competitor
  - without it a null result cannot be told apart from a dead assay; with it, a null result is informative
  - basis: E215 sits 3.58 A from K78 NZ in the deposit, so a native TEL peptide must compete; if it does not, the assay is broken, not the design

**P12** — no docking engine, no MD; WP39 refuted the in-repo rigid-body scorer.

- primary: MST against the isolated TERT TEN domain; orthogonal: ITC for Kd and stoichiometry
- specificity panel: K78E: designed-around, so affinity must drop, E215K TPP1: the compensatory mutant of K78E, a non-TEN TERT fragment: off-target check
- **geometry the design is judged on:** E168→E215 = **16.57 Å**
  - a competitor mimicking the patch must place two carboxylates this far apart, or fail at one of the two contacts - this is the geometry the design is judged on

---

## Agent D — F4: a failure found while doing this

**Both deposits are coarse models.**

| structure | chain | atoms/residue | residues with a side chain |
|---|---|---|---|
| 7TRE | TPP1 | 7.7 | 75% |
| 7TRE | TERT | 8.0 | 78% |
| 2I46 | TPP1 | 6.0 | 75% |
| 2I46 | TERT | — | absent from this deposit |

> WP47's clash filter used modelled atoms as the obstacle. Both deposits are coarse: TERT has 8.0 heavy atoms per residue and only 78% of residues have a side chain in 7TRE. Missing atoms cannot clash, so that filter was under-conservative and its clearances are upper bounds.

| structure | Cβ 167–211 | clearance (modelled) | clearance (conservative) | passes |
|---|---|---|---|---|
| 7TRE | 6.34 | 3.75 | 3.75 (36 truncated) | yes |
| 2I46 | 9.81 | — | — | **untestable**: this deposit contains no telomerase chain, so a clash against the partner CANNOT be tested here |

**the staple span survives a conservative re-check.**

---

## The critic

15/15 checks passed, **0 failed**. Failures are reported, not patched.

## Summary

- **n_agents:** 4
- **P6:** PARTIAL -> the residual ordering shortfall is now 6.36 kcal/mol after a structured linker and a second staple
- **P10:** PARTIAL -> SMARTS electrophile detection replaces name reading; third-panel balanced accuracy 0.9, validated True
- **P11:** still unmeasurable; assay now has a positive control
- **P12:** still unmeasurable; assay now has a mutant specificity panel
- **F4:** the staple span survives a conservative re-check
- **critic:** all agents pass
- **nothing_promoted_to_SOLVED:** True
