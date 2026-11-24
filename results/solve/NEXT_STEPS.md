<!-- generated, do not edit by hand: python scripts/solve/wp49_master_record.py -->

# Next steps — what to order, what to run, what would change the answer

    python scripts/solve/wp49_master_record.py

Nothing here has been made. This is the materials and protocol list, written
so the work can leave the computer.

## 1. Peptides to order

All solid-phase synthesis. Crosslinked residues are α-methyl, α-alkenyl
non-natural amino acids (Fmoc-S5 / Fmoc-R8 series), incorporated by standard
SPPS and closed by ring-closing metathesis on resin.

| id | sequence | length | staples | purpose |
|---|---|---|---|---|
| GRAFT-LINEAR | `DWEEKEFGGGGGFSLLPTE` | 19 | 0 | the unstapled control |
| WP47-STAPLED | `DXEEKEFGGGGGFSXLPTE` | 19 | 1 | W167-L211 staple, glycine linker (WP45's, kept as control) |
| WP48-NEW | `DXEEKEFPAPAPFSXLPTE` | 19 | 1 | W167-L211 staple + structured PAPAP linker |
| SCRAMBLE | same composition, order shuffled | 19 | 0 | specificity control |
| NATIVE-TEL | TPP1 166–215 as expressed | 50 | 0 | **positive control** |

The positive control is not optional. Without it a null result cannot be
distinguished from a dead assay.

## 2. Assays, in order

### 2.1 Binding (P12)

- **primary:** MST or SPR against the TERT TEN domain (7TRE chain A TERT is present in the deposit) - label-free Kd
- **orthogonal:** ITC for the compounds that give a Kd, enthalpy and stoichiometry
- **specificity panel:**
  - K78E TERT: WP44/WP45 designed around K78, so the design must lose affinity against this mutant. A design that does not is binding elsewhere.
  - E215K TPP1: the compensatory mutant of K78E; a design that mimics E215 should be competed by it

### 2.2 Function (P11)

- readout: TRAP or direct primer extension, in the presence and absence of POT1-TPP1
- n = 93 per arm at a 20 % effect (388 at 10 %)
- acceptance: a competitor is a hit only if it reduces POT1-TPP1-stimulated processivity by >30 % without reducing basal activity >10 %

## 3. The decision rule, declared before the data

| outcome | interpretation | next move |
|---|---|---|
| Kd measurable, K78E abolishes it | site-specific binding | proceed to function |
| Kd measurable, K78E does not matter | binding elsewhere | discard the design, not the target |
| no Kd at any concentration | ordering cost wins | needs the scaffold (below) |
| function changes without binding | artefact | reassess the assay |

## 4. The one thing computation cannot do

**Find the presenting scaffold.** 6.36 kcal/mol of ordering remains unpaid even
after the linker redesign and the staple. Closing it needs a scaffold that
pre-organises the two segments, and finding one needs:

1. PDB-wide discontinuous-motif search (`Multigraft Match`, Azoitei 2011)
2. loop modelling of the connecting segments (MODELLER or Rosetta
   `LoopModel`/`NextGenerationKIC`)
3. scoring with Rosetta `FastRelax` + DOPE, or a comparable energy function

None of the three is available in this sandbox. That is a hard boundary, not
an unfinished task.

## 5. Cheap things to do next, in this sandbox

| task | why it is worth it | cost |
|---|---|---|
| fetch a higher-resolution TERT–TPP1 model if deposited | the current deposits are 8 atoms/residue | one fetch |
| test the graft against 7TRC's fuller TPP1 | an independent copy of the same segment | one measurement |
| extend the fold-back linker candidates | PAPAP was the best of nine, not of all | search |
| curate an electrophile library | the remaining filter misses are vocabulary, not logic | curation |
