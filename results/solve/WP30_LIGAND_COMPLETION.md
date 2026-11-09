# WP30 — TELIN-1's missing part: located, rebuilt, and tested

The request was to find the part a molecule is missing and add it from a
molecule that works. TELIN-1 turns out to be a textbook case: **it is not one
molecule.** Everything below is measured from the file, not recalled.

    python scripts/solve/wp30_ligand_completion.py
    -> results/solve/wp30_ligand_completion.csv
    -> data/ligands/telin1_repaired.pdb  +  .sdf   (a coordinate repair, not a claim)

## 1. The defect, measured from the file

| property | value |
|---|---|
| heavy atoms | 36 (declared `C26 F N2 O6 S1` — matches) |
| connected components | **2** — sizes **29** and **7** |
| bonds between them | **none** |
| closest approach between fragments | **2.744 Å** |
| shortest CONECT bond | **1.208 Å**, in a six-membered ring |
| longest CONECT bond | **1.715 Å** (C–C) |
| ring angles | **158.9°** and **95.8°** (sp³: 109.5°, sp²: 120°) |

A CONECT list is an assertion of bonding. This one asserts nothing across the
gap, so the file describes two molecules 2.744 Å apart and calls them one.

## 2. The core is not merely broken — it is lens-shaped

All six ring bonds are present, and they cycle twice through one extra atom:

| bond | length | normal |
|---|---|---|
| C1–N1 | 1.393 Å | C–N 1.33–1.47 |
| N1–C2 | 1.393 Å | C–N 1.33–1.47 |
| C2–N2 | 1.393 Å | C–N 1.33–1.47 |
| N2–C3 | 1.575 Å | C–N 1.33–1.47 |
| C3–C4 | 1.463 Å | C–C 1.50–1.55 |
| **C4–C1** | **1.208 Å** | C–C 1.50–1.55 |

Three C–N bonds in a 3:3 alternating pattern, all at exactly 1.393 Å, cannot be
right — the sequence charges the ring twice over. And one ring C–C bond is
1.208 Å, *shorter than the C=O two atoms away*, in a ring that is **flat**
(all torsions < 8°) with interior angles of **158.9°** and **95.8°**. That is a
lens, not a hexagon.

The distortion is **local**: the exocyclic `=S` is 1.847 Å (normal 1.60–1.70)
and the C=O is 1.221 Å (normal), so this is not a scale or unit error.

## 3. Why the arm cannot simply be attached where it sits

The 7-atom piece is pentanoic acid (`C5H10O2`). Its nearest atom across the gap
is its **terminal methyl carbon**. The 29-atom piece's nearest atom is a
**terminal thiocarbonyl sulfur** (`C=S`). A methyl carbon bonded to a
thiocarbonyl sulfur is not a chemistry that exists.

And there is no place to move it to. The core scaffold is
**hexahydropyrimidine-2-thione-4,6-dione**: C3 carries `=S`, N1 and N2, and each
nitrogen already bears an aryl ring. Every valence on that scaffold is used.

> **So the defect is not a bad coordinate. It is a bad molecule.** The arm was
> drawn as if it had somewhere to attach, and it does not.

## 4. Where the arm attaches — the declared formula decides

The header's own composition is the clue, and it is exact:

```
core     C21 H19 F N2 O4 S      MW 414.46    29 heavy atoms
+ arm    C5  H9  O2             (loses one H at the junction)
=        C26 H27 F N2 O6 S      MW 514.58    36 heavy atoms   <- the file's own
```

Reproduced to the last atom. **The intended molecule is "core + one C5 carboxylic
arm, bonded through an aryl C–H"** — the two aryl rings are the only positions
with a free valence.

| attachment | formula | MW | heavy |
|---|---|---|---|
| *para* to N | C26H27FN2O6S | 514.58 | 36 |
| *meta* to N | C26H27FN2O6S | 514.58 | 36 |
| *ortho* to N | C26H27FN2O6S | 514.58 | 36 |

All three match. A composition cannot distinguish them, so **the positional
ambiguity is recorded as open, not guessed.** What is settled is the
connectivity: *nothing needs to be added except one bond.*

## 5. The part worth borrowing — and it is a subtraction

TELIN-1's claim is protein–protein interface disruption through a carboxylate.
The one **approved** molecule built to do exactly that is **lifitegrast**
(FDA 2016, LFA-1/ICAM-1 antagonist, [1](https://drugs.ncats.io/drug/038E5L962W)). Its design rules:

- **one** carboxylate in the binding epitope, not several;
- it sits on a **short flexible linker** off an aromatic core;
- an **electron-poor aryl** and an amide complete the epitope.

TELIN-1 presents **two carboxylates plus a thione** with no linker control
between them. The borrowed part is therefore *"keep exactly one anchor
carboxylate on the arm, and delete the competing one"*. That is a subtraction,
and it is the whole lesson.

## 6. Tested — and the repair does **not** produce a good drug

| molecule | MW | logP | HBD | HBA | TPSA | COOH | Lipinski |
|---|---|---|---|---|---|---|---|
| TELIN-1 core only | 414.46 | 4.28 | 1 | 4 | 77.92 | 1 | pass |
| **TLN-1 repaired (2 COOH)** | **514.58** | **5.08** | 2 | 5 | 115.22 | 2 | **2 fail** |
| TLN-2 one anchor | 470.57 | 5.38 | 1 | 4 | 77.92 | 1 | 1 fail |
| **TLN-3 anchor + short arm** | **442.51** | **4.60** | 1 | 4 | 77.92 | 1 | **pass** |
| TLN-4 no alkyl tail | 414.46 | 3.96 | 1 | 4 | 77.92 | 1 | pass |
| lifitegrast (approved) | 615.49 | 4.77 | 2 | 6 | 133.99 | 1 | **1 fail** |

Two results, both worth stating plainly:

1. **The straight formula repair yields a chemically valid but non-drug-like
   molecule** — MW 514.58 > 500 and logP 5.08 > 5. That is a finding about the
   original design, not a failure of the repair.
2. **Applying the borrowed rule fixes both flags.** TLN-3 (one anchor, C3 arm)
   is MW 442.51, logP 4.60, and passes all four rules.

**But do not over-read that table.** Lifitegrast — an approved drug — itself
fails `MW > 500`. Lipinski is a permeability heuristic for oral small molecules;
it is not evidence of binding, and passing it does not make TLN-3 a candidate.
TLN-1 (514 Da) and lifitegrast (615 Da) are both outside the rule and both were
designed for a protein surface rather than a deep pocket, which is exactly where
the rule is known to be weakest.

## 7. Tested — is the rebuilt 3D structure sane?

| | original | repaired |
|---|---|---|
| connected components | 2 (29 + 7) | **1** |
| atoms | 36 | 36 |
| bond-length range | 1.208 – 1.847 Å | **1.217 – 1.676 Å** |
| six-membered ring angles | 158.9° present | **none** |

The conformer is an **embedded MMFF minimum**, not a binding pose. It exists so
downstream code has one connected molecule instead of two.

## 8. The bigger missing part is not in the ligand at all

The complex file docks TELIN-1 into a structure containing **TERT, TERC and DNA
— and no TPP1.** Chain names present: `A`, `B`, `D`, `L`. There is no TPP1
chain.

> **A molecule cannot disrupt a TERT–TPP1 interface that is not in the model.**
> The docked pose was never testing the claimed mechanism — it cannot have been.

That is the real missing part, and it is not something a ligand edit can fix.
The design is missing a **defined epitope**. Recovering it needs the TERT–TPP1
structure, and RCSB is TLS-blocked from this sandbox (`files.rcsb.org`,
`models.rcsb.org` and the PDBe mirror all fail to negotiate SSL; `github.com`
returns 200 in the same breath, so the block is specific, not general).

## 9. Scoreboard

| | state |
|---|---|
| one connected covalent species | **fixed** |
| composition matches the declared header | **fixed** |
| 1.208 Å ring bond, 1.715 Å C–C, 158.9°/95.8° angles | **fixed** in the rebuild |
| which aryl position carries the arm | **open** — three candidates tie |
| docked pose | **open** — 160 clashes < 2.8 Å, ligand C19 0.17 Å from a backbone N |
| defined TPP1 epitope | **fixed in WP31** — TPP1 was on disk (7TRE chain C, 2I46 chain A); this row was wrong |
| any affinity (Kd, IC50) | **none, and not obtainable here** |

**Quarantine unchanged.** TELIN-1 remains an unvalidated in-silico proposal with
zero empirical evidence. A repaired formula is not an activity claim, and a
passing Lipinski filter is not a result.

## 10. Tests

25 checks, 0 failed — including: the original really is two fragments
(29 + 7), the gap really is ~2.7 Å, the impossible bonds really exist, the core
really is 29 heavy atoms, **all three repairs reproduce `C26H27FN2O6S` exactly**,
each repair is a single component, the repair adds exactly one carboxylate,
the written PDB is **one** component of 36 atoms with no bond below 1.20 Å, and
lifitegrast really has exactly one carboxylic acid while repaired TELIN-1 has
two.
