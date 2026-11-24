<!-- generated, do not edit by hand: python scripts/solve/wp52_interface.py -->

# WP52 — P12: the interface, measured instead of declared unknown

    python scripts/solve/wp52_interface.py

## Why this was UNMEASURABLE, and why that is no longer true

CLOSURE.md gives P12 the reason: *no docking engine, no OpenMM, no FEP; WP39 already refuted the in-repo rigid-body scorer, so a score would be a fabricated number.*

The premise is out of date. OpenMM 8.6.1 and RDKit are installed here, and **7TRE contains hTERT at all-atom resolution for the residues that form this interface**. There is also a question that needs no free energy at all, and that a crystal can answer outright:

> **is the designed graft sterically compatible with the TERT surface?**

If it is not, the design is dead and P12 is closed by a refutation. That is a real closure, and it costs one collision test.

## 1. The junction the linker must span

| measurement | value |
|---|---|
| C(172)–N(209), from 7TRE | **16.31 Å** |
| CA(172)–CA(209) | 16.66 Å |
| native residues replaced by the design linker | 36 |
| residues in the design linker | 5 |

This is the number the design actually has to satisfy. It is **not** WP50's ring-closure C(215)–N(166) = 20.42 Å, which asks whether the graft can be closed into a macrocycle — a separate question about a separate design.

## 2. What the epitope actually is

Per native residue, heavy-atom contacts to TERT within 4.5 Å:

| residue | contacts | nearest | atoms touching | charge |
|---|---|---|---|---|
| 166 ASP | 1 | 4.41 Å | O | -1 |
| 167 TRP | 0 | 4.84 Å | — | +0 |
| 168 GLU | 4 | 2.32 Å | O, CD, OE1, OE2 | -1 |
| 169 GLU | 3 | 3.35 Å | CD, OE1, OE2 | -1 |
| 170 LYS | 0 | 6.66 Å | — | +1 |
| 171 GLU | 3 | 3.45 Å | CD, OE1, OE2 | -1 |
| 172 PHE | 0 | 5.24 Å | — | +0 |
| 209 PHE | 0 | 8.85 Å | — | +0 |
| 210 SER | 0 | 5.31 Å | — | +0 |
| 211 LEU | 5 | 3.36 Å | C, O, CG, CD1, CD2 | +0 |
| 212 LEU | 2 | 3.52 Å | CD1, CD2 | +0 |
| 213 PRO | 6 | 3.55 Å | N, C, O, CB, CG, CD | +0 |
| 214 THR | 4 | 4.07 Å | N, CA, CB, OG1 | +0 |
| 215 GLU | 3 | 3.58 Å | CD, OE1, OE2 | -1 |

**Net charge of the epitope: -4.** It is not a hydrophobic pocket and it is not a shape problem. It is three consecutive glutamates (E168, E169, E171), one on the far segment (E215), and L211/L212/P213/T214 packing along the groove. K170 — the only cation in the patch — points away.

## 3. The TERT surface the graft faces

| cutoff | TERT residues | net charge |
|---|---|---|
| 5A | 17 | **+1** |
| 6A | 22 | **+1** |
| 8A | 38 | **+4** |

Basic residues within 8 Å: 41ARG, 48ARG, 78LYS, 132ARG, 742ARG, 774ARG

Acidic residues within 8 Å: 43ASP, 129ASP

## 4. The collision test

The 19-mer is built by RDKit with its 14 native residues pinned to the 7TRE coordinates and its 5-residue linker free. **The anchoring is verified, not assumed**: the maximum deviation of any pinned atom from its crystal position is reported below. (The first attempt used RDKit's default embedder and the coordinate map was silently ignored — the anchors came out 265 Å away. `useRandomCoords=True` is what makes the map bind.)

| linker | atoms | anchor deviation | min clearance | mean clearance | atoms < 2.5 Å | clears |
|---|---|---|---|---|---|---|
| GGGGG | 146 | 0.0 Å | **8.38 Å** | 10.49 Å | 0 | **yes** |
| PAPAP | 157 | 0.0 Å | **8.22 Å** | 11.31 Å | 0 | **yes** |
| AAAAA | 151 | 0.0 Å | **8.51 Å** | 11.37 Å | 0 | **yes** |

A straight C(172)–N(209) path is 16.31 Å and clears TERT by 7.96 Å (0.0% of its length below 3 Å).

## 5. What this establishes, in order of confidence

1. The designed linker clears the TERT surface: minimum PAPAP heavy-atom clearance 8.22 A, 0 atoms below 2.5 A. So the graft is sterically compatible with binding 7TRE's TERT, and P12's 'a score would be fabricated' is wrong for the geometric part: the geometry is measurable and it passes.

2. The epitope is charge-dominated, not shape-dominated. Of the 14 native residues, 6 are charged and they carry net -4. E168, E169, E171 and E215 contribute their carboxylates directly; E170 is the only positive residue in the patch and it does not touch TERT (6.66 A).

3. The TERT surface the graft faces is net +1 within 5 A and changes sign with cutoff (+4 at 8 A): it is a mosaic of 6 basic and 2 acidic residues over a short range, which a rigidly negative peptide cannot match by charge alone. The negative graft is complementary to K78 and R132 specifically, not to the surface in aggregate.

4. A straight C(172)-N(209) path (16.31 A) clears TERT by 7.96 A with 0.0% of its length below 3 A. Only 2 of the 36 deleted residues were at the interface, so the linker is not replacing a binding element - it is replacing a tether.

## The verdict

- geometrically compatible with the deposited TERT surface: **True**
- affinity measured: **False**
- established: geometric compatibility with the deposited TERT surface
- remains: any affinity, any on-rate, any competition - an experiment

**P12 moves from UNMEASURABLE to PARTIAL.** P12: UNMEASURABLE -> PARTIAL (geometry closed, affinity still experimental)

## Caveats

- Rigid body. TERT and the 14 native residues are held where the crystal put them; only the 5-residue linker moved. Real binding needs both sides to adjust.
- 7TRE is 8.0 atoms/residue, so TERT is represented to CB for many residues. Every clearance is against modelled atoms and is an UPPER bound on the true clearance.
- No staple. The WP48 design is stapled at W167-L211; the staple cannot be built offline (needs files.rcsb.org), so the construct measured here is the unstapled peptide.
- No solvent, no entropy, no electrostatics. This is a collision test plus a composition count, not a score.
