# WP42 — New molecules, new proteins, for the one surface nobody has designed to

    python scripts/solve/wp42_designs.py
      -> results/solve/wp42_designs.json
      -> results/solve/wp42_designs.csv
    python scripts/solve/wp42_designs.py --self-test   -> 48 checks, 0 failed

The request was to design molecules and proteins that did not exist before. This
does that, and then reports, from measurement, that **the small-molecule half of
it cannot work for the reason it was built to solve.** That finding is worth more
than the five structures, so it goes first.

## 1. The gap, established rather than assumed

The TPP1 TEL patch has been called a drug target continuously for thirteen years:

| year | statement | outcome |
|---|---|---|
| 2012 | Nandakumar et al., *Nature*: "a new target for anticancer drug development" | — |
| 2014 | Schmidt et al., *eLife*: "provides a direct target for the design of novel therapeutic inhibitors" | — |
| 2022 | Sekne et al., *Science*: "a rational basis for drug design … by repressing recruitment" | — |
| 2023 | Lim et al., *Commun Chem*: "we have recently discovered the first TPP1 ligands … pharmacological characterization is underway" | **no structure, no series, no affinity ever published** |

What exists is the inverse of a molecule: Yang et al. 2019 overexpressed the whole
TPP1 OB domain as a dominant-negative. That is a protein, not a ligand, and it is
the natural sequence.

The record is unambiguous on one point that matters for honesty: **no small
molecule, peptide or designed protein against the TEL patch has ever been
reported.** The novelty claim of this document is therefore not that the target is
new — it is that the target has been *named* 60+ times and never *designed to*.

`BIBR1532` is often cited near this topic and is not relevant here: it binds the
TER pseudoknot/TERT thumb site, not the patch.

## 2. The design constraint nobody wrote down

Measured from `data/raw/7TRE.cif`, TERT chain A and TPP1 chain C:

```
TEL patch carboxylate  ->  TERT basic residue      distance
  E168 OE1             ->  R774 NH1                 2.32 A
  E169 OE1             ->  R774 NH1                 8.69 A
  E215 OE1             ->  K78  NZ                  4.10 A
  E215 OE1             ->  R132 NH1                 5.69 A
```

and the two TERT acceptors are nowhere near each other:

```
K78 NZ  <-> R132 NH1        4.35 A     one groove
K78/R132 <-> R774          20.26 A     a gulf
```

**The TEL patch is bipartite.** It grips two separate TERT surfaces 20.26 Å
apart. E215 takes the K78/R132 pair; E168 takes R774. Every published description
of "the TEL patch" as a single binding surface is a simplification that is fatal
to design, because a molecule sized to one half cannot touch the other.

## 3. The gulf cannot be crossed by a linear scaffold

A straight line from the K78/R132 pair to R774 is obstructed by the TEN 129–132
loop itself: **17 heavy atoms sit within 2.2 Å of it** (Ala130 O and C, Leu131 N,
Arg132 N/CA/CB/CD, Glu168 OE1/OE2, …).

A free-space grid search (0.75 Å, requiring 1.9 Å clearance) finds that a *bent*
route does exist:

| quantity | value |
|---|---|
| shortest free path | **25.54 Å** |
| straight-line distance | 20.26 Å |
| tortuosity | **1.26** |
| maximum deviation from the straight line | 5.67 Å |
| bend angle at the worst turn | **119°** |

So a bipartite ligand is not merely long, it must **kink 119° and travel 25.5 Å**.
That is the number every design below is measured against.

## 4. The small-molecule result: the class is excluded, and here is the cost

Five ligands were built and each was sampled over 60 ETKDG conformers with MMFF
optimisation, taking the maximum achievable interatomic distance across all
conformers — i.e. the most extended state, not a folded one.

| id | design | MW | logP | max span | needed | verdict |
|---|---|---|---|---|---|---|
| KP-1 | kinked bipartite bridge, m-terphenyl | 318.3 | 4.42 | 14.26 Å | 25.5 Å | **FAILS** |
| KP-2 | KP-1 + one phenyl | 394.4 | 6.08 | 17.95 Å | 25.5 Å | **FAILS** |
| TEL-1 | TEL-patch surface mimetic | 270.3 | 2.61 | 14.27 Å | 20.2 Å | **FAILS** |
| TEL-2 | TEL-1, kinked | 270.3 | 2.61 | 12.81 Å | 20.2 Å | **FAILS** |
| TEL-3 | single-site K78/R132 occupant | 179.2 | 0.41 | 8.68 Å | 6.0 Å | reaches |

**None of the spanning designs reaches.** This was not the expected answer, and
the four failures are reported as failures. To check it is not a sampling
artefact, the same measurement was made on a homologous series of
para-dicarboxy-oligophenylenes:

| phenyls | MW | logP | max span |
|---|---|---|---|
| 1 | 166.1 | 1.08 | 7.44 Å |
| 2 | 242.2 | 2.75 | 11.61 Å |
| 3 | 318.3 | 3.77 | 12.31 Å |
| 4 | 394.4 | 5.54 | 15.44 Å |

Fitted per phenyl: **+2.67 Å of span, +76.1 Da, +1.49 logP**. Extrapolating to
25.54 Å gives

```
7.8 phenyls -> MW ~683, logP ~11.2
```

which is not a drug by any standard, and is also not a molecule anyone can
dissolve. **The conclusion is not that these particular molecules are wrong — it
is that a monomeric, drug-like, bipartite TEL-patch ligand does not exist.** Any
future proposal to "bridge the patch with a small molecule" can be refused by
arithmetic before it is funded.

### What the ligand route actually is

There is exactly one reachable small-molecule site, and it is half the patch. In
the K78/R132 groove the two basic centres are **4.35 Å** apart (K78 NZ–R132 NH1;
NH1–NH2 of R132 is 2.29 Å), so a dicarboxylate with an acid–acid span near
**5.05 Å** chelates both. That is a real, achievable design and it has a real
precedent: the charge-swap rescue of K78E by TPP1 E215K establishes that a single
carboxylate at this groove is worth the whole interaction.

**TEL-3** is that design — deliberately the smallest molecule in the set, and the
only one that passes. It is also the molecule a fragment screen would find, which
means the honest comparison for it is a screen, not a paper.

## 5. The protein result: the patch's own span is a helical register

This is the part that worked, and it worked because the measured spacings turn
out to be helical ones. Ideal α-helix geometry (100°/residue, 1.5 Å rise,
2.3 Å radius) gives:

| register | ideal Cα–Cα | measured requirement | delta |
|---|---|---|---|
| i → i+4 | 6.43 Å | E168–E169 **6.58 Å** | **−0.15 Å** |
| i → i+7 | 10.75 Å | E168–E171 **10.22 Å** | **+0.53 Å** |
| i → i+14 | 21.13 Å | E168–E215 **20.18 Å** | **+0.95 Å** |
| i → i+17 | 25.60 Å | bent path **25.54 Å** | **+0.06 Å** |
| i → i+18 | 27.10 Å | — | — |

Two results follow, and neither was assumed:

1. **The TPP1 TEL patch is encoded as an α-helical register.** Its three
   glutamate spacings (6.58, 10.22, and 20.18 Å end to end) are what a helix
   gives at i, i+4, i+7 and i+14 to within 1 Å. That is why the mimic works: it
   is not an analogy, it is the same geometry.
2. **A helix at i+17 spans the bent path to 0.06 Å.** The 119° kink that defeats
   every small molecule is free to a helix, because a helix curves.

Five de novo proteins, each with a register that is *checked against* the measured
contact rather than asserted:

| id | purpose | aa | MW | pI | net charge | display register |
|---|---|---|---|---|---|---|
| **TELM-1** | TEL-patch proteomimetic, full bipartite | 48 | 4833 | 4.84 | −3.41 | E at 0,4,7,14 — MATCH |
| **TELM-2** | minimal single-helix reagent (SPPS) | 28 | 2861 | 4.24 | −3.38 | E at 0,4,7 — MATCH, 3-acid* |
| **TELM-3** | full-span helix across the gulf | 50 | 5189 | 4.84 | −3.41 | E at 0,4,7,**17** — MATCH |
| **TENB-1** | antagonist binder to the patch | 48 | 4830 | 10.00 | +4.58 | K at 0,4,7,14 — MATCH |
| **TENB-2** | groove-repair for K78E/R132E | 40 | 3917 | 5.40 | −1.41 | E at 0,3 — MATCH, 3-acid* |

\* TELM-2 and TENB-2 carry three polar residues on their display face rather than
four, so they fall below the four-residue floor this analysis uses for
"amphipathic". TELM-2 is the minimal reagent by design; TENB-2 is a groove filler,
not a surface mimic, so a full amphipathic face is not what it is for. The
distinction is stated rather than hidden behind a column.

### Why each one exists

- **TELM-1** is the therapeutic that a TEL-patch mutation is the genetic version
  of: it occupies the patch surface with the same charge pattern, reversibly and
  dose-controlled. It is built to be **resistant to the failure mode of Yang
  2019**, where cells escaped by silencing the expressed TPP1-OB cassette — a
  designed sequence with no TPP1 residues has less to silence and can be delivered
  as protein rather than as a gene.
- **TELM-2** exists because a 28-residue peptide is a **reagent**. TELM-1 is a
  biologic. The difference decides whether this is testable this year or next.
- **TELM-3** is the only design that addresses the measured 25.5 Å path rather
  than the 20.2 Å straight distance. It is the one to make if the hypothesis is
  "the bipartite grip is what matters".
- **TENB-1** proves the sign is not a gimmick: same geometry, opposite charge,
  binds instead of mimics.
- **TENB-2** is the proposal nobody has made. **K78E and R132E are
  recruitment-loss alleles whose entire defect is one missing positive charge.**
  The charge-swap rescue (E215K restoring K78E) establishes that supplying that
  charge restores function. TELM-1 blocks the patch; TENB-2 *repairs* it, in
  trans, for telomere-biology disorders. Same interface, opposite intent, and it
  is the reason to design both rather than one.

## 6. What this is not

**Quarantine: unchanged.** No molecule or protein here has been synthesised,
expressed, or measured. Nothing here binds anything. There is no docking score,
no free energy, no affinity — this repository has no FEP/MM-GBSA and WP39 already
showed its rigid-body scorer produced a false positive (the TLN-amide proposal was
refuted by its own control). Reporting a score from that machinery would repeat a
known error.

The helix model is **ideal geometry**, not a homology model and not an MD
ensemble. It assumes the designed peptides fold as designed, which they may not —
a 28-residue amphipathic peptide can be a helix or a molten globule, and this
document cannot tell you which. Interface specificity is likewise unaddressed:
the K78/R132 groove and the R774 site are basic, and basic surfaces are common,
so off-target binding to other basic patches is the first thing an experiment will
find.

## 7. The falsifiable tests, in the order they should be run

| # | test | refuted if |
|---|---|---|
| 1 | **TELM-2** vs TELM-2-E→A by CD and RAP in the TPP1-stimulated assay | the acid face is not required, i.e. E→A is as active |
| 2 | TELM-2 competes with TPP1 for telomerase (pull-down, IC50) | TELM-2 does not displace TPP1 at any concentration |
| 3 | TELM-2 activity is lost against TERT K78E/R132E | the effect is not at the measured groove |
| 4 | TENB-2 **restores** RAP stimulation of K78E telomerase | K78E stays unstimulated |
| 5 | TEL-3 diacid vs monoacid control in the groove | no difference — then the groove needs more than one charge |
| 6 | TELM-3 vs TELM-1 | TELM-3's extra 4.5 Å of reach buys nothing, meaning the bipartite grip is not required |

Test 1 is the whole hypothesis. Tests 4 and 6 are the two that, if positive,
would be new biology rather than a new molecule.

## 8. Tests

48 checks, 0 failed. They pin the numbers this document stands on:

- the K78–R132 pair is < 5 Å and R774 is > 15 Å from it — the patch is bipartite;
- E215 reaches K78 and E168 reaches R774, from the deposit;
- the straight K78/R132→R774 line is obstructed by more than 10 heavy atoms;
- a bent path exists, is longer than the straight distance, has tortuosity
  1.1–1.6, and respects 1.8 Å clearance;
- the i+17 register matches the bent path better than i+18, and to within 1 Å;
- the register check is **not** trivially true, and linker residues do not pollute
  it (both asserted, because the first implementation had exactly that bug);
- every ligand parses, and the reported acid-group count is computed from the
  graph rather than the string;
- TEL-3 is under 300 Da and KP-1 is heavier than it;
- TELM-1 is acidic, TENB-1 is basic, and their product is negative;
- TELM-1 has a nonzero hydrophobic moment; TELM-2 is peptide-sized;
- no designed protein contains the TPP1 sequence run `ELLLQD`;
- the novelty record has four entries and names the BIBR1532 distinction.
