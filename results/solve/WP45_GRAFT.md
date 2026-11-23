# WP45 — Solving the site WP44 could not, and the mistake that made it look unsolvable

    python scripts/solve/wp45_graft.py
      -> results/solve/wp45_graft.json
      -> results/solve/wp45_graft.csv
    python scripts/solve/wp45_graft.py --self-test   -> 30 checks, 0 failed

WP44 concluded that the TEL triad was "not helix-solvable at all". **That
conclusion was correct and its premise was wrong.** The site was never
unsolvable; the wrong scaffold was being fitted to it.

## 1. The mistake, as a measurement

WP44 modelled the triad as an ideal α-helix. The deposit says the segment is not
one:

| check | value | α-helix expects |
|---|---|---|
| ψ(E169) | **+24.4°** | −60 to −45° |
| CA168–CA171, 7TRE | **7.68 Å** | 5.05 Å (i+3) |
| CA168–CA171, 2I46 | **7.70 Å** | 5.05 Å (i+3) |
| 3-10 helix at i+3 | 6.10 Å | still not 7.68 |

A **positive ψ** breaks the helix, and two independent crystals agree on the
7.68/7.70 Å spacing to 0.02 Å. The segment between E168 and E171 is a **kinked
turn**, and fitting an ideal helix to it produced a lattice with too few degrees
of freedom — which is why 7.31 Å "could not be reached".

### Why that error was possible, and what it cost

The error was a **model substitution**: an idealised helix stood in for the
deposit, and once it was in place every geometry downstream was computed against
the substitute rather than the structure. The self-test could not catch it
because the self-test checked *consistency with the helix model*, not
*consistency with the deposit*. A test that only re-derives its own assumption
will pass forever.

What it cost: WP44 reported a real site as impossible, and reported four
"clean" modules whose registers were fitted to a geometry the protein does not
have.

## 2. The method that works

The literature's answer to exactly this problem is **computation-guided backbone
grafting of a discontinuous motif** (Azoitei et al., *Science* 334:373, 2011) —
keep the functional segments as the crystal shows them and solve for what
connects them. Their Multigraft search is on scaffolds; the same problem is
solved here as geometry, which is exact at this level of modelling.

```
segment A   DWEEKEF   166-172   verbatim from the deposit   (E168, E169, E171)
segment B   FSLLPTE   209-215   verbatim from the deposit   (E215)
linker      must span the junction
```

Measured junction options:

| junction | CA–CA | linker residues needed |
|---|---|---|
| CA171–CA210 | **14.38 Å** | 5 |
| CA172–CA210 | 14.72 Å | 5 |
| CA171–CA209 | 15.96 Å | 5 |
| CA172–CA209 | 16.66 Å | 5 |

Conservative in both directions: **5 residues** using the tightest junction and
the measured 3.5 Å/residue extension.

**The result: 3 of 3 contacts preserved**, because two of them come from segment A
*and* segment B is carried in its crystallographic position relative to them. The
triangle 6.58 / 10.22 / 7.31 Å is reproduced to 0.02 Å — not by fitting it, but by
not throwing the structure away in the first place.

| | WP42 helix | WP44 module | **WP45 graft** |
|---|---|---|---|
| triad | 2 of 3 | 2 of 3 | **3 of 3** |
| K78/R132 groove | no (20 Å, needs a kink) | separate module | **yes — E215 at 4.10 Å** |
| idealised scaffold | yes | yes | **no** |
| payload fraction | — | — | **21 %** |

`DWEEKEF-GGGGG-FSLLPTE` · 19 aa · 2055 Da · charge −4.44.

## 3. What is solved, and what is not

**This is the part that must not be blurred.** WP44's unsolvable problem was
*geometric*, and it is solved. The graft raises a **different** problem — whether
a 19-residue peptide can hold that geometry — and that one is **not** solved.

Both cheap covalent routes were measured, and both fail:

| route | requirement | measured | verdict |
|---|---|---|---|
| disulfide across the segments | Cβ–Cβ 3.5–5.0 Å | closest **6.21 Å** (W167/S210) | **fails** |
| head-to-tail cyclization | CA–CA ≤ ~6 Å | closest **14.28 Å** | **fails** |

So an unconstrained linear graft has nothing holding it. Pricing the alternative
— let the target itself hold it (route a) — settles it:

```
freezing 20 flexible residues     12 – 20 kcal/mol
measured contact energy          -5.4 kcal/mol  (4 salt bridges)
                                 -> cannot pay
```

**An unconstrained linear graft cannot pay for its own ordering.** The peptide
would spend more on entropy than the interface returns.

### Therefore: the graft is a granter, not a molecule

`DWEEKEF-GGGGG-FSLLPTE` is not proposed as a drug. It is the **exact geometric
specification** that any presenting scaffold must satisfy, and it is testable on
its own as a weak competitor. What is missing — and stated as missing — is the
rigid scaffold that holds the two segments 14.4 Å apart at the right angles.

**That scaffold search cannot be run here**: no scaffold structures offline, no
Rosetta, no loop-modelling engine. It needs the Multigraft-class search the
literature uses, which is what this document says to run next.

## 4. The three routes, and what each costs

| route | status here | what it needs |
|---|---|---|
| **(a) target-templated** — linear graft, target holds it | **priced, insufficient** (12–20 vs 5.4 kcal/mol) | nothing new; test it as a weak competitor |
| **(b) covalent constraint** — disulfide or cyclization | **measured, infeasible** (6.21 Å / 14.28 Å) | would need a designed cross-link spanning 6 Å |
| **(c) scaffold grafting** — the literature's method | **specified, not buildable offline** | PDB scaffold search + loop modelling |

Route (c) is the one that solves it, and it is the one that cannot be done in this
sandbox. Saying so is the difference between a plan and a claim.

## 5. Correction record

This is the second error of this class in the programme, and it is recorded rather
than quietly fixed:

- **WP44** asserted the triad was not helix-solvable. The assertion is retained as
  true; the **reason** given for it was a property of the wrong model.
- The helix register tables in WP42/WP44 for the **triad** (i+4 → 6.43 Å, i+7 →
  10.75 Å) should be read as fits to *idealised* geometry, not as descriptions of
  TPP1. The K78/R132 numbers (i+3 → 5.05 Å vs 4.35 Å measured) are unaffected:
  that site is genuinely close to helical.
- Anywhere this repository says the TEL patch is "helical", the correct statement
  is: **one segment of it (167–171) is a kinked turn; the residues are not on a
  single helix.**

## 6. Quarantine

**Unchanged.** Nothing synthesised, nothing measured. The energetic numbers are
screened-Coulomb and entropic estimates with stated constants, not free energies.
No structure was built — the graft's geometry is inherited from the crystal, which
is the one thing here that is not a model.

## 7. Tests

30 checks, 0 failed. The two that matter most are the ones that would have caught
the original error:

- **ψ(E169) is asserted NOT α-helical and positive**, and CA168–CA171 is asserted
  too long for both an α-helix *and* a 3-10 helix;
- **two independent deposits are asserted to agree** on the 7.68/7.70 Å anomaly,
  so the finding cannot be a single-crystal artefact;
- segment A is asserted to equal the deposited sequence residue for residue, and
  to still carry all three glutamates;
- the chosen linker is asserted to span the measured junction, to be uncharged,
  and to use only flexible residues;
- **no covalent constraint is asserted to be available**, the closest disulfide
  pair and the head-to-tail distances are asserted to fail, and the unconstrained
  route is asserted unable to pay its entropic cost;
- the graft is asserted to be minimal — segments plus linker and nothing else.
