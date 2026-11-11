# WP33 — Fixing it for good: four defects, each with the test that now blocks it

WP32 reported three bug fixes, a control test, and a design conclusion. Four of
its claims do not survive scrutiny. Each is retired below with a measurement,
and each has a test that now fails if the defect returns.

    python scripts/solve/wp33_valid_scoring.py
    -> results/solve/wp33_valid_scores.csv

## Defect 1 — the control was confounded by size

| patch (WP32) | residues | atoms | net q | q/atom |
|---|---|---|---|---|
| interface | 188 | **377** | −17.89 | −0.0475 |
| control | 112 | **224** | −1.87 | −0.0083 |

Raw total energies were compared across patches differing by a factor of
**1.68 in atom count**. So "all four ligands prefer the control" measured *size*
as much as chemistry, and was not attributable to charge.

**Fix.** A greedy outward walk from the interface centroid builds a control
matched on **atom count** (−10 %) rather than residue count: **26 interface
residues / 205 atoms vs 27 control residues / 186 atoms, ratio 0.91**, ≥ 25 Å
away. Energies are also reported **per ligand heavy atom**, removing the
trivial "more atoms, more contacts" trend. The achieved ratio is printed, and if
it misses tolerance the run says so instead of hiding it.

## Defect 2 — the prose contradicted the numbers

WP32 wrote that the interface is acidic **and** that this agrees with "the TEL
patch is basic". Those are opposite statements about the same surface.

The measurement decides:

| TEL-patch residue in the interface | formal charge |
|---|---|
| ARG92 | +1 |
| ASP166 | −1 |
| GLU168 | −1 |
| GLU169 | −1 |
| LYS170 | +1 |
| ARG180 | +1 |
| GLU215 | −1 |
| **net** | **−1** |

The TEL patch is **ambipolar with a net of −1**, not uniformly basic. And the
−17.9 is not the patch — it is the whole 47-residue contact surface, most of
whose negative charge sits on TERT (ASP129, GLU800) and on TPP1's far face
(GLU215).

**Withdrawn:** "the interface is acidic, so the anchor should be basic." That
was an overstatement built on a total that was never the patch's own charge.

## Defect 3 — five of six "variants" were different molecules

WP30's TLN-2/3/4 each dropped an aryl ring to shed the core carboxylate,
changing scaffold, size and aromatic surface **at once**. Only TLN-5 was a clean
isosteric swap, so the six-way comparison was never evidence about anchor
charge.

**Fix.** A series where only the proton — and therefore the charge — differs:

| variant | formula | heavy | rings | q |
|---|---|---|---|---|
| acid / carboxylate :: neutral | C23H23FN2O4S | 31 | 3 | 0 |
| acid / carboxylate :: charged | C23H22FN2O4S− | 31 | 3 | −1 |
| amine / ammonium :: neutral | C22H24FN3O2S | 29 | 3 | 0 |
| amine / ammonium :: charged | C22H25FN3O2S+ | 29 | 3 | +1 |
| guanidine / guanidinium :: neutral | C23H26FN5O2S | 32 | 3 | 0 |
| guanidine / guanidinium :: charged | C23H27FN5O2S+ | 32 | 3 | +1 |
| amide (neutral control) :: neutral | C23H24FN3O3S | 31 | 3 | 0 |

**Within each pair the heavy-atom count is identical** — asserted by the test.
A five-member "same heavy atoms" series is *not* achievable (a carboxylate and a
guanidinium differ by nature), and claiming one would have repeated the defect.

**A second bug on the way:** the first version of this pass scored everything
through the pH 7.4 ioniser, so both rows of every pair collapsed to the same
charged species and the comparison was **empty** — visibly, as identical
numbers. The neutral rows are now scored as **explicit protonation states**
(`ionise=False`), a controlled perturbation, not a claim about what exists in a
cell.

## Defect 4 — the hydrogen-bond test held the wrong thing fixed

WP32 compared a 3.2 Å collinear contact against a 2.83 Å bent one: distance
*and* angle both changed, and 2.83 Å sits near the 2.6 Å lower bound.

**Fix.** Hold **D···A at exactly 3.0 Å** and vary only the angle, with a probe
whose realised angle is verified:

| constructed | realised D–H···A | D···A | counted |
|---|---|---|---|
| 160° | 160.0° | 3.00 Å | **yes** |
| 150° | 150.0° | 3.00 Å | **yes** |
| 130° | 130.0° | 3.00 Å | **yes** |
| 125° | 125.0° | 3.00 Å | **yes** |
| 121° | 121.0° | 3.00 Å | **yes** |
| 119° | 119.0° | 3.00 Å | **no** |
| 115° | 115.0° | 3.00 Å | **no** |
| 100° | 100.0° | 3.00 Å | **no** |

The 120° boundary holds from both sides at constant distance. Writing this test
also exposed a bug in *the test*: the naive construction (rotating the acceptor
about the donor) requested 125° and actually produced **106°**, so the first
version was checking an angle it was not creating. The probe now solves
`L² − 2L·cos β + 1 = d²` and the realised angle is asserted.

## What the corrected comparison actually says

| variant | E/atom (interface) | E/atom (control) | Δ | prefers |
|---|---|---|---|---|
| **amide (neutral)** | **−0.2965** | −0.0977 | **−0.1988** | **INTERFACE** |
| amine :: charged (+1) | −0.1123 | −0.1204 | +0.0081 | control |
| acid :: neutral (0) | −0.0563 | −0.1298 | +0.0735 | control |
| guanidinium :: charged (+1) | −0.0295 | −0.1633 | +0.1338 | control |
| guanidinium :: neutral (0) | +0.0272 | −0.1218 | +0.1490 | control |
| amine :: neutral (0) | −0.0038 | −0.1560 | +0.1522 | control |
| acid :: charged (−1) | +0.0029 | −0.2085 | +0.2114 | control |

**The neutral amide — no charge at all — is the only variant that prefers the
interface, and it wins by a wide margin.**

Charge does not order the series monotonically:

- adding a proton to the acid (0 → −1): **Δ +0.138, worse**
- adding a proton to the amine (0 → +1): **Δ −0.144, better**

That is the fingerprint of an **ambipolar** patch. A uniformly acidic patch
would punish both; a uniformly basic one would reward both. The interface holds
both ARG92/LYS170/ARG180 and ASP166/GLU168/GLU169/GLU215, so the correct
prediction is that charge sign matters **locally** and no single global charge is
right — which is what the data show. And it is why "acid *or* base" could never
have found the amide: that framing assumed one sign for the whole surface.

**Caveat that must travel with this.** The amide result is one rigid placement
protocol, one dielectric, no solvent, no entropy, no flexibility. A neutral
anchor avoiding both charge mismatches is *chemically coherent*, and coherence
is not evidence. It is a hypothesis with a mechanism behind it, which is more
than the previous variants had.

## Scoreboard

| | state |
|---|---|
| control matched on atom count | **fixed** — 0.91 ratio, was 0.59 |
| energies size-normalised | **fixed** — per ligand heavy atom |
| TEL patch charge characterised | **fixed** — ambipolar, net −1 |
| "anchor must be basic" | **WITHDRAWN** |
| variant series isolates charge | **fixed** — matched protonation pairs |
| H-bond angle test | **fixed** — constant distance, verified angle |
| a variant prefers the interface | **yes — the neutral amide** |
| affinity | **none, and none obtainable here** |

**Quarantine unchanged.** No Kd, no IC50, no free energy. A relative per-atom
energy is not a binding measurement, and the ranking above is a ranking of a
surrogate.

## Tests

36 checks, 0 failed. The ones that matter: every matched pair has an identical
heavy-atom count; the pH 7.4 pass really ionises the neutral acid to −1; the
control is matched within 25 % on atoms, shares no residue with the interface,
and is ≥ 25 Å away; the TEL patch is ambipolar with a net near −1; the H-bond
probe **realises the angle it claims** and honours the 120° boundary at constant
distance; the WC32 clash fix is still repulsive with electrostatics off at
overlap; and the acid still ionises to −1.
