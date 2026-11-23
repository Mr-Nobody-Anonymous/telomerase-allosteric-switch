# WP44 — One module per sub-problem, merged, and the merge stress-tested

    python scripts/solve/wp44_modular_designs.py
      -> results/solve/wp44_modular_designs.json
      -> results/solve/wp44_modular_designs.csv
    python scripts/solve/wp44_modular_designs.py --self-test   -> 39 checks, 0 failed

WP42 built multifunctional helices: one 48-residue chain that had to span 20 Å,
read a carboxylate triad *and* present a hydrophobic face. That is fragile by
construction. The problem is split here into the four sites that were actually
measured, **one module per site**, each held to two rules — and then merged, with
the merge tested for stability.

The rules are machine-checked, not aspirational:

> **FULL** — a module must cover *every* measured contact of its own site.
> **ONLY** — a module must register on *no other* site.

## 1. The four sites, measured from 7TRE

| site | anchors | separation | kind |
|---|---|---|---|
| **S1** K78/R132 groove | K78 NZ, R132 NH1 | **4.35 Å** | electrostatic |
| **S2** R774 | R774 NH1 | single | electrostatic |
| **S3** TEL triad E168/E169/E171 | 3 carboxylates | 6.58, 10.22, **7.31** Å | electrostatic |
| **S4** TEN 129–132 loop | A130/L131/R132 | 7.56, 7.33, 7.23 Å | hydrophobic |

## 2. The five modules — one site each

| id | site | job | register | verdict |
|---|---|---|---|---|
| **MD1-GROOVE** | S1 | block the groove | E at i, i+3 → 5.05 Å vs 4.35 Å measured | **ONE-SITE-FULL** |
| **MD2-R774** | S2 | engage R774 | single E anchor | **ONE-SITE-FULL** |
| **MD3-PATCH** | S3 | read the patch | K at i, i+4, i+7 | **SITE-UNSOLVABLE** |
| **MD4-ANCHOR** | S4 | greasy anchor | hydrophobic face, no charge | **ONE-SITE-FULL** |
| **MD5-REPAIR** | S1 | supply the lost charge | E at i, i+3 | **ONE-SITE-FULL** |

### Finding 1 — one site cannot be solved by a helix at all

S3 has three pairwise separations, and they are not independent: 6.58, 10.22 and
7.31 Å. A helix offers 6.20 (i+4) and 10.53 (i+7), which match two of them to
0.4 Å. The third, **E169–E171 at 7.31 Å, falls between i+3 (5.05) and i+4 (6.20)
and cannot be reached at any register** — the closest achievable is 6.20, an
error of 1.11 Å.

A single α-helix therefore **cannot** present TPP1's carboxylate triad. Two of
three is the geometric maximum. TPP1 does it because the residues are on a loop
and a β-hairpin, not one helix — which is exactly the kind of thing a design has
to admit rather than paper over.

### Finding 2 — solving one site can force you onto another

MD3's register on S3 is K at i, i+4, i+7. But i+4 → **6.2 Å**, and S1's K78–R132
separation is **4.35 Å**… and i+3 → **5.05 Å**, inside tolerance. So the module
built to read the patch *also* fits the groove.

That is recorded as `ONE-SITE-FULL-WITH-RISK`, not as a pass and not as a
failure: whether those two lysines land on the groove is a docking question, and
no docking can be run here. It is a risk to design around, and precisely the
thing the ONLY rule exists to catch.

### Finding 3 — two modules can want the same site with opposite intent

MD1 blocks the S1 groove; MD5 repairs it. The merge engine **refuses** the
combination rather than silently building a molecule that asks one groove to be
simultaneously occupied and vacant:

> *"MD1-GROOVE and MD5-REPAIR both target K78/R132 groove with opposite intent;
> one groove cannot be both blocked and repaired."*

## 3. The merges, and whether fusion is stable

| merge | modules | linker | length | verdict |
|---|---|---|---|---|
| **MRG-MIN** | MD3 + MD1 | rigid | 50 aa | STABLE |
| **MRG-BLOCK** | MD3 + MD1 + MD2 | rigid | 76 aa | STABLE |
| **MRG-SAMESIGN** | MD1 + MD2 + MD4 | rigid | 76 aa | STABLE |
| **MRG-REPAIR** | MD5 + MD4 | rigid | 50 aa | STABLE |
| **MRG-FLEX** | MD3 + MD1 | **flexible** | 50 aa | **COLLAPSE-RISK** |
| MRG-ILLEGAL | MD1 + MD5 | — | — | **REFUSED** |

Four stable, one collapse, one refused.

### The stability result that matters

`MRG-MIN` and `MRG-FLEX` are **the same two modules**. The only difference is the
linker, and it decides the answer:

- **rigid `EAAAK`** holds the helices ~12 Å apart → no attraction
- **flexible `GSGSG`** lets them collapse to ~5 Å → **−6.04 kcal/mol** of
  intramolecular attraction between a basic module and an acidic one

Each module is fine alone. Together, with the wrong linker, they fold onto each
other rather than onto the target. **That is the failure mode modular design
introduces, and it is only visible after merging.**

### Two models disagree, and the disagreement is printed

Biopython's instability index calls **every one of these constructs "predicted
unstable"** (48.9 → 102.9, threshold 40) while the electrostatic and geometric
model calls them STABLE. The index is calibrated on globular proteins with a
dipeptide scale and flags Ala/Glu-rich designed helices systematically, so it is
reported as **advisory and does not decide the verdict** — but it is printed on
every construct, because a silent disagreement between two models is how one of
them later gets quoted out of context.

## 4. Which verdicts are real and which are a choice

The aggregation rule has thresholds, and a verdict that flips inside a plausible
range of a constant is not a finding. So the thresholds were swept — 12
combinations of R ≥ {0.30, 0.35, 0.40, 0.45} against hydrophobic fraction ≥
{0.45, 0.50, 0.55}:

| construct | flagged at how many of 12 settings |
|---|---|
| **MD4-ANCHOR** | **1 of 12** — threshold-sensitive, not a result |
| **MD3-PATCH** | **3 of 12** — threshold-sensitive, not a result |
| MD1, MD2, MD5, and all four stable merges | **0 of 12** — clean |

MD4 *must* present a greasy patch to cover S4, and the patch that covers S4 is
close to the patch that aggregates. That tension is a property of the site, not of
this design; the honest statement is "MD4 sits on the boundary", not "MD4 is
unstable".

## 5. Why "fully solves only one problem" is only partly achievable

Reported plainly, because it is the direct answer to the request:

- **4 of 5 modules** solve their site fully and touch nothing else.
- **1 site (S3) is not helix-solvable at all** — a geometric fact about TPP1's
  triad, not a weakness of the design.
- **2 of the modules that do solve their site carry a register risk** of fitting a
  second site (MD3 onto S1).
- **1 merge is refused** for asking one groove to do opposite jobs.
- **1 merge collapses** purely because of the linker, not the modules.

The modular architecture did what it was supposed to: it localised every failure
to a specific, named cause. That is worth more than a merged molecule that passes
a battery while hiding four separate compromises inside one chain.

## 6. Quarantine

**Unchanged.** Nothing here has been synthesised or measured. The energy terms are
**screened-Coulomb heuristics with stated constants** (ε = 40, λ_D = 10 Å) — they
support *relative* comparisons between constructs scored identically, and are not
binding or folding free energies. No docking was run and no pose is claimed.

## 7. Tests

39 checks, 0 failed. They pin the findings, not just the arithmetic:

- every site's anchors resolve from the deposit, and the tolerance is one named
  constant (1.0 Å) used everywhere;
- every module either fully covers its site **or** is declared unsolvable with a
  named contact — asserting all five are covered would hide Finding 1;
- E169–E171 at 7.31 Å is asserted **not** reachable, which is the Finding 1 datum;
- MD3's coincidental i+3 fit to the K78/R132 groove is asserted to be a warning,
  not a pass;
- the illegal merge is refused, the refusal names the conflict, and it has no
  sequence;
- the flexible linker collapses the same module pair that the rigid one holds
  apart, and the collapse is **electrostatic** (energy < 0);
- merge integrity equals "carries no unsolvable module", which is a real check
  rather than a hopeful one;
- the threshold sweep runs over ≥ 9 combinations and at least one verdict is
  exposed as threshold-sensitive;
- a purely greasy helix is flagged as a risk at every setting.
