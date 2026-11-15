# WP39 — The two things WP37 called impossible, done

    python scripts/solve/wp39_executable.py
    python scripts/solve/wp39_executable.py --self-test      -> 32 checks, 0 failed

WP37 closed with a table of six items it said could not be closed here. **Two of
those were wrong**, and both were wrong for the same reason: they were classified
as blocked before the literature or the code on disk had been checked.

| WP37 said | actually |
|---|---|
| `c > 1` extensions/division — "nothing establishes it" | **measured**, per cell line, since 2011 |
| TEL-patch affinity — "no FEP/MM-GBSA, so not possible" | no *affinity* is possible; the **refutation test was runnable and was never run** |

The second one is the more serious mistake. WP37 wrote a refutation criterion —
a specific number, a specific mutant, a specific negative control — and then
never evaluated it. **A criterion you never evaluate is not a prediction.**

---

## 1. `c` is measured, and the news is bad for the gate

`c` is the fraction of telomere ends extended per cell cycle. It is not a free
parameter; it has been measured repeatedly.

| line | `c` | bp/PD at R½ = 8 | ≥50? | ≥90? | source |
|---|---|---|---|---|---|
| **H1299** (adenocarcinoma) | **1.00** | **48.0** | **NO** | NO | Zhao 2011, *Mol Cell* 43:251, Fig 6E–F |
| HeLa (cervical) | 0.70 | 33.6 | NO | NO | Zhao 2011, Fig 6C–D |
| hESC | 0.70 | 33.6 | NO | NO | Zhao 2011, Fig S5 |
| HeLa, pH 7.4 | 0.77 | 37.0 | NO | NO | *NAR* 44:8086 (2016), Fig 2B–C |
| HeLa, pH 6.8 | 0.48 | 23.0 | NO | NO | *NAR* 44:8086 (2016) |
| *S. cerevisiae* | 0.07 | 3.4 | NO | NO | Teixeira 2004, STEX |

**This closes the gate harder, not softer.** WP35/WP37 recorded `c > 1.04` as the
one escape hatch from a CLOSED gate. The best measured value anywhere is
`c = 1.00` (H1299, where 100 % of lagging-daughter ends are extended every
cycle) and it still lands **2 bp/PD short** of the 50 floor.

| requirement | needs `c` ≥ | reachable at measured `c`? |
|---|---|---|
| 50 bp/PD floor | 1.04 | **no** — max measured is 1.00 |
| 90 bp/PD (HUVEC-like) | 1.88 | **no** |

Equivalently, the gate now depends only on R½, which would need **8.33** to clear
the floor and **15.0** for HUVEC, against **8.0** measured. Note that `c` cannot
exceed 1 by much and still mean anything — it is a fraction of ends.

**What follows:** the escape hatch is not an assumption any more, and it is
slightly *airtight*. Anyone reopening this gate must supply an R½ above 8.33,
not an argument about `c`.

---

## 2. The docking criterion, evaluated — and refuted

WP37 proposed: neutralise the ambipolar TEL patch (net −1 → 0) and the neutral
amide should **lose** its interface preference, while a hydrophobic control
should not. Threshold −0.15 per heavy atom. It was never run.

It is run here with WP33's size-normalised rigid-body scorer.

**The patch is the seven named residues, not the interface.** `WP33` is explicit
that the interface total (−17.89) is the *whole* 47-residue TERT+TPP1 contact
surface, not the TEL patch. Scoring against it would re-introduce exactly the
confound WP33 retired — and the first draft of this script did exactly that. The
patch is restricted to ARG92, ASP166, GLU168, GLU169, LYS170, ARG180, GLU215:
**7 residues, 66 atoms, formal charge −1.00**.

| ligand | WT | mutant | Δ (mut−wt) | prefers |
|---|---|---|---|---|
| neutral amide (TLN-amide anchor) | −0.322 | −0.368 | **−0.046** | neutralised |
| carboxylic acid | −0.335 | −0.258 | +0.077 | WT |
| carboxylate | −0.269 | −0.258 | +0.011 | WT |
| amine / ammonium | −0.236 | −0.352 | −0.116 | neutralised |

WP37's criterion, applied literally:

| condition | result |
|---|---|
| amide loses preference when neutralised | **False** |
| effect exceeds 0.15 per heavy atom | **False** |
| amide moves more than the controls | **False** (0.046 vs 0.068) |

> **VERDICT: the proposal as written is REFUTED by its own scorer.**

The amide is **not** charge-sensitive at the patch. It sits 0.046 apart across a
full neutralisation of a −1 patch, which is *less* than the controls move. The
charge hypothesis behind the neutral anchor does not survive its own test.

**This is not evidence that nothing binds.** It is a rigid-body, distance-
dependent-dielectric, no-desolvation ranking over 4 000 random placements. It has
no entropy, no induced fit, no explicit water. Entirely possible that the real
answer differs — but the criterion WP37 wrote said this test would settle it, and
on this scorer it says no.

### Two quantities, one symbol (again)

The scorer's "net charge" for the patch is **−2.54**, not the formal **−1.00**,
because `wp32.residue_charge()` sums formal side-chain charges together with
partial polar charges on backbone N/O. Both are reported separately here. This is
the same defect class as C7 and the `N_INCORP` conflation — a single name over
two quantities — caught this time before it propagated.

---

## 3. What is actually still missing

| # | missing | status |
|---|---|---|
| 1 | `n`, measured | **design exists** (two-chase, WP37); needs wet lab |
| 2 | a measured second-enzyme fold | no published k_off fold for any second system |
| 3 | P2 run as a two-chase | needs wet lab |
| 4 | **any affinity at the TEL patch** | **still none** — and the docking proposal that was meant to substitute is now refuted |
| 5 | ~~`c > 1`~~ | **RESOLVED — measured, max 1.00, gate closes harder** |
| 6 | Gate 3 per-band densitometry | not deposited; provably unreachable |

One item moved from "blocked" to "closed", and one moved from "proposed" to
"refuted". Net: the programme is smaller than it was, which is the correct
direction.

---

## 4. Tests

32 checks, 0 failed. They pin:

- the six measured `c` values, each with a source, each a fraction ≤ 1;
- `c = 1.00` gives exactly 48 bp/PD and **misses** the 50 floor — the test asserts
  the miss, not a pass;
- `c_required(50) = 1.042 > 1`, i.e. unreachable at any measured `c`;
- yeast at 0.07 as the contrast that makes the point;
- the mutant adds/removes no atom and moves no atom — charge is the only change;
- the patch is exactly 7 residues, formal −1, scorer total −2.54, and the test
  asserts those two numbers **differ**;
- the verdict logic can return *refuted* in three distinct ways (wrong direction,
  too small, less than control) — a guard that can only say yes is not a guard;
- no function in the file computes a dissociation constant or reports kcal/mol;
- the scan for those forbidden tokens is scoped to production code, because the
  first two versions of the check matched **their own literal**.

## 5. Provenance

Real measured inputs: `data/external/tesmer2023/*.csv` plus the literature `c`
values above. The docking scores are from WP33's scorer in this repository.
No affinity, Kd, or IC₅₀ is reported anywhere, because none was measured.
