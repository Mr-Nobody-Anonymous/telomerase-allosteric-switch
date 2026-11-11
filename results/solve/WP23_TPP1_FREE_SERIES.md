# WP23 — The TPP1-free matched series (gap 7, closed)

Every structural comparison in this repository was **TPP1-bound vs TPP1-bound**
(7TRE / 7TRF / 7QXA / 7QXS / 9SHY). Such a design cannot measure the structural
*change* on TPP1 binding at all — with no apo state, the question is never asked.

The apo state was already on disk. `docs/ARCHITECTURE.md` calls **7TRD** the
"apo-unbound" state; chain A has 934 modelled residues spanning 7–1132 and the
deposition contains no TPP1 chain.

    python scripts/solve/wp23_tpp1_free_series.py
    -> results/solve/wp23_tpp1_free_series.csv

## The pair is genuinely matched

| | 7TRE (bound) | 7TRD (apo) | common |
|---|---|---|---|
| TERT chain A | 934 CA | 934 CA | **934 (100 %)** |
| TEN 1–196 | 160 | 160 | 160 |
| TRBD 326–540 | 188 | 188 | 188 |
| RT 541–935 | 385 | 385 | 385 |
| CTE 936–1132 | 197 | 197 | 197 |

Residue sequences are identical position-by-position, so author numbering is
aligned here — unlike the cross-entry cases the repo warns about.

## Result — and it is negative

| domain | n | mean disp (RT-core frame) | max |
|---|---|---|---|
| TEN | 160 | **0.435 Å** | 1.173 Å |
| TRBD | 188 | 0.318 Å | 0.735 Å |
| RT | 385 | 0.261 Å | 1.089 Å |
| IFD-TRAP | 31 | 0.304 Å | 0.820 Å |
| CTE | 197 | 0.529 Å | 0.969 Å |

Frame self-RMSD 0.292 Å; the global fit is 0.356 Å over all 934 CA.

**The honest reading is not "TEN moves when TPP1 binds."** It is:

> No displacement larger than the coordinate uncertainty is detectable between
> the TPP1-bound and TPP1-free TERT models. The largest of any kind is
> **0.53 Å**.

So any account in which TPP1 binding throws a *large-scale* mechanical switch in
TERT is **not supported by these coordinates**. That is consistent with a local,
sub-Å rearrangement at the TEN–TRAP–TPP1 surface — which is what the WP1/WP4
geometry already implied, now measured rather than assumed.

## Two caveats, both material

1. **7TRD and 7TRF share one TERT model.** Their chain-A RMSD is **0.0005 Å** —
   they are the same coordinates. There is really *one* TPP1-free TERT model,
   from the same group as 7TRE, so part of even this 0.4 Å may be refinement
   variation rather than a TPP1-induced change. This also disposes of the
   earlier "7TRD/7TRF are two apo structures" framing.
2. At 3.2–3.3 Å nominal resolution, a 0.4 Å domain motion is at the edge of what
   these maps distinguish. It cannot be *attributed* to TPP1 with confidence.

## What this establishes

The matched series now exists and is computed, where previously the absence of a
large motion was not even testable. It converts an untested assumption into a
measured null result, and it sharpens the mechanism question: TPP1 does not act
by moving TERT's domains by a resolvable amount.

## Note on the implementation

The first version of this script used a transposed Kabsch rotation
(`U·D·Vᵀ` instead of `Vᵀ·D·Uᵀ`) and reported **3.55 Å** where the correct value
is **0.356 Å**. It looked like a number, so it survived review. It is now pinned
by `tests/test_wp23_series.py`, which checks two *exact* cases — a structure
against itself, and a random rigid move — both of which must return zero.
