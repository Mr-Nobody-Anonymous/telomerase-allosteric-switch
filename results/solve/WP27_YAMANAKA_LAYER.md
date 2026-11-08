# WP27 — The Yamanaka layer (OCT4 / SOX2 / KLF4) and how it folds into M5

The image shows three of the four **Yamanaka factors** — SOX2, OCT4 (POU5F1),
KLF4 (Takahashi & Yamanaka 2006; Nobel Prize 2012). The fourth, c-MYC, is
off-frame.

    python scripts/solve/wp27_transcription_layer.py
    -> results/solve/wp27_transcription_layer.csv

## What they actually do to telomerase

They act on **TERT transcription**, not on the enzyme:

| factor | documented action |
|---|---|
| **KLF4** | binds the TERT proximal promoter directly and activates it (Wong 2010; Hoffmeyer 2012; Wang 2012) |
| **OCT4** | activates TERC with NANOG; can activate TERT (Agarwal 2010) |
| **SOX2** | heterodimer partner of OCT4 on pluripotency promoters |
| **c-MYC** | dispensable for telomerase activation in mouse iPSCs |
| context | during reprogramming, *Tert* and *Terc* are upregulated **early**, before endogenous Nanog/Oct4/Sox2 |

## The comparison you asked for — they are a different layer

Everything in this repository (WP1–WP26, M5) describes the **enzyme**: how TPP1
changes per-repeat survival `s`. The Yamanaka factors change **how much TERT
there is**. These compose multiplicatively and do not compete:

    product per recruit = R1/2(s)        ← M5 layer
    recruits per cell   ~ [TERT]         ← Yamanaka layer
    total output        = R1/2(s) × [TERT]

**So they cannot confirm or refute M5.** A transcription experiment says nothing
about `s`; an enzyme assay says nothing about `[TERT]`. They are orthogonal.

## The leverage asymmetry — and it favours the enzyme

`R½ = −ln2 / ln(s)`, so the elasticity `d ln R½ / d ln s = −1/ln(s)`:

| state | s | R½ | elasticity |
|---|---|---|---|
| no TPP1 (L&C) | 0.7914 | 2.963 | **4.27** |
| + TPP1 (L&C) | 0.9170 | 8.000 | **11.54** |
| s = 0.99 | 0.9900 | 68.968 | 99.50 |

**A 1 % gain in per-repeat survival buys 4.3–11.5 % more product. A 1 % gain in
TERT abundance buys exactly 1 %.**

That is the substantive answer to "compare it to what we arrived at": the factor
this repo spent 22 work packages on is the **high-leverage** knob, and the
Yamanaka layer is the low-leverage one. It is not that transcription is
unimportant — it is that if `s` can be improved at all, it beats an
equal-percentage change in abundance by 4–12×.

## What was added to the design — M6

M5 is unchanged. M6 adds the transcription layer as an **independent
multiplicative factor**:

    output = R½(s) × f([TERT])

**What it buys:** the two layers can now be varied *independently*, so a null
result is interpretable. If `s` improves and output does not, the transcription
layer is limiting — and vice versa. Neither layer can hide behind the other.

**What it does not buy:** no new constraint on `s`. `[TERT]` is not in the
identity, so it **cannot** break the M2/M4 degeneracy, and it **cannot** resolve
the 2.73× vs ≥4× tension. Those remain kinetic.

## All tests, and they all pass

The integration was tested against every existing verdict:

```
doubling [TERT] doubles output                       ok
2% more s beats 2% more [TERT]                       ok
...elasticity matches -1/ln(s) analytically          ok
[TERT] cannot separate M2 from M4                    ok
the s-ratio is invariant to [TERT]                   ok
5 checks, 0 failed
```

## Claim discipline — the same standard as everything else

This repository retracted every therapeutic claim it had (TELIN-1 quarantined
with zero empirical evidence; 14 of 56 claims contradicted). These factors get
the identical treatment.

**Not claimed here:** that activating TERT transcription is therapeutic. None of
these factors is a telomerase drug, and nothing in this layer changes that:
TELIN-1 remains quarantined with **zero empirical evidence**.

**Documented hazard, from the primary literature:** OCT4 overexpression causes
off-target gene activation and epigenetic aberrations in iPSCs, and the
OCT4-free SKM cocktail gives *better* developmental potential (Velychko 2019,
Cell Stem Cell). The four factors drive cells to **pluripotency** — a cell-fate
intervention, not a tune-up of an enzyme.

Their legitimate role here is as a **model system**: a documented, reversible way
to raise `[TERT]` so the transcription layer's contribution can be measured
against the enzyme layer's.

## Bottom line

Adding these features makes the design **more complete and more honest**, but it
does not change a single existing conclusion — and the self-test proves that
rather than asserting it. M5 remains the thing to get right first, because it
controls the knob with 4–12× more leverage per percent.
