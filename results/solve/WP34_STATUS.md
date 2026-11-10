# WP34 — What is not achieved, how to achieve it, the best candidate, and the odds

> **SUPERSEDED IN PART (WP35).** Three claims below do not survive: P2 is not a
> re-analysis (it needs the chase +/- a DNA trap), the "88 % explained" figure
> compares a probability to the anomaly it is meant to explain and is negative
> in the +PT arm, and the maintenance gate is **UNDETERMINED** rather than
> CLOSED — the two physically derived stall couplings give delta* ≈ 57 bp/PD,
> above the 50 floor. See `WP35_P2_IDENTIFIABILITY.md`. The candidate ranking,
> the pipeline arithmetic and the quarantine status are unchanged.

    python scripts/solve/wp34_status.py
    -> results/solve/wp34_status.csv

## One new result first

The 1.46× tension (M5 predicts 2.73×, L&C measure ≥4×) has been carried as an
open discrepancy. **It is 88 % explained by a quantity L&C already report.**

```
required retry multiplier     = 4.00 / 2.73                    = 1.465
required retry probability    = 1 − 1/1.465                    = 0.318
stalled fraction of survivors = (0.7914 − 0.570)/0.7914        = 0.280
agreement                     = 0.280 / 0.318                   = 0.88
```

The stalled population is not a fudge factor — L&C write that *"stalled and
dissociated complexes could contribute to the untranslocated population."* So
the prediction is now **sharp**: the retry probability should measure **≈ 0.32**,
and **P2 measures exactly that**. That is a much stronger claim than "stalling
probably explains it", and it is falsifiable in days with no new reagent.

---

## 1. What is not achieved

| # | not achieved | blocker | how to achieve it |
|---|---|---|---|
| 1 | **M5 not confirmed** | only candidate that is physical *and* falsifiable, but untested independently | **P2 DNA-trap re-analysis** of existing gels — no new reagent |
| 2 | **4× vs 2.73× unresolved** | 88 % explained, not measured | same P2 run reports the retry probability, predicted ≈ 0.32 |
| 3 | M2 vs M4 degeneracy | `φ = k_tr/(k_tr+k_off_tr)` is 2 unknowns to 1 observation | **not load-bearing** — `φ` left the product, `s` never sees `k_tr` |
| 4 | **Maintenance unreachable** | δ* = 48.0 bp/PD at 2.73×, below the 50 bp/PD floor | measure the **true** k_off ratio, not the ≥4× bound |
| 5 | No compound binds the TEL patch | best number is a rigid-body surrogate | MM-GBSA or FEP at the patch |
| 6 | TELIN-1 has zero affinity data | never assayed | synthesise + SPR/ITC (it is now buildable) |
| 7 | Gate 3 unreachable | depositors released a threshold fraction, not per-band counts | request per-band densitometry |
| 8 | 9SHZ/9SI0 not ingested | TLS-blocked from this sandbox | one run elsewhere; parser is written and tested |
| 9 | Gate 4 / Gate 5 | dNTP chase not run; δ is a judgement | wet lab / a decision the repository cannot make |
| 10 | Y176A mutant panel | designed, never run | wet lab: processivity + TPP1 stimulation |

---

## 2. Best candidates — three different questions

### (a) Best mechanism: **M5**, "two-state, stalls retry"

The only candidate that is **both physical and falsifiable**. M1, M2, M3 and M4
all break `s₀ ≤ φ₀` on L&C's own no-TPP1 data, so every earlier model is
unphysical on the data it came from. M5 predicts a number (2.73×) it was **not
fitted to**, and derives δ* = **47.98 bp/PD** — matching L&C's independently
observed *"eight repeats (48 nt)"*.

### (b) Best molecule: **TLN-amide** (neutral amide anchor)

```
scaffold: O=C1C(=S)N(Ar)C(=O)N(Ar')C1CC(C)C   with a –CH2CH2–C(=O)NH2 arm
MW 441.53   logP 4.00   TPSA 83.71   HBD 1   HBA 4   rotB 7   heavy 31
Lipinski: PASS all four
```

The **only one of seven variants** that prefers the interface over a
size-matched control patch (Δ −0.199 per ligand heavy atom).

**Why it is plausible:** the TEL patch is **ambipolar** — net −1, with ARG92,
LYS170, ARG180 (+1) and ASP166, GLU168, GLU169, GLU215 (−1). A **neutral** anchor
avoids both charge mismatches, which is exactly why an "acid *or* base" search
could never have found it. The data support this: protonating the acid makes it
**worse** (+0.138), protonating the amine makes it **better** (−0.144) — the
signature of an ambipolar surface, not a uniformly charged one.

**What it is not:** a drug. It has never been made, never been assayed, and its
ranking comes from one rigid protocol with no solvent, entropy or flexibility.

### (c) Best next experiment: **P2, the DNA-trap re-analysis**

Re-analysis of data that already exists. It confirms M5 near 2.73×, falsifies
it if the ratio stays ≥4×, **and** now reports the retry probability.

---

## 3. Success rate — the arithmetic

### Mechanism

| stage | P(pass \| prior) |
|---|---|
| M5 derives the observed 48 nt exactly | 90 % |
| the 4× gap is retry (88 % explained, one free parameter) | 65 % |
| `s = (1+a)⁻⁶` survives an independent enzyme | 60 % |
| **M5 survives as the correct mechanism** | **≈ 35 %** |

These are **my priors, not frequencies** — stated so they can be argued with.
The weakest link is the third, and it is also the cheapest to test.

### The hard gate: maintenance

| | value |
|---|---|
| δ* at the measured ≥4× | **70.08 bp/PD** |
| lowest measured human attrition | **50 bp/PD** |
| verdict at 2.73× | **CLOSED** (48.0 < 50) |
| verdict at 4× | **open, but insufficient** (70.1 < 90) |
| needed for a 90 bp/PD cell | **5.14×** — have ≥4× |
| needed for a 50 bp/PD cell | **2.85×** — have ≥4× |

A 90 bp/PD cell needs **1.29×** the measured ratio: the true value would have to
sit **29 % above the reported lower bound**. Not absurd — the bound is a bound —
but not established either.

### Drug pipeline, conditional on the gate opening

| stage | P(pass) | cumulative |
|---|---|---|
| hit at the TEL patch (FEP/MM-GBSA) | 30 % | 30 % |
| IC50 < 100 µM in vitro | 40 % | 12 % |
| selectivity counterscreen | 45 % | 5.4 % |
| cellular recruitment readout | 50 % | 2.7 % |
| permeability / metabolic stability | 35 % | 0.95 % |
| in vivo PD + safety (no oncogenesis) | 30 % | 0.28 % |
| phase I → approval | 8 % | **0.023 %** |

Combined with the mechanism estimate: **P(approved therapy) ≈ 0.008 %,
conditional on the gate opening.**

**Unconditional it is zero** — because the gate is *closed*, not merely
improbable. 48.0 bp/PD cannot outrun 50 at any [TERT]; coverage saturates at 1
and there is nothing above it.

### What is *not* near zero

> **P(a publishable, falsifiable mechanism result) ≈ 35 %.**

That is the realistic value of this programme, and it is high for what it cost.
The repository's honest output is a mechanism with one cheap decisive test
attached, not a drug.

---

## 4. Order of operations that maximises the odds

| # | step | cost | buys |
|---|---|---|---|
| 1 | **P2 DNA-trap re-analysis** | days, no reagent | confirms/falsifies M5; predicts retry ≈ 0.32 |
| 2 | **true k_off ratio** (not the bound) | one kinetics run | decides whether 5.14× is reachable |
| 3 | synthesise **TLN-amide** + SPR | weeks | the first real affinity number, ever |
| 4 | MM-GBSA/FEP at the patch | compute | ranks analogues before synthesis |
| 5 | Y176A mutant panel | weeks | tests the anchor site independently |
| 6 | cell maintenance assay | months | the only test of the real objective |

**Steps 1 and 2 are cheap, and either can end the programme honestly.** Doing 3
before 2 would be synthesising against an unknown target.

---

## 5. Tests

18 checks, 0 failed — the retry multiplier is 1.465 and the needed probability
0.318; the stalled fraction is 0.280 and explains 80–95 % of the gap; δ* is 48.0
at 2.73× and 70.1 at 4×; **2.73× cannot reach the 50 floor but 4× clears it**;
4× still cannot reach a 90 bp/PD cell; both required ratios round-trip through
`a = a₀/2^(1/δ)`; the drug pipeline product is between 0.01 % and 1 %; and
TLN-amide is 31 heavy atoms and passes all four Lipinski rules.
