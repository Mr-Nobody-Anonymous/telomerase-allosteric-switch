# WP29 — The actually-right objective: maintenance, not output

WP27 optimised *product per cell*. WP28 optimised the same thing with a
saturating recruitment term. **Both optimised the wrong quantity**, and setting
the target correctly does more to the design than either of them did.

A cell is not maintained because it makes a lot of repeats. It is maintained
because the repeats it makes per division **exceed the repeats it loses**:

```
gain per division = 6 · c(x) · R½(s)   bp       (6 nt per repeat)
loss per division = δ                  bp
net               = gain − loss        ← has a SIGN
```

A cell that gains 48 bp and loses 90 bp per division **still senesces**. Output
is not maintenance; the difference is, and the difference has a sign.

    python scripts/solve/wp29_maintenance.py
    -> results/solve/wp29_maintenance.csv

## 1. The whole design collapses to one exact identity

Setting `net = 0`:

```
δ* = 6 · R½(s) · c(x)   ≤   6 · R½(s)              since c ≤ 1
   = 6 · ln2 / (6 ln(1+a))                         since R½ = ln2/(6 ln(1+a))
δ* = 1 / log₂(1 + a)
```

**`δ*` is a property of the enzyme alone.** Coverage drops out of the
derivation and survives only as a multiplier ≤ 1 — so full coverage is the best
case, and nothing above it exists.

| enzyme state | a | R½ | δ* at full coverage |
|---|---|---|---|
| no TPP1 (L&C) | 0.03976 | 2.963 | **17.78 bp** |
| M5 prediction (2.73×) | 0.01455 | 7.997 | **47.98 bp** |
| L&C measured (≥4×) | 0.00994 | 11.680 | 70.08 bp |
| Wang RRL (4.48×) | 0.00887 | 13.075 | 78.45 bp |

## 2. M5's own number was already the ceiling — and it matches

L&C write: *"half the chains elongated by eight repeats (48 nt)"*. M5's 2.73×
**derives** `δ* = 1/log₂(1.01455) = 47.98 bp`. Eight repeats × 6 nt = **48 nt**.

That agreement is the strongest single number in the whole programme. It is
also the ceiling: *48 is not 90.*

## 3. Every real human attrition rate exceeds that ceiling

Attrition is 50–100 bp/PD in human fibroblasts, 90 bp/PD measured in HUVEC, up
to 200 bp/PD across the literature ([5](https://pmc.ncbi.nlm.nih.gov/articles/PMC2195566/), [4](https://pmc.ncbi.nlm.nih.gov/articles/PMC6496299/)).

| δ (bp/PD) | source | coverage needed |
|---|---|---|
| 50 | fibroblasts, lower bound | 1.04 |
| 60 | fibroblasts, typical | 1.25 |
| 100 | fibroblasts, upper bound | 2.08 |
| 90 | HUVEC (measured 89–90) | 1.88 |
| 200 | literature maximum | 4.17 |

Coverage **cannot exceed 1** — it is the fraction of telomeres carrying an
elongating complex. At 50 bp/PD the design would need 1.04 telomerases per
telomere. **There is no [TERT] at which this works**, which is precisely why no
amount of the Yamanaka layer fixes it. This turns WP28's "sizing problem" into
an impossibility result.

## 4. The one lever with no ceiling

`δ*` is monotone in the off-rate ratio, and *that* is not capped:

| target δ | a required | off-rate fold | have we got it |
|---|---|---|---|
| 50 bp/PD | 0.01396 | **2.85×** | no: 2.73–4× |
| 90 bp/PD | 0.00773 | **5.14×** | no: 2.73–4× |
| 100 bp/PD | 0.00696 | 5.72× | no |
| 200 bp/PD | 0.00347 | 11.45× | no |

**The programme reduces to one number: the `k_off` ratio, target ≥ 5.14× for a
90 bp/PD cell.** M5 predicts 2.73×; L&C measure ≥4×. Both are short, and the
gap between "what we have" and "what would work" is now a single measurable
quantity rather than a list of gaps.

## 5. The transcription layer is a dead lever, not a weak one

`a` enters only through coverage, so

```
d ln(1+a) / d ln(x) = (ln2/δ) · c · (1−c)   ≤   ln2/(4δ)
```

| δ | rate per e-fold of [TERT] | [TERT] fold needed for 4× |
|---|---|---|
| 50 | 0.347 % | ~400 |
| 90 | 0.193 % | ~720 |
| 100 | 0.173 % | ~800 |

Buying a 4× off-rate reduction by abundance alone needs **~1e174-fold**. This is
the quantitative form of WP28's point: for *this* objective abundance is not
merely low-leverage, it is **inert**.

## 6. Lifespan arithmetic — and a check on it

HUVEC-like: `L₀ = 11.5 kb`, `L_crit = 5 kb`, `δ = 90 bp/PD`.

| condition | gain/PD | net/PD | PD to senescence |
|---|---|---|---|
| untreated | 0.00 | 90.00 | **72** |
| +TPP1, full coverage (2.73×) | 47.98 | 42.02 | 155 |
| at L&C measured 4×, full coverage | 70.08 | 19.92 | 326 |
| at 5.14×, full coverage | 90.01 | −0.01 | **IMMORTAL** |

The untreated figure (~72 PD) **reproduces the observed HUVEC limit (~69–70)**.
That is a check on the arithmetic, not a result — and it is why the 155 and 326
are worth taking seriously. TPP1–TEN extends lifespan ~2.2× and **still
senesces**. Only past 5.14× does the sign of the net change.

## 7. The improved design — M7

M5/M6′ are unchanged as *descriptions*. What changes is the objective:

```
objective:  net = 6 · c(x) · R½(s) − δ        maximise, need net ≥ 0
target:     δ* = 1/log₂(1+a) ≥ δ_cell
only free lever: a = k_off_ext / k_pol        target ≥ 5.14× for δ = 90
```

**What this buys:** a single falsifiable design target with a sign, replacing
"increase output". It also explains why the previous framing kept producing
bounded, unsatisfying answers — it was optimising a quantity that does not
decide the outcome.

**What it does not buy:** any new constraint on the M2/M4 degeneracy, or any
resolution of the 2.73× vs ≥4× tension. Both are untouched. What it does is
tell you the tension *matters less than it looked*: neither 2.73× nor 4× is
sufficient, so the 1.46× disagreement is not the difference between success and
failure — **the whole range is on the failing side.**

## 8. Still not a therapy

Immortalising a somatic cell is not treating anything; it is the canonical
first step of oncogenesis. Nothing here changes the OCT4 hazard from WP27, and
**no therapeutic claim is made or supported.** This is a statement about what a
design would have to achieve and which quantity to optimise.

## 9. Tests

```
delta* == 1/log2(1+a)                                    ok
delta* == 6*R1/2 from s=(1+a)^-6                         ok
2.73x predicts 48.0 bp (L&C observed 48)                 ok
no [TERT] immortalises a 90 bp/PD cell at 2.73x          ok
required fold reproduces delta* = target                 ok
net > 0 iff a < a_needed                                 ok
lifespan None iff net <= 0                               ok
untreated HUVEC lifespan in 65-75 PD                     ok
substitution rate < 1% per e-fold                        ok
gain is linear in coverage (WP28 c(x) enters once)       ok
10 checks, 0 failed
```

Check 3 is the substantive one: the model independently reproduces the 48 nt
L&C observed. Check 8 guards the arithmetic that everything else rests on.
