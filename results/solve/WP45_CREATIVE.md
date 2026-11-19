# WP45 — The creative route: L&C's own two numbers over-determine τ

    python scripts/solve/wp45_creative.py
    python scripts/solve/wp45_creative.py --self-test    -> 24 checks, 0 failed

Every work package since WP34 treated the stall share τ as a free parameter that
only a **new** experiment could fix. That was the wrong frame. L&C 2010 published
**two independent per-repeat quantities from the same preparation**, and nobody
had set them against each other:

```
s_obs(−PT) = 0.7914    R½ = 2.963, from the product distribution
φ(−PT)     = 0.570     translocation efficiency, Fig 2F
```

Under the stall renewal both are functions of τ, so the pair **fixes it**:

```
s_obs  = s_true / (1 − τ(1 − s_true))       stall inflates survival
s_true = φ                                   translocation is the un-retried success

⇒  τ₀ = (s_obs − φ) / (s_obs (1 − φ)) = 0.6506
```

**No new experiment.** τ was "unmeasured" only because the two published numbers
were never cross-referenced.

---

## 1. It refutes M5 — by a factor of 2.4

With τ fixed, M5 makes a **held-out** prediction. `a₀` comes from φ(−PT),
`a₁ = a₀/2.7335` is M5's own fold, and the stall re-inflates the survival. The
+PT data are then read only to judge the result.

| | |
|---|---|
| M5 predicts (+PT fold) | **2.734** |
| implied by L&C's stored anchor | **6.654** |
| L&C stated measured | **≥ 4.00** |

M5 **under-predicts by 2.4×**. That is not a near miss at a definitional edge —
at the stall share L&C's own two numbers force, M5's fold is less than half what
the data show.

This is the strongest evidence *against* M5 the programme has produced, and it
comes from published data with **no new assumption beyond `s_true = φ`**.

---

## 2. The gate has been evaluated with a contradictory pair of assumptions

This is the more damaging finding. The stored `R½(+PT) = 8.0` is reproduced only
if the +PT stall is assumed to **vanish** (`s_true_pt → s_obs_pt`):

| reading | R½(+PT) | δ* | ≥50 |
|---|---|---|---|
| TRUE (stall at τ₀) | 1.314 | 7.88 | no |
| APPARENT | 1.333 | 8.00 | no |
| **STORED 8.0** | 8.000 | 48.00 | **YES** |

True and apparent agree *closely* here, because `τ(+PT)` is tiny (0.015). So the
discrepancy is not between those two — it is against the **stored anchor**, and
it is **6.0×**.

The programme has therefore been computing the gate with:

- `τ(−PT)` **large enough to move the fold** (0.65), and
- `τ(+PT)` **simultaneously small enough not to move the survival it is compared
  against** (→ 0).

Those two cannot both hold under one coupling rule. **That is C7 again, one level
up**: a single number carrying two readings, and choosing between them decides
the entire margin in dispute.

---

## 3. What this does and does not settle

**Settled:** M5 is refuted on L&C's own published numbers under the renewal model
with `s_true = φ`. And the stored `8.0` is not a model output — it is a
consequence of assuming no +PT stall, which was never chosen explicitly.

**Not settled:** whether `s_true = φ` is right. It is a **modelling step**, and
it is the one load-bearing assumption here. φ is L&C's window-limited
translocation efficiency; WP22 established that φ is not a per-repeat survival
probability in the naive sense, which is exactly what makes this renewal reading
possible. But it is a reading, not a theorem.

**If it is wrong**, the whole of §1 and §2 falls with it. There is no way to
check it from this repository.

---

## 4. What definitely does not move

| | |
|---|---|
| A 90 bp/PD HUVEC-like cell | δ* = 7.9 bp/PD in true terms — **short by 82** |
| Any affinity | still none; no Kd, no IC₅₀, no ΔG |
| A drug | still none |
| TELIN-1 | still quarantined |

The 90 bp/PD ceiling has survived every work package since WP29 and survives this
one. On the true reading it is not merely closed — it is *further away* than
believed.

---

## 5. Sensitivity

L&C report a **lower** bound (≥4×), so larger values are admissible:

| observed | 3.5× | 4.0× | 5.0× | 6.0× |
|---|---|---|---|---|
| τ₀ | 0.245 | 0.353 | 0.504 | 0.602 |
| implied fold | 3.74 | 4.26 | 4.97 | 5.72 |

Every admissible value strengthens the refutation in §1 rather than weakening it.

---

## 6. Errors made and caught in this work package

Recorded because they are the point, not an embarrassment:

1. **The inversion had the wrong sign.** My first `τ` formula returned **−1.86**
   for these data. Caught by a forward-substitution check, which is now a test.
2. **A docstring asserted `τ₀ = 0.2605` without computing it.** It was a guess
   written as a result. The real value is 0.6506 — and the guess would have
   reversed the conclusion.
3. **I reported a `6.1×` true-vs-apparent gap** before checking: at `τ(+PT) =
   0.015` those two readings almost coincide. The real gap is against the *stored
   anchor*, not between true and apparent.
4. **A patch deleted `main()` and five functions.** Caught by `NameError`, not by
   a test — the self-test does not import the module.

Each was caught by executing something, never by re-reading.

---

## 7. Where this leaves the programme

The honest ledger after WP45:

| row | value |
|---|---|
| M5 | **refuted** on L&C's published pair (2.4× under-prediction) |
| The gate | evaluated with inconsistent τ; the stored 8.0 assumes no +PT stall |
| The 2.73× fold | M5's own output, fed back — as established in WP40 |
| τ | **0.6506**, derived, not free |
| Affinity / drug / `n` | unchanged: none / none / not load-bearing |

The programme is now in the position of having refuted its only surviving
mechanism using data it already had. That is a better position than having 35 %
confidence in something untested — but it is not a result anyone would frame.

**The single load-bearing assumption is `s_true = φ`.** Everything in §1 and §2
rests on it, one experiment could test it, and this repository cannot run that
experiment.
