# WP28 — Improving M6: the operating point decides which knob binds

WP27 answered *"compare the Yamanaka factors to what we have"* with a number:
**1 % more `s` buys 4.3–11.5 % more product, 1 % more [TERT] buys 1 %.** That
number is right and it is still the wrong decision rule. Two defects, both of
which change what the design should do:

1. **`4.3–11.5` is not a range.** It is two local slopes evaluated at two
   particular operating points. The elasticity is `ε = −1/ln(s)`, which
   **diverges** as `s → 1`. At `s = 0.99` it is 100, at `s = 0.999` it is 1000.
2. **M6 assumed `output ∝ [TERT]`.** Telomerase must recruit to a telomere
   first, and occupancy saturates. The linear term is the `x ≪ 1` limit of a
   saturating one.

Fixing both is WP28. The conclusion flips from *"the enzyme knob always wins"*
to something narrower and more useful: **per percent the enzyme knob always
wins; in total range the winner depends on where you start.**

    python scripts/solve/wp28_operating_point.py
    -> results/solve/wp28_operating_point.csv

## 1. The elasticity is not bounded

| s | R½ | ε = −1/ln s |
|---|---|---|
| 0.7914 | 2.963 | 4.27 |
| 0.9170 | 8.000 | 11.54 |
| 0.9500 | 13.513 | 19.50 |
| 0.9800 | 34.310 | 49.50 |
| 0.9900 | 68.968 | **99.50** |
| 0.9990 | 692.801 | **999.50** |

WP27 listed 4.27 and 11.54 and called them a range, then used 0.99 as a third
row without noticing it breaks the range. Quoting the anchors as bounds
**understates** the enzyme layer and implies its advantage is capped. It isn't.

## 2. M6 is the linear limit of a saturating model — M6′

    c(x)   = x / (1 + x)                    x = [TERT] / K_d   (coverage)
    output = N_tel · c(x) · R½(s)           M6′
    ε_x    = 1 / (1 + x)   (exactly = 1 − c, verified to 1e-15)

| x | coverage | ε_x | output at s₀ |
|---|---|---|---|
| 0.01 | 0.0099 | 0.990 | 0.029 |
| 0.10 | 0.0909 | 0.909 | 0.269 |
| 1.00 | 0.5000 | 0.500 | 1.482 |
| 10.0 | 0.9091 | 0.091 | 2.694 |
| 100 | 0.9901 | 0.010 | 2.934 |

Two regimes that behave nothing alike: `x ≪ 1` linear (every molecule helps),
`x ≫ 1` flat (every telomere already has one). **M6 was not wrong, it was
unconditional** — and an unconditional model cannot tell a design which regime
it is in. That is the defect that actually mattered.

## 3. The corrected comparison

**Per percent, at every operating point** — this part of WP27 survives intact:

- enzyme `ε_s ≥ 4.27` (minimum, at the +TPP1 anchor)
- transcription `ε_x ≤ 1.00` (maximum, at `x → 0`)

The enzyme knob wins per percent, everywhere, always.

**In total range, from a stated start** — this is new and it does not always
favour the enzyme knob:

| start | transcription span to full coverage |
|---|---|
| x₀ = 0.01 | 4.6 ln units |
| x₀ = 0.10 | 2.3 ln units |
| x₀ = 1.00 | 0.69 ln units |
| enzyme, measured | **0.99 ln units** (2.70×) |

The enzyme figure is a **measured** range — the 2.70× M5 predicts and L&C
measure, expressed as a median-length ratio. The transcription spans are
**geometric** — how far you *can* travel, not how far anything measured *has*
travelled. They are not the same kind of number, and separating them is the
point: a TERT-negative somatic cell has 4.6 ln units of headroom on the
transcription knob and 0.99 on the enzyme knob. **The knobs don't have a fixed
ranking; the operating point has a ranking.**

## 4. The transcription knob's price rises without limit

To match the enzyme layer's measured 2.70×:

| start x₀ | coverage | fold [TERT] required |
|---|---|---|
| 0.001 | 0.0010 | 2.70 |
| 0.010 | 0.0099 | 2.75 |
| 0.100 | 0.0909 | 3.25 |
| 0.500 | 0.3333 | **18.00** |
| 1.000 | 0.5000 | **unreachable** |
| 2.000 | 0.6667 | **unreachable** |

From near-zero TERT it costs ~2.7× — the linear regime, where the buy is
almost free. At **one-third** coverage it already costs **18×** of [TERT], and
at half coverage it **cannot be bought at any [TERT]**: matching 2.70× would
require coverage 2.70 > 1. The price does not double or quadruple, it diverges.
The enzyme knob has no such penalty. That, not the local elasticity, is the real
argument for prioritising it.

## 5. Why 26 work packages never touched this layer

Every assay in this repository is **in vitro**: purified enzyme, fixed DNA,
fixed [TERT]. In the model that is `∂f/∂x = 0` *by construction* — the assay
cannot move `x`, so it has **zero** sensitivity to the transcription layer
however long it runs. That is not an oversight in the experiments; it is a
structural property of the assay, and it makes the blind spot legitimate:
WP1–WP26 measure `s` correctly and say nothing about `x` **because they
cannot**. Any test of the transcription layer has to be a **cell** experiment
with TERT abundance as the manipulated variable.

## 6. Revised recommendation (supersedes WP27's bottom line)

> **Superseded in turn by WP29.** This section still optimises *output*. WP29
> shows output is the wrong objective — a cell is maintained only if gain
> exceeds loss per division — and that at the measured enzyme kinetics **no**
> [TERT] reaches the maintenance threshold. Read WP29 for the objective that
> replaces this one; the analysis below stands as arithmetic.

1. Stop quoting "4–12×". Quote `ε = −1/ln(s)` **with the operating point
   attached**, and state that it diverges as `s → 1`.
2. Replace **M6** with **M6′**; M6 is its `x ≪ 1` limit and remains valid there.
3. **Name the operating point before comparing knobs.** From near-zero [TERT]
   transcription is the binding constraint; from saturating coverage it is dead
   and the enzyme layer is all there is.
4. The in vitro programme is still the right first investment — not because its
   per-percent leverage is 4–12×, but because it is **the only knob this assay
   can see**, and its price does not rise with the operating point the way
   transcription's does.
5. **Still no therapeutic claim.** Reaching coverage *is* immortalisation, and
   the OCT4 hazard from WP27 is unchanged.

## 7. Tests

```
eps_x == 1 - coverage analytically                     ok
enzyme elasticity > transcription elasticity at all x  ok
M6' -> M6 (linear) as x -> 0                           ok
output saturated at high x (<0.1% for 10x TERT)        ok
fold_to_match hits its target ratio                    ok
unreachable target returns None                        ok
output monotone in s and in x                          ok
enzyme range == R1/2 ratio (2.69997)                   ok
8 checks, 0 failed
```

Check 3 is the important one for the existing model: **M6′ reduces exactly to
M6** in the regime WP27 used, so nothing WP27 computed is retracted. Only its
*scope* was too wide, and check 2 pins the one inequality that was always true.
