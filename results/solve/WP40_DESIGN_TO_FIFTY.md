# WP40 — Designing past 50 %: the odds and the floor

    python scripts/solve/wp40_design_to_fifty.py
    python scripts/solve/wp40_design_to_fifty.py --self-test    -> 43 checks, 0 failed

Two 50 % targets were requested and both are addressed. The work turned out to be
less about hitting a number than about noticing that **the number was being read
off the wrong input.**

| target | before | after this analysis |
|---|---|---|
| **P(M5) > 50 %** | 0.35 (0.70 × 0.50) | **reachable** — needs only the P2 term at 0.71, i.e. **1.43× in odds** |
| **δ\* > 50 bp/PD** | "CLOSED at 48.0" | **OPEN at 70.1** on L&C's *measurement*; the 48 was M5's own prediction |

---

## 1. The anchor was M5's output, not L&C's measurement

Every gate statement since WP29 has used `R½ = 8.0` as "the +PT measurement".

```
8.0 / 2.963 = 2.7000×      ←  M5's PREDICTED fold
```

`WP29` labels that row *"M5 prediction (2.73×)"*. L&C's actual **measured**
statement is **≥4×**, which is `R½ ≥ 11.85`.

| anchor | R½ | bp/PD | provenance |
|---|---|---|---|
| no TPP1 (L&C) | 2.963 | 17.8 | **measured** |
| "8.0, +PT" | 8.00 | 47.4 | **M5's prediction, fed back** |
| L&C measured ≥4× | 11.85 | **70.1** | **measured lower bound** |

**So "the gate is CLOSED at 48" is a statement about M5's prediction.** On the
measurement the gate is open at 70.1 — past the 50 floor and within reach of a
90 bp/PD HUVEC-like requirement.

This inverts what the closed gate means. It is not evidence against the
*therapeutic* target; it is evidence against **M5**: the mechanism predicts less
processivity than has already been observed.

This is the same defect class as C7 — one symbol, two quantities — but here the
consequence is the opposite of before. Previously a conflation made the gate look
**open**; now one has made it look **closed**. Both were wrong for the same
reason, and the fix is the same: name which quantity you mean.

---

## 2. The gate sits at τ ≥ 0.042, not at a choice of 0 or 0.18

`δ*` is not constant — it is a function of the stall share τ, which is the very
parameter P2 exists to measure.

| τ | true fold | R½ | bp/PD | ≥50 | ≥90 |
|---|---|---|---|---|---|
| 0.000 | 2.7335 | 8.00 | 48.0 | no | no |
| **0.042** | 2.848 | 8.36 | **50.2** | **YES** | no |
| 0.100 | 3.007 | 8.79 | 52.8 | YES | no |
| 0.1833 | 3.235 | 9.46 | 56.7 | YES | no |
| 0.350 | 3.690 | 10.78 | 64.7 | YES | no |
| 0.500 | 4.100 | 11.97 | 71.8 | YES | no |

WP35's "OPEN at 57.3" (τ = 0.1833) and WP37's "CLOSED at 48.0" (τ = 0) are the
**same function evaluated at two modelling choices**. Neither was a measurement.

The honest statement: **the gate is open for any stall share above 4.2 %.** That
is a far weaker thing to need than 18 %, and it is testable.

A 90 bp/PD HUVEC-like requirement needs τ ≥ 0.88, which is beyond any plausible
stall share — **that ceiling does not move**, and this document does not pretend
otherwise.

---

## 3. E2 needs half of all complexes to stall

E2 (the branch where the stored ≥4× is stall-inflated down to 2.73×) requires

```
(1 + τ₋)/(1 + τ₊) = 4.00 / 2.7335 = 1.4633
```

| τ(−PT) | τ(+PT) required | verdict |
|---|---|---|
| 0.10 | −0.248 | **impossible** |
| 0.1833 | −0.191 | **impossible** |
| 0.30 | −0.112 | **impossible** |
| 0.40 | −0.043 | **impossible** |
| **0.4633** | 0.000 | boundary |
| 0.50 | +0.025 | possible |

E2 is **arithmetically excluded** below τ(−PT) = 0.4633. The repository's own
τ₀ = 0.1833 is far below that — so **E2 as written contradicts the stall model it
is supposed to rescue.** Either nearly half of all untrapped complexes stall, or
E2 is dead and L&C's 4× is real.

---

## 4. A test M5 had already been run against

Zhao et al. 2011 measured **~60 nt added per telomere end per cell cycle *in
vivo*** (H1299 and HeLa, under maintenance). That number depends on no chase, no
trap, and no model.

| τ | M5 predicts | vs 60 nt |
|---|---|---|
| 0.00 | 48.0 | 0.80 |
| 0.1833 | 56.7 | 0.945 |
| 0.50 | 71.8 | 1.20 |

M5 lands within a generous ±30 % band for τ ∈ [0, 0.88], and its τ = 0 value is
already inside it. **This does not confirm M5** — the band is wide and the
prediction was not pre-registered against this datum. But it is an independent
measurement M5 does not fail, which is more than C6 can say for Wang.

### A mistake I made here, and why it matters

My first draft multiplied this datum in as a further `× 0.80`, which *lowered*
P(M5) to 0.45. That is wrong in kind. A factor means "M5 must **also** pass
this" — but M5 has largely already passed it. A hurdle already cleared is not a
further hurdle, and no likelihood ratio is available for it, so the honest
treatment is to leave it **out** of the product and report it as a datum.

Multiplying it in would have been the same error class as **C1**: counting
something as evidence without asking what it could have falsified.

---

## 5. What moves P(M5) past 50 %

```
P(M5) = P(physical on L&C) × P(P2 confirms) = 0.70 × 0.50 = 0.35
need 1.857× in odds to reach 0.50
```

The lever is the **P2 term**. Holding 0.70 fixed:

| P(P2) | P(M5) |
|---|---|
| 0.50 (current) | 0.350 |
| 0.60 | 0.420 |
| **0.714** | **0.500** ← target |
| 0.80 | 0.560 |
| 0.85 | 0.595 |
| 0.20 (P2 fails) | 0.140 |

**The target needs the two-chase result to be only 1.43× more likely under M5
than under its alternatives.** A decisive experiment clearing 1.4× in odds is an
ordinary ask. This is the smallest sufficient design, not a hopeful one.

---

## 6. The design, ordered by information per unit cost

| # | measurement | cost | moves |
|---|---|---|---|
| 1 | **Two-chase ± DNA trap** (WP37 inversion) | 1 assay | τ exactly, fold exactly |
| 2 | **R½ measured directly at +PT** | 1 gel | the anchor, 8.0 → measured |
| 3 | ***In vivo* extension per end per cycle ± TPP1** | 2 lines | the gate vs a model-free datum |
| 4 | Stall-share sweep to τ ≈ 0.5 | expansion of 1 | E2 from excluded to testable |

Each with its falsifier:

**1. Two-chase.** *De-confounds:* `n` and `τ` cancel in `m/c`, so no model
assumption survives. *Falsified if:* a trap **raises** processivity (`m < c`) —
that refutes stalling before any fitting.

**2. Direct R½ at +PT.** *De-confounds:* removes M5's own prediction from its own
evidence — the defect in §1. *Falsified if:* R½ at +PT comes out **below 8.0**,
which kills M5's fold outright.

**3. In vivo extension ± TPP1.** *De-confounds:* no chase, no trap, no model; the
baseline is already 60 nt. *Falsified if:* TPP1 knockdown **does not lower** that
figure — a processivity mechanism predicts it must.

**4. Stall sweep to τ ≈ 0.5.** *De-confounds:* settles whether ≥4× is
stall-inflated at all (§3). *Falsified if:* τ(−PT) measures **well below 0.4**,
killing E2 outright.

---

## 7. What this does not change

- **No affinity exists** for TLN-amide, and WP39 refuted the docking proposal
  that was meant to substitute for one.
- **The HUVEC 90 bp/PD ceiling does not move** — it needs τ ≥ 0.88 or R½ ≥ 15.0.
- **`c` cannot be pushed past 1.00**; it is a fraction of ends (WP39).
- **TELIN-1 stays quarantined.** Nothing here touches it.
- **P(M5) is still a personal prior times a hypothesis.** Reaching 0.56 by
  assuming P2 clears at 0.80 is a *plan*, not a belief about the world.

The one thing that genuinely improved: the gate is open for any stall share above
**4.2 %**, and the reason it looked closed was that the input was M5's own
prediction.
