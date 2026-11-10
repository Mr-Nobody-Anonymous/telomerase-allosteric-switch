# WP35 — P2 implemented: the sign that was wrong, and what it does to the gate

    python scripts/solve/wp35_p2_identifiability.py
    -> results/solve/wp35_p2.csv

    python scripts/solve/wp35_p2_identifiability.py --self-test
    -> 54 checks, 0 failed

WP34 ended with a six-step order of operations whose step 1 was "P2 DNA-trap
re-analysis — days, no reagent". This work package is that step, built and run.
Building it broke the step: P2 is not a re-analysis, it costs more than "days",
and the arithmetic WP34 attached to it does not hold. Fixing that arithmetic
moved the maintenance gate, which is the only number in WP34 that mattered.

**Quarantine unchanged.** Nothing here is a therapeutic claim.

> **Follow-on (WP36).** The audit that followed this found that the 48 nt
> "agreement" WP34 relied on is an *identity* (`δ* = n·R½`), and that four other
> results share the same error class. `P(M5)` is **≈ 21 %**, not 35 % — both of
> the factors that made it look well-supported were circular. See
> `WP36_CIRCULARITY.md`.

---

## 1. Three defects in WP34's step 1

### (a) P2 is not a re-analysis. The repository had already proved it could not be.

WP34, WP25 and WP22 all describe P2 as re-analysis of existing gels. It is not:

| candidate data source | what it actually is |
|---|---|
| `data/external/tesmer2023/` | **one scalar per condition** — fraction of signal at ≥16 repeats |
| `data/benchmarks/rap_independent_benchmark.csv` | **generated**, `np.random.seed(42)`, intensities synthesised from chosen k_on/k_off/R_max |

WP24 already proved the first is hopeless: two product-length distributions
that agree at ≥16 repeats to ten decimal places differ **4.98×** in the fraction
ending at ≤3 repeats. A per-repeat survival cannot be recovered from a
threshold fraction at any level of care. And "confirming" a model against the
second would be recovering the numbers used to write it.

**So P2 is a new experiment**: the existing pulse-chase run ± a DNA trap.

### (b) The stall correction had the wrong sign — and this was in my own first draft

WP34/WP25 treat stalling as inflating the un-trapped chase. That is right. But
they apply it to the chase **alone**, and the product distribution is a
stall-inflated observation too. A retry means the complex keeps elongating, so
the observed per-repeat survival *overstates* the true one:

```
s_obs  = s_true / (1 − τ(1 − s_true))          retries inflate survival
s_true = s_obs (1 − τ) / (1 − τ s_obs)         exact inverse

κ_dist(τ)  = a(s_true0(τ)) / a(s_true1(τ))     distribution fold   — RISES with τ
κ_chase(τ) = 4.00 (1 − τ)                      trap dilutes chase  — FALLS with τ
```

| τ | s_true₀ | s_true₁ | κ_dist | κ_chase | δ* (bp/PD) |
|---|---|---|---|---|---|
| 0.000 | 0.7914 | 0.9170 | 2.7335 | 4.0000 | 48.00 |
| 0.050 | 0.7828 | 0.9169 | 2.8605 | 3.8000 | 50.21 |
| 0.100 | 0.7735 | 0.9168 | 3.0002 | 3.6000 | 52.65 |
| 0.150 | 0.7633 | 0.9167 | 3.1548 | 3.4000 | 55.35 |
| **0.183** | **0.7560** | **0.9167** | **3.2670** | **3.2670** | **57.30** |
| 0.250 | 0.7400 | 0.9166 | 3.5190 | 3.0000 | 61.69 |
| 0.300 | 0.7265 | 0.9165 | 3.7356 | 2.8000 | 65.47 |
| 0.400 | 0.6948 | 0.9163 | 4.2631 | 2.4000 | 74.67 |

The two readouts **converge from opposite sides on one τ**. That is what makes
P2 identifiable, and it is the opposite of what WP34 described.

My first draft had this running forward instead of inverted. It produced a fold
that barely moved (2.7335 → 2.9121) and an identity "δ\* = 1 + τ" that was
simply **false**. The self-test caught both before either reached a document.
That is the third time in this programme a hard-coded verdict contradicted the
numbers above it; the difference is that this time a check was in place.

### (c) "88 % explained" does not survive, on either of its two terms

WP34 compared `1 − 2.73/4.00 = 0.318` with `(0.7914 − 0.570)/0.7914 = 0.280`.

1. **0.280 is a function of (s, φ).** Those are the two *different* quantities
   whose mismatch is the anomaly M5 exists to explain. `(s₀−φ₀)/s₀` **is** that
   anomaly, rescaled — it cannot explain itself. In the +PT arm the same
   expression is **−0.0796**, and a fraction of a population cannot be negative.
2. **0.318 is not a physical retry probability.** τ is defined by the renewal
   model above, and the crossing is at **τ₀ = 0.1833**.

What replaces it is stronger, not weaker: stalling **reconciles** the two
observations at one τ with no freedom left over, and that τ is measurable.
"Reconciled at a measurable τ" beats "88 % explained by the gap itself".

---

## 2. The assumption that decides the sign — derived, not chosen

A retry can only follow a **failure**, and failures per repeat are `1−φ`:
**0.430** without POT1-TPP1, **0.010** with it. Saturating POT1-TPP1 leaves
almost nothing to retry, so the retry probability should scale with the failure
rate:

```
τ₁ = τ₀ · (1−φ_PT)/(1−φ_noPT) = τ₀ · 0.0233        [default coupling]
```

Three couplings are solved and reported so that no single one is quoted alone:

| coupling | τ₀ | fold | δ\* | vs 50 bp/PD | vs 90 bp/PD |
|---|---|---|---|---|---|
| **failure-scaled** (default) | 0.1833 | 3.2670 | **57.30** | OPEN | closed |
| none (M5 as written) | 0.1815 | 3.2741 | 57.43 | OPEN | closed |
| same τ both arms | 0.3345 | 2.6621 | 46.76 | closed | closed |

Only "same τ in both arms" keeps the gate shut, and the φ data contradict it.
But the spread from **46.8 to 57.4 bp/PD is driven by the coupling choice, not
by data** — so reporting one number without the other two would be exactly the
error WP33 caught in WP32.

---

## 3. What this does to the maintenance gate

WP34: *"gate CLOSED at 48.0 < 50, therefore P(approved drug) = 0 %, because the
gate is closed, not merely improbable."*

**48.0 bp/PD is δ\* at τ₀ = 0 — and τ₀ = 0 is precisely the value P2 exists to
test.** That is a conditional presented as a closure. The honest status:

| state | fold | δ\* | gate vs 50 | gate vs 90 |
|---|---|---|---|---|
| τ₀ = 0 (M5 as written) | 2.73 | 48.0 | **closed** | closed |
| τ₀ = 0.18 (P2, both derived couplings) | 3.27 | **57.3** | **OPEN** | closed |
| fold ≥ 5.14 (needs τ₀ ≳ 0.51) | 5.14 | 90.0 | open | open |

**Gate status: UNDETERMINED, decided by one measurement.** Not closed.

So P(approved drug) is not 0 %. It is WP34's own pipeline product, conditional
on the gate opening — which WP34 computed as **≈ 0.008 %**. Small, and *not a
finding of zero*. The distinction matters because "0 %" told the reader to stop,
and the correct instruction is "go and do the one measurement".

What does **not** change: **δ\* = 57.3 is still below a 90 bp/PD cell.** Even the
optimistic branch cannot maintain a HUVEC-like attrition rate. The maintenance
objective remains out of reach at the reconciled kinetics, and only a genuine
5.14× change in the true off-rate — step 2 of WP34's order — could reach it.

---

## 4. P2, pre-registered in both directions

| trap chase returns | τ | true fold | δ\* | verdict |
|---|---|---|---|---|
| **4.00×** | 0.000 | 2.7335 | 48.00 | τ = 0 — M5 unchanged, gate closed |
| 3.60× | 0.100 | 3.0002 | 52.65 | gate open |
| 3.40× | 0.150 | 3.1548 | 55.35 | gate open |
| **3.27×** | 0.182 | 3.2644 | 57.26 | **pre-registered prediction** |
| 3.00× | 0.250 | 3.5190 | 61.69 | gate open |
| 2.80× | 0.300 | 3.7356 | 65.47 | gate open |

Falsification: a trap chase **≥ 3.99×** means no detectable stall contribution,
M5 keeps 2.73×, and the gate is closed on the evidence. **Both outcomes are
publishable and neither is a rescue.**

---

## 5. What P2 costs — and why "days" was wrong

The primary contrast is **3.27× against 4.00×, a ratio of 1.224**, measured
against the only experiment-level CV this repository holds from a **real** gel —
Tesmer Fig 3B WT, 3 independent experiments (0.5088 / 0.5547 / 0.3929),
**CV = 17.2 %**. (`data/benchmarks/*.csv` is generated and cannot be used for a
power calculation about a real gel.)

| contrast | ratio | n per arm |
|---|---|---|
| 3.27× prediction vs 4.00× bound | 1.224 | **12** |
| 2.73× (τ = 0) vs 4.00× bound | 1.463 | 4 |
| 3.27× vs 2.73× (direct τ estimate) | 1.195 | 15 |

| pairing ρ | n per arm | CV tolerated |
|---|---|---|
| 0.0 (independent arms) | 12 | 0.177 |
| 0.5 | 6 | 0.177 |
| **0.8 (within-reaction)** | **3** | **0.198** |
| 0.9 | 2 | 0.229 |

L&C report 3 per arm. So P2 is **one extra paired experiment** (ρ = 0.8) or a
**campaign of twelve** (ρ = 0) — and which one it is depends on a *design
choice*, not on anything in this script. "Days, no reagent" was only supportable
if you silently assumed ρ = 0.8. At n = 3 the 95 % prediction interval is
**[2.48, 4.30]×**, which straddles 4.00 — so even the cheap version is marginal,
and that is worth knowing before committing to it.

---

## 6. Status against WP34's scoreboard

| WP34 said | WP35 says |
|---|---|
| P2 = re-analysis, days, no reagent | P2 = chase ± trap; n = 12, or n = 3 if paired |
| stalling inflates the chase only | stalling inflates **both** readouts |
| 4× gap 88 % explained by 0.280 | that quantity **is** the gap, rescaled; negative in +PT |
| retry probability needed ≈ 0.32 | **τ₀ = 0.1833**; 0.318 was not a τ |
| M5 prediction stays 2.73× | reconciled fold **3.267×** (2.662× under 'same') |
| gate **CLOSED** at 48.0 < 50 | **UNDETERMINED**: 57.3 / 46.8 bp/PD |
| P(approved drug) = 0 % | 0 % was the assumption τ₀ = 0, not a result |
| P(M5) ≈ 35 % | **≈ 21 %** (WP36 — both former support terms were circular; see below) |

The probability did not move *within WP35* — the gap became a **bound with a
measurable parameter** rather than a fit, which is stronger evidence. But the
follow-on audit (WP36) then removed both of WP34's scoring terms as circular, and
the honest value fell to **≈ 21 %**. See `WP36_CIRCULARITY.md`.

**A is the weakest link**, and A — `s = (1+a)⁻⁶` — is the one thing no amount of
gating, scoring or re-analysis in this repository can test.

---

## 7. What is still not achieved

| # | not achieved | why | next |
|---|---|---|---|
| 1 | **P2 not run** | needs the pulse-chase ± trap at n = 3 paired | wet lab, pre-registered above |
| 2 | gate undetermined | decided by P2 | same run |
| 3 | 90 bp/PD unreachable | needs 5.14×; reconciled is 3.27× | true k_off ratio, step 2 |
| 4 | identity A untested | `s = (1+a)⁻⁶` never tested on another enzyme | Wang 2007 data exists; not yet run |
| 5 | no compound binds the TEL patch | rigid-body surrogate only | MM-GBSA / FEP |
| 6 | TLN-amide never made | no synthesis, no Kd | synthesise + SPR (WP30/33) |
| 7 | Gate 3 blocked | per-band densitometry absent | request from depositors |
| 8 | 9SHZ/9SI0 not ingested | TLS-blocked here | one run with network access |

---

## 8. Tests

`--self-test` runs **54 checks, 0 failed**. They pin, among others:

- the inverse correction round-trips exactly (`s_obs → s_true → s_obs`);
- `κ_dist` **increases** with τ under the derived couplings — the guard against
  silently reverting to the wrong-sign version that trapped the first draft;
- `κ_chase` falls by exactly 4.0 per unit τ, and the crossing is unique on
  (0, 0.6) for all three couplings;
- all three couplings reconcile, at three different τ, and "same" is the
  pessimistic one;
- δ\* is 48.00 at τ = 0, 57.30 at the default crossing, 46.76 under "same";
- the 90 bp/PD cell needs 5.1426× and **needs a τ above the crossing**;
- the real assay CV is 17.2 % and `n = 12` independent / `n = 3` paired;
- `cv_tolerated` is the exact inverse of `n_per_arm`;
- **WP34's 0.280 is reproducible** and its +PT analogue is **negative** — both
  asserted, so the distinction between the two scales cannot be lost.
