# Preregistration — A2: is the retry memoryless?

**Registered before execution.**

- **Test ID:** A2
- **Gap class:** G2 (measurement) — issue [#16](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/16)
- **Informational value:** EIG **0.3295 bits**, rank 3 of 5 (WP46)
- **Decides:** the validity of the **stall renewal itself** — the load-bearing
  structure of every work package since WP34

---

## 1. The claim under test

WP44 introduced, and every later package used, the stall renewal:

```
s_obs = s_true / (1 − τ (1 − s_true))
```

This is **itself geometric only if the retry is memoryless** — that is, if a
molecule that stalls and retries has the *same* success probability as one that
never stalled. Nothing in the programme has checked this.

If the retry is instead **state-dependent** (a molecule that just stalled is more
likely to stall again, or has been given a fresh processive run), then the renewal
is wrong, and **every `R½` derived since WP34 is an effective value masquerading as
a true one.**

## 2. Prediction, stated numerically in advance

Under a stall-inducing perturbation, `s_obs` as a function of `s_true` must follow
the renewal curve:

| `s_true` | predicted `s_obs` at τ = 0.1833 |
|---|---|
| 0.600 | 0.6475 |
| 0.700 | 0.7407 |
| 0.800 | 0.8304 |
| 0.917 | 0.9312 |
| 0.950 | 0.9588 |

**If memoryless:** the observed points lie on this curve, and `τ` fitted from them
is *independent of* `s_true`.
**If not:** the curve saturates — observed `s_obs` falls below prediction at low
`s_true` and the fitted τ drifts with `s_true`.

**The discriminating statistic: the constancy of fitted τ across the `s_true`
range.**

## 3. Design

- **System:** a telomerase preparation where `s_true` can be varied
  independently — e.g. primer length or sequence series, or a processivity mutant
  panel — **crossed with** a stall-inducing lesion.
- **Grid:** ≥4 levels of `s_true` × 2 levels of stall (present/absent), n ≥ 3.
- The no-stall arm gives `s_true` directly; the stall arm gives `s_obs`.
- **Both arms in the same session** — day-to-day variance would otherwise be read
  as drift.

## 4. Stopping rule

Fixed grid, n = 3, no optional continuation. Any dropped well is reported.

## 5. Analysis, fixed in advance

1. Fit τ at each `s_true` level by inverting the renewal.
2. **Primary statistic:** the slope of τ against `s_true`.
3. **PASS (memoryless):** |slope| < 0.15 and CoV(τ) < 0.20.
4. **FAIL:** slope outside that, or observed `s_obs` below prediction by > 2 SE at
   two or more points.
5. **AMBIGUOUS:** reported as ambiguous. Not a pass.

## 6. Outcomes

| outcome | consequence |
|---|---|
| **PASS** | The renewal stands. `R½` values may be called true values. |
| **FAIL** | **The stall extension is wrong.** Every `R½` since WP34 becomes an *effective* value, the true/apparent distinction in WP45 §3 changes meaning, and R13 must be re-derived — the 6× discrepancy may move in either direction. Reported as a structural failure, not a parameter change. |
| **AMBIGUOUS** | The renewal stays **unvalidated**, and no result resting on it may be called a finding. |

## 7. Why this ranks third and not first

A1 (EIG 0.4547) outranks A2 (0.3295) because A1 decides whether the *refutation is
real*, while A2 decides whether the *vocabulary* is real. If A1 fails, R12/R13 are
withdrawn and A2's answer applies to a framework nobody is using.

**But A2 is the cheaper of the two to run in practice**, and if only one experiment
is ever funded, the ordering above is the argument for running A1.

## 8. Independent replication slot

One full grid at a **second site**, with the analysis script run by someone who
did not write this document. Booked before the first run.

If none is available, that limitation is stated **up front** — the reproduction is
then "same-lab only" and every claim resting on A2 carries that qualifier. It is
not a detail to discover at submission.

## 9. What this cannot settle

- Does not test **A1**.
- Does not distinguish **which** mechanism causes the stall — only whether the
  renewal describes it.

---

*Scanned by `tests/test_open_findings.py`; removing the prediction table breaks the
build.*
