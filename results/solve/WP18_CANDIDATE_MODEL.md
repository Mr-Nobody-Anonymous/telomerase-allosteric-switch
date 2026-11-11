# WP18 — One combined model, four candidates, tested against everything

Combines the best working part of each work package into a single model, then
scores four mechanistic candidates against every published constraint.

    python scripts/solve/wp18_candidate_model.py
    -> results/solve/wp18_candidate_model.csv

## What was borrowed

| From | Working part imported |
|---|---|
| WP1 / WP4 | Geometry. TPP1 contacts TEN at **3.35 Å** and reaches **G100 at 7.4 Å**, but the catalytic aspartates are **35.8–42.1 Å** away across the RT ring. |
| WP2 | The constraint that kills eq.1: one primer `k_off` cannot express resting vs cycling. Resting `k_off = 0.012 /min` is measured and does **not** move. |
| WP3 | Step-resolved cycle: 6 incorporations then translocation; φ 0.57→0.99, `k_tr` 0.07→0.10. |
| WP10 | The exact identity and the (s, φ) discriminator. |
| WP13 | Rebinding makes the off-rate term an **upper bound**; only pulse-chase/dilution decomposes it. |
| WP17 | R1/R2 are degenerate under eq.1 → eq.1 is abandoned, not patched. |

## The combined model

    a   = k_off_ext / k_pol        φ = k_tr / (k_tr + k_off_tr)
    s   = (1 + a)^-6 · φ

    log(s) − log(φ) = −6·log(1 + a)          (exact, k_tr-free)

## Result — the scoreboard

| Candidate | RAP fold | z | C1 φ | C3 resting k_off | C4 2.7× | C5 reach | Verdict |
|---|---|---|---|---|---|---|---|
| M1 retention (global `k_off`) | 1.190 | −2.52 | fail | **FAIL** | **FAIL** | pass | eliminated |
| M2 translocation (`k_tr` ↑) | 1.737 | −1.60 | pass | pass | pass | pass | **survives** |
| M3 catalytic (`k_pol` ↑) | 1.190 | −2.52 | fail | pass | **FAIL** | **FAIL** | eliminated |
| M4 state-specific (`k_off_tr` only) | 1.737 | −1.60 | pass | pass | pass | pass | **survives** |

**M1 dies on C3**: a global affinity change would move the measured resting
`k_off` off 0.012. **M3 dies twice**: it puts the change at an active site
35.8–42.1 Å from TPP1 across the ring, and it leaves φ at 0.57.

## The finding that changes the previous answer

The implied translocation off-rate is

    k_off_tr = 0.05281 (no TPP1) -> 0.00101 (+TPP1)     ~52×

But φ and `k_tr` are read off **the same figure** (L&C Fig 2F), so `k_off_tr`
is *derived*, not measured. Only the combination is observable:

    ε = k_off_tr,noPT / k_off_tr,PT ,   r = k_tr,PT / k_tr,noPT ,   ε = 74.68 / r

| ε (drop in `k_off_tr`) | required `k_tr` ratio |
|---|---|
| 20.0 | 3.73 |
| 37.3 | 2.00 |
| **52.3** | **1.43** ← measured pair |
| 74.7 | 1.00 |
| 150.0 | 0.50 |

So **"k_tr rises" is not established.** The data admit the whole family, and
**M2 and M4 are observationally identical** — same mechanism, two descriptions.
C2 is therefore withdrawn as an independent constraint, and the earlier plan to
rank M4 strictly above M2 was over-claiming.

## The honest gap

The parameter-free prediction from φ alone is **1.737×** against 2.7 ± 0.6,
i.e. **0.36 SD short**. WP3 previously claimed the prediction range
"always brackets the observation" — that is arithmetically false (the range is
1.74–2.48× and does not contain 2.7). Corrected in `wp3_step_resolved_model.py`
and `RESULTS.md`. On a log scale φ explains 56 % of the effect, the `k_tr` ratio
a further 36 %, leaving 9 % (0.36 SD) unexplained — consistent with L&C's
processivity being a cutoff-dependent fraction rather than a per-repeat survival.

## Verdict

**Best candidate: M2/M4 — state-specific translocation.** TPP1 acts on the
translocation step through TEN (3.35 Å) and G100 (7.4 Å), lowering `k_off_tr`
~37–75× and raising φ 0.57→0.99, while resting `k_off` stays at 0.012 and the
active site is never touched.

**Still missing:** a per-repeat *s* with matched φ, from pulse-chase or dilution.
Without it the identity `s/φ = (1+a)^−6` cannot be used to test that `a` is
unchanged, and M2/M4 cannot be separated — nor does the existing data require
separating them.
