# WP19 — Closing gaps with measurements that already exist

WP18 left one item open: an absolute per-repeat survival *s* with a matched φ,
to test the identity `s/φ = (1+a)^−6` against an independently measured
off-rate. **That measurement exists.** WP19 encodes it and, in doing so,
retracts one of the repository's headline claims.

    python scripts/solve/wp19_literature_closure.py
    -> results/solve/wp19_literature_closure.csv

## What was found, verified against primary text

**1. The cycling off-rate *has* been measured.** Latrick & Cech 2010
(EMBO J 29:924, PMC2837173), pulse-chase, exact wording:

> "With actively cycling enzyme, POT1–TPP1 clearly slowed the dissociation of
> the ensemble of extending primers … the time required for half-dissociation
> was increased at least four-fold at 16° or 30 °C."

The claim *"nobody has measured primer k_off ± TPP1 for cycling enzyme"* is
**false and is retracted.** README and `MECHANISM_ASSESSMENT.md` are corrected.

**2. k_pol does not change.** Same paper: "we observed no increase k_pol in the
presence of POT1–TPP1." **M3 (catalytic) is now excluded by measurement**, not
only by the 35.8 Å reach argument.

**3. WP18 overstated the degeneracy.** φ (A = 0.57 → 0.99) and k_tr
(k = 0.07 → 0.10 /min) are fitted as *separate* parameters of one curve. WP18
claimed "only the combination is observable" — wrong, k_tr is reported. They
are **correlated through one fit**, not degenerate. C2 is reinstated.

**4. The absolute anchor.** Wang et al 2007 (Nature 445:559) report R₁/₂ — the
number of repeats before half the chains dissociate — as **0.78 → 3.3**. Using
their own relation R₁/₂ = −ln2/ln(s):

| | no POT1-TPP1 | + POT1-TPP1 |
|---|---|---|
| R₁/₂ (repeats) | 0.78 | 3.3 |
| **s = 2^(−1/R₁/₂)** | **0.4112** | **0.8105** |
| measured s ratio | | **1.971×** |

Cross-check from L&C's own text (≥ half reach 8 repeats with PT, < 10 % without):
s ≥ 0.917 and s ≤ 0.750. Same direction, same order — consistent, different
primers, so treated as a bound.

## The result that matters

Inverting the identity on real numbers (`a = (s/φ)^(−1/6) − 1`):

| | no POT1-TPP1 | + POT1-TPP1 |
|---|---|---|
| a = k_off_ext / k_pol | 0.0559 | 0.0339 |

With k_pol **measured unchanged**, that is an extension-phase off-rate fall of
**1.65×**, while the resting k_off stays at 0.012 /min.

But the pulse-chase measures the cycling t½ as **≥ 4×** slower. Those disagree:

    pulse-chase ...................... >= 4.00x
    identity + absolute R_1/2 ........    1.65x
    disagreement .....................    2.42x

**This is the quantitative signature of rebinding**, exactly as WP13 predicted:
a primer that dissociates can reassociate before the next extension, inflating
the apparent bound fraction. R₁/₂ counts only committed extensions, so it is
rebinding-free.

### New, falsifiable prediction

> Under dilution or a DNA trap (rebinding abolished), the pulse-chase t½
> fold-change should fall from **≥ 4× toward ~1.7×**, a factor of ~2.4.
> **If it stays ≥ 4× under trap, the model is falsified.**

This is the first claim in the repository that a single existing assay type
could kill.

## M1 revisited

TPP1 tightens POT1–DNA binding 3.5× (K_d 26 → 7.4 nM, Wang 2007), so global
affinity *does* change — on the **DNA/POT1 axis**. It does **not** change on
TERT–primer `k_off`, which L&C measured as unchanged. M1 in its strong form
(global TERT k_off) stays dead; a DNA-affinity component is real and retained.
The mechanism remains **state-specific**: it acts on cycling enzyme.

## Circular-frame bug fixed in `src/`

`compute_domain_rmsds()` superposed on the catalytic core 326–1132 while
reporting TRBD, RT, IFD-TRAP and CTE — each 100 % inside that range. The frame
is now the invariant RT core. The published table changed materially:

| domain | old (circular) | new (corrected) | change |
|---|---|---|---|
| TEN | 0.4341 | **0.4719** | +8.7 % |
| TRBD | 0.2941 | 0.3459 | +17.6 % |
| RT | 0.3378 | 0.2963 | −12.3 % |
| IFD-TRAP | 0.4268 | 0.3429 | −19.7 % |
| **CTE** | 0.3746 | **0.5606** | **+49.6 %** |

The CTE result is the substantive one: **CTE moves almost twice as much as
previously reported**, and it is the domain physically adjacent to TPP1. TEN,
TRBD, IFD-TRAP and CTE are now 0 % inside the fit frame.

## What is still genuinely open

| Item | Status |
|---|---|
| M2 vs M4 — does `k_tr` itself rise? | k_tr = 0.10 /min is reported, so M2 is supported; A and k are fitted from one curve and their covariance is unpublished, so M4 is not excluded. Needs the raw time-course, not a new experiment. |
| Rebinding correction | quantified by running the pulse-chase under DNA trap/dilution (WP19 prediction). |
| RAP-state superposition 9SHZ / 9SI0 | coordinates are retrievable one residue at a time from models.rcsb.org; not yet ingested. |
| 10 vs 8.4 repeats | dropped — 8 × 6 nt = 48 nt exactly, so L&C's count is self-consistent. |
