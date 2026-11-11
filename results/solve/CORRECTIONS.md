# Corrections to the first-pass answer

This file exists because one claim I made in `RESULTS.md` was wrong, and because outside input
(balch 2026, the repository's own November status, and a reviewer's objection) exposed it. It is
kept separate from `RESULTS.md` rather than quietly edited into it.

---

## C1. "A no-primer incorporation assay is the decisive experiment" — WITHDRAWN

**What I claimed.** That single-turnover nucleotide incorporation ± TPP1 *in the absence of
primer* would separate retention from allostery, because retention has no route to act with no
substrate present.

**Why it is wrong.** Telomerase is a ribonucleoprotein: nucleotide incorporation requires a
primer base-paired to the RNA template. With no primer there is no substrate for the chemical
step, so there is no k_chem to measure. The assay cannot be run. `wp4_decisive_experiment`'s
decision rule was built on an impossible measurement.

This matters because it was my headline recommendation, and following it would waste an
experimental programme.

**What replaces it.** Primer is necessary, which means DNA occupancy is always present at the
active site, so *occupancy* can never be held fixed by removing the primer. It must be held fixed
experimentally instead. Two measurements do that:

1. **Translocation at matched DNA residence time.** Measure the per-repeat survival *s* and the
   primer dissociation rate in the same samples. The retention model's single claim is that
   `s = k_tr / (k_tr + k_dissoc)` with `k_dissoc` and not `k_tr` moving. If TPP1 changes `k_tr` at
   matched residence time, retention is incomplete. This distinction *is* measurable within the
   published design — Latrick & Cech already report φ = 0.57 → 0.99 and k_tr = 0.07 → 0.10 min⁻¹,
   a 1.43× change in k_tr. My WP3 folded that factor in as one of two admissible weights instead of
   testing whether it has an independent origin. That was the weakest step in the analysis.
2. **Single-turnover incorporation at matched, saturating primer occupancy**, using the reaction's
   own chemistry as the readout, with the residence-time control of (1) measured in parallel. A
   k_chem change at matched occupancy and matched *s* is the allosteric signature.

**Corrected decision rule.**

| outcome | reading |
|---|---|
| *s* changes with `k_dissoc` change and no change in `k_tr` at matched residence | retention-only |
| *s* changes with a `k_tr` change at matched DNA residence time | translocation mechanism, distinct from retention |
| *k_chem* changes at matched occupancy and matched *s* | catalytic allostery |
| all three move | combined; report the decomposition, do not pick one |

Note that row 2 is the one the data already hint at (k_tr 0.07 → 0.10 min⁻¹, 1.43×), and it is the
one my first-pass analysis absorbed into a fitting weight rather than treating as a finding.

---

## C2. "Retention" was used for two different mechanisms — SPLIT

WP1–WP8 used "retention" for both

- **(R1)** TPP1 changing the enzyme's affinity for DNA, and
- **(R2)** TPP1 changing how the TEN domain engages DNA (a conformational effect that *produces*
  retention).

They are different mechanisms with different structural requirements, and only R1 was tested. R1
is excluded: TPP1 is 28.7 Å from the DNA and makes no DNA contact. R2 is not excluded — and the
structure argues for it (WP9):

```
residue   d(TPP1)   d(DNA)    published requirement
A46        3.35     36.46     closest TPP1 contact on TEN
F101       5.81     26.77     FGF loop
G100       7.38     24.34     Zaug 2010: G100V abolishes POT1-TPP1 stimulation
F99        8.81     18.66     FGF loop
Y176      19.21     11.11     Liu 2022: DNA anchor (RAP reduced by Y176A/Q177A)
```

G100 is 7.4 Å from TPP1 but **24.3 Å from the DNA**. It cannot be gripping DNA, yet its mutation
removes TPP1's effect entirely. The functional element for *stimulation* is on the TPP1-proximal
face of TEN; the DNA anchor (PLYQ) is on the opposite face, 19.2 Å from TPP1.

So TPP1 does not hold the DNA, and does not reach the active site. It acts on TEN — and TEN is a
**tethered** domain: residues 197–321 are unmodelled in all three TERT entries, i.e. TEN hangs off
the TERT ring on a 124-residue linker, with only 36 of 160 Cα atoms having a ring Cα within 10 Å.
Its orientation relative to the ring is a hinge degree of freedom.

This is the branch the repository's binary "retention vs allostery" framing cannot represent, and
it is the mechanism the 2026 structures describe.

---

## C3. WP3's 2.48× is a consistency check, not a confirmation — DOWNGRADED

The reviewer's objection is correct in substance even though its stated reason is not. The
prediction is built from Latrick & Cech's own measured φ and k_tr. Agreement with an effect
measured on the same system is therefore **calibration**, not independent confirmation: it shows
the retention model is *sufficient and self-consistent*, not that it is *true*.

Two of the four inputs are also not independent parameters but ratio statistics that partly define
the output. The honest status of WP3 is:

- it **falsifies** eq. 1 (which cannot express unchanged-resting/4×-slower-cycling k_off at all),
- it **bounds** the retention contribution to within the observed effect,
- it does **not** demonstrate retention, because a model calibrated on the same numbers cannot.

The independent-information content of my analysis is in WP2, WP4, WP7, WP8 and WP9 — the
algebra, the geometry and the entry audit — not in the 2.48×.

---

## C4. "One bench measurement settles it" — WRONG

I said the remaining question was one experiment. It is not. The 2026 structures (9SHY/9SHZ/9SI0)
show TEN and IFD-TRAP moving through the RAP cycle *with TPP1 present in all three states*, so they
cannot separate TPP1-driven from cycle-driven motion. And the paper itself reports that the
locked-nucleotide substrates used for the structural work supported minimal RAP activity, so the
deposited states are not the productive cycle.

Resolving TPP1's contribution needs a TPP1-free matched series at the same catalytic states —
which does not exist — or the matched-residence-time measurement of C1. Both are real experimental
programmes, not "one measurement".

---

## C6. "A ≈ 0.99 positively excludes extension-phase retention" — TOO STRONG

**What I claimed (WP11/WP12).** That a baseline channel-A value near 1 is a *positive exclusion* of
extension-phase retention rather than a null result.

**Why it is wrong.** A near-1 value means the assay has little sensitivity to extension-phase
dissociation *in that regime*. It does not establish that the effect is absent. Quantified at 1%
precision on the A-ratio:

| A_ctrl | a = k_off/k_pol | R if k_off falls 4× | upper bound on fold-change from a null | useful? |
|---|---|---|---|---|
| 0.99 | 0.0017 | 1.008 | **99-fold** | no |
| 0.97 | 0.0051 | 1.023 | 1.49-fold | marginal |
| 0.90 | 0.0177 | 1.082 | 1.11-fold | yes |
| 0.70 | 0.0613 | 1.304 | 1.03-fold | yes |

At A = 0.99 a null is consistent with k_off_ext falling 99-fold. **The bound is vacuous.**

**Corrected statement.** A near-1 reading is *uninformative*, not negative. A defensible negative
requires a **narrow confidence interval** at an A value where the assay had power. Non-significance
alone is never a negative result — the repository's preregistration says the same thing
(equivalence testing, not "no significant difference").

**Also corrected:** resting k_off must not be used as a stand-in for k_off_ext. At k_pol = 0.1 the
resting value predicts a 1.65× effect where the extension-phase value gives 1.14× — a design sized
on the resting number is underpowered by an order of magnitude in replicates.

---

## C5. What was right

Recorded so the record is even-handed:

- **WP2's algebra.** `R = k_off²/(k_off + k_on[P])`, n² scaling, and the finding that at
  [P]/K_d = 100 the affinity channel carries 0.99% of the loss budget. Exact, independent of any
  parameter choice.
- **WP2's falsification of eq. 1.** One k_off cannot be both unchanged at rest and ≥4× slower while
  cycling. The model class fails regardless of parameters.
- **WP7.** 7TRD and 7TRF share a TERT model with bit-identical B-factors. The repository has one
  independent apo structure, not two.
- **WP8.** The 7-mer's 3′ end is on the catalytic aspartates and its 5′ end on the anchor strand;
  the exit path runs into Tyr176 at 6.8 Å.
- **WP9.** G100 (required for stimulation, 24 Å from DNA) and Y176 (DNA anchor, 19 Å from TPP1) are
  on opposite faces of TEN. This is the strongest new constraint produced here.
- **The mutant-panel null (WP6).** Nothing correlates; the sign flips on one flagged point.

---

## Summary of the corrected position

| claim | first pass | corrected |
|---|---|---|
| decisive experiment | no-primer single turnover | **impossible**; use translocation at matched residence time |
| retention | "demonstrated" | **R1 not excluded (it is assumed away in the fit) but also not unsupported: the resting `k_off` control (0.012 /min, unchanged by saturating POT1-TPP1) argues against it; R2 supported but not demonstrated** — see WP17; a primer titration cannot separate them |
| 2.48× prediction | quantitative agreement | **calibration, not confirmation** |
| what TPP1 acts on | DNA grip | **TEN orientation, via the FGF-loop face, on a hinged linker** |
| remaining work | one experiment | **a matched TPP1-free structural series + residence-time-matched kinetics** |
