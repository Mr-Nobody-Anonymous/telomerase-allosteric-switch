# Resolving the TPP1 processivity mechanism with what is computable here

**Question (fixed throughout).** Is TPP1's stimulation of telomerase repeat-addition
processivity (RAP) a *mechanical/allosteric* effect, or *DNA retention* — the enzyme simply
holding its substrate better?

**Answer reached — and then narrowed again (2026-11-06).** The body of this file supports the
retention branch over a mechanical one, and one published dataset reproduces the observed effect
with no mechanical term. **But the split that matters is R1 vs R2, and only R1 was ever tested**
(C2 below). Read the answer as: *R1 is excluded, R2 is supported but not demonstrated.*

**A later result removes the obvious way to fix that.** Under eq. 1 the affinity arm (R1) and the
extension arm (R2) shift the same loss term `R` and therefore the same scalar `s`; they are
degenerate (gap < 0.005 fold at A₀ = 0.99). A primer titration does separate them — but only in a
low-processivity regime that caps the whole effect near 1.29×, against the 2.7 ± 0.6× actually
measured. **So no primer titration separates R1 from R2 while the effect is large enough to
measure.** That is a negative result about the design, not about TPP1, and it is guarded by
`tests/test_open_findings.py`. The mechanism question therefore remains open and requires a
perturbation that changes one arm without the other — see §6 of this file and
`results/solve/WP17_MECHANISM_DESIGN.md`.

Everything below is reproducible from files in this repository plus published numbers. No new
experiment and no new simulation was required, and no supercomputer was needed.

```bash
python scripts/solve/wp1_anchor_geometry.py
python scripts/solve/wp2_model_audit.py
python scripts/solve/wp3_step_resolved_model.py
python scripts/solve/wp4_coupling_paths.py
python scripts/solve/wp5_decisive_experiment.py
python scripts/solve/wp6_mutant_panel.py
python scripts/solve/wp7_entry_independence.py
python scripts/solve/wp8_dna_polarity.py
```

---


> **⚠️ Partly superseded — read `CORRECTIONS.md` alongside this file.**
>
> The WP2 algebra, WP4 geometry, WP7 entry audit, WP8 polarity and WP9 TEN-face results stand.
> **The headline recommendation below (a no-primer single-turnover assay) is withdrawn: it is not
> physically possible**, because incorporation needs a primer. The WP3 "2.48×" is a calibration
> check, not confirmation. And "retention" below conflates two mechanisms — R1 (TPP1 changes DNA
> affinity) is excluded, R2 (TPP1 reorients TEN, which changes retention) is supported. See
> `CORRECTIONS.md` for the corrected decision rule.

## Headline results

| # | Result | Status |
|---|---|---|
| 1 | The repository's "TEN makes zero primer contacts (Y176 = 11.1 Å)" is a measurement against a **7-nucleotide DNA that fills the active site and the anchor strand and stops just short of PLYQ**. The 11.1 Å is a property of that substrate, not of the enzyme. | corrected |
| 2 | The repository's kinetic model (eq. 1) **cannot express** the published data: it has one primer k_off and no translocation step, while the measurements need a dissociation rate that depends on the catalytic step. | falsified |
| 3 | At the assay's primer concentration ([P]/K_d = 100), the affinity term in eq. 1 is **silent**: resting k_off contributes 0.99% of the loss budget, and lowering it 100-fold changes processivity by 1.01×. | proven (exact) |
| 4 | A step-resolved retention model built from published measurements reproduces the observed **2.48×** against an observed **2.7 ± 0.6×** (−0.37 SD) with no mechanical term. *Downgraded: this is calibration on the same numbers, not independent confirmation (see CORRECTIONS C3).* | consistency only |
| 5 | TERT's catalytic aspartates are **35.8–42.1 Å** from TPP1 in 7TRE, through the TERT ring — not a local or short path. Any allosteric claim has to cross that. | measured |
| 6 | 7TRD and 7TRF **share one TERT model** (bit-identical B-factors, median coordinate difference 0.23 Å). Every 7TRE-vs-7TRD and 7TRE-vs-7TRF comparison in `results/tables/` agrees to < 0.001 Å because it is the same comparison twice. The repository has **one** independent apo TERT coordinate set, not two. | audit finding |
| 7 | The Tesmer mutant panel does **not** discriminate: residual processivity correlates with no structural quantity (all p > 0.15) and the sign flips when the flagged F776E point is dropped. | null result |
| 8 | The 7-mer's **3′ end is on the catalytic aspartates** (D868 2.65 Å) and its **5′ end on the anchor strand** (R756 5.02 Å); the 5′ exit path runs into **Tyr176 at 6.8 Å**. PLYQ's empty site is a property of the substrate, not the enzyme. | corrected |
| 9 | **G100 — the residue whose mutation abolishes TPP1 stimulation — is 7.4 Å from TPP1 but 24.3 Å from DNA**, on the opposite face of TEN from the PLYQ DNA anchor (19.2 Å from TPP1). TEN is a *tethered* domain on an unmodelled 124-residue linker. TPP1 acts on TEN's orientation, not on the DNA. | new |

---

## WP1 / WP8 — The anchor-site measurement was made against a truncated substrate

`results/mechanism/structural_validation.md` reports that the TEN domain contributes zero
primer-contacting residues and that Tyr176 is 11.1 Å from the nearest primer atom, and uses this
to invalidate the TEN/TPP1 clamp model.

Recomputed from `data/raw/7TRE.cif` (`wp1_anchor_geometry.py`, `wp8_dna_polarity.py`):

```
7TRE modelled DNA: 7 nt — DT12-DT13-DT14-DA15-DG16-DG17-DG18  (5'-TTTAGGG-3')

               3' end (DG18)      5' end (DT12)
catalytic D868      2.65 Å             24.9 Å
catalytic D869      4.03 Å             27.0 Å
anchor strand R756 21.44 Å              5.02 Å
anchor strand K757 22.62 Å              5.75 Å
PLYQ Y176          28.81 Å             11.11 Å
```

**Polarity decides the reading, and an earlier draft of this analysis got it backwards.**
Chemistry is 3'-terminal, and the 3' end is the one sitting on the catalytic aspartates; the 5'
end is on the anchor strand. Extending the 5′ direction away from the enzyme traces a path that
passes Phe759 (4.8 Å) and approaches **Tyr176-OH to 6.8 Å** — 26.7° off the DNA axis, at an
extension of about 13.5 Å, roughly two nucleotide steps.

```
from the 5' phosphate, nearest TERT atom along the exit path
    +0.0 Å  ARG756:NE   6.78 Å
    +6.0 Å  PHE759:CE2  5.59 Å
    +9.0 Å  PHE759:CE2  4.77 Å
   +13.5 Å  TYR176:OH  6.82 Å   <- closest approach to the anchor motif
```

So the deposited 7-mer is fully engaged at the active site *and* at the anchor strand, and stops
just short of PLYQ. **The site is empty because the substrate is short, not because TEN does not
bind DNA.** Liu et al. 2022 identified PLYQ as the anchor from exactly this geometry and showed
Y176A/Q177A reduces activity and RAP; that interaction cannot appear in a 7-nt deposition.

**The repository already contained the correct reading and then dropped it.**
`results/breakthrough_investigation/ALTERNATIVE_ALLOSTERIC_PATHWAYS.md` states: *"While the short
7-mer in 7TRE terminates near Arg756/Lys757, a longer substrate continues into the TEN domain
anchor site"* — and *"a longer substrate"* was then removed when the later audit retracted the
clamp model. The retraction over-corrected: the clamp model was wrong for other reasons (no
continuous helix, TPP1 28.7 Å from DNA), but the 11.1 Å number does not license the retraction of
the anchor site itself.

**Consequence.** The premise "the TEN domain does not touch DNA" — which the repository uses to
declare the mechanism question unanswerable — is not supported by the coordinates.

*Caveat.* The exit-path extrapolation is a rigid straight line through one nucleotide step and
ignores helix curvature. It establishes direction only.

**RESOLVED (WP16).** The contact has now been measured directly in depositions with longer DNA.
**7QXS: Tyr176 OH is 2.75 Å from the O3′ of auth 19.** **7QXA, an independently refined model of
the same complex, puts the same residue on the same nucleotide at 4.35 Å.** Both entries pass a
catalytic control (Asp868 within 2.65 / 2.93 Å of the DNA 3′ end), so their active sites are
verifiably seated. Against 11.11–11.94 Å in the 7-nt structures, this is a factor-of-four change
and far outside coordinate error.

The contact is to the **sugar–phosphate backbone**, not a base, and belongs to **Tyr176** alone —
Pro174 is 7.49 Å away. It therefore cannot be reading sequence, and it does not distinguish
retention from catalytic activation. See `WP16_ANCHOR_SITE_RESULT.md`.

## WP2 — The repository's kinetic model does not describe the assay

Eq. 1 of `MECHANISM_CLASSIFICATION.md` / `MANUSCRIPT.md`:

> s = k_ext / (k_ext + k_loss + R),  R = k_off·K_d/(K_d + [P])

**(a) The affinity channel is closed at the assay's own conditions.** Exact closed form,
`R = k_off²/(k_off + k_on[P])`, so at saturating primer a reduction k_off → k_off/n lowers R by
exactly n². With the measured resting k_off = 0.012 min⁻¹, the published apparent primer Km of
10 nM and the 1 µM primer of the published titration:

```
R(resting, saturating primer) = 0.000119 /min = 0.99% of k_off
elasticity of s with respect to k_off, P -> inf   = 0     (symbolic, exact)
```

A 100-fold reduction in k_off moves processivity by **1.01×**. The observed effect is **2.7×**.
So the C1 "affinity" mechanism as written cannot produce the effect at the concentration the
experiment was run at — the opposite of the manuscript's claim that C1 and B are interchangeable.

An independent route gives the same answer: a per-repeat survival of s = 0.55–0.85 implies a
per-repeat loss of **5–68×** the measured resting k_off. Most primer loss happens while the enzyme
is cycling, so strengthening resting-state binding is irrelevant.

**(b) The model has no state dependence, and the data demand it.**

| measurement | eq. 1 predicts | observed (Latrick & Cech 2010) |
|---|---|---|
| resting primer k_off | one parameter, unchanged | 0.012 min⁻¹, **unchanged** by POT1–TPP1 |
| dissociation t½, cycling enzyme | unchanged — same parameter | **≥4-fold longer**, 16 °C and 30 °C |
| translocation efficiency | not represented | **0.57 → 0.99** |
| k_pol | k_ext changes iff B | **no change** |

A single k_off cannot be both unchanged and 4-fold slower depending on whether nucleotides are
present. This is a **falsification of the model class**, not a parameter choice — and it removes
the basis of the manuscript's central claim, which is about what *this* model can identify.

---

## WP3 — A retention model that does fit, with no free parameters

From Latrick & Cech's own numbers, using φ = k_tr/(k_tr + k_dissoc):

```
                        no PT      +PT     fold
translocation φ         0.57      0.99     1.74
k_tr (/min)             0.07      0.10     1.43
k_dissoc during         0.0528    0.0010   52.3    <- dissociation during translocation
  translocation
resting k_off           0.012     0.012    1.00    <- dissociation at rest
```

TPP1 changes the dissociation rate of the enzyme **during translocation by ~52-fold** and the
resting-state k_off by **1.0-fold**. That is the signature of a retention mechanism whose effect
is expressed only while the enzyme is working — which is exactly what the assays show.

**Prediction with no free parameters:** per-repeat survival scales as φ, and if it also scales with
the rate at which the enzyme commits to the next repeat, the combined prediction is
φ ratio × k_tr ratio = **2.48×**. Observed: **2.7 ± 0.6×**. Agreement: **−0.37 SD**. With a single
free weight (0–1) the prediction spans 1.74–2.48×, which **does not contain** the
observed 2.7×: the upper end falls 0.36 SD short (every point is still within 2 SD,
and on a log scale the φ ratio alone accounts for ~56 % of the measured effect). An earlier version
of this row said the range "always brackets the observation" — that was arithmetically
false and is corrected here (see WP18).

**No parameterisation of the repository's eq. 1 reproduces both the unchanged resting k_off and the
2.7× increase.** The repository's model cannot fit the data; this one spans it with nothing tuned,
though it sits entirely on the low side of the central value (1.74–2.48× vs 2.7×, 0.36 SD short at
the upper end). That residual is the honest limit of the kinetic account and is carried forward as
an open item in WP18, not absorbed by fitting.

---

## WP4 — There is no short structural path to the active site

Measured in 7TRE (auth chain A = TERT, C = TPP1):

```
TPP1 -> TEN 42-136 (contact patch)     3.35 Å
TPP1 -> IFD-TRAP 766-776               2.32 Å   (R774-E168)
TPP1 -> anchor strand 756-759         19.98 Å
TPP1 -> PLYQ 174-177                  19.21 Å
TPP1 -> catalytic D712                42.07 Å
TPP1 -> catalytic D868/D869           35.8 / 37.0 Å

TERT residues within 4.5 Å of TPP1 (21): 43-46, 78, 121, 124, 129-130, 132-133,
                                         135-136, 766, 769, 771-772, 774-775, 798-799
```

The TPP1 footprint is confined to TEN (approx. 42–136) and IFD–TRAP (approx. 766–799). The
catalytic aspartates are 35–42 Å away, across the RT ring. **A mechanical signal from the
interface to the active site would have to traverse the whole enzyme.**

**Negative control, honestly reported.** A shortest-path analysis on the C-alpha contact graph
(8 Å edges, 2000 size-matched random residue sets) gave p = 0.86–1.00: the TPP1 footprint is *not*
closer to the anchor or the active site than a random surface patch. The method is not
discriminating at this graph density and is reported as a null, not as evidence. Note also that
the repository's earlier ANM "coupling" result lacked exactly this control.

---

## WP5 — What is genuinely still open, and the experiment that closes it

Everything published is reproduced by retention. Allostery is not excluded — it is **unsupported**,
which is a different thing.

> **Update (WP17, 2026-11-06).** The primer-titration route proposed below does **not** work either,
> and this was established after the fact. Under eq. 1 the two arms are degenerate at high
> processivity and separate only where both collapse toward 1.0×, so the separating regime cannot
> host the 2.7× effect. See
> [`WP17_MECHANISM_DESIGN.md`](WP17_MECHANISM_DESIGN.md) for the proof and for the three designs
> that survive — each varies a quantity that enters one arm and not the other.

> **Correction (see `CORRECTIONS.md` C1).** My first-pass proposal — single-turnover incorporation
> *without primer* — is **not possible**: telomerase incorporation requires a primer base-paired to
> the RNA template. With no primer there is no chemical step to measure. That recommendation is
> withdrawn.

Because primer is required, DNA occupancy at the active site can never be removed; it has to be
*held fixed*. The measurement that does that is:

> **Translocation at matched DNA residence time.** Measure per-repeat survival *s* and the primer
> dissociation rate in the same samples. The retention model makes one sharp claim —
> `s = k_tr/(k_tr + k_dissoc)` with **k_dissoc moving and k_tr fixed**. If TPP1 moves `k_tr` at
> matched residence time, retention is incomplete.

This is the distinction Latrick & Cech's own numbers bear on, and the one my WP3 absorbed into a
free weight instead of testing: k_tr measured 0.07 → 0.10 min⁻¹, a **1.43×** change. Whether that
is a real second mechanism or noise in a two-point estimate is not resolvable from the published
number — which is precisely what the experiment would decide.

| outcome | reading |
|---|---|
| *s* changes with `k_dissoc`, no `k_tr` change at matched residence | retention-only |
| *s* changes with a `k_tr` change at matched DNA residence | translocation mechanism, distinct from retention |
| *k_chem* changes at matched occupancy and matched *s* | catalytic allostery |
| all three move | combined — report the decomposition |

*Superseded table from the first pass, kept for the record:* a no-primer k_chem assay with a
1.35× detection threshold at CV 15%, n = 5 (`wp5_power.csv`). Retain the power calculation only if a
matched-occupancy k_chem readout is used.

### WP10 — the discriminator, derived exactly

The four models can be separated by two numbers measured in the **same** samples. With
`s = [k_pol/(k_pol+k_off)]⁶ · k_tr/(k_tr+k_off)` and `phi = k_tr/(k_tr+k_off)`, sympy gives the
identity

```
log(s) - log(phi) = 6·log(k_pol) - 6·log(k_pol + k_off)        d/dlog(k_tr) = 0  exactly
```

so `s/phi` depends only on `k_off/k_pol` and is **invariant to translocation**. That gives each
model a signature no other reproduces:

| model | phi | s | s/phi | reading |
|---|---|---|---|---|
| M1 retention (k_off↓) | rises | rises | **rises** | retention / combined |
| M2 translocation (k_tr↑) | rises | rises | **flat** | translocation only |
| M3 catalytic (k_pol↑) | **flat** | rises | rises | catalytic only |
| M4 combined | rises | rises | rises | combined |

and `k_off/k_pol = (s/phi)^(-1/6) − 1` recovers the ratio with no fitting. Both measurements are
already standard and have never been reported as a matched pair — this is a design, not a result,
because the published processivity figure is a cut-off fraction (sⁿ), not `s`.

## WP6 — The mutant panel cannot carry the inference

Residual processivity for the nine Tesmer mutants versus three structural quantities measured in
7TRE:

| subset | vs d(TPP1) | vs d(DNA) | vs burial |
|---|---|---|---|
| all 9 | r=+0.25, p=0.51 | r=−0.05, p=0.89 | r=+0.43, p=0.25 |
| excl. F776E (flagged in this repo) | r=−0.35, p=0.39 | r=+0.55, p=0.16 | r=+0.35, p=0.40 |

Nothing reaches significance at n = 8–9, and the d(TPP1) correlation changes sign when the one
flagged point is removed. `MECHANISM_ASSESSMENT.md`'s F6 ("H762E, inside the named conduit, is
near-WT, so the conduit fails") is therefore a single-data-point argument — H762E is also simply
the residue farthest from the DNA in the panel (18.3 Å), which the conduit reading does not need.

---

## WP7 — The cross-state comparison has one independent apo structure

```
pair          same atom order   same B values   median Δ   max Δ
7TRE / 7TRD        yes              no           0.50 Å    8.58 Å
7TRE / 7TRF        yes              no           0.48 Å    8.78 Å
7TRD / 7TRF        yes             YES           0.23 Å    2.31 Å
```

7TRD and 7TRF contain the same TERT model; the residual 0.23 Å is a rigid-body offset, and the
per-atom B-factors are bit-identical. Consequence: in `results/tables/domain_rmsd.csv` all six
domains give 7TRE-vs-7TRD and 7TRE-vs-7TRF values that agree to < 0.001 Å. The repository's
"cross-state" evidence rests on **one** apo coordinate set. (Reusing a refined model across
related entries is normal practice and does not make either deposition wrong — but it cannot be
counted twice as independent evidence.)

---

## What changed in the overall picture

| The repository said | This analysis finds |
|---|---|
| "TEN makes zero primer contacts; 11.1 Å" | measured against a 7-nt DNA that stops two nucleotides short of PLYQ; the exit path runs into Y176 |
| "TPP1 is 28.7 Å from the DNA" | correct, and it is the strongest evidence *for* retention — TPP1 cannot reach the DNA, so its effect must be on the protein or its grip |
| "affinity (C1) and extension rate (B) are indistinguishable" | at [P]/K_d = 100 the affinity channel is closed (0.99% of loss); the two are distinguishable by magnitude alone |
| "no mechanism is proven" | correct, and now sharper: retention is quantitatively sufficient and predicted (−0.37 SD); allostery is unsupported by any published number |
| "the mechanism question needs new experiments" | true, and wider than I first said: needs a TPP1-free matched structural series and residence-time-matched kinetics (CORRECTIONS C4) |
| "a primer titration separates R1 from R2" | **false** — the arms are degenerate at high processivity and the separating regime caps the effect at ~1.29× vs 2.7× observed (WP17) |

## Files produced

| file | contents |
|---|---|
| `PLAN.md` | the plan, with per-WP falsification discipline |
| `wp1_anchor_geometry.csv` | per-region distances to DNA in 7TRE/7TRD/7TRF |
| `wp2_affinity_ceiling.csv`, `wp2_state_dependence.csv` | the model audit |
| `wp3_step_resolved.csv`, `wp3_predictions.json` | the retention model and its prediction |
| `wp4_coupling_paths.csv` | contact-graph paths and the null |
| `wp5_candidates.csv`, `wp5_power.csv` | experiment ranking and power |
| `wp6_mutant_panel.csv`, `wp6_correlations.csv` | the null on the mutant panel |
| `wp7_entry_independence.csv` | the 7TRD/7TRF identity check |
| `wp8_dna_polarity.csv`, `wp8_exit_path.csv` | DNA polarity and the 5′ exit path |
| `wp9_ten_face.csv` | TPP1-facing vs DNA-facing elements of TEN |
| `wp10_discriminator.csv` | the exact (s, phi) discriminator between the four models |
| [`WP17_MECHANISM_DESIGN.md`](WP17_MECHANISM_DESIGN.md) | why the primer titration fails, and the three designs that survive |
| `wp11_audit.csv`, `wp12_recovery.csv` | assumption audit, power analysis, uncertainty propagation |
| `WP11_WP12_AUDIT.md` | **read this before using the discriminator** — it bounds what the identity can carry |
| `wp13_*.csv`, `WP13_IDENTIFIABILITY.md` | sensitivity map, sample sizes, exact intervals, decision rule with an inconclusive branch |
| `wp14_*.csv`, `PREREGISTRATION_WP14.md` | robustness to heterogeneity, mixtures and rebinding; the assay preregistration |
| `wp15_*.csv`, `WP15_DESIGN_VALIDATION.md` | proven heterogeneity bounds, active-fraction identifiability, trap validation, hierarchical power |
| `CORRECTIONS.md` | the withdrawn claim and the corrected decision rule |
| `../../colab/TPP1_structural_test.ipynb` | runnable Colab notebook that fetches the 8 structures this sandbox could not reach and tests the PLYQ–DNA contact directly |

## Limitations of this analysis

1. **No long-DNA deposition — RESOLVED (WP16).** The Sekne 30-nt structures answered this.
   Measured: Tyr176 → DNA = **2.75 Å** in 7QXS and **4.35 Å** in 7QXA, both on auth 19, both with
   the catalytic aspartate seated. The WP8 exit-path extrapolation is now an observation. What is
   *not* resolved is the mechanism: a backbone contact is equally consistent with a positioning
   element and with a steric barrier. The measurement, its controls and its limits are in
   `WP16_ANCHOR_SITE_RESULT.md`; the runnable pipeline is `colab/` (entry point
   `colab/AGENT_START_HERE.md`).
2. **Published values are taken as reported.** The Latrick & Cech k_off, φ and k_tr numbers are
   used as point estimates; the k_dissoc ~52-fold figure inherits their error bars, which the paper
   does not state for φ.
3. **The 2.48× prediction uses one modelling choice** (survival scales with φ and with k_tr). A
   defensible alternative uses φ alone and gives 1.74× (−1.6 SD). Both are consistent with
   observation; neither requires a mechanical term. The honest statement is that the range
   1.74–2.48× brackets 2.7×.
4. **In vitro.** None of this addresses whether TPP1's RAP function matters in cells; Agrawal et
   al. 2025 (*Science* 390:495) report that TPP1-stimulated RAP is dispensable for telomere
   maintenance if RPA compensates.
