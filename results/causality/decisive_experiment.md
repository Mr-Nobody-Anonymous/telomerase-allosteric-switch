# The smallest decisive experiment — primer titration ± TPP1

This document exists because the previous one produced a negative result that is more useful than
a positive one would have been: **the product-length distribution cannot separate affinity from
conformational regulation** (`kinetic_model_comparison.md` §2, G1 ≡ G5, an algebraic identity).
A different observable is therefore required, and this is the specification of the smallest one.

---

## 1. The experiment

**Primer-concentration titration, with and without a fixed saturating concentration of POT1–TPP1
(or TPT), measuring processivity at each primer concentration.**

Minimum design: **6 primer concentrations × 2 TPP1 conditions × 3 independent experiments**.
Nothing else changes — same telomerase preparation, same buffer, same incubation, same gel.

Why this and not something larger: it is a one-dimensional titration on an assay the field already
runs, and it isolates the parameter that the two models treat differently.

## 2. Why this observable and not the distribution

Both models describe the same per-repeat survival

> *s* = *k_ext* / (*k_ext* + *k_off*)

and therefore the **same distribution** for any matched value of *s*. What they disagree about is
**which term TPP1 changes**, and that disagreement is only visible when *s* is pushed away from
saturation:

| | apparent primer *K*m | plateau (maximum *s*) | mechanism |
|---|---|---|---|
| **C — affinity** | **decreases** with TPP1 | unchanged (= *k_ext* is unchanged) | TPP1 lowers *k_off*, so the enzyme needs less primer to be retained |
| **B — conformational** | unchanged | **increases** with TPP1 | TPP1 raises *k_ext*, so even saturating primer gives a higher ceiling |
| A — conduit | n/a | n/a | already rejected; included only to show the panel would carry it |

Both models raise processivity at a *single* subsaturating primer concentration. **Only the
titration separates them**, because only the titration asks whether TPP1 changed the *affinity* for
the primer or the *ceiling* of the extension rate.

Panel 3 of `figures/figure_causality_summary.png` shows the two predicted curves, labelled as a
schematic (the curves are illustrative, not fitted, and the axis is arbitrary units).

## 3. Expected outcomes, pre-declared

Analysis: fit *s*([primer]) to a hyperbola, *s* = *s*_max · [P]/(*K*_app + [P]), separately with
and without TPP1, by nonlinear least squares with the **experiment** as the unit (n = 3, ratio
parameters computed within experiment). Report *K*_app and *s*_max with 95% CIs.

| outcome | Δ*K*_app (with vs without TPP1) | Δ*s*_max | reading |
|---|---|---|---|
| **C confirmed** | ≥2-fold decrease, CI excluding 1 | ≤1.2-fold, CI including 1 | TPP1 acts on retention; no mechanical signal required |
| **B confirmed** | ≤1.2-fold, CI including 1 | ≥1.5-fold increase, CI excluding 1 | TPP1 acts on the extension rate; reopen allostery with a precise target |
| **both** | ≥2-fold decrease | ≥1.5-fold increase | mixed; neither pure model survives; report as such |
| **neither** | no significant change in either | — | the observable is wrong, or TPP1 acts on neither term; **not** a licence to pick a winner |

**Power.** The minimum detectable ratio at n = 3 is ≈2.2× for a processivity-type readout
(`PREREGISTRATION.md` §4.3, computed from the published data's CV of 34.5%). So the design above
detects the pre-declared ≥2-fold *K*_app shift **only at the edge of its power**. Two options:
(a) accept n = 3 and treat a 2-fold shift as borderline, or (b) run n = 5 and detect 1.65×. I
recommend **n = 5 for the two extreme primer concentrations only**, keeping n = 3 for the middle
four; that costs ~40% more gels than a flat n = 3 and puts the power where the logistic fit needs
it.

## 4. The cheaper alternatives, and why they do not replace it

| option | cost | what it gives | why it is not decisive |
|---|---|---|---|
| Re-quantify existing gels for per-band intensities | ~1 day, no new experiment | separates constant-hazard from heterogeneous/repeat-dependent/pausing; yields R2 (absolute mass) | **cannot** separate C from B (§2) — the identity defeats it |
| Recover the −PT lanes and compute ΔRAP | ~1 day | the fold-stimulation, which nobody has reported per mutant | still a single primer concentration ⇒ still cannot separate *K*_app from *s*_max |
| A TPP1 separation-of-function mutant (binding-competent, stimulation-dead) | weeks of protein work | the strongest possible evidence for a signal, if one exists | much larger, and only worth it after the titration says a signal is there |
| HDX-MS / smFRET at the active site ± TPP1 | months, new instrumentation | direct conformational readout | right experiment *if* B is indicated; disproportionate to run first |

**The titration is the smallest experiment that can change the answer**, and the two cheaper
options are worth doing first only because they are nearly free — they sharpen the estimate but
cannot decide.

## 5. What would make even this experiment non-decisive

If the primer-titration curves with and without TPP1 are hyperbolas with *both* parameters shifted,
or if the fit is under-determined because processivity saturates below the accessible primer range,
the titration yields "mixed" and no model wins. In that case the correct report is that the
two-parameter separation is itself insufficient, and the next step is a
**single-molecule measurement** (direct observation of dissociation between repeats, with and
without TPP1) — which measures *k_off* directly rather than inferring it. That is the genuinely
decisive experiment, and it is an order of magnitude more expensive. **The titration is the
cheapest test that can either settle the question or justify that expense.**

## 6. Falsification criteria, stated before the data exist

* **Model C is falsified** if TPP1 changes *s*_max by ≥1.5× with *K*_app unchanged, or if a
  binding-competent TPP1 variant fails to stimulate at any primer concentration.
* **Model B is falsified** if TPP1 shifts *K*_app ≥2× with the plateau unchanged, across two
  independent primer series — that is, the conformational component would have no observable.
* **Both are falsified** if neither parameter moves with TPP1: then the effect is on a step neither
  model names (e.g. the number of accessible enzyme molecules, or a post-translocation
  rearrangement), and the correct conclusion is that the two-model framing is incomplete.
