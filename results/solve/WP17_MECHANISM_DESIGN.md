# WP17 — Why the primer titration does not settle the mechanism, and what would

**Status:** the negative result below is established. The positive design has not been run.
**Date:** 2026-11-06
**Guards:** `tests/test_open_findings.py::test_*` (mechanism section)

---

## 1. The question, stated so it can fail

TPP1 raises telomerase repeat-addition processivity by ~2.7 ± 0.6×. Two mechanisms remain:

- **R1 — affinity.** TPP1 changes how tightly the enzyme holds the primer. It acts on `K_d`.
- **R2 — retention/extension.** TPP1 changes how fast or how persistently the enzyme extends,
  independently of primer affinity. It acts on the catalytic step.

The repository previously reported "R1 excluded, R2 supported". That is not what the evidence
says. In `results/solve/wp3_predictions.json` the fit *sets*

```
"resting_koff_no_PT": 0.012
"resting_koff_PT":    0.012      <- identical; R1 is zero by construction
```

R1 was **assumed away inside the model**, not excluded by a measurement. The structural argument
(TPP1 is 28.7 Å from DNA and makes no DNA contact) excludes a *direct* TPP1–DNA interaction; it
does not exclude TPP1 changing the enzyme's primer affinity allosterically.

## 2. The model, and the degeneracy

Per-repeat survival under the repository's eq. 1:

```
theta = [P]/(K_d + [P])                  primer occupancy
r     = k_off / k_pol                    resting ratio
s     = 1 / (1 + r*(1 - theta))          per-repeat survival
```

Both arms act on `s` and nothing else:

| arm | TPP1 does | effect on `s` |
|---|---|---|
| R1 | lowers `K_d` by `n` | `theta -> P·n/(1+P·n)` |
| R2 | raises `k_pol` by `a` | `r -> r/a` |

At A₀ = 0.99, `r ≈ 0.0101`, and `1−theta` is already ~0.0099. Both corrections land on a term that
is already negligible, so **the two arms agree to 0.0005 fold**. This is the repository's own
thesis (`paper/MANUSCRIPT.md`, result 1) re-derived for R1/R2 specifically.

## 3. The titration separates them — in the wrong regime

Gap between the arms, as a function of baseline processivity and primer occupancy:

| A₀ | best \[P\]/K_d | R1 | R2 | gap | R2/R1 |
|---|---|---|---|---|---|
| 0.99 | 0.01 | 1.0004 | 1.0075 | 0.0071 | 1.007 |
| 0.90 | 0.01 | 1.0038 | 1.0803 | 0.0765 | 1.076 |
| 0.70 | 0.01 | 1.0115 | 1.2877 | 0.2762 | 1.273 |
| 0.50 | 0.01 | 1.0193 | 1.5952 | 0.5759 | 1.565 |

The ratio grows as A₀ falls — but only because **both arms collapse toward 1.0×**. The maximum
TPP1 effect the model can produce at A₀ = 0.70 is **1.29×**; at `[P]/K_d = 1` it is **1.15×**. The
effect being explained is **2.7×**.

> **The regime that separates R1 from R2 is a regime in which the phenomenon has largely
> disappeared.** Designing the titration is therefore not a fix.

At A₀ = 0.70, `[P]/K_d = 1`, the gap is 0.022 fold-units — about **1.9%** of the 1.15× effect
being measured. No densitometry resolves that.

## 4. What was actually verified in the structures

Independent of the kinetics, the TEN anchor geometry was re-measured from deposited coordinates
(`results/solve/WP16_ANCHOR_SITE_RESULT.md`, now CONTACT in all three long-DNA entries):

| entry | modelled DNA | Tyr176 → closest DNA |
|---|---|---|
| 7QXS | 22 of 30 nt | **2.75 Å** (OH···O3′ of nt19) |
| 7QXA | 12 of 30 nt | **4.35 Å** (OH···C5′ of nt19) |
| 9SHY | 10 of 34 nt | **2.78 Å** (OH···OP1 of nt25) |
| 7TRE/7TRD/7TRF | 7 nt | 11.11 / 11.94 / 11.94 Å |

This settles **geometry, not mechanism.** The contact is to the sugar–phosphate backbone, so it
cannot read sequence, and it is equally consistent with a positioning element and with a steric
barrier. It does not distinguish R1 from R2.

**Superposition caveat (fixed).** `compare_states()` superposed on RT 541–935 while reporting
IFD-TRAP 755–785, which lies entirely inside that range — the fit absorbed the motion it was
measuring (~16% attenuation, and ~0.5 Å of invented motion in the "rigid" ring). The frame is now
the invariant RT core (RT minus IFD-TRAP). The RAP-state superposition itself is **still
uncomputed**: 9SHZ and 9SI0 coordinates have never been retrieved, so no per-domain displacement
for those states exists in this repository.

## 5. The experiment that would settle it

The degeneracy is broken by anything that changes one arm **without** the other. In the model the
two arms sit in different factors, so the discriminating observable is not the size of the effect
but its **dependence on a variable that enters only one arm**.

**Design A — vary `k_pol` independently, hold `K_d` fixed.**
Use a substrate that changes extension rate without changing primer affinity — e.g. a primer
with a matched-length 3′ mismatch-adjacent base, or a Mn²⁺/Mg²⁺ substitution that lowers `k_pol`
with a measured, independent `K_d`. Predictions:

- **R2** (acts on `k_pol`): the TPP1 effect **moves** when `k_pol` is varied.
- **R1** (acts on `K_d`): the effect is **invariant** to `k_pol`.

This is the correct axis. It is the one variable the two arms disagree about.

**Design B — measure `K_d` directly, in the presence and absence of TPP1.**
A direct primer-binding measurement (filter binding or anisotropy) is the single observable that
adds rank to the processivity design (README, identifiability result 2). If TPP1 does not change
`K_d`, R1 is dead by measurement rather than by assumption, and R2 is all that remains.

**Design C — separate the effect from the anchor.**
`Y176A` removes the DNA anchor (2.75 Å contact) and reduces RAP stimulation. If anchor loss
abolishes stimulation, the effect requires DNA engagement. If stimulation survives anchor loss,
it does not, and R1 is dead. This is a perturbation, not a titration, and it does not depend on
the A₀ problem above.

**Power.** At measured A₀ = 0.99 the effect itself is hard to resolve; `wp15_hierarchical_power.csv`
puts it at 7 purifications for a 5% effect. Designs A–C do not require resolving a *gap* — they
require resolving whether a **2.7× effect moves**, which is a far easier measurement.

## 6. What this document does not claim

- It does not show R2 is true. It shows R1 has not been excluded and that the proposed
  titration cannot exclude it.
- It does not produce a RAP-state superposition. 9SHZ/9SI0 were never pulled.
- It does not convert the Tyr176 contact into a mechanism. A 2.8 Å backbone contact is
  compatible with several mechanisms, including a purely steric one.
- The Colab notebook that motivated this (`Untitled0.ipynb`) computed cells 24–26 but its saved
  outputs are absent; the numbers above were re-derived and are guarded by tests, not quoted.
