# Telomerase Allosteric Switch — Full Programme Report

**WP27 → WP45** · branch `arena/3d3bb488-telomerase-allosteric-switch`
Generated from the ledger and the executed scripts, not from prose memory.
Every number below is produced by a script in `scripts/solve/` that can be re-run.

---

## 0. The one-paragraph version

This programme set out to design a small molecule that switches telomerase off
allosterically, and to design it against a maintenance floor. Forty-five work
packages later the honest position is this: **the floor is reachable only for a
cell that needs about 48 bp per population doubling, it is not reachable for a
demanding cell, no molecule has been shown to bind anything, and the one
mechanism that survived every earlier audit — M5 — has now been refuted on data
the programme already had.** The work package that did the refuting (WP45) is
also the one that found a way to measure the stall share without new bench work.
That is the shape of the result: a sharper negative and a cheap measurement, not
a drug.

If you read one section, read **§6 (What is refuted)** and **§8 (What would
change the answer)**.

---

## 1. What was asked, and what was delivered

| ask | delivered | where |
|---|---|---|
| Design a switch to beat 50 % odds | Design built; P fell from 0.70 to **0.35**, then M5 was **refuted** | WP40, WP45 |
| Cover both referents (>50 % P **and** 50 bp/PD) | Both tracked; the floor is met only at c ≈ 0.7–1.0 | WP31, §5 |
| Explore freely, label assumptions | Every load-bearing assumption is now named in §7 | throughout |
| Try to solve all six blocked items | Six attempted; **none** converted to a pass; one new model-free pass found *en route* | WP44 |
| Solve it a more creative way | **Yes — and it refuted M5** | WP45 |

---

## 2. The model in one page

The whole programme rests on one geometric survival model. Telomerase adds
copies of `TTAGGG` (6 nt) to the 3′ overhang in processive repeats.

- `s` = probability of adding the **next** repeat
- `s^m` = probability of adding `m` repeats
- `R½` = mean repeats added = `ln2 / (−6 ln s)`
- `δ` = nucleotides per telomere per event = `6 · R½`
- `c` = fraction of telomeres actually extended in a cell cycle (**measured**, see §5)

Processivity is then converted to a maintenance rate by multiplying by `c`. That
is the entire model. It has no free shape parameters — only `s` and `c`, and `c`
is measured.

**The stall extension (WP44).** If a stalled molecule is allowed to retry, the
*observed* survival inflates:

```
s_obs = s_true / (1 − τ (1 − s_true))      τ = stall share
```

Everything derived since WP34 lives inside this equation. **It is the load-bearing
structure of the second half of the programme**, and its validity rests on one
unchecked premise: that the retry is memoryless (§7).

---

## 3. What was built (and survives)

These are real, tested artefacts. The full suite is **652 passed**.

| component | what it is | tests |
|---|---|---|
| `wp30`–`wp33` | the survival model, δ\*, the `c` calibration | 94 |
| `wp35_p2_identifiability` | identifiability of the P2 chase design | 54 |
| `wp37_missing` | the missing-measurement analysis | 28 |
| `wp38_archive` | the retraction register (**21 rows**, C/R/X namespaces) | 63 |
| `wp39` | affinity surrogate — *built, then refuted* | 39 |
| `wp40_design_to_fifty` | the design-to-fifty calculation and the gate | 45 |
| `wp41_status_ledger` | the generated status ledger | 44 |
| `wp44_attempt_all` | all six blocked items, attempted | 31 |
| `wp45_creative` | τ from published numbers; the refutation | 24 |

The register in WP38 is the programme's most useful single artefact: it is a
machine-checked list of every claim that was withdrawn, with a **predicate** for
each — what observation would have falsified it — plus document-scanning
negative controls that fail the build if a superseded claim is still stated
anywhere without a banner.

---

## 4. Every retraction, in one table

21 retractions. Nine are **circular inference**, the programme's dominant failure
mode: a quantity derived from an input and then used to "predict" that input.
Every one of them passed its arithmetic checks.

| id | what was claimed | why it fell |
|---|---|---|
| C1–C7 | assorted circular derivations | derived → re-used as evidence |
| R1–R7 | early fold/affinity claims | same family |
| R8 | the `frac → s` inversion | wrong sign convention, 2.2× ambiguous |
| R9 | a stale gate claim left in prose | prose fix without a banner |
| **R10** | the two-chase design measures τ | **it is an identity** — `m/c` returns τ by construction, error 8e-17 |
| **R11** | `n` (repeats per binding event) must be measured | **`n` is never load-bearing** — 7.4 % spread over n = 1…64 |
| **R12** | M5 is live at P ≈ 0.35 | **refuted 2.4× on L&C's published pair** |
| **R13** | δ\* = 48 bp/PD at the anchor decides open/closed | the anchor assumes the +PT stall *vanishes* while τ(−PT) is large enough to move the fold |

**R10 and R11 are worth dwelling on.** Both were the programme's stated reasons
for needing new experiments. Both dissolved when actually executed: the two-chase
"measurement" was algebraically the same number it claimed to measure, and the
elusive parameter `n` cancels out of every quantity anyone wanted it for.

---

## 5. The one thing that is genuinely measured: `c`

Unlike almost everything else in this report, `c` is not an assumption. It comes
from published single-telomere work.

| system | `c` | δ (bp/PD) | source |
|---|---|---|---|
| H1299 (telomerase-high) | 1.00 | **48.0** | Zhao 2009 |
| HeLa | 0.70 | **33.6** | Zhao 2009/2011 |
| HeLa, pH 7.4 | 0.77 | 37.0 | — |
| HeLa, pH 6.8 | 0.48 | 23.0 | — |
| hESC | 0.70 | 33.6 | — |
| *S. cerevisiae* | **0.07** | **3.4** | Teixeira 2004 |

The required value to hold a 50 bp/PD floor is `c_required = 1.042`. **Only the
most telomerase-positive cell type sits at the line.** Yeast is 15× short — the
model reproduces that, which is a point in its favour.

---

## 6. What is refuted

### 6.1 M5 — the mechanism (WP45)

M5 was the last surviving candidate. Its core claim: the allosteric switch works
by slowing dissociation, so `a = k_off/k_pol` falls by ~2.73×.

Latrick & Cech published **two independent per-repeat numbers** from the same
preparation and nobody had set them against each other:

```
s_obs(−PT) = 0.7914     φ(−PT) = 0.570
```

Both are functions of τ, so the pair fixes it: **τ₀ = 0.6506**. With τ fixed, M5
makes a held-out prediction:

| | |
|---|---|
| M5 predicts (+PT fold) | 2.734 |
| implied by L&C's own data | 6.654 |
| L&C stated measured | ≥ 4.00 |

**M5 under-predicts by 2.4×.** This is the strongest evidence against M5 the
programme produced, and it required no new experiment.

### 6.2 The gate number (WP45)

The stored `R½(+PT) = 8.0` is recovered **only** if the +PT stall is assumed to
vanish. Applying τ consistently gives `δ* = 7.9 bp/PD`, six times lower.

The programme has been computing the gate with τ(−PT) large enough to move the
fold and τ(+PT) small enough not to move the survival. Those cannot both hold.

### 6.3 The affinity surrogate (WP39)

Built, executed, **refuted**: amide Δ = −0.046 against a −0.15 threshold, and
wrongly signed. Zero calibration ligands exist. **Affinity: none.**

### 6.4 Gate 3 (WP44)

Not invertible. A single chain fraction maps to 2.3× different `R½` values
depending on convention (0.9529 → 2.395 vs 0.8945 → 1.036).

---

## 7. Every load-bearing assumption, named

This is the section that decides whether anything above is trustworthy.

| # | assumption | carries | testable? |
|---|---|---|---|
| A1 | **`s_true = φ`** | the entire WP45 refutation | **yes — one experiment** |
| A2 | the retry is **memoryless** | the whole stall extension | yes |
| A3 | `δ = 6·R½` with mean-field `c` | all δ values | partially |
| A4 | L&C's ≥4× is a lower bound, not a point | §6.1 magnitude | yes |
| A5 | Tesmer's IFD panel shares a convention across mutants | WP42 reciprocity | yes |
| A6 | TLN-amide is the active species, not racemic | the drug branch | yes |

**A1 is the one that matters.** Every word of §6.1 and §6.2 rests on reading φ as
the un-retried per-repeat success. WP22 established that φ is *window-limited*, so
this reading is possible — and that is precisely why it is a modelling step rather
than a theorem. If A1 is wrong, the refutation falls.

I am labelling it as an assumption, not a finding, exactly as instructed.

---

## 8. What would change the answer

In order of information per unit cost:

1. **Test A1.** Diffusive-access or competitor-titration experiment separating
   translocation efficiency from repeat-extension efficiency. Resolves whether M5
   is refuted or restored.
2. **Test A2.** Repeat the chase with a stall-inducing lesion and check whether
   `s_obs` follows the renewal or saturates. Bounds τ independently.
3. **`k_off` for a second enzyme pair.** Still blocked: the Tesmer panel has
   fractions only — no Kd, no k_off anywhere in the deposited data (§9).
4. **Any binding datum at all** for TLN-amide. Zero exist. This is the single
   largest hole in the drug branch.

---

## 9. What is still blocked, and why "blocked" is now earned

Each of these has an executed failed attempt behind it (WP44), not an assertion:

| item | verdict | evidence |
|---|---|---|
| measure `n` | **never load-bearing** | 7.4 % spread, n = 1…64 |
| second-enzyme fold | **blocked** | panel is fractions only, no kinetics |
| P2 two-chase | **design is vacuous** | identity, error 8e-17 |
| affinity | **refuted** | Δ = −0.046, wrongly signed |
| Gate 3 | **unreachable** | 2.3× span from one number |
| Y176A panel | **blocked** | BrdU/CsCl counts repeats, not attempts |

---

## 10. The bonus pass

WP44 turned up one **model-free, published** discriminator nobody had used: under
partial inhibition the observed fraction falls (1.00 → 0.50 → 0.30) while the
**length stays constant at ~60 nt**. That is the signature of a **processive**
enzyme, and it is independent of the model. Small, but it is a genuine pass.

---

## 11. The bottom line

**Probability that the current design achieves a switch with >50 % odds: below
0.35, and the mechanism that carried that 0.35 is now refuted.** The number is not
zero only because A1 is unvalidated.

**The maintenance floor is met for c ≈ 0.7–1.0 (33–48 bp/PD) and not for 50+
bp/PD cells at the true reading.**

**No molecule in this report has been shown to bind its target.** There is no
drug here. There is a model, a tested gate, a machine-checked retraction register,
21 withdrawn claims, and one sharp, cheap, falsifiable refutation of the
programme's own leading hypothesis.

That is an honest negative result. It is not what was hoped for, and it is worth
more than a 0.70 that was never tested.

---

## 12. Reproducing any of it

```bash
python -m venv .venv && . .venv/bin/activate
pip install numpy scipy pandas matplotlib biopython gemmi sympy requests pytest rdkit
pip install -e .
pytest -q                                   # 652 passed
python scripts/solve/wp45_creative.py       # the refutation
python scripts/solve/wp41_status_ledger.py  # regenerate the ledger
python scripts/solve/wp38_archive.py        # the retraction register
```

Every script has `--self-test`. The register's document-scanning controls will
fail the build if a superseded claim reappears anywhere without a banner.
