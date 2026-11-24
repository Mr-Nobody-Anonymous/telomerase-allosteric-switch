# The Resolution Plan — what is missing, and how to close it

**WP46** · generated from `scripts/solve/wp46_the_plan.py`
Re-run any number here: `python scripts/solve/wp46_the_plan.py`

> **What this document is.** A plan, not a finding. Every probability in it is a
> **labelled assumption** and is meant to be argued with. What is *derived* is the
> **structure**: how many independent routes exist, which missing measurement is
> decisive, and which of the gaps this repository can close by itself.

---

## 1. The thing that was actually missing, and nobody had named

The programme tracked one number for 45 work packages:

```
P(M5) = P(allostery is real) × P(the P2 chase confirms) = 0.70 × 0.50 = 0.35
```

You asked for the design to be better than 50 %. The programme reported 0.35 and
then spent WP34–WP45 trying to lift it — first by chasing the chase, then by
refuting M5.

**Both were measuring the wrong quantity.** `P(M5)` is the probability of *one
mechanism*. What you asked about is the probability that *a switch exists at all*.
Those are different numbers that happened to share a label, and the difference is
the whole reason 46 work packages felt like they kept arriving at bad news.

Executed:

| quantity | value |
|---|---|
| `P(M5)` alone — the old number | **0.350** |
| `P(any of M1–M5)`, routes treated as independent | 0.702 ← **wrong**, see §2 |
| `P(any of M1–M6)`, with the shared premise | **0.541** |

Refuting M5 removed **one route out of six**. It did not remove the outcome.

---

## 2. Why the structure matters more than any single prior

The six routes are **not independent**. Five of them share one premise: that a
druggable pocket coupling to processivity exists at all. If that premise is false,
M1–M5 all fail together. The old `0.70 × 0.50` silently assumed independence, which
is why it looked far more informative than it was.

```
P(any) = 1 − [1 − premise × P(any allosteric)] × [1 − P(M6)]
```

| premise | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 | 0.80 |
|---|---|---|---|---|---|---|
| **P(any)** | 0.420 | 0.461 | **0.501** | 0.541 | 0.581 | 0.621 |

The decision crosses 0.5 at **premise = 0.498** (computed, not eyeballed).

### The lever ordering — and the result I did not expect

I expected the shared premise to dominate. Per unit of input it does not:

| lever, +0.10 | Δ P(any) |
|---|---|
| **M6 — non-allosteric catalytic** | **+0.0656** |
| premise | +0.0402 |
| M1 | +0.0230 |
| M4 | +0.0224 |
| M3 | +0.0218 |
| M2 | +0.0211 |
| M5 | +0.0183 |

**M6 is the largest single lever in the entire model — and M6 is the route the
programme excluded as out of scope.** It dominates because it is the only route
that does not get discounted by the allosteric premise.

That is the most consequential single line in this programme's planning: the
**scoping decision**, made early and never revisited, is worth more than any
experiment since.

---

## 3. What is missing — five classes of gap

This is the direct answer to your question. It is not a list of experiments; it is
a list of **kinds** of gap, and only two kinds are closeable here.

**Live tracking:** every gap below is an open GitHub issue.

| gap | issue |
|---|---|
| G1 reference standard | [#14](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/14) |
| G2 measurement | [#16](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/16) |
| G3 formal (label collision) | [#15](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/15) |
| G4 material | [#17](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/17) |
| M6 scoping decision | [#18](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/18) |

| # | kind | what is missing | closes it | simulable here |
|---|---|---|---|---|
| **G1** | **reference standard** | a positive control — *any* compound in this repo shown to bind its target, or *any* assay with a known-answer ligand | nothing downstream is interpretable without it | **no** |
| **G2** | **measurement** | the A1 test (`s_true = φ`) and the A2 test (memoryless retry) | decides whether the WP45 refutation stands | **no** |
| **G3** | **formal** | the probability was the wrong quantity | re-express as P(any of M1–M6) — **done in WP46** | **yes** |
| **G4** | **material** | a compound series with SAR; one quarantined proposal is not a series | no med-chem iteration is possible | **no** |
| **G5** | **process** | preregistration against a real assay, a stopping rule, an independent replication slot | prevents the circular-inference class (9 of 21 retractions) | **yes** |

**G1 is the one that should worry you most.** There is no positive control anywhere.
Affinity: none. Drug: none. Every "no" in the ledger is a *no result*, not a
*measured zero* — the difference between "the assay didn't fire" and "the compound
doesn't bind." Without a known-answer ligand, that distinction cannot be made, and
no amount of analysis closes it.

---

## 4. The order to do things in — and why the programme's order was different

Ranked by **expected information gain** (bits of decision-uncertainty removed):

| rank | test | EIG (bits) | cost | what it settles |
|---|---|---|---|---|
| **1** | **A1** — is `s_true = φ`? | **0.4547** | 1 experiment | whether the WP45 refutation stands |
| 2 | **B1** — any binding datum for TLN-amide | 0.3391 | 1 SPR run | the whole drug branch |
| 3 | A2 — memoryless retry | 0.3295 | 1 experiment | the stall renewal itself |
| 4 | B2 — `k_off`, second enzyme pair | 0.1893 | 1 stopped-flow series | M5/M3 separation |
| 5 | C1 — expert prior elicitation | 0.0661 | 3 interviews | the individual priors |

*I first ranked these by "how often the result flips the decision" and that metric
was **degenerate**: at P ≈ 0.54 every moderately informative test flips about half
the time, so all five scored ≈0.50 and the ranking was noise. EIG is the right
measure. The bad metric is left in the script's output as a record.*

**A1 is first, and it is first because WP45 made it load-bearing.** Running WP45
before A1 was the correct order anyway — it was cheap and it fixed what the
expensive test should be.

---

## 5. The plan

### Track 1 — close G3 and G5 here (this repository can do it)

| step | deliverable | status |
|---|---|---|
| 1.1 | replace `P(M5)` with `P(any)` everywhere; retract the label collision | **done, WP46** |
| 1.2 | mechanism inventory M1–M6 with labelled priors, re-elicitable | **done, WP46** |
| 1.3 | VOI ranking so effort goes where information is | **done, WP46** |
| 1.4 | track every gap as a GitHub issue with an owner and a close condition | **done — issues #14–#18** |
| 1.5 | preregistration for A1 and A2: prediction, stopping rule, replication slot | **done — `PREREG_A1.md`, `PREREG_A2.md`** |
| 1.6 | document scan to catch any recurrence of the `P(M5)` framing | **done — caught two live offenders, both fixed** |

### Track 2 — close G2 (needs a bench; can be specified here)

| step | experiment | falsifies |
|---|---|---|
| 2.1 | **A1**: separate translocation efficiency from extension efficiency (competitor titration or diffusive-access) | if `s_true ≠ φ`, R12 and R13 fall — **fully specified in `PREREG_A1.md`** |
| 2.2 | **A2**: chase under a stall-inducing lesion; does `s_obs` follow the renewal or saturate? | if saturating, the stall extension is wrong — **fully specified in `PREREG_A2.md`** |
| 2.3 | **B1**: SPR or ITC on TLN-amide, **with a known-answer ligand on the same plate** | if no binding, the drug branch closes |
| 2.4 | **B2**: `k_off` on a second enzyme pair with kinetics, not fractions | separates M5 from M3 |

**Rule for Track 2: no step starts before its predecessor is executed.** The
programme's habit of specifying six parallel experiments is what produced six
blocked items.

### Track 3 — close G1 and G4 (needs chemistry; the real bottleneck)

| step | what | why it blocks everything | status |
|---|---|---|---|
| 3.1 | obtain **any** validated telomerase inhibitor as a positive control | the first measured "yes" in the whole programme | **specified — `ASSAY_DESIGN_G1_G4.md` §2** |
| 3.2 | build a series, not a molecule: ≥10 analogues around one scaffold | G4 — no SAR without a series | **specified — §3** |
| 3.3 | state the target product profile *before* synthesis | stops the goalposts moving | **fixed — §4** |

The assay design also fixes the **analysis pipeline and blinding rule** (§6) so
the data cannot be analysed after the fact in a way that favours the hypothesis.
None of §2–§4 can be executed here; that is the honest boundary.

### Track 4 — the scoping decision, revisited (the biggest lever)

**M6 is the largest per-unit lever in the model.** The programme excluded
non-allosteric catalytic inhibition as "not a switch." That is a *definitional*
exclusion, not an evidential one, and it is the single most expensive choice made
here.

Simulated rather than waited for (`scoping_options()` in the script). Both numbers
are reported for every option, because conflating them is R14 over again:

| option | P(allosteric) | P(outcome) |
|---|---|---|
| **A** — outcome matters, M6 prior → 0.45 | 0.3443 | **0.6393** |
| **B** — mechanism matters, M6 out | **0.3443** | 0.3443 |
| **B** at premise 0.30 | 0.1721 | 0.1721 |
| **C** — both, tracked separately | 0.3443 | **0.5410** |

**Recommended: C.** At premise 0.60, choosing B costs 0.54 → 0.34. That cost
should be visible in advance rather than discovered in a ledger row later — and
if the premise is only 0.30, option B is 0.1721.

**The choice is still yours** — but it is now priced, so it can be made on the
numbers instead of on a definition.

> **Correction, logged.** The first version of this section quoted 0.36 for the
> allosteric figure and 0.420 for "option B". Both were wrong: 0.3443 is the
> allosteric value, and 0.420 is the *outcome* figure including M6 — so it was
> not option B at all. The error was R14 re-committed inside the document that
> named R14. Corrected on issue #18 and pinned in the self-test.

---

## 6. What "solved once and for all" would actually require

Being blunt, because the alternative is a 47th work package that re-derives a
number from a number:

1. **A positive control (G1).** Without it every negative is ambiguous.
2. **A1 and A2 executed (G2).** These decide whether the current mechanism story
   survives at all.
3. **A compound series with measured binding (G4).** Currently: zero compounds with
   any binding datum.
4. **A preregistered assay with a stopping rule (G5).** Otherwise the 22nd
   retraction is already written.
5. **A decision on M6 (§4).** Worth more than items 1–4 combined.

**Items 1, 2 and 3 cannot be done in this repository.** Not because of effort —
because there is no wet lab and no chemical supplier. A repository can do algebra,
audit its own claims, and rank experiments by value. It cannot produce a binding
measurement, and every honest path to "solved" runs through one.

---

## 7. What has genuinely been closed

So that this does not read as pure negation:

- **The model** — one geometric survival model, `c` measured, no free shape
  parameters. Reproduces the yeast-vs-human 15× gap.
- **The identifiers** — `n` proved not load-bearing (R11); the two-chase design
  proved vacuous (R10). Two experiments were deleted from the plan for free.
- **The register** — 21 retractions, machine-checked, with a falsification
  predicate attached to each and document-scanning controls that break the build
  if a superseded claim reappears.
- **The refutation** — M5 refuted on published data, 2.4× under-prediction (R12).
- **The gate defect** — the stored anchor assumes the +PT stall vanishes while
  τ(−PT) moves the fold; 6× discrepancy (R13).
- **The framing** — `P(M5)` → `P(any)`, with a VOI ranking (WP46).

Six real closures. None of them is a drug.

---

## 8. Reproduce

```bash
python scripts/solve/wp46_the_plan.py            # every number in §1-§4
python scripts/solve/wp46_the_plan.py --self-test # 31 checks
python scripts/solve/wp45_creative.py            # the refutation
python scripts/solve/wp41_status_ledger.py       # the ledger
pytest -q                                        # full suite
```
