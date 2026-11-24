# Track 3 — Closing G1 (reference standard) and G4 (material)

**The bottleneck.** Tracks 1 and 2 are specification work. Track 3 is the only
track that produces a *measurement*, and until it runs, every negative in this
repository is ambiguous — see issue [#14](https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch/issues/14).

This document specifies what would have to be bought, made, and measured. It does
**not** claim any of it has been done.

---

## 1. Why G1 blocks everything, stated precisely

WP39 built an affinity surrogate, ran it, and refuted it (amide Δ = −0.046 against
a −0.15 threshold, and wrongly signed). With **zero calibration ligands**, two very
different explanations fit the data identically:

| explanation | what it means |
|---|---|
| the surrogate is bad | fix the model, try again |
| the compounds do not bind | the whole drug branch is closed |

**The experiment cannot tell these apart.** That is what "no positive control"
means, and it is why the ledger's every "none" is a *no result* rather than a
*measured zero*.

## 2. Step 3.1 — obtain a positive control

The requirement is not "a telomerase inhibitor". It is **a compound with a
published, quantified, reproducible binding or inhibition datum against a defined
telomerase preparation**, run on the same plate as anything we test.

Candidate classes, in order of how well-documented the known answer is:

| # | class | why it is a good control | known answer to expect |
|---|---|---|---|
| 1 | **catalytic-site nucleoside analogue** (e.g. a BIBR1532-class or a nucleoside analogue with published IC₅₀) | quantified IC₅₀ in the literature against a named enzyme prep; mechanism understood | a specific IC₅₀ within ~2× of published |
| 2 | **oligonucleotide template antagonist** (e.g. GRN163L-class) | dose-dependent, published ~50 % inhibition at a stated concentration | inhibition rising monotonically with dose |
| 3 | **POT1-TPP1 protein** as a *processivity* control | not a small molecule, but it is the one reagent whose effect on `R½` is measured (0.78 → 3.3) | `R½` rises ~4× |

**Class 3 is the one to run first**, because the programme's entire model is about
processivity: POT1-TPP1 is the only reagent in the literature whose effect on `R½`
is already quantified, so it doubles as a *model* control and a *reagent* control.

## 3. Step 3.2 — build a series, not a molecule

The programme has **one** in-silico proposal (TELIN-1, quarantined, two open
sub-items: which aryl position, and a 160-clash pose). One quarantined molecule is
not a series, and without a series there is no SAR and no iteration.

Minimum viable series:

- **≥10 analogues** around one scaffold
- **≥1 negative control** structurally matched but predicted inactive
- **the active stereoisomer resolved** for at least one member (TLN-amide is
  currently racemic — WP42 — so it is not even established which species is active)
- every member carries a **computed** property (logP, polar surface area, a docking
  score) *and* a **measured** datum, or it is not a series, it is a list

## 4. Step 3.3 — the target product profile, written before synthesis

Goalposts move when the goal is not written down. Fixed now:

| property | acceptable | unacceptable |
|---|---|---|
| measured binding | K_d or IC₅₀ stated with a positive control on the same plate | "computational affinity" of any kind |
| selectivity | effect on a non-telomerase control enzyme stated | unstated |
| mechanism | stated, with the assay that distinguishes it from the alternative | "consistent with" |
| dose | a dose–response with ≥5 points | single-concentration data |
| stereochemistry | defined | racemic, with the active species unidentified |

**Any member failing a row is not a candidate.** The programme's history is of
candidates that were really hypotheses; this table is the line between them.

## 5. What running it costs, honestly

| item | what it needs |
|---|---|
| a positive control | commercial purchase or collaboration; **not available in this repo** |
| a series | synthetic chemistry, ≥ weeks of bench time |
| binding measurement | SPR / ITC / a TR-FRET binding assay, plus a known-answer ligand |
| the processivity control | recombinant POT1-TPP1 |

**None of this can be executed here.** There is no wet lab and no chemical
supplier. What *can* be delivered from here is the specification above, the
analysis pipeline that would consume the data, and the preregistered pass/fail
criteria — all three of which are now in the repository.

## 6. The analysis pipeline that would consume the data

Written now so that the data cannot be analysed after the fact in a way that
favours the hypothesis:

1. **Positive-control gate.** If the control does not return its published value
   within 2×, the whole plate is void. No results are read from a void plate.
2. **Blinding.** Compound identity hidden from the analyst until the primary
   statistic is computed.
3. **Primary statistic fixed in advance** per assay (IC₅₀ with 95 % CI; K_d with
   a stated model).
4. **Negative control must read inactive.** If it is active, the series is
   confounded and the run is void.
5. **Any compound whose only evidence is a computed score is reported as a
   hypothesis, never as a candidate.**

---

## 7. Status

- [x] positive-control requirement specified
- [x] control classes ranked, with the model control identified
- [x] series specification written
- [x] target product profile fixed
- [x] analysis pipeline and blinding rule written
- [ ] **positive control obtained** — needs purchase or collaboration (G1 open)
- [ ] **series synthesised** — needs a chemistry partner (G4 open)
- [ ] **any binding datum measured** — zero exist today

The last three are the whole distance between this repository and a result.
