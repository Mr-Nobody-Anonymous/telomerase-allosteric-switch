# WP16 — The anchor-site question, answered

**Status:** the measurement that could not be made in this repository has now been made.
**Date:** 2026-11-06
**Question:** does the TEN-domain anchor motif PLYQ (residues 174–177) physically contact
the DNA once the substrate is long enough to reach it?

**Short answer: yes — but only for one of the four residues, and only against the
sugar–phosphate backbone.**

---

## 1. Why the question was open

Every structure already on disk (7TRE, 7TRD, 7TRF) carries a **7-nucleotide** DNA. The
5′ end of that substrate stops two nucleotides short of PLYQ. The previously reported
figure — Tyr176 at **11.11 Å** from DNA — therefore describes a site the substrate
**never reached**. Extrapolating the 5′ exit direction gave ~6.8 Å, but that was a
straight line through one nucleotide step: an inference, not an observation.

The repository's own preregistration records the gap (`results/causality/PREREGISTRATION.md`, F1):
product-length and distribution evidence "cannot be retrieved here".

## 2. Why the first Colab runner would not have answered it

Two defects in `colab/tpp1_colab.py` as first shipped, both fatal in the same direction:

1. **A ≥15-nucleotide gate decided whether the test had run.** It required at least one
   entry with ≥15 *modelled* nucleotides. But 7QXA models 12 and 9SHY models 10, so
   both were discarded outright — even though 9SHY's catalytic aspartate sits **2.78 Å**
   from its DNA 3′ end, which proves its active site is correctly assembled and its
   coordinates are usable. The gate threw away good structures and then reported
   "TEST DID NOT RUN".
2. **Gating on size instead of on a control.** A structure cannot answer "is PLYQ near
   the DNA" merely because it is big; it can answer only if its active site is verifiably
   seated on the substrate. Size is a proxy. The aspartate is the measurement.

The corrected rule, implemented in `colab/verdict.py`:

- **the control decides usability** — catalytic Asp868 must be within 6.0 Å of the DNA
  3′ end, otherwise the entry cannot speak to the question;
- **a contact needs no coverage** — if an atom is within H-bonding distance, it is;
- **a non-contact needs coverage** — "no contact" is only admissible if the nucleotide
  that would have to reach the anchor is actually modelled, otherwise the result is
  `NOT RESOLVED`, never `NO CONTACT`.

## 3. What was measured

Coordinates taken from the deposited entries through the PDBe model server, restricted
to the residues needed for the measurement and stored in `colab/fixtures/` so the numbers
below are reproducible offline (`python colab/verdict.py --selftest`).

| entry | modelled DNA | Asp868 → DNA 3′ end | **Tyr176 → closest DNA** | atom pair | call |
|---|---|---|---|---|---|
| **7QXS** (Sekne 2022) | 22 of 30 nt (auth 7–16, 19–30) | **2.65 Å** | **2.75 Å** | OH ··· O3′ of auth 19 | **CONTACT** |
| **7QXA** (Sekne 2022, no POT1 side chains) | 12 of 30 nt (auth 19–30) | **2.93 Å** | **4.35 Å** | OH ··· C5′ of auth 19 | **CONTACT** |
| 9SHY (Balch 2026, initiation) | 10 of 34 nt (auth 25–34) | 2.78 Å | 30.23 Å | nearest modelled nt is auth 33 | NOT RESOLVED |
| 7TRE / 7TRD / 7TRF (7-nt reference) | 7 nt | — | 11.11 / 11.94 / 11.94 Å | — | out of reach |

Two independently refined models place the contact on **the same nucleotide** (auth 19),
which is the part of this result that does not depend on refinement detail.

Per-nucleotide profile for 7QXS, Tyr176 OH to the closest heavy atom of each modelled
nucleotide (auth numbering):

```
nt      16    19    20    21    22    23
d (Å)  11.75  2.75  3.50  9.02 13.79 19.00
```

The chain runs from 11.8 Å at nt 16, dips to **2.75 Å at nt 19**, and recedes again.
Nucleotides 17–18 are unresolved in the model, so the approach happens across a gap in
the deposited coordinates.

Pro174 is **7.49 Å** from nt 19 — so of the four PLYQ residues, only the tyrosine reaches
the DNA. The motif name is useful as a label; the contact belongs to Y176.

## 4. What this does and does not establish

**Established.** The geometric premise of a TEN anchor site is correct. There is a real,
short, chemically sensible contact between Tyr176 and the DNA backbone once the substrate
is long enough. The 11.11 Å figure was an artefact of a 7-mer that never reached the site.
This is no longer an extrapolation.

**Not established — and this is the important half:**

- **It is a backbone contact, not a base contact.** The nearest atom is the O3′ of a
  thymidine. Y176 touches the sugar–phosphate chain, so it cannot be reading sequence.
  That is equally consistent with a positioning element and with a steric barrier.
- **A contact is not a mechanism.** Nothing here distinguishes DNA retention from
  catalytic or allosteric acceleration. The identity `log(s/φ) = −6·log(1 + k_off/k_pol)`
  constrains a *ratio*, and a structural contact does not supply the missing kinetic term.
- **The resolution does not support three-significant-figure precision.** These are 3.2–3.8 Å
  cryo-EM maps. The exact atom pair and the 2.75 Å value carry roughly ±1 Å uncertainty.
  What survives any plausible error is the ordering: ~3 Å in the long-substrate structure
  versus ~11 Å in the 7-mers. That is a factor of four, far outside coordinate error.
- **9SHY does not contradict this.** It reports Y176 ~30 Å from DNA, but only auth 25–34
  are modelled and the anchor-contacting nucleotide (19) is absent. Under the corrected
  rule that is `NOT RESOLVED`. It would have been recorded as `NO CONTACT` under the old
  ≥15-nt gate — the precise failure mode this project has already committed three times:
  reading a missing measurement as a negative one.
- **Structural biology cannot settle function.** The discriminating experiment is still
  the one preregistered in `results/solve/PREREGISTRATION_WP14.md` and gated by the
  stopping rules in `results/solve/WP15_DESIGN_VALIDATION.md`. A structure cannot report
  whether an anchor site buys processivity; only kinetics can.

## 5. Bearing on the standing question

TPP1-stimulated RAP: mechanical/allosteric, or DNA retention? Unchanged. R1 stays excluded,
R3 stays unsupported, R2 stays the leading candidate, and the mechanism stays experimentally
unresolved. What changes is that one previously unmeasured premise — "is there a TEN site in
contact with the DNA at all" — now has an answer, and it is yes.

The result raises the priority of one specific experiment: a **PLYQ cavity mutant**
(most cleanly Y176A) tested for processivity and for TPP1 stimulation. If the anchor site
is functionally load-bearing, removing it should cost processivity while leaving recruitment
intact. That is falsifiable by the same readouts WP15 already specifies, and it does not
require new structural work.

## 6. Reproduce

```bash
bash colab/run_colab.sh            # self-test, then download + analysis + verdict
python colab/verdict.py --selftest # the geometry and rule checks, offline
python colab/test_colab.py         # validates every claim in colab/AGENT_START_HERE.md
```

`colab/fixtures/` holds only the residues needed to reproduce the table in §3, with
provenance recorded in each file. They are subsets of deposited coordinates — not models
and not reconstructions.
