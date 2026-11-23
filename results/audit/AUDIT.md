# No-lab audit: four analyses, and what they do and do not establish

**Status:** complete and reproducible. `python results/audit/audit.py` regenerates
every number below in about two seconds; `--selftest` checks six invariants.
Output: [`audit_results.json`](audit_results.json).

**Scope.** None of this is a result about telomerase biology. It is an audit of
what four existing bodies of evidence can support: deposited coordinates, a
deposited mutant panel, the variance in that panel, and our own mixture
estimator. All four need no laboratory, no GPU and no network.

---

## A. The K757 claim is not supported by the deposited coordinates

Read directly from `data/raw/7TRE.cif`, `7TRD.cif`, `7TRF.cif` (chain A, author
numbering):

| measurement | 7TRE | 7TRD | 7TRF |
|---|---|---|---|
| Tyr176 OH – Gln794 OE1 | **2.34 Å** | 4.37 Å | 4.37 Å |
| Lys757 NZ – Gln794 OE1 | 11.44 Å | 11.90 Å | 11.90 Å |
| Lys757 NZ – Tyr176 OH | 9.91 Å | 10.30 Å | 10.30 Å |
| nearest atom to Lys757 NZ | Val755 O, 4.36 Å | Val755 O, 4.06 Å | Val755 O, 4.06 Å |

Lys757 does not contact the Tyr176–Gln794 pair: 9.9–11.9 Å is not a contact, and
the nearest atom to the Lys757 side-chain nitrogen is a backbone carbonyl at
4.1–4.4 Å, which is not a hydrogen bond. Any model that couples Lys757 to this
pair is coupling residues that the structure does not place in contact.

**A 2.34 Å O···O separation is shorter than any hydrogen bond (2.6–3.1 Å).** At
3.5 Å resolution a side-chain pair cannot be placed that precisely, so this is a
property of the deposited model, not a measurement of the protein. The modelled
junction files inherit it unchanged (`junction_bound_wt.pdb` = 2.34 Å,
`junction_apo_wt.pdb` = 4.37 Å); the 3.5 and 5.3 Å values quoted elsewhere are
the *means of the first frames of the randomised replicas*, not the deposited
geometry.

**The two models disagree with each other, which locates the uncertainty.** The
chemically comparable Tyr176 OH – Gln794 NE2 distance is 3.40 Å in 7TRE and 3.80 Å in
7TRD/7TRF — a 0.39 Å spread. The OE1 distance moves 2.03 Å between the same deposits,
and the per-atom B-factor for Tyr176 OH is **40.3 Å² in 7TRE against 130.2 Å² in
7TRD/7TRF**. So what differs is the Gln794 side-chain rotamer, and at that B-factor it is
not determinate. See `structure_audit()` and `salt_bridge_audit()`.

**The salt bridge is unresolvable from one deposit.** Of the six pairings of Arg774
NH1/NH2/NE with Glu168 OE1/OE2, only one is short (NH1–OE1, 2.32 Å); the other five are
**3.85–5.95 Å**. The TPP1 partner is modelled in only one of the three deposits (7TRD has
no TPP1 chain; 7TRF models it as chains E/F), so the contact cannot be checked across
deposits at all. That is why `CLM-005` was reclassified from "uncertain" to "untested":
a single 3.5 Å model cannot establish a specific salt bridge.

**What this does not say.** It does not say the interface is unchanged by TPP1,
and it does not say the density fails to support these side chains. Both need
EMD-26087 inspected in ChimeraX at a workstation. Three separate readings of
the coordinate file now agree (two using gemmi, one reading the mmCIF columns as
plain text), which rules out a parser bug but is not an independent test of the
density.

---

## B. Re-analysis of the deposited mutant panel

Paired log-ratio to WT within each experiment, Holm-Bonferroni across mutant
tests, 95% intervals. Processivity uses the raw Fig 3B values, so WT variance is
included; activity is already WT-normalised per experiment and is tested as a
one-sample log ratio.

### Processivity

| mutant | ratio to WT | 95% CI | *p* | Holm | class |
|---|---|---|---|---|---|
| K760E | 0.54 | [0.51, 0.59] | 0.001 | **survives** | reduced (strong) |
| R756E | 0.58 | [0.42, 0.81] | 0.020 | no | reduced (uncorrected only) |
| K757E | 0.41 | [0.22, 0.75] | 0.024 | no | reduced (uncorrected only) |
| L786E | 0.34 | [0.16, 0.74] | 0.027 | no | reduced (uncorrected only) |
| F759E | 0.21 | [0.06, 0.80] | 0.038 | no | reduced (uncorrected only) |
| K749E | 0.54 | [0.22, 1.30] | 0.094 | no | not resolvable |
| D788R | 0.61 | [0.29, 1.27] | 0.100 | no | not resolvable |
| H762E | 0.84 | [0.58, 1.21] | 0.179 | no | not resolvable |
| F776E | — | — | — | — | **unusable** |

Nine mutants were tested. **One survives correction.** F776E has a replicate of
−0.338: a fraction of signal cannot be negative, so that column cannot be used
as deposited.

### Sensitivity: the pairing assumption

Pairing is inferred from row order and cannot be verified from the file. Dropping
it (Welch on logs) changes the picture: R756E falls to *p* = 0.074, K760E becomes
marginal [0.30, 0.99], while K757E strengthens to *p* = 0.006. **The conclusion is
fragile to the one assumption the deposited file cannot confirm.**

### Activity

Every interval spans 1 (0.03–5.18, 0.02–4.48, …), because within-mutant CVs run
56–137% at df = 2.

**This is a statement about the deposited summary numbers, not about the authors'
work.** They hold the raw gels, which may have smaller variance. Wide intervals
mean "cannot tell", not "no effect" — and independent data (Sekne 2022, Fig 4E–G)
do report a significant activity reduction for K757A/F759A.

---

## C. The experiment was never powered to answer its own question

Independent preparations are the unit; bands within a gel are not replicates.
80% power, α = 0.05:

| endpoint | CV | n for 1.5× | n for 2× | n for 3× | detectable at n = 3 |
|---|---|---|---|---|---|
| processivity, panel median | 35% | **8** | 4 | 3 | only **≥ 2.75×** |
| processivity, K757E | 8% | 3 | 3 | 2 | ≥ 1.30× |
| activity, panel median | 74% | **23** | 9 | 5 | nothing up to 6× |
| activity, D788R | 137% | >40 | 19 | 9 | nothing up to 6× |

With Holm across nine tests: 13 preparations for a 1.5× processivity change, 38
for a 1.5× activity change.

**The "activity preserved" arm cannot be made from n = 3 at all.** To claim it,
an equivalence test is needed, and its cost is:

| CV | exclude a 33% loss | exclude a 19% loss |
|---|---|---|
| 35% | 8 preparations | 25 |
| 74% | 28 preparations | 93 |

This is why a `p > 0.05` criterion was never viable for that arm: the correct
design is one to two orders of magnitude larger than anything proposed.

---

## D. The recommended estimator recovers *s* but cannot detect the mixture

Two-species mixture, true majority *s* = 0.60, lognormal noise, bands 1–80,
detection floor at 10⁻³ of peak. Tail-slope estimator over three read windows:

| φ | *s*₂ | window 3–10 | window 10–30 | window 30–60 |
|---|---|---|---|---|
| 0.005 | 0.95 | 0.602 | 0.632 | empty |
| 0.010 | 0.95 | 0.605 | 0.680 | empty |
| 0.020 | 0.95 | 0.611 | 0.787 | empty |
| **0.050** | **0.95** | **0.627** | **0.909** | **0.947** |

The same distribution, same noise, returns 0.63 or 0.95 depending only on which
bands are read. But the effect needs φ ≳ 2%: at φ ≤ 1% the bias never exceeds
0.08 and the far window is empty. **The manuscript's worked example sits at the
top of the observable range** (φ = 5%, *s*₂ = 0.95), not at a "rare" contaminant.

A joint fit with the mixture weight free recovers the majority *s* to within
±0.02 across all conditions — against a mid-window tail slope off by +0.31. But
it does **not** identify φ: at φ = 1% the fit returns φ̂ between 0.19 and 0.31.
The correct statement is that the full-distribution fit is robust to a mixture it
cannot see, **not** that it estimates the mixture.

---

## What survives, in one paragraph

Lys757 is not in contact with the Tyr176–Gln794 pair in the deposited models, and
the 2.34 Å contact at their heart is below hydrogen-bond distance at 3.5 Å
resolution. The deposited mutant panel supports one strongly reduced processivity
mutant after correction (K760E), four before correction, and no activity effect at
all. The panel could never have detected less than a 2.75-fold processivity change
or any activity change, and an equivalence claim needs 8–93 preparations depending
on the margin. The tail-slope estimator is genuinely window-dependent, but only
when a minority species exceeds about 2%. **No mechanism is established, nothing
here is evidence about how TPP1 raises processivity, and the one experiment that
could settle it remains unrun.**
