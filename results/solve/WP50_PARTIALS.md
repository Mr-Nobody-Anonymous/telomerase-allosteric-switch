<!-- generated, do not edit by hand: python scripts/solve/wp50_partials.py -->

# WP50 — the three PARTIAL items, attacked by other means

    python scripts/solve/wp50_partials.py

Four work packages tried to close these three the same way: add a
constraint, add a filter term, add a protocol. This one tries different
means, and two of the three attempts fail. They are reported as failing.

---

## P6 — three routes that were never tried

Starting point: **6.36 kcal/mol unpaid** after the WP48
construct. Every prior attempt ADDED something. These try subtraction, a
different ring, and a different payer.

### Route 1 — delete the flank that has no function

WP45 kept `DWEEKEF` and `FSLLPTE` verbatim. Only E168, E169, E171, E215
carry payload, and only W167/L211 hold the staple. The rest is flank:

| residue | function |
|---|---|
| D166, K170, S210, L212, P213, T214 | **none** — neither payload nor anchor |
| F172, F209 | junction position only |

Searching every segment pair that keeps **all four payload residues AND both
staple anchors**, the smallest construct is a 15-mer:

| | segment A | segment B | linker | total | junction |
|---|---|---|---|---|---|
| WP45 | 166–172 | 209–215 | 5 | 19 | 16.66 Å |
| **WP50** | **167-172** | **211-215** | 4 | **15** | **13.48 Å** |

Dropped: D166, F209, S210.

> K170 and L212/P213/T214 look like dead flank by function, but K170 sits BETWEEN E169 and E171 and L212-P214 sit between L211 and E215. Removing an internal residue changes the spacing that defines the triangle, so it is not flank at all - it is the strut. Only TERMINAL residues can be trimmed without moving a payload. That is why 3 come off, not 6.

**What this does and does not buy.** The fraction of the molecule held by
one staple rises from 0.105 to 0.133. It does **not** buy a kcal/mol figure:

> Trimming removes FLEXIBLE residues, which removes cost, but the cost model is not linear in length - WP45's 12-20 kcal/mol was priced for a 20-residue flexible pool, not per residue. So the trim is reported as a FACT (fewer movable residues, all functions preserved) and NOT as a fabricated kcal/mol saving. What it does buy is a smaller construct the same staple can rigidify: a 6-torsion staple over 15 residues constrains a larger FRACTION of the molecule than over 19.

### Route 2 — macrocyclisation (and why it is now closed for good)

| pair | distance |
|---|---|
| **C(215) → N(166)** — the construct's actual free ends | **20.42 Å** |
| CA(215) → CA(166) | 19.1 Å |
| C(172) → N(209) — the junction WP45 tested | 16.31 Å |

A peptide bond reaches ~7 Å. The closure distance is **20.42 Å**, so it would need **4 extra residues** — which turns the
macrocycle back into the linear construct it was meant to replace.

> head-to-tail cyclisation is INFEASIBLE: 20.42 A against a ~7 A peptide-bond reach, so it would need 4 extra linker residues - which converts a macrocycle back into the linear construct it was meant to replace

> This is a genuine negative result of a kind WP45 did not have: not 'head-to-tail is too far' but 'the closure distance exceeds the reach by a factor of three, so no variant of it works'. Ring closure is closed off for good.

### Route 3 — pay the entropy with binding instead of folding

Every prior attempt assumed the peptide must **pre-fold, then bind**, so the
ordering cost is paid by the peptide alone: 12–20 kcal/mol against 5.41
available. A covalent tether changes *who pays*.

| engineered position | distance to anchor | nearest peptide Cβ |
|---|---|---|
| D129→Cys | 6.65 Å | P213 at 4.57 Å |
| T128→Cys | 7.12 Å | P213 at 9.0 Å |
| S121→Cys | 7.2 Å | E215 at 10.29 Å |
| V81→Cys | 7.03 Å | E215 at 12.0 Å |

**Specification:** a non-natural electrophilic residue (chloroacetamide or acrylamide on a cysteine-targeting side chain) at the peptide position nearest P213, paired with a D129Cys mutation. Both are standard; neither is modelled here..

| term | kcal/mol |
|---|---|
| measured contact energy | 5.41 |
| shortfall before a tether | 6.59 |
| tether gain (freezing one rotor pair) | +4.0 |
| **residual** | **2.36** |
| still short? | **True** |

> A tether does not create binding energy. What it removes is the translational/rotational entropy of the second component, worth roughly 3-5 kcal/mol, and it converts the order-disorder cost from a folding penalty into a capture penalty. Priced at the conventional ~4 kcal/mol for freezing one rotor pair, it closes part, not all, of the 6.36 kcal/mol.

> This route is a SPECIFICATION, not a demonstration. It states which mutation and which peptide position, and prices the entropy it saves. It does not show the tethered construct folds.

**P6 verdict: STILL PARTIAL - three new routes measured. Ring closure is now closed off for good (20.42 A vs ~7 A reach). Trimming gives a 15-residue construct with every function kept. The covalent tether is the only route that addresses WHY the cost exists, and it leaves 2.36 kcal/mol unpaid.**

---

## P7 — the isolation criterion was the wrong criterion

The blocker was that MD3 was built on a helix that is really a kinked turn.
Re-fitting registers would just repeat the error. The question the
instruction actually poses is: **does the module solve its own site and
touch no other?**

| module | presents | sites engaged | verdict |
|---|---|---|---|
| MD3-PATCH (WP44, helix) | 3 carboxylates | S1_K78_R132 | **ONE-SITE-FULL** |
| MD3-GRAFT (new, grafted) | 3 carboxylates | S3_TEL_TRIAD | **ONE-SITE-FULL** |
| MD1-GROOVE (WP44) | 2 carboxylates | S1_K78_R132 | **ONE-SITE-FULL** |
| MD4-ANCHOR (WP44) | 0 carboxylates, hydrophobic | S4_TEN_LOOP | **ONE-SITE-FULL** |

**The WP44 helix module does not reach S3 at all. Its carboxylate spacings [6.43, 5.08, 10.75] A miss the measured triad triangle (6.58, 10.22, 7.31) A by more than the 1.0 A tolerance, and its shortest pair matches the 4.35 A K78-R132 separation instead. So the helix module was not merely suboptimal for S3 - it was pointed at a DIFFERENT site. The graft's spacings are the measured triangle itself, so it reaches S3 and nothing else.**

> By fixing the SPACING, not the count. The helix module's shortest carboxylate pair is 5.08 A - within 1.0 A of the 4.35 A K78-R132 separation, so it reaches S1. The graft's shortest pair is 6.58 A, 2.23 A away from S1, so it does not reach it. The discriminator is measured spacing, and only the graft has the measured triangle.

> P7 is answered, and the answer is worse for WP44 than expected: MD3 was never pointed at S3 - it was aimed at S1. It is closable by adding a SPACING test to the isolation criterion, which every module now passes with its correct site named - but the WP44 triad module is reclassified, not repaired.

**P7 verdict: PARTIAL -> the isolation criterion was the wrong criterion. Scored on spacing, the WP44 triad module reaches S1_K78_R132 and the graft reaches S3_TEL_TRIAD, so the helix was aimed at a different site, not merely mis-registered. Every module now has exactly one named site - which is the property the one-site-full instruction actually requires.**

---

## P10b — change the question, not the threshold

Lipinski rejects atorvastatin at MW 559. Rather than refit, two changes of
question:

### (1) A filter is triage, so measure enrichment

| rule | survives (good/bad) | bad survival | good survival | enrichment |
|---|---|---|---|---|
| `lipinski_strict` | 4/5 : 6/6 | 1.0 | 0.8 | 0.0 |
| `reactivity_only` | 5/5 : 0/6 | 0.0 | 1.0 | **unbounded** |
| `reactivity_and_<=1_lipinski` | 4/5 : 0/6 | 0.0 | 0.8 | 5.0 |

**`reactivity_only` rejects every bad molecule and accepts every approved
drug** — unbounded enrichment, zero false rejections. The Lipinski objection
vanishes, not because the threshold moved, but because a Lipinski violation
stops being a rejection.

### (2) Exposure — an alert is a hazard, not a risk

| molecule | electrophile | dose (mg) | mmol | rejected at 0.10 mmol |
|---|---|---|---|---|
| OK-amoxicillin | — | 500 | 1.368 | no |
| OK-atorvastatin | — | 20 | 0.036 | no |
| OK-warfarin | — | 5 | 0.016 | no |
| OK-nalidixic-acid | — | 1000 | 3.973 | no |
| OK-vancomycin-frag | — | 125 | 0.644 | no |
| BAD-chloroacetamide | alkyl_halide | 500 | 4.113 | yes |
| BAD-benzyl-bromide | alkyl_halide | 500 | 2.923 | yes |
| BAD-acrylonitrile | michael_acceptor | 50 | 0.942 | yes |
| BAD-ethyl-mesylate | alkyl_halide | 500 | 4.027 | yes |
| BAD-phenyl-isocyanate | isocyanate_isothiocyanate | 100 | 0.839 | yes |
| BAD-formaldehyde | aldehyde | 10 | 0.333 | yes |

| dose threshold (mmol) | correct |
|---|---|
| 0.01 | 10/11 |
| 0.03 | 10/11 |
| 0.05 | 10/11 |
| 0.1 | 10/11 |
| 0.2 | 10/11 |
| 0.5 | 9/11 |
| 1.0 | 7/11 |

> Dose co-factoring is not invented here: a published analysis of structural alerts for drug-induced autoimmune risk reports false positives falling from 12% to 4% and PPV rising from 25% to 42% when a >=100 mg daily-dose cut is added to a benzene-substituent alert. The mechanism is that an alert describes a HAZARD, and hazard needs exposure to become RISK.

> The doses are order-of-magnitude clinical figures entered by hand to show the calculation. They are not the paper's data and are labelled as illustrative wherever they appear.

> Switching the metric to triage solves P10b's stated problem without moving a single threshold: 'reactivity only' rejects every bad molecule and accepts every approved drug, so the Lipinski objection disappears because a Lipinski violation stops being a rejection. Whether that is CORRECT depends on whether the gate is meant to be an oral-drug filter or a specificity filter, and that is a decision, not a measurement.

**P10b verdict: PARTIAL -> framed as triage, the Lipinski objection disappears without moving a threshold; framed as an oral-drug filter, it stands. The dose-aware sweep is reported whole.**

---

## Summary

- **P6:** STILL PARTIAL - three new routes measured. Ring closure is now closed off for good (20.42 A vs ~7 A reach). Trimming gives a 15-residue construct with every function kept. The covalent tether is the only route that addresses WHY the cost exists, and it leaves 2.36 kcal/mol unpaid.
- **P7:** PARTIAL -> the isolation criterion was the wrong criterion. Scored on spacing, the WP44 triad module reaches S1_K78_R132 and the graft reaches S3_TEL_TRIAD, so the helix was aimed at a different site, not merely mis-registered. Every module now has exactly one named site - which is the property the one-site-full instruction actually requires.
- **P10b:** PARTIAL -> framed as triage, the Lipinski objection disappears without moving a threshold; framed as an oral-drug filter, it stands. The dose-aware sweep is reported whole.

Nothing here was promoted to `SOLVED`. Two routes failed outright (ring
closure, and the tether closing the last 4 kcal/mol), and both failures are
results: the first closes a whole family of designs permanently, the second
locates exactly where the remaining cost lives.
