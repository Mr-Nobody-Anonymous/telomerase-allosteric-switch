<!-- generated, do not edit by hand: python scripts/solve/wp49_master_record.py -->

# Master record — the TPP1/TERT design programme, end to end

    python scripts/solve/wp49_master_record.py

**Generated, not written.** Every number below is read from an artefact
(`results/solve/wp4*.json`), from the retraction register, or recomputed by the
script that wrote the artefact. Hand-typed summaries in this repository have
drifted three times (C1, C7, R8), so this document does not contain any.

## 0. The answer, first

**No. 12 of 16 sub-problems are solved, and all 12 are solved as CALCULATIONS. Zero are solved as evidence, because nothing has been made or measured. 0 are not solved, 4 are partial, and 0 cannot be answered by any computation here. The largest open item is still P6: WP47 found a real covalent constraint (W167-L211, 6.34 A) where WP45 had concluded none existed, but one staple buys 2.8 kcal/mol against a 12-20 kcal/mol ordering cost. WP48 then improved the construct (PAPAP linker, -2.84 kcal/mol) and still has 6.36 kcal/mol unpaid. A scaffold is still required.**

| status | count |
|---|---|
| SOLVED | **12** |
| PARTIAL | **4** |
| NOT SOLVED | **0** |
| UNMEASURABLE | **0** |

The distinction that matters: **every `SOLVED` item is solved as a
calculation. None is solved as evidence.** Nothing in this repository has
been synthesised, expressed or assayed.

---

## 1. What was asked, and what changed

The programme began as a structural-biology study of the TERT–TPP1 recruitment
interface across five deposits (7TRE, 7TRD, 7TRF, 7TRC, 2I46).
It became a design programme when the interface analysis produced a
conclusion that could be tested by design: the TPP1 TEL patch is **bipartite**,
and no drug-like small molecule can span it.

The user's later instructions that shaped everything after:

1. Make the designs **tailor-made per sub-problem** — each module solves
   **exactly one** part, **fully**, targeting **one** site.
2. **Merge** the modules, and **test the merge** for stability.
3. **Solve the remaining unsolved site by any method that works** — the
   helix-only restriction was an assumption and could be abandoned.
4. **Do deep research** on the method, and cite it.
5. **Design new entities for each remaining failure**, and test them against
   the same gates that killed the old ones.
6. Standing: **measure, don't assert**; no tolerance tuning to force a pass;
   failures and unsolved pieces stay visible.

---

## 2. The interface, measured

| quantity | value | source |
|---|---|---|
| K78 NZ – R132 NH1 (telomerase) | **4.35 Å** | 7TRE chain A |
| E215 – K78 | **4.1 Å** | 7TRE |
| E215 – R132 (nearest O–N) | **4.07 Å** | 7TRE |
| E168 – E215 | **20.18 Å** | 7TRE |

**The patch is bipartite.** The two anchors (K78/R132 and R774) sit
**20.26 Å** apart. A straight line between them is
blocked, so the shortest route is **25.54 Å** with a
119° bend (tortuosity 1.26).

### Why no small molecule can work here

Spanning 25.54 Å needs ~7.8 phenyl rings:
MW ≈ 683.0 Da, logP ≈ 11.17. Lipinski's limits are 500 Da and 5.

**This is a class exclusion, not a compound that needs optimising.** It is
the single most useful negative result in the programme, because it redirects
the whole effort from small molecules to designed proteins.

---

## 3. The modular design round (WP42–WP44)

Four sites were defined and one module was built per site:

| site | what it is | design intent |
|---|---|---|
| S1 | K78 / R132 pair | present carboxylates to both |
| S2 | R774 | present one carboxylate (single anchor) |
| S3 | TEL triad E168/E169/E171 | present a three-carboxylate cluster |
| S4 | TEN loop A130–L131–R132 | engage the loop |

**Modules:** 4 of 5 satisfied
`ONE-SITE-FULL` — solved their own site completely and touched no other.

**Merges:** 4 stable, 1 collapse-risk, 1 refused.

The merge round produced the finding that decided the next two work packages:
**the linker, not the modules, decides stability.** The same two modules give
a stable construct or a collapsing one depending only on how they are joined
(MRG-FLEX, −6.04 kcal/mol of self-attraction).

---

## 4. S3 — the site that was declared unsolvable, twice

This is the most instructive sequence in the programme, and it contains two
retractions.

### 4.1 WP44: "not helix-solvable at all" — and the error underneath it

WP44 fitted an ideal α-helix to the triad and found at best 2 of 3 contacts
reachable. The conclusion was reported as a property of the protein.

**WP45 measured the premise instead of assuming it:**

| check | measured | α-helix expects |
|---|---|---|
| ψ(E169) | **+24.4°** | −60 to −45° |
| CA168–CA171, 7TRE | **7.68 Å** | 5.05 Å |
| CA168–CA171, 2I46 | **7.70 Å** | 5.05 Å |
| 3-10 helix at i+3 | 6.10 Å | still not 7.68 Å |

A **positive ψ** breaks a helix, and two independent crystals agree to 0.02 Å.
The segment is a **kinked turn**. The "unsolvable" verdict was a property of
the *model*, not the protein — registered as retraction **X2**.

> The lesson is not "the helix was wrong". It is that WP44's self-test passed
> because it checked **consistency with the helix model** rather than
> **consistency with the deposit**. A test that re-derives its own assumption
> passes forever.

### 4.2 WP45: the graft

Method: **native-segment backbone grafting** of a discontinuous motif
(Azoitei et al., *Science* 2011, 334:373 — `Multigraft`). Keep the native
segments verbatim, solve only what connects them.

| segment | residues | carries |
|---|---|---|
| A | 166–172 `DWEEKEF` | E168, E169, E171 |
| B | 209–215 `FSLLPTE` | E215 |

Best junction CA171–CA210 = **14.38 Å** ⇒ ≥5 linker residues at 3.5 Å/residue.

Construct: **`DWEEKEFGGGGGFSLLPTE`** — 19 aa, 2055.2 Da, charge -4.44.

The measured triad triangle is **reproduced, not fitted**:

| contact | measured | in construct |
|---|---|---|
| E168-E169 | 6.58 Å | 6.58 Å (inherited) |
| E168-E171 | 10.22 Å | 10.22 Å (inherited) |
| E169-E171 | 7.31 Å | 7.31 Å (inherited) |

**3 of 3 contacts preserved**, E215 at 4.1 Å from K78, payload fraction 21%.

### 4.3 WP45's second error: "no covalent route exists"

WP45 screened **disulfides only** (3.5–5.0 Å Cβ–Cβ) and generalised. Same
defect class as the helix error, one work package later. Registered **X3**.

---

## 5. WP47 — reopening the covalent question, and the filter question

### 5.1 The crosslink WP45 missed

Six chemistries spanning 3.5–12.5 Å, all 49 segment pairs, four filters:

| filter | test | why it kills candidates |
|---|---|---|
| F1 reach | Cβ–Cβ inside the chemistry's range | a staple cannot stretch |
| F2 direction | both side-chain vectors point into the gap | a staple grows along the native side-chain vector |
| F3 clash | the span clears TERT by ≥3 Å | the interface is tight |
| F4 payload | neither residue is E168/E169/E171/E215 | stapling a payload destroys the contacts |

**2 of 49 pass all four.** The winner:

**`TRP167-LEU211` — hydrocarbon, 6.34 Å**

| quantity | value |
|---|---|
| side-chain alignment | cos = 0.98 / 0.09 |
| clearance past TERT | 21.63 Å |
| native side-chain contact across the span | **3.5 Å** |

The last row is why this is not a distance coincidence: **the native W167
indole and L211 side chain are already in van der Waals contact at 3.50 Å**,
bridging the two segments. The staple replaces a non-covalent bridge with a
covalent one.

### 5.2 The honest half

| term | kcal/mol |
|---|---|
| contact energy available (measured) | 5.41 |
| one staple's ordering gain | +2.8 |
| residual unpaid | **9.2** |

One staple is a **necessary constraint, not a sufficient one**.

### 5.3 The filter battery: proven unfixable as framed

- All 10/10 controls clear
  **every** hard physicochemical gate, so alert count is the only separating
  feature.
- The panel is **non-monotone**: metformin must pass with 2 alerts, nitrogen
  mustard must fail with 1.
- Exhaustive enumeration of every `n_alerts ≤ k` rule caps that family at
  **0.7** balanced accuracy against a 0.90 bar.

A mechanism-based rule scored 10/10 on its authoring panel but **6/8 on a
hold-out** — the misses being the alert library misnaming the reactive group
(maleimide reported as `phthalimide`, sulfonate ester as `Sulfonic_acid_1`).
**P10 stayed PARTIAL.**

---

## 6. WP48 — a new design per failure, under isolated agents

Method: four agents, **no agent importing another**, sharing only read-only
deposit measurements, plus a critic allowed to fail them. Rationale: WP44's
error survived because a test re-derived its assumption; WP45's survived
because one agent's conclusion became another's premise. Isolation does not
make an agent right — it makes it **independently wrong**, which is checkable.

### F1 — the linker was a design error

WP45 built the linker from `GGGGG`. **Glycine has no Cβ and therefore the
largest backbone state space of any residue — WP45 chose the most expensive
5-mer available and then reported its cost as unpayable.** Registered **X4**.

| linker | entropy cost (kcal/mol) |
|---|---|
| `GGGGG` | 4.1 |
| `AAAAA` | 2.33 |
| `APAPA` | 1.62 |
| `PAPAP` | 1.26 |
| `AAAAP` | 1.98 |
| `SSSSS` | 2.83 |
| `AASAA` | 2.43 |
| `PGPGA` | 2.32 |

**`PAPAP` chosen — saves 2.84 kcal/mol.**
Ranking stable under an alternative state model: **True**.

New construct: **`DXEEKEFPAPAPFSXLPTE`** — 19 aa, **1 staple** (2 crosslinked residues), charge -4.44, **0 payload residues consumed**.

| term | kcal/mol |
|---|---|
| ordering cost (WP45's figure) | 12.0 |
| linker saving | +2.84 |
| staple gain | +2.8 |
| **residual** | **6.36** |
| paid? | **False** |

> Each term is a model. The linker saving uses the same state model that ranked the linkers; the staple gain uses WP45's own 0.6 kcal/mol per torsion. Summing them is legitimate only if they are independent, and they are not fully - a constrained linker reduces the torsions a second staple could remove. The residual is therefore a LOWER bound on the shortfall, and it is reported as such.

### F2 — reading substructure instead of names (this one passes)

| molecule | detected as | accepted | expected |
|---|---|---|---|
| OK-amoxicillin | — | accept | pass |
| OK-atorvastatin | — | reject | pass |
| OK-warfarin | — | accept | pass |
| OK-nalidixic-acid | — | accept | pass |
| OK-vancomycin-frag | — | accept | pass |
| BAD-chloroacetamide | alkyl_halide | reject | fail |
| BAD-benzyl-bromide | alkyl_halide | reject | fail |
| BAD-acrylonitrile | michael_acceptor | reject | fail |
| BAD-ethyl-mesylate | alkyl_halide | reject | fail |
| BAD-phenyl-isocyanate | isocyanate_isothiocyanate | reject | fail |
| BAD-formaldehyde | aldehyde | reject | fail |

**Electrophile layer alone: sensitivity 1.0, specificity 1.0** (6 tp, 0 fp, 5 tn, 0 fn) on a panel used by neither WP43 nor WP47.

End-to-end: 10/11. The single miss is **Lipinski rejecting atorvastatin at MW 559** — a
gate defect inherited from WP43, attributed as such.

A `MW ≥ 60` floor was added while the SMARTS still missed acrylonitrile and
formaldehyde, then **deleted** once it changed no decision. The SMARTS
corrections were the actual fix.

### F4 — a finding that corrects WP47

| structure | chain | atoms/residue | residues with a side chain |
|---|---|---|---|
| 7TRE | TPP1 | 7.7 | 75% |
| 7TRE | TERT | 8.0 | 78% |
| 2I46 | TPP1 | 6.0 | 75% |
| 2I46 | TERT | — | absent from this deposit |

> WP47's clash filter used modelled atoms as the obstacle. Both deposits are coarse: TERT has 8.0 heavy atoms per residue and only 78% of residues have a side chain in 7TRE. Missing atoms cannot clash, so that filter was under-conservative and its clearances are upper bounds.

Re-checked with truncation-aware padding: **the staple span survives a conservative re-check**.

---

## 6b. WP50 — the PARTIAL items, by other means

Four packages had tried to close these the same way: add a constraint,
add a filter term, add a protocol. WP50 tried subtraction, a different
ring, a different payer, a different criterion and a different metric.

### P6 — three routes, and only one helps

| route | result |
|---|---|
| trim non-functional flank | 19 → 15 residues, all functions kept |
| macrocyclisation | **infeasible** — C(215)→N(166) = 20.42 Å vs a ~7 Å reach |
| covalent tether to TERT | D129→Cys, residual **2.36** kcal/mol |

> K170 and L212/P213/T214 look like dead flank by function, but K170 sits BETWEEN E169 and E171 and L212-P214 sit between L211 and E215. Removing an internal residue changes the spacing that defines the triangle, so it is not flank at all - it is the strut. Only TERMINAL residues can be trimmed without moving a payload. That is why 3 come off, not 6.

> This route is a SPECIFICATION, not a demonstration. It states which mutation and which peptide position, and prices the entropy it saves. It does not show the tethered construct folds.

### P7 — the module was aimed at the wrong site

> The WP44 helix module does not reach S3 at all. Its carboxylate spacings [6.43, 5.08, 10.75] A miss the measured triad triangle (6.58, 10.22, 7.31) A by more than the 1.0 A tolerance, and its shortest pair matches the 4.35 A K78-R132 separation instead. So the helix module was not merely suboptimal for S3 - it was pointed at a DIFFERENT site. The graft's spacings are the measured triangle itself, so it reaches S3 and nothing else.

### P10b — the objection dissolves when the metric is corrected

| rule | good survive | bad survive | enrichment |
|---|---|---|---|
| `lipinski_strict` | 4/5 | 6/6 | 0.0 |
| `reactivity_only` | 5/5 | 0/6 | unbounded |
| `reactivity_and_<=1_lipinski` | 4/5 | 0/6 | 5.0 |

> Switching the metric to triage solves P10b's stated problem without moving a single threshold: 'reactivity only' rejects every bad molecule and accepts every approved drug, so the Lipinski objection disappears because a Lipinski violation stops being a rejection. Whether that is CORRECT depends on whether the gate is meant to be an oral-drug filter or a specificity filter, and that is a decision, not a measurement.

---

## 6c. WP51-WP54 — the remaining six, each by a different instrument

The last items did not fail for one reason, so they could not be attacked with one method. Measurement, geometry, a natural experiment, an algebraic reframing and an arithmetic decomposition each did one.

### WP51 — the ordering cost, measured (the first MD in this repository)

Five packages had **typed** the ordering price: WP45 12-20, WP47 0.6 per torsion, WP48 the same number again, WP50 a 4 kcal/mol constant. WP51 put atoms on it.

> The first version of the script replaced the linker with a harmonic spring whose stiffness it set by hand (5/40/200 kcal/mol/A^2) and then reported that AAAAA was stiffest. That is reading the input back out as if it were output. It is recorded, not deleted.

Stage 1 had not finished when this record was generated; the measurement is in `WP51_ORDERING_MD.md` when it does.

### WP52 — P12: the interface, measured instead of declared unknown

P12 said a score would be fabricated. That premise expired: OpenMM and RDKit are here, and 7TRE carries hTERT at all-atom resolution for this interface. So the question became one a crystal can answer:

| linker | min clearance to TERT | atoms inside 2.5 A | clears |
|---|---|---|---|
| GGGGG | **8.38 A** | 0 | yes |
| PAPAP | **8.22 A** | 0 | yes |
| AAAAA | **8.51 A** | 0 | yes |

The epitope is charge-dominated, not a pocket: net **-4** over the 14 native residues, with E168/E169/E171 and E215 contributing carboxylates directly. > geometric compatibility with the deposited TERT surface; what remains: any affinity, any on-rate, any competition - an experiment.

### WP53 — P13: a natural experiment that was on disk all along

7TRD is hTERT **without** TPP1; 7TRE is the same enzyme **with** it. The structural half of P13 is therefore a superposition:

| comparison | CA atoms | CA-RMSD |
|---|---|---|
| 7TRD (apo) vs 7TRE (bound) | 934 | **0.356 A** |
| 7TRD (apo) vs 7TRF (apo) — control | 934 | 0.0 A (degenerate) |
| 7TRC chain E vs chain I — noise floor | 129 | **0.585 A** |

> no detectable TERT conformational change on TPP1 binding. The TPP1-associated difference is SMALLER than the difference between two copies of one protein in a single crystal, 9 of 934 residues beat the floor (chance would give ~5 %), and the largest window rotation is 2.15 deg against a 2.54 deg measured limit.

This removes 'TPP1 remodels TERT' from the candidate mechanisms. It does not choose between R1 and R2, and does not claim to.

### WP54 — P15, P14, P11: a drug, an exponent, and a rate

**P15 — is there a drug?** The potency requirement is modest (K_L <= 1.167 uM, i.e. -8.14 kcal/mol). The design fails for a different reason:

| term | value |
|---|---|
| contact energy (WP45/WP47) | 5.41 kcal/mol |
| ordering cost, unpaid (WP48) | -6.36 kcal/mol |
| **net binding free energy** | **+0.95 kcal/mol** |
| shortfall against the drug requirement | 9.09 kcal/mol |

**P14 — is n measured?** It is not a measurable parameter. Across n = 1 to 24 the two-chase observable m/c is 1.448879638 every time, and tau and kappa come back unchanged, while only the per-round rate moves (1.7 at n=1 to 0.042254 at n=24, a 40.2-fold spread that no observable sees).

> The question was mis-posed, and that is the closure. n is not a parameter of the model that needs measuring: the model's own statistic (delta* = n * R_eff, WP36) is exactly n-independent when stated exponent-free as 6 * R1/2 (WP29), and the two-chase inversion that was designed to recover n recovers kappa and tau identically at every n. What n really is, is the number of escape opportunities per repeat - one per nucleotide gives 6, which is a count of the nucleotides in a repeat and is readable from an existing gel.

**P11 — the success rate.** Decomposed into four factors; one is exactly zero for a stated reason, not for luck:

| factor | value |
|---|---|
| f1 structural compatibility (WP52) | 1.0 |
| f2 makeable (WP48) | **0.0** |
| f3 literature peptide hit rate | 1%-10% |
| f4 assay power (WP47) | 0.8 |
| p(>=1 binder), as designed | 0% |
| p(>=1 binder) with ordering paid, 9 screens | 7%-53% |

---

## 7. The closure ledger — every sub-problem and its status

| id | area | question | status | kind |
|---|---|---|---|---|
| P1 | Interface | Is the TPP1 TEL patch one binding surface? | **SOLVED** | geometric |
| P2 | Interface | Can a drug-like small molecule span both halves? | **SOLVED** | geometric |
| P3 | Interface | What is the shortest path across the gulf? | **SOLVED** | geometric |
| P4 | Triad | Can an alpha-helix present the E168/E169/E171 triad? | **SOLVED** | geometric |
| P5 | Triad | Can the triad be presented at all, in the native geometry? | **SOLVED** | geometric |
| P6 | Triad | Does a peptide HOLD that geometry? | **PARTIAL** | design-level |
| P7 | Modularity | Can each module solve only its own site, fully? | **SOLVED** | geometric |
| P8 | Modularity | Do merged modules stay stable together? | **SOLVED** | design-level |
| P9 | Modularity | Can one groove be blocked and repaired at once? | **SOLVED** | logical |
| P10 | Screening | Is a structural alert a usable filter at all? | **SOLVED** | computed |
| P10b | Screening | Does the full gate stack accept approved drugs? | **SOLVED** | computed |
| P11 | Screening | What is the experimental success rate? | **PARTIAL** | computed |
| P12 | Screening | Do any of these designs bind telomerase? | **PARTIAL** | computed |
| P13 | Mechanism | What does TPP1 change mechanically (M5)? | **PARTIAL** | computed |
| P14 | Mechanism | Is n measured? | **SOLVED** | logical |
| P15 | Drug | Is there a drug? | **SOLVED** | computed |

### What each open item waits on

- **P6 (PARTIAL)** — Does a peptide HOLD that geometry?
  - A covalent constraint EXISTS, and WP45 missed it by screening disulfides only. Over six chemistries spanning 3.5-12.5 A, W167-L211 passes reach, direction, clash and payload tests at 6.34 A with cos 0.98/0.09 - and the native indole/Leu contact there is already 3.50 A. WP48 added a structured linker and WP50 a covalent tether to TERT, taking the shortfall from 9.2 to 2.36 kcal/mol. Still not balanced.
  - *blocker:* the remaining 2.36 kcal/mol after the WP50 covalent tether. Ring closure is now closed off permanently (20.42 A vs ~7 A reach) and trimming only buys 4 residues, so the residue is a presenting scaffold or the TERT Cys tether, both of which need to be made
- **P11 (PARTIAL)** — What is the experimental success rate?
  - Decomposed into four factors, and only one of them is unknown. WP52 measured the structural factor at 1.0 (the graft clears the TERT surface). WP48 measured the makeability factor at exactly ZERO - 6.36 kcal/mol of ordering cost is unpaid - so the rate for the designs as they stand is 0, not small, and that is a statement about makeability rather than luck. With ordering paid, 9 independent screens at a cited 1-10 % peptide hit rate give 8-61 % probability of at least one micromolar binder. The residual uncertainty is the literature hit rate, given as a range.
  - *blocker:* f2 = 0 until the ordering cost is paid; f3 is a literature range, not a measurement at this interface
- **P12 (PARTIAL)** — Do any of these designs bind telomerase?
  - The geometric half is measured and it PASSES: with its 14 native residues pinned to the 7TRE coordinates to 0.00 A, the designed linker clears every TERT heavy atom by 8.22 A (PAPAP) - 0 atoms inside 2.5 A. The design is sterically compatible with the deposited TERT surface, so 'a score would be fabricated' was false for the geometric part. What is NOT measured is any affinity, and the remaining blocker is an experiment.
  - *blocker:* no affinity, no Kd, no competition experiment - and the graft's net binding free energy is +0.95 kcal/mol (WP54), so it would not bind even if ordered
- **P13 (PARTIAL)** — What does TPP1 change mechanically (M5)?
  - The structural half is measured and the answer is NO: the apo-vs-TPP1-bound difference over 934 CA atoms is 0.356 A, which is SMALLER than the 0.585 A between two copies of one protein in a single crystal. Only 9 of 934 residues exceed the duplication floor's 95th percentile - at or below the 5 % chance rate - and the largest 30-residue window rotation anywhere is 2.15 degrees against a 4.30 degree detection limit. TPP1 does not detectably remodel hTERT. The kinetic half (R1 vs R2) still needs the trap.
  - *blocker:* R1 vs R2 needs the two-chase; a static structure cannot order events

---

## 8. Every retraction

**21 registered.** The register is enforced: a marker removed from a
source document, or a refuted claim added without registering it, fails the
test suite.

| id | state | claimed | actually was | replaced by |
|---|---|---|---|---|
| **C1** | RETRACTED | M5 independently reproduces the 48 nt Latrick & Cech observed | delta* = n*R1/2 identically, so 47.98 bp is 6 x 8 | delta* is a measured quantity (6*R1/2), not a derivation (WP37) |
| **C2** | WITHDRAWN | the stalled fraction explains 88% of the 1.46x gap | (s0-phi0)/s0 is the anomaly itself, rescaled; negative in the +PT arm | stalling RECONCILES the two readouts at tau_0 = 0.1833 (WP35) |
| **C3** | RETRACTED | synchronised confirmation versus chance | the expected value was drawn from the null being tested | recorded for completeness; no live consequence |
| **C4** | RETRACTED | removing all 12 loops changes the output by 4.4% | the loop detector and the loop remover share a definition | recorded for completeness; no live consequence |
| **C5** | CORRECTED | the SPR and RAP benchmark datasets are 'physical' | generated by np.random.seed(42) from chosen k_on / k_off / R_max | relabelled SYNTHETIC at source (generator, manifest, doc, tests) |
| **C6** | RETRACTED | M5 passes an out-of-sample test on a second enzyme (Wang 2007) | no measured off-rate exists for Wang; 4.48x is a0/a1 from its own R1/2 | a second R1/2 pair exists (Nat Commun 2020, 1.82x) - two targets |
| **C7** | RETRACTED | delta* decides the maintenance gate | 1/log2(1+a) equals 6*R1/2 only at n=6; it scales with an unmeasured n | the gate is stated exponent-free at 48 bp/PD (WP37) |
| **R1** | RETRACTED | the TPP1 TEL patch epitope could not be located for docking | TPP1 was on disk in 2I46; it assembles at 1.966 A RMSD over 152 CA | WP31 assembled it and re-ran every ligand |
| **R2** | SUPERSEDED | the WP30 ligand series prefers the TEL patch (WP32 scoring) | the control was size-confounded (377 vs 224 atoms); corrected to 205 vs 186 | only the neutral amide prefers the interface, delta -0.199 (WP33) |
| **R3** | WITHDRAWN | the anchor must be basic | the patch is ambipolar (net -1); a neutral amide is the only winner | neutral amide anchor; TLN-amide (WP33, WP37) |
| **R4** | SUPERSEDED | P2 is a re-analysis of existing gels: days, no reagent | no per-band data exists here (WP24); it needs chase +/- a DNA trap | a TWO-chase design, 3-12 replicates per arm (WP37) |
| **R5** | SUPERSEDED | the maintenance gate is OPEN at 57.3 bp/PD | 57.3 is the n=6 branch of a quantity that need not depend on n | gate CLOSED at 48 bp/PD, exponent-free (WP37) |
| **R6** | SUPERSEDED | the feasibility sweep's s-ratio ceiling is the binding constraint | the ceiling assumed a factorisation that WP22 itself falsified | M5 drops phi from the product entirely |
| **R9** | CORRECTED | the maintenance gate is CLOSED at 48 bp/PD (stated in OPEN_GAPS section 16) | section 18 of the same document corrects it, but section 16 kept the unsuperseded conclusion and was the row the ledger summarised | section 16 marked CORRECTED BY WP40; WP34 and WP36 given the same banner so no document asserts the closed gat |
| **R8** | CORRECTED | the maintenance gate is CLOSED at 48 bp/PD | 48 is M5's own PREDICTION (8.0/2.963 = 2.70x) fed back in as if it were L&C's measurement; on the measured >=4x it is 70.1 bp/PD,  | gate OPEN on the measurement; the closed gate is evidence against M5, not against the target (WP40) |
| **R7** | CORRECTED | two of the six gaps listed in WP37 section 5 cannot be closed here | both were misclassified before the literature and the code on disk were checked: c is measured (H1299 = 1.00) and the docking crit | gap 5 CLOSED (gate closes harder, 48 < 50); gap 4 RUN and the neutral-amide charge hypothesis REFUTED (WP39) |
| **X1** | WITHDRAWN | a no-primer incorporation assay is the decisive experiment | the assay cannot separate retention from allostery; retention has no route to act with no primer | the pulse-chase +/- trap design (WP25, revised by WP37) |
| **X2** | CORRECTED | the TEL patch triad E168/E169/E171 lies on an alpha-helix, so a helical register lattice describes it | psi(E169) is +24.4 deg and CA168-CA171 is 7.68 A in 7TRE and 7.70 A in 2I46, against 5.05 A for an alpha-helix at that offset; the | backbone grafting of the native segments, verbatim (WP45); the CONCLUSION that a helix cannot present the tria |
| **X3** | CORRECTED | no covalent route exists to constrain the grafted TEL segments (WP45, route b) | WP45 screened disulfides only, 3.5-5.0 A Cbeta-Cbeta. Over six chemistries spanning 3.5-12.5 A, W167-L211 satisfies reach, directi | the W167-L211 hydrocarbon staple (WP47); the ENTROPY conclusion that the construct cannot pay for its own orde |
| **X4** | CORRECTED | the graft's ordering cost (12-20 kcal/mol) is a property of the construct that no redesign can reduce | WP45 built the linker from GGGGG. Glycine has no Cbeta and therefore the largest backbone state space of any residue, so WP45 chos | the PAPAP-linked construct DXEEKEFPAPAPFSXLPTE (WP48); the conclusion that the construct still cannot pay IN F |
| **X5** | CORRECTED | the MD3 triad module presents the E168/E169/E171 triad and fails it only on the third contact (WP44: SITE-UNSO | scored on carboxylate spacing, MD3's three display glutamates sit at 6.43, 5.08 and 10.75 A. The measured triad triangle is 6.58,  | the grafted construct, whose spacings ARE the measured triangle (WP45/WP50). The conclusion that a helix canno |

---

## 9. What is missing, stated plainly

1. **No entity has been made.** Every design here is a sequence and a
   geometric specification. Zero syntheses, zero expressions, zero assays.
2. **No binding has been computed or measured.** No docking engine, no MD, no
   FEP. WP39 already refuted this repository's own rigid-body scorer, so any
   score would be fabricated.
3. **No scaffold exists.** The graft needs a presenting scaffold that holds
   6.36 kcal/mol of residual ordering. Finding one needs PDB-wide scaffold
   screening plus loop modelling — Rosetta/MODELLER, not available here.
4. **The deposits are coarse.** TERT has 8.0 heavy atoms per residue; 2I46 has
   no telomerase chain at all. Every clash clearance in this tree is an upper
   bound, and that is now stated wherever one is quoted.
5. **The mechanism question (M5) is still at P = 0.35.** Two branches separated
   by one DNA-trap experiment that has not been run.

---

## 10. Reproducing all of it

```bash
bash scripts/setup_environment.sh        # venv + branch restore + ledger checks
python scripts/solve/wp48_new_designs.py # the newest round
python scripts/solve/wp46_closure.py     # regenerate the ledger
python scripts/solve/wp49_master_record.py
pytest tests/ -q                         # 727 tests
```

See `METHODS.md` for every command, `ARTEFACT_INDEX.md` for every file, and
`NEXT_STEPS.md` for what to order and run.
