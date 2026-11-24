<!-- generated, do not edit by hand: python scripts/solve/wp49_master_record.py -->

# Master record — the TPP1/TERT design programme, end to end

    python scripts/solve/wp49_master_record.py

**Generated, not written.** Every number below is read from an artefact
(`results/solve/wp4*.json`), from the retraction register, or recomputed by the
script that wrote the artefact. Hand-typed summaries in this repository have
drifted three times (C1, C7, R8), so this document does not contain any.

## 0. The answer, first

**No. 8 of 16 sub-problems are solved, and all 8 are solved as CALCULATIONS. Zero are solved as evidence, because nothing has been made or measured. 3 are not solved, 3 are partial, and 2 cannot be answered by any computation here. The largest open item is still P6: WP47 found a real covalent constraint (W167-L211, 6.34 A) where WP45 had concluded none existed, but one staple buys 2.8 kcal/mol against a 12-20 kcal/mol ordering cost. WP48 then improved the construct (PAPAP linker, -2.84 kcal/mol) and still has 6.36 kcal/mol unpaid. A scaffold is still required.**

| status | count |
|---|---|
| SOLVED | **8** |
| PARTIAL | **3** |
| NOT SOLVED | **3** |
| UNMEASURABLE | **2** |

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

## 7. The closure ledger — every sub-problem and its status

| id | area | question | status | kind |
|---|---|---|---|---|
| P1 | Interface | Is the TPP1 TEL patch one binding surface? | **SOLVED** | geometric |
| P2 | Interface | Can a drug-like small molecule span both halves? | **SOLVED** | geometric |
| P3 | Interface | What is the shortest path across the gulf? | **SOLVED** | geometric |
| P4 | Triad | Can an alpha-helix present the E168/E169/E171 triad? | **SOLVED** | geometric |
| P5 | Triad | Can the triad be presented at all, in the native geometry? | **SOLVED** | geometric |
| P6 | Triad | Does a peptide HOLD that geometry? | **PARTIAL** | design-level |
| P7 | Modularity | Can each module solve only its own site, fully? | **PARTIAL** | geometric |
| P8 | Modularity | Do merged modules stay stable together? | **SOLVED** | design-level |
| P9 | Modularity | Can one groove be blocked and repaired at once? | **SOLVED** | logical |
| P10 | Screening | Is a structural alert a usable filter at all? | **SOLVED** | computed |
| P10b | Screening | Does the full gate stack accept approved drugs? | **PARTIAL** | computed |
| P11 | Screening | What is the experimental success rate? | **UNMEASURABLE** | wet-lab |
| P12 | Screening | Do any of these designs bind telomerase? | **UNMEASURABLE** | wet-lab |
| P13 | Mechanism | What does TPP1 change mechanically (M5)? | **NOT SOLVED** | wet-lab |
| P14 | Mechanism | Is n measured? | **NOT SOLVED** | wet-lab |
| P15 | Drug | Is there a drug? | **NOT SOLVED** | design-level |

### What each open item waits on

- **P6 (PARTIAL)** — Does a peptide HOLD that geometry?
  - A covalent constraint EXISTS, and WP45 missed it by screening disulfides only. Over six chemistries spanning 3.5-12.5 A, W167-L211 passes reach, direction, clash and payload tests at 6.34 A with cos 0.98/0.09 - and the native indole/Leu contact there is already 3.50 A. But the staple buys only 2.8 kcal/mol of ordering against a 12-20 kcal/mol cost: residual 9.2 kcal/mol unpaid. Constraint found, energy still not balanced.
  - *blocker:* the remaining 6.36 kcal/mol even after a structured PAPAP linker and a second staple attempt - neither closes it, so a presenting scaffold is still required (PDB screening + loop modelling, unavailable here)
- **P7 (PARTIAL)** — Can each module solve only its own site, fully?
  - 4 of 5 do. One site (S3) was declared unsolvable on a helix model that was wrong - corrected by P4/P5 - so the module result is retained but its premise is superseded.
  - *blocker:* the MD3 module needs rebuilding on the graft (P5), not on a helix
- **P10b (PARTIAL)** — Does the full gate stack accept approved drugs?
  - No. 10/11 end-to-end on the third panel; the one miss is Lipinski rejecting atorvastatin (MW 559, logP > 5). That is a pre-existing gate defect inherited from WP43, not a reactivity failure - it is attributed, not blamed on the layer that was fixed.
  - *blocker:* Lipinski false rejection of large approved drugs; needs a rule that does not reject a marketed molecule at MW 559
- **P11 (UNMEASURABLE)** — What is the experimental success rate?
  - Unknown. N = 0 designs have been made, so the rate has no denominator.
  - *blocker:* synthesis and assay; no amount of computation supplies this
- **P12 (UNMEASURABLE)** — Do any of these designs bind telomerase?
  - Unknown - and it is reported as unknown rather than scored.
  - *blocker:* no docking engine, no OpenMM, no FEP; WP39 already refuted the in-repo rigid-body scorer, so a score would be a fabricated number
- **P13 (NOT SOLVED)** — What does TPP1 change mechanically (M5)?
  - M5 carries P = 0.35. Two branches separated by one DNA-trap experiment.
  - *blocker:* the P2 two-chase experiment has not been run
- **P14 (NOT SOLVED)** — Is n measured?
  - No. The two-chase inversion is exact and unrun.
  - *blocker:* wet lab
- **P15 (NOT SOLVED)** — Is there a drug?
  - No, and the small-molecule route for the bipartite patch is excluded (P2). The protein route is unmeasured (P11, P12).
  - *blocker:* P6 for the graft; P11/P12 for everything else

---

## 8. Every retraction

**20 registered.** The register is enforced: a marker removed from a
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
