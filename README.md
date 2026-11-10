# Telomerase Allosteric Switch (`telomerase-allosteric-switch`)

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Test Suite](https://img.shields.io/badge/pytest-587%20passed-success.svg)](tests/)
[![Code Hygiene](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Simulation Engine](https://img.shields.io/badge/OpenMM-8.6.1-orange.svg)](https://openmm.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Manuscript](https://img.shields.io/badge/manuscript-identifiability%20analysis-purple.svg)](paper/MANUSCRIPT.md)
[![Status](https://img.shields.io/badge/Status-Mechanism%20not%20proven-model%20M5-yellow.svg)](results/solve/WP22_RAP_MODEL_V2.md)

**A structural and experimental-evidence audit of how TPP1 stimulates human telomerase repeat-addition processivity — including the record of which claims in this repository failed and why.**

## Where the open questions stand

A single audited ledger of every gap, each checked against the tree rather than
recalled: [`results/solve/OPEN_GAPS.md`](results/solve/OPEN_GAPS.md).

Four of the seven gaps that were open are now closed, and the three that remain
are blocked for stated reasons rather than neglected:

- **Gate 3 (densitometry)** — the repository's data cannot satisfy it, and this is
  now a *proof*: two product-length distributions can share the recorded ≥16-repeat
  fraction exactly while differing 4.98× in shape. See
  [`WP24`](results/solve/WP24_GATE3_STATUS.md).
- **The 2.73× vs ≥4× tension** — reduced to two branches separated by a single
  DNA-trap experiment with distinct numeric predictions. See
  [`WP25`](results/solve/WP25_RESOLVING_TESTS.md).
- **Simulation-only closure** — quantified. The repo's MD is 3.0 ps on a 97-residue
  fragment, ~8 orders of magnitude short, and omits the active site. The structural
  half is closable (alchemical free energy first); the kinetic half is not closable by
  any amount of compute. See [`WP26`](results/solve/WP26_SIMULATION_READINESS.md).
- **9SHZ / 9SI0** — parser written and verified offline; the sandbox blocks TLS.
  One networked run of `scripts/solve/fetch_rap_fixture.py` finishes it.

Newly closed: the **TPP1-free matched series**. 7TRD is apo and matches 7TRE at
934/934 residues; the largest domain displacement on TPP1 binding is **0.53 Å**,
i.e. at the coordinate uncertainty — a *measured null*, not an untested
assumption. See [`WP23`](results/solve/WP23_TPP1_FREE_SERIES.md).

## Manuscript

**What Processivity Distributions Cannot Tell Us About TPP1** — an identifiability analysis of
telomerase repeat-addition processivity, with a corrected treatment of heterogeneous populations.

**Beamlak Sintayehu** · 2026-11-01 · [`paper/MANUSCRIPT.md`](paper/MANUSCRIPT.md) · [figures](paper/figures/README.md)

The paper asks what a processivity measurement can *identify* about TPP1's mechanism, and answers:
less than is commonly assumed. It has four moving parts —

1. **An extension-rate mechanism and a retention mechanism are the same function of the
   observable.** Both act only through the per-repeat survival probability *s*, so under matched
   *s* they produce identical distributions. A distribution measured at one primer concentration
   carries no mechanistic information beyond *s*.
2. **No absolute kinetic rate is recoverable from processivity**, at any precision: every
   observable is invariant under scaling all rates by a common factor (symbolically, exactly zero).
   A processivity design has structural rank 3 of 4; direct primer binding is the only independent
   measurement that adds rank.
   > **Correction (WP19).** "Direct primer binding has never been measured" was wrong. Latrick &
   > Cech 2010 measured it twice by pulse-chase: the *resting* `k_off` = 0.012 /min and is
   > **unchanged** by saturating POT1-TPP1, while the **cycling** complex's half-dissociation time
   > increases **at least four-fold**. They also report `k_pol` unchanged. The open item is
   > therefore narrower and sharper: the cycling measurement needs a DNA trap or dilution to
   > remove rebinding (WP13), and it is that corrected number which would add rank.
   > **Correction (WP22).** Two further errors in how that literature was used, now fixed.
   > (i) The **2.7× is a ratio of R₁/₂ — a median length** — not of per-repeat survival `s`.
   > Latrick & Cech give "half the chains elongated by eight repeats (48 nt)", so R₁/₂ = 8 and
   > the `s` ratio is **1.159**, not 2.700. The two differ by **2.33×**, and WP3/WP18/WP20 all
   > compared them directly. (ii) Wang 2007's absolute anchor comes from RRL-reconstituted
   > enzyme while L&C's 2.7× comes from HEK "super-telomerase" — **different enzymes**, which
   > the model had multiplied into one equation. Scored separately, the repository's proposed
   > form `s = (1+a)⁻⁶·φ` is **falsified**: on L&C's own no-POT1-TPP1 numbers
   > `s₀/φ₀ = 0.7914/0.57 = 1.388 > 1` forces a **negative** rate constant. `φ` is a
   > window-limited efficiency, not a survival probability, and must not multiply into `s`.
   > The corrected model **M5** — `s = (1+a)⁻⁶`, with stalled complexes retrying — is physical
   > and predicts a **2.73×** off-rate fall against L&C's independently measured **≥4×**.
   > See [`results/solve/WP22_RAP_MODEL_V2.md`](results/solve/WP22_RAP_MODEL_V2.md).
3. **Affinity is not one hypothesis, and this was wrong in the earlier draft.** Lowering *k*_off
   and raising *k*_on at matched *K*_d are *not* equivalent. The loss term rewrites exactly as
   *R* = *k*_off²/(*k*_off + [P]·*k*_on), so the two differ by an exact factor λ — which is
   algebraically confounded with *k*_loss, the parameter a mechanism claim would rest on.
4. **The tail-slope estimator is window-dependent, and that disqualifies it.** For a two-species
   mixture the adjacent-band ratio has a closed form; it is bounded by the most extreme
   subpopulation, it converges to the **longest-lived** species rather than the majority, and a
   minority takes over beyond a computable length. In a worked case a true *s* = 0.60 reads as 0.61
   or 0.95 depending only on the band window — a **12-fold** change in the implied *k*_loss/*k*_ext
   from identical data. Every mechanism claim must therefore come from a full-distribution fit with
   the mixture weight free.

**The distinction the paper exists to protect:** a mathematically exact result within a model is not
an experimentally established biological fact. Everything above is of the first kind. **No mechanism
is demonstrated here**, and the next real step is real per-band densitometry of an existing gel —
not another simulation.

Reproduce the paper's every number:

```bash
python scripts/mechanism_classification_checks.py   # ~2 s, writes JSON + CSV
python paper/make_figures.py                        # all six figures
pytest tests/test_paper_manuscript.py -q            # 29 guards: fails if a number or claim drifts
```

> ## ⚠️ Status (2026-11-01): audit complete, identifiability established, no mechanism proven
>
> The question this repository exists to answer is: **does TPP1 stimulate processivity through a
> local structural conduit, through broad conformational allostery, or through
> affinity/anchoring that requires no mechanical propagation?** The audit that answered as much as
> the published data allow is in [`results/mechanism/`](results/mechanism/MECHANISM_ASSESSMENT.md):
>
> | Work package | Deliverable | Result |
> |---|---|---|
> | WP1 structural validation | [`structural_validation.md`](results/mechanism/structural_validation.md) | Recomputed from raw mmCIF; **retracts the "one continuous IFD-TRAP helix" claim** |
> | WP2 experimental evidence | [`experimental_evidence.md`](results/mechanism/experimental_evidence.md) | Recruitment is separable from processivity; a mechanical signal is not |
> | WP3 model comparison | [`model_comparison.md`](results/mechanism/model_comparison.md) | Affinity/anchoring (C) best supported; conduit (A) contradicted; allostery (B) untested |
> | WP4 decision matrix | [`experimental_decision_matrix.md`](results/mechanism/experimental_decision_matrix.md) | Two cheap experiments would settle it |
> | WP5 assessment | [`MECHANISM_ASSESSMENT.md`](results/mechanism/MECHANISM_ASSESSMENT.md) | Findings ranked by confidence, with falsification criteria |
>
> ### WP16 — the TEN anchor site has now been measured (2026-11-06)
>
> The long-standing structural gap is closed. The Sekne 2022 depositions with a **30-nt** DNA put
> **Tyr176 within 2.75 Å (7QXS), 4.35 Å (7QXA) and 2.78 Å (9SHY)** of the DNA — a backbone contact —
> against **11.11–11.94 Å** in this repository's 7-nucleotide structures. All three pass a
> catalytic control (Asp868 within 2.65 / 2.93 / 2.78 Å of the DNA 3′ end), so their active sites
> are verifiably seated. Two independently refined models agree on the same residue and the same
> nucleotide.
>
> **9SHY was previously reported as `NOT RESOLVED` at 30.23 Å.** That was a bug, not a result: the
> old code measured to the *furthest* modelled nucleotide instead of the nearest, and read the
> `NOT RESOLVED` branch while sitting 2.78 Å from a contact. The row now reads CONTACT, and
> `tests/test_open_findings.py` fails if a single row ever reports a contact distance and
> `NOT RESOLVED` at the same time again.
>
> **This settles geometry, not mechanism.** The contact is to the sugar–phosphate backbone, so it
> cannot be reading sequence, and it is equally consistent with a positioning element and with a
> steric barrier. It does not distinguish DNA retention from catalytic activation, and it does not
> change the standing mechanism verdict. See
> [`results/solve/WP16_ANCHOR_SITE_RESULT.md`](results/solve/WP16_ANCHOR_SITE_RESULT.md).
>
> ### WP17 — R1 was never tested, and the obvious experiment does not fix that (2026-11-06)
>
> **"R1 excluded, R2 supported" was too strong.** The fit *sets* the resting primer off-rate equal
> with and without TPP1 (`resting_koff_no_PT == resting_koff_PT == 0.012` in
> `wp3_predictions.json`) — R1 is assumed away *inside the model* — and the one control that exists points the other
> way: the resting primer `k_off` is **measured at 0.012 /min and is unchanged by saturating
> POT1-TPP1**, which is evidence *against* a global affinity change but does not by itself
> exclude R1 during cycling.
> The structural argument only rules out a direct TPP1–DNA contact.
>
> The natural fix — titrate primer until R1 and R2 diverge — **fails**, and it fails in a way that
> is worth stating plainly. Under eq. 1 both arms act on the same scalar `s`, so they are
> degenerate at A₀ = 0.99 (gap < 0.0005 fold). They do separate as processivity falls, but only
> because *both* effects collapse toward 1.0×:
>
> | A₀ | best \[P\]/K_d | R1 | R2 | R2/R1 |
> |---|---|---|---|---|
> | 0.99 | 0.01 | 1.0004 | 1.0075 | 1.007 |
> | 0.90 | 0.01 | 1.0038 | 1.0803 | 1.076 |
> | 0.70 | 0.01 | 1.0115 | 1.2877 | 1.273 |
>
> At A₀ = 0.70 the model's **largest possible** TPP1 effect is 1.29×. The effect to be explained
> is **2.7 ± 0.6×**. The separating regime is not the regime where the phenomenon lives.
>
> **What would settle it** is a variable that enters one arm and not the other — vary `k_pol` at
> fixed `K_d` (R1 predicts the effect is invariant; R2 predicts it moves), measure `K_d` directly
> ± TPP1, or use `Y176A` to remove the anchor without touching the catalytic step. Designs, power
> and falsification criteria are in
> [`results/solve/WP17_MECHANISM_DESIGN.md`](results/solve/WP17_MECHANISM_DESIGN.md).
>
> **Also fixed.** `compare_states()` superposed on RT 541–935 while reporting IFD-TRAP 755–785,
> which sits entirely inside that range — the fit absorbed the motion it was measuring (~16 %
> attenuation, ~0.5 Å of invented motion in the "rigid" ring). The frame is now the invariant RT
> core and IFD-TRAP is defined as **755–785** in both `colab/` and `src/` (they previously
> disagreed: 745–805 vs 755–785).
>
> **RAP states (9SHZ / 9SI0) — status changed.** These are no longer unobtainable. Both entries
> are served by `models.rcsb.org`, and `scripts/solve/fetch_rap_fixture.py` now fetches their CA
> atoms in the same fixture format as 7TRE/7TRD/7TRF/2I46/9SHY. The parser is tested offline
> against a live-verified single-atom pull (Y176 CA = 133.669, 153.363, 143.095 Å,
> `tests/test_rap_fixture_fetch.py`). The remaining step is one network-enabled run of that
> script; the sandbox this was developed in blocks TLS egress from Python, so **the fixtures are
> fetched but not yet ingested**, and no per-domain displacement for those states is claimed here.
>
> Reproduce, offline, in one command:
>
> ```bash
> bash colab/run_colab.sh --selftest   # three suites, no internet, ~5 s
> ```
>
> The pipeline an agent should run is documented in
> [`colab/AGENT_START_HERE.md`](colab/AGENT_START_HERE.md); `colab/test_colab.py` asserts that
> document still matches the code.

> **Next phase (open): affinity–processivity causality.** Does primer retention quantitatively
> account for the TPP1-dependent processivity increase, without a mechanical signal? The
> [preregistration](results/causality/PREREGISTRATION.md) is frozen before any fitting, with its
> source-data provenance and SHA-256 in [`data/external/tesmer2023/`](data/external/tesmer2023/).
> Two constraints were found at the gate: the numeric source data contain **no product-length
> distribution** (one scalar per replicate), and **n = 3** with CVs of 35–90% can only resolve
> ~2.2× differences in processivity. Neither is a reason to fit anyway.
>
> **That phase is now closed, in the negative.** Under matched effective survival the
> extension-rate model and the affinity/retention model are **algebraically identical** — no
> processivity observable distinguishes them — and the published numbers support only a coarse
> scalar comparison (n = 3, CVs 35–90%, ~2.2× resolution). The result is recorded in
> [`MECHANISM_CLASSIFICATION.md`](results/causality/MECHANISM_CLASSIFICATION.md), which classifies
> every parameter at four levels (molecular event / microscopic rate / effective rate / observable)
> so that no effective parameter is renamed into a mechanism.
>
> The governing question is now narrower and belongs to experiment, not to another simulation:
> **does TPP1 change how tightly telomerase retains the primer, how quickly it extends the primer,
> or both?** The classification document states which measurements are needed, what each would
> refute, and which parameters are provably unrecoverable from processivity data at any precision.
> Nothing in that phase upgrades a mechanistic claim: **no mechanism is proven**, and the progress
> it records is methodological. Writing the document also caught three errors in its own first draft —
> a false claim that the two affinity variants are equivalent at matched *K*_d, a sign error in the
> crossing formula, and an unreproducible saturation bias — all now guarded by tests, not by trust.
>
> **That analysis is now written up as a paper** — [`paper/MANUSCRIPT.md`](paper/MANUSCRIPT.md) —
> which also proves three results the classification document only sketched: the mixture ratio
> *r*(*n*) is bounded by min(*s*,*s*₂) and max(*s*,*s*₂); the tail slope converges to the
> **longest-lived** species; and a minority contaminant takes over the signal beyond a computable
> length *n*\* = 1 + ln(1/κ)/ln(ρ). The consequence is the sharpest thing in the repository: for a
> true *s* = 0.60 the same distribution reads as 0.61 or 0.95 depending only on which band window is
> used, a 12-fold change in the implied *k*_loss/*k*_ext. **The next step is densitometry of an
> existing gel, not another model.**
>
> **Four things this repository keeps separate** (they were previously merged, which is how a
> distance gradient became a "mechanism"):
>
> 1. The original **TEN–TPP1–TRAP bipartite clamp** hypothesis — refuted on its structural
>    premise: the TEN domain makes **zero** primer contacts (Tyr176 is 11.1 Å away) and TPP1 is
>    28.7 Å from the DNA.
> 2. The **structural audit** that invalidated it — stands, and is now extended: 752–65 is a
>    β-strand, 766–67 a single turn, 768–71 a loop with P771, 772–80 helix α25. There is no
>    single helix between the DNA anchor and the TPP1 footprint.
> 3. The **conduit hypothesis (760–766)** built on the audit — not supported: only **L766** of
>    that window is at the TPP1 interface (3.6 Å), and **H762E** inside it retains 86% of WT
>    processivity. Retired.
> 4. The **ANM coupling result** — **inconclusive, not evidence against allostery.** Its "specific"
>    effect (−0.041) equals the short-range positive control (+0.041) and is +0.34σ from a
>    distance-matched null. It is retired as evidence.
>
> Two facts from the published literature carry the assessment: (i) adding the TPT complex
> **rescues** — rather than fails to stimulate — telomerase mutants with a damaged DNA anchor
> (Sekne et al. 2022), the signature of affinity compensation; (ii) processivity-defective
> IFD-TRAP mutants are still **recruited to telomeres normally** (Tesmer et al. 2023), so
> recruitment is a separate function. **What remains untested is whether any mechanical component
> exists at all**, because no published perturbation separates it from DNA affinity.
>
> The **pytest badge (522 passed)** is verified on a clean clone with `pip install -e .`; the
> previously advertised "342 passed" and then "357 passed" were both reachable only after fixing
> fixtures, and are superseded.
>
> **Historical record, retained deliberately:** the older banner below this section, and sections
> 1–5, describe claims that have since been retracted (K757A "breakthrough", 25 Å allosteric
> switch, MD-based mechanism, TELIN-1 ligand). They are kept because deleting a negative result is
> worse than publishing it. [`results/verification/CLAIM_VERIFICATION_REPORT.md`](results/verification/CLAIM_VERIFICATION_REPORT.md)
> records 38 verified / 14 contradicted / 2 untestable / 2 partial of 56 headline claims.

---

## 1. Scientific Overview & Central Discovery Hypothesis *(historical — retracted, see the status banner)*

For over two decades, the structural mechanism by which the shelterin component **POT1–TPP1** stimulates human telomerase **Repeat Addition Processivity (RAP)** while tethering the enzyme to chromosome ends has remained an unsolved problem in telomere biology. While the TPP1 TEL-patch is known to recruit telomerase to telomeres, how recruitment translates across distance into enhanced primer retention during translocation has lacked structural explanation.

This repository hosts the computational discovery, all-atom molecular dynamics validation, and pre-registered biochemical testing protocols for the **TEN–TPP1–TRAP Bipartite Clamp Mechanism**:

```
           [ Telomere Tethering: TPP1 TEL-Patch ]
                         │
                         ▼
             TERT TEN Domain (Tyr176)
                         │  (Inter-domain cleft compaction)
                         ▼
             TERT TRAP Domain (Gln794 ── Lys757)
                         │
                         ▼
    [ DNA Primer Exit Channel Retention Rheostat ]
```

### The Primary Candidate: TERT `Lys757A` — *superseded 2026-10-30, see §0 banner*
> **Prior-art caution.** TERT Lys757 is part of the **DNA anchor site**, and the K757A/F759A double
> mutant has already been reconstituted and assayed: it *reduces* both activity and RAP (Sekne et
> al., *Science* 2022, Fig. 4E–G). The K757E/F759E mutants are likewise defective in activity and
> processivity while still localising to telomeres (Tesmer et al., *Life Sci. Alliance* 2023).
> The prediction below ("preserves basal catalysis") is therefore expected to be false, and this
> section is retained only as a record. The corrected, falsifiable claim is in
> [`results/verification/corrected_breakthrough_claim.md`](results/verification/corrected_breakthrough_claim.md).

We propose that **TERT Lys757** operates as the critical mechanical hinge of an allosteric lever:
- **Separation-of-Function Hypothesis**: Mutation `Lys757A` selectively uncouples TPP1-stimulated processivity while preserving high-affinity TPP1 recruitment and basal catalytic polymerization.
- **Stereochemical Rescue Control**: Mutation `Lys757R` preserves positive charge and sidechain reach, predicting full rescue of processivity stimulation ($\chi_{\text{RAP}} \ge 2.0$).
- **Therapeutic Implication**: Tuning processivity without shutting down basal catalysis creates a potential therapeutic window to starve rapidly dividing cancer cells of telomere length while sparing quiescent tissue stem cells.

---

## 2. Computational Audit & Molecular Dynamics Findings

All-atom molecular dynamics simulations were conducted in **OpenMM 8.6.1** (AMBER `ff14SB` force field, OBC2 implicit solvent, Langevin Middle Integrator, $T = 300\text{ K}$) across independent replicas ($N = 3$) for both the **TPP1-Bound** (`PDB 7TRE`) and **Apo** (`PDB 7TRD`) states.

To prevent pseudoreplication, independent simulation replicas were used as the statistical unit ($df = 4$):

| Biophysical Metric | TPP1-Bound WT ($N = 3$) | Apo-Unbound WT ($N = 3$) | Statistical Test ($df = 4$) | Finding & Status |
| :--- | :--- | :--- | :--- | :--- |
| **Cleft Distance** (`Tyr176:OH` $\leftrightarrow$ `Gln794:OE1`) | **$3.820 \pm 0.633\text{ \AA}$** | **$5.620 \pm 0.196\text{ \AA}$** | Two-sample $t = -4.705$, **$p = 0.00928$** (Welch $p = 0.0299$) | **Statistically Significant Cleft Compaction** ($p < 0.01$) |
| **H-Bond Occupancy** ($\le 3.2\text{ \AA}$) | **$32.2\% \pm 28.7\%$** | **$0.0\% \pm 0.0\%$** | Two-sample $t = 1.946$, **$p = 0.1236$** | **Non-Significant Latch Difference** ($p > 0.05$) |
| **Cα RMSD Plateau** | $1.32 \pm 0.13\text{ \AA}$ | $1.31 \pm 0.06\text{ \AA}$ | Stable structural plateau | Rapid convergence post-relaxation |

### Critical Scientific Caveats & Epistemic Boundaries
1. **Simulation Duration (3.0 ps)**: The 1,500-step simulations with $\Delta t = 2.0\text{ fs}$ capture **$3.0\text{ ps}$** of physical time ($2.4\text{ ps}$ production post-burnin). This timescale samples local atomic relaxation but is several orders of magnitude below biological turnover ($\approx 0.5\text{--}2.0\text{ s}$ per repeat addition).
2. **Cleft Compaction vs Permanent Latch**: TPP1 binding is dynamically associated with cleft compaction, but **a persistent Tyr176–Gln794 hydrogen-bond latch has not been established in silico**. In Replicate 103, `Gln794` sampled an open sidechain rotamer ($4.48\text{ \AA}$), reducing occupancy to $0\%$ in that trajectory window.
3. **Primacy of Wet-Lab Falsification**: Computational models only generate testable hypotheses; definitive proof requires physical biochemical evidence.

---

## 3. Structural Models & Telin-1 Inhibitor Clash Falsification

In addition to apo and bound holoenzyme reconstructions, this repository archives the small-molecule inhibitor models investigated during the campaign:

1. **TELIN-1 Small Molecule Inhibitor** ([`data/ligands/telin1_inhibitor.pdb`](data/ligands/telin1_inhibitor.pdb)):
   - 36 heavy atoms ($\text{C}_{26}\text{N}_2\text{O}_6\text{SF}$).
   - Adversarial coordinate audit revealed chemical anomalies (disconnected fragments, stretched bonds, planar distorted geometry).
2. **TELIN-1 Docked Complex** ([`data/ligands/7tre_telin1_complex.pdb`](data/ligands/7tre_telin1_complex.pdb)):
   - Docked into cryo-EM structure 7TRE.
   - Comprehensive clash audit ([`scripts/verify_telin1.py`](scripts/verify_telin1.py)) revealed:
     * **160 steric clashes $< 2.80\text{ \AA}$**
     * **21 severe clashes $< 1.50\text{ \AA}$**
     * **8 extreme overlaps $< 1.00\text{ \AA}$**
     * **Minimum interatomic distance of $0.17\text{ \AA}$** directly penetrating the backbone nitrogen of TERT Ala130.
   - Complete absence of TPP1 (Chain C) in the model proved that static docking attempts severely disrupted recruitment rather than modulating processivity, directly motivating the pivot toward allosteric targeting of `Lys757`.
   - Verified by automated test suite: [`tests/test_ligand_telin1.py`](tests/test_ligand_telin1.py).

---

## 4. Sequence-Level Verification & Construct Design

Canonical sequence audit against **NCBI `NM_198253.3`** and **UniProt `O14746`** (3,399 bp CDS, 1,132 aa):
- **Residues 753–761**: $\text{Gly753 - His754 - Val755 - Arg756 - \textbf{Lys757} - Ala758 - Phe759 - Lys760 - Ser761}$ (`G-H-V-R-K-A-F-K-S`).
- **Codon 757**: Wild-type is **`AAG`**.
- **Mutagenesis**:
  - `Lys757A`: `AAG` $\to$ `GCG` (creates novel diagnostic **BstUI** `CG^CG` and **HinP1I** `G^CGC` sites in amplicons).
  - `Lys757R`: `AAG` $\to$ `CGG` (creates novel diagnostic **MspI / HpaII** `C^CGG` sites).
- Full cloning and sequencing specifications: [`CONSTRUCT_DESIGN_WT_K757A_K757R.md`](results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md).

---

## 5. Pre-Registered Breakthrough Evidence Criteria

To eliminate confirmation bias, the experimental decision framework is cryptographically frozen:

```
+─────────────────────────────────────────────────────────────────────────────────────────────────────────+
│                                  PRE-REGISTERED WET-LAB DECISION CRITERIA                               │
+─────────────────────────────────────────────────────────────────────────────────────────────────────────+
│ 1. RECRUITMENT AFFINITY (Biacore 8K SPR):                                                               │
│    - Kd(Lys757A) <= 300 nM (equivalent to WT reference: 150 ± 25 nM; p > 0.05).                        │
│ 2. BASAL CATALYTIC ACTIVITY (Primer Extension):                                                         │
│    - k_cat(Lys757A) >= 85% of WT without TPP1 (p > 0.05 vs WT; p < 0.001 vs catalytic dead Asp712A).    │
│ 3. PROCESSIVITY SELECTIVITY & RESCUE (12% Urea-PAGE RAP Gel):                                           │
│    - chi_RAP(Lys757A) <= 1.8 fold stimulation (statistically impaired vs WT ~2.4 fold; p < 0.001).     │
│    - chi_RAP(Lys757R) >= 2.0 fold stimulation (statistically rescued; p > 0.05 vs WT).                 │
│ 4. MECHANISTIC CONFIRMATION (Pulse-Chase Primer Retention):                                             │
│    - Complex half-life t_1/2(Lys757A + TPP1) <= 3.5 min (primer dissociation not suppressed).          │
+─────────────────────────────────────────────────────────────────────────────────────────────────────────+
```

---

## 6. Repository Directory Architecture

```
.
├── paper/                               # CURRENT deliverable: the manuscript
│   ├── MANUSCRIPT.md                    # full paper — identifiability of TPP1's mechanism
│   ├── make_figures.py                  # regenerates all six figures from the same equations
│   ├── figures/                         # fig1-fig6, PNG (300 dpi) + PDF, with legends index
│   └── COVER_LETTER.md                  # covering-note template (NOT sent; addresses unverified)
│
├── results/causality/                   # the analysis the paper is built on
│   ├── MECHANISM_CLASSIFICATION.md      # four-level parameter classification + release gates
│   ├── classification_checks.json       # every quoted number, machine-readable
│   └── PREREGISTRATION.md               # frozen before any fitting, with source provenance
│
├── results/mechanism/                   # the earlier mechanism audit (WP1-WP5)
│   ├── structural_validation.md         # WP1: interface, numbering, contacts, 2ry structure, retractions
│   ├── experimental_evidence.md         # WP2: what published experiments do and do not establish
│   ├── model_comparison.md              # WP3: three models scored against the evidence
│   ├── experimental_decision_matrix.md  # WP4: experiments ranked by discriminating power
│   └── MECHANISM_ASSESSMENT.md          # WP5: findings ranked by confidence + falsification criteria
│
├── docs/                                # Technical architecture, methodology, & SOP manuals
│   ├── ARCHITECTURE.md                  # Comprehensive software design, pipeline flow, & data contracts
│   ├── METHODOLOGY.md                   # Biophysical derivations, OpenMM setup, & kinetic equations
│   ├── EXPERIMENTAL_SOP_MANUAL.md       # Unified wet-lab operational manual (Expi293F, SPR, RAP)
│   └── TELOMERASE_UNSOLVED_FRONTIER.md  # Scientific foundation report on telomerase recruitment
│
├── src/
│   └── telotrade/                       # Core structural bioinformatics engine
│       ├── analysis.py                  # Structure alignment, Cα RMSD, and residue mapping
│       ├── cross_state.py               # Conformational transition analysis (7TRE vs 7TRD/7TRF/2I46)
│       ├── interface_geometry.py        # Buried Surface Area (BSA) and contact geometry
│       ├── metrics.py                   # Quantitative structural metric calculation
│       ├── sequence_mapping.py          # Sequence alignment & UniProt cross-referencing
│       └── validate.py                  # SHA-256 data manifest and coordinate validator
│
├── scripts/                             # Automated analysis, simulation, & audit scripts
│   ├── spr_kinetic_analysis.py          # 1:1 Langmuir/MTL kinetic fitting engine (Biacore 8K)
│   ├── rap_gel_quantification.py        # Phosphorimager 12% Urea-PAGE densitometry engine
│   ├── recalculate_trajectory_statistics.py # Replica-level raw MD statistics recalculator
│   ├── audit_openmm_equilibration.py    # OpenMM 1,500-step burn-in equilibration audit
│   ├── export_raw_trajectories.py       # Raw trajectory timeseries generator and exporter
│   ├── generate_independent_benchmarks.py # Independent SPR and RAP benchmark dataset generator
│   ├── verify_telin1.py                 # Coordinate clash and chemical integrity validator for Telin-1
│   ├── structural_validation.py         # WP1: independent gemmi/numpy recomputation of the TPP1 interface
│   ├── allosteric_coupling.py           # ANM engine (result: inconclusive; retired as evidence)
│   ├── mechanism_classification_checks.py # symbolic + complex-step verification of every claim
│   └── generate_code_freeze_manifest.py # Computes SHA-256 manifest and creates freeze document
│
├── tests/                               # Comprehensive test suite (522 passing tests)
│   ├── e2e/                             # End-to-end pipeline execution tests
│   ├── test_adversarial_b2_challenger.py# Adversarial superimposition and SVD stress tests
│   ├── test_adversarial_m1.py ... m5.py # Data integrity, corrupted CIF, and figure generation tests
│   ├── test_experimental_scripts.py     # 10 unit tests for SPR, RAP, and MD analysis pipelines
│   ├── test_ligand_telin1.py            # 2 unit tests verifying Telin-1 clashes and composition
│   ├── test_mechanism_classification.py # 25 guards on the classification document
│   ├── test_paper_manuscript.py         # 29 guards: manuscript numbers must match the code
│   └── test_breakthrough_investigation.py # Invariant tests for allosteric models
│
├── data/
│   ├── raw/                             # Verified wwPDB cryo-EM & crystal structures (7TRE, 7TRD, 7TRF, 2I46)
│   ├── ligands/                         # Small-molecule inhibitor and docked complex PDB models
│   ├── structures/                      # Curated 3D coordinate models (junctions and complexes)
│   ├── benchmarks/                      # Physical independent benchmark datasets (SPR & RAP)
│   ├── trajectories/                    # Raw step-by-step MD timeseries CSVs (6 replicas, 3.0 ps)
│   └── data_manifest.csv                # Authoritative SHA-256 checksum registry
│
├── results/
│   ├── breakthrough_investigation/      # Computational breakthrough analysis & manuscript draft
│   │   ├── MANUSCRIPT_FOR_PUBLICATION.md# Complete scientific manuscript draft
│   │   └── ...                          # Structural reports, rankings, and figure legends
│   └── experimental_validation/         # Wet-lab operational packages & protocols
│       ├── PRE_REGISTRATION_CODE_FREEZE.md # Cryptographic SHA-256 freeze manifest (33 components)
│       ├── RECONSIDERED_POWER_ANALYSIS.md  # Statistical power plan (N=4-5, CV=15-25%)
│       ├── LABORATORY_EXECUTION_CHECKLIST.md # Phases 1-7 operational checklist
│       ├── CONSTRUCT_DESIGN_WT_K757A_K757R.md # Verified primers & restriction audit
│       ├── RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md # Expi293F expression SOP
│       ├── SPR_BINDING_DATA_PACKAGE.md     # Biacore 8K multi-cycle protocol
│       ├── QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md # 12% Urea-PAGE densitometry SOP
│       └── MOLECULAR_DYNAMICS_TRAJECTORY_REPORT.md # All-atom MD replica analysis report
│
├── pyproject.toml                       # Python package configuration
├── requirements.txt                     # Pinned project dependencies
├── pytest.ini                           # Test runner configuration
├── LICENSE                              # MIT License
└── README.md                            # Comprehensive project documentation
```

---

## 7. Installation & Quickstart

### Prerequisites
- Python 3.10 or higher
- Git and C++ compiler (for Gemmi and OpenMM)

### Installation
```bash
# Clone the repository
git clone https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch.git
cd telomerase-allosteric-switch

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .
```

> **`freesasa` is optional.** It is a C extension with no wheel on several platforms and
> needs CPython headers (`python3-dev`) to build. The code does not require it:
> `interface_geometry.py` imports it inside a `try/except` and falls back to Biopython's
> `ShrakeRupley`, and the reference SASA values are cross-checked against that fallback.
> `pip install -e .` therefore succeeds without it. If you do want the fast C path:
> `pip install -e ".[freesasa]"` (add `--no-build-isolation` and install your platform's
> Python headers if the build fails).
>
> `openmm` and `pdbfixer` in `requirements.txt` are only used by the MD-preparation scripts
> (`scripts/prepare_md_systems.py`); neither the library nor the test suite imports them.
> They can be skipped on a lean install.

### Reproduce Key Analyses

```bash
# 0. The manuscript: verify every number, then rebuild the figures
python scripts/mechanism_classification_checks.py   # ~2 s -> classification_checks.json/.csv
python paper/make_figures.py                        # fig1-fig6, PNG + PDF

# 1. Run full project test suite
pytest

# 2. Run independent SPR kinetic fitting engine on benchmark data
python scripts/spr_kinetic_analysis.py --input data/benchmarks/spr_independent_benchmark.csv

# 3. Run quantitative RAP gel densitometry engine on benchmark data
python scripts/rap_gel_quantification.py --input data/benchmarks/rap_independent_benchmark.csv

# 4. Recalculate replica-level statistics directly from raw MD trajectories
python scripts/recalculate_trajectory_statistics.py

# 5. Run Telin-1 inhibitor clash verification script
python scripts/verify_telin1.py

# 6. Run end-to-end structural bioinformatics pipeline from scratch
python run_pipeline.py
```

---

## 8. Quality Assurance & Verification

The codebase maintains strict engineering standards:
- **Test Suite**: full suite passing across 24 test modules; the manuscript's own 29 guards
  are in `tests/test_paper_manuscript.py`, and the classification document's 25 in
  `tests/test_mechanism_classification.py`. Both fail if a quoted number stops matching
  the code.
- **Linting & Formatting**: `ruff check src` passes clean, and so do the four modules added for
  the manuscript (`paper/make_figures.py`, `scripts/mechanism_classification_checks.py`, and
  the two guard modules). The older `tests/` and `scripts/` files carry a number of
  pre-existing warnings (unused variables, missing `check=` on `subprocess.run`) that are
  listed rather than silenced; the previous "zero errors" claim here was not accurate.
- **Pre-Registration Integrity**: 33 components locked with SHA-256 hashes in [`PRE_REGISTRATION_CODE_FREEZE.md`](results/experimental_validation/PRE_REGISTRATION_CODE_FREEZE.md).

---

## 9. License & Citation

Distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.

If you build upon this computational model or experimental protocol, please cite:
```bibtex
@article{telomerase_allosteric_switch_2026,
  title={An Allosteric Lever in Human Telomerase Reverse Transcriptase Uncouples Repeat Addition Processivity from Shelterin Recruitment},
  author={Telotrade Investigation Team},
  year={2026},
  journal={Pre-Registration Archive},
  note={GitHub: Mr-Nobody-Anonymous/telomerase-allosteric-switch}
}
```
