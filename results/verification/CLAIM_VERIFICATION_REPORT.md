# Claim-by-Claim Verification of the "Telomerase Allosteric Switch" Project

**Date:** 2026-10-30
**Scope:** every quantitative and mechanistic claim in `README.md`, `report.md`,
`results/breakthrough_investigation/FINAL_RESEARCH_REPORT.md`, `results/audit/evidence_ledger.csv`,
the TELIN-1 ligand files, and the MD trajectory data — re-derived from the deposited coordinates
in `data/raw/`, not from the project's own output tables.
**Reproduce:** `python scripts/verify_all_claims.py`
→ `results/verification/claim_recheck.csv` (56 claims, machine-readable)
→ `results/verification/verification_log.txt` (every intermediate number)

---

## 0. Bottom line

| # | Finding | Consequence |
|---|---|---|
| 1 | **The geometry engine is sound.** BSA (856.82 Å²), the contact networks (9/46/123 atom contacts; 5/21/32 residue pairs), the tier partition (5/16/11), the core RMSDs (0.338 / 0.356 / 0.434 Å), the TPP1 conformational shifts (0.512 / 4.293 / 3.809 / 3.320 Å) and the R774–E168 salt bridge (2.320 Å) all reproduce to ≤0.1 % from scratch. | 38 of 56 claims verified. These are real measurements. |
| 2 | **14 claims fail on re-derivation**, mostly where atom pairs were mislabelled or two different numbers were quoted for the same quantity. | See §2, §4.1. |
| 3 | **The "state-dependent switch" is a snapshot difference dressed as a mechanism.** Y176:OH–Q794:OE1 is genuinely 2.338 Å in 7TRE and 4.372 Å in 7TRD, but the B-factor argument used to support it is invalid: 7TRE is globally ADP-scaled 5× lower than 7TRD (mean 16.5 vs 85.9 Å²). | §4.2. The distance is a fact; "TPP1 quenches motion at the DNA exit pore" is not supported. |
| 4 | **The MD cannot support anything about mechanism.** 2.36 ps of Langevin dynamics on a **97-residue four-fragment model with no hTR and no DNA**, and the file labelled `junction_bound_k757a_fixed.pdb` is **byte-identical to the WT file** (MD5 `1b284b2f…`), i.e. no K757A system was ever built. | §4.4. Every MD-derived claim collapses. |
| 5 | **TELIN-1 does not exist as a molecule.** 36 heavy atoms in two disconnected fragments (one is an unbonded pentanoic-acid tail), zero literature, zero measured affinity, and the "docked complex" interpenetrates TERT down to 0.17 Å with 160 clashes < 2.8 Å. | §6. It is not a lead compound; it is a rendering artifact. |
| 6 | **The headline "breakthrough" is not novel and its predicted phenotype is contradicted by published biochemistry.** Sekne et al. (Science 2022) already made and assayed TERT **K757A/F759A** and reported reduced activity **and** reduced RAP; Tesmer et al. (Life Sci. Alliance 2023) made **K757E/F759E** and found them defective in activity and processivity while still localising to telomeres. The repository's criterion "k_cat(K757A) ≥ 85 % of WT" is therefore almost certainly false. | §7. The claim must be rewritten, not repaired. |
| 7 | **The "validated" SPR/RAP data are synthetic.** `scripts/generate_independent_benchmarks.py` draws them with `np.random.seed(42)` from hand-chosen k_on/k_off, then the fitting pipeline recovers them and stamps "CRITERION 2A PASSED". | §4.5. This is circular validation and must never be presented as evidence. |

**The single defensible breakthrough claim is not the one currently in the repository.** It is
stated, with its falsification criteria, in
[`results/verification/corrected_breakthrough_claim.md`](corrected_breakthrough_claim.md).

---

## 1. Method

* **Independence.** `scripts/verify_all_claims.py` imports nothing from `telotrade`. CIF parsing
  uses gemmi; superpositions use a hand-written Kabsch SVD; SASA uses a hand-written
  Shrake–Rupley (Fibonacci sphere, 960 points, 1.4 Å probe); contacts use an independent
  `scipy.cKDTree` query. Agreement with the repository is therefore genuine cross-validation,
  not a re-run of the same code.
* **Verdict definitions.**
  * `VERIFIED` — reproduced within stated tolerance.
  * `MINOR_DRIFT` / `PARTIALLY_VERIFIED` — right order, wrong detail (usually a convention or a
    reporting-level mismatch).
  * `CONTRADICTED` — the number, the atom pair, or the stated basis is wrong.
  * `UNTESTED` — cannot be reproduced from what ships (missing dependency, missing computation).
* **Tolerances.** Distances ±0.02 Å; RMSDs ±0.005 Å; B-factors ±3 Å²; BSA ±2 % (sphere point-set
  differences between SR implementations); integer contact counts exact.
* **Conventions are named, not assumed.** Whole-residue nearest-atom distances are distinguished
  from side-chain-only distances throughout (§2 rows D757-*, D756-RNA). This matters: the
  project's "Lys757 is 5.75 Å from the primer" is true for the residue's **backbone carbonyl O**
  and false for the side-chain NZ that a K757A mutation deletes (6.31 Å).

---

## 2. Master table (56 claims)

All numbers below are from `results/verification/claim_recheck.csv`. "Repo" = the value as printed
in the repository; "recomputed" = independent re-derivation.

### 2.1 Data provenance and structure metadata

| Claim | Repo | Recomputed | Verdict |
|---|---|---|---|
| SHA-256 of `7TRE.cif` (FINAL_RESEARCH_REPORT) | `2c681322533afea5a1f1…` | `2c681322533afea566f7…` | **CONTRADICTED** |
| SHA-256 of `7TRD.cif` | `faf6d1aa6401d3ddcfd2…` | `faf6d1aa6401d3dd0f66…` | **CONTRADICTED** |
| SHA-256 of `7TRF.cif` | `e62419d5c2aac19d28c7…` | `e62419d5c2aac19d9f6c…` | **CONTRADICTED** |
| SHA-256 of `2I46.cif` | `6f70cb1e6a92881bbfaf…` | `4f4402ce3e9aadfa5145…` | **CONTRADICTED** |
| 2I46 resolution | 1.85 Å (report.md) / 1.50 Å (final report) / 2.70 Å (manifest) | **2.70 Å** (`_reflns.d_resolution_high`) | **CONTRADICTED** (two of three values) |
| 2I46 chain A range | 90–241 | **90–301** (203 res; chain B 90–287) | **CONTRADICTED** |
| 2I46 chain B DE-loop 168–171 unmodelled | missing | A has 168–171; B has none | VERIFIED |
| hTERT 753–761 = GHVRKAFKS | GHVRKAFKS | GHVRKAFKS | VERIFIED |

### 2.2 Interface geometry (7TRE, TERT auth chain A vs TPP1 auth chain C)

| Claim | Repo | Recomputed | Verdict |
|---|---|---|---|
| Interface BSA (Bondi, /2 convention) | 856.82 Å² | 857.46 Å² (Δ +0.07 %) | VERIFIED |
| ΔSASA total | 1713.65 Å² | 1714.91 Å² | VERIFIED |
| Interface BSA (ProtOr) | 923.16 Å² | not reproducible offline | UNTESTED |
| BSA sensitivity (probe + radius set) | ~250 Å² (~30 %) | 410 Å² (39.5 %) over probe 1.2–1.6 Å × Bondi/Chothia | PARTIALLY_VERIFIED |
| Heavy-atom contacts ≤ 3.5 Å | 9 (5 residue pairs) | 9 (5) | VERIFIED |
| ≤ 4.0 Å | 46 (21) | 46 (21) | VERIFIED |
| ≤ 4.5 Å | 123 (32) | 123 (32) | VERIFIED |
| Tier partition core/intermediate/peripheral | 5 / 16 / 11 | 5 / 16 / 11 | VERIFIED |
| Arg774:NH1–Glu168:OE1 | 2.32 Å | 2.320 Å | VERIFIED |
| Lys78:NZ–Glu215:OE2 | **2.89 Å** | 3.584 Å (closest pair in the residue pair = NZ–OE2 3.584 Å) | **CONTRADICTED** |
| Leu769:CD2–Leu183:CD1 | **3.52 Å** | 5.532 Å (correct pair is CD2–CD2 = 3.665 Å) | **CONTRADICTED** |
| Pro771:CG–Leu212:CD1 | **3.86 Å** | 5.630 Å (correct pair is CG–CD2 = 3.524 Å) | **CONTRADICTED** |

> The three contradictions are **atom-pair labelling errors**, not measurement errors. In each case
> the residue pair is a genuine contact and the ledger's alternative number is the right one
> (`evidence_ledger.csv` CLM-006/008/009 are correct; `report.md` §6.1 is not). This is exactly
> the error class that a reader cannot detect without re-opening the coordinates — which is why
> the independent recomputation matters.

### 2.3 Cross-state comparison

| Claim | Repo | Recomputed | Verdict |
|---|---|---|---|
| Core Cα RMSD 7TRE vs 7TRD (326–1132, 770 Cα) | 0.338 Å | 0.3377 Å | VERIFIED |
| Global Cα RMSD (934 Cα) | 0.356 Å | 0.3557 Å | VERIFIED |
| TEN excursion after core fit (160 Cα) | 0.434 Å | 0.4341 Å | VERIFIED |
| Core Cα RMSD 7TRE vs 7TRF | 0.338 Å | 0.3377 Å | VERIFIED (but see next row) |
| 7TRD and 7TRF are independent unbound models | treated as such | superpose at **0.0005 Å**; identical B-factors | **CONTRADICTED** |
| TPP1 OB-fold RMSD 7TRE vs 2I46 (152 Cα, sequence-verified) | 1.966 Å | 1.966 Å | VERIFIED |
| Leu183 displacement | 0.512 Å | 0.512 Å | VERIFIED |
| Lys170 displacement | 4.293 Å | 4.293 Å | VERIFIED |
| Glu169 / Glu171 displacement | 3.809 / 3.320 Å | 3.809 / 3.320 Å | VERIFIED |
| Glu168 displacement | **2.748 Å** (report.md §8) | 1.253 Å (ledger CLM-018 is right) | **CONTRADICTED** |
| Leu212 displacement | 1.040 Å | 1.040 Å | VERIFIED |
| DE-loop B-factors apo vs holo | >100 / ~53 Å² | 88.2 / 53.6 Å² | VERIFIED (weakly) |

> **7TRD ≡ 7TRF.** All 7,488 TERT heavy atoms superpose at median 0.0006 Å (max 2.32 Å at two
> side chains, residues 125/127) and the mean B-factor difference is 0.0000 Å². The two "unbound"
> entries are one model reported twice, so "7TRD and 7TRF agree with 7TRE at 0.338 Å" is a single
> measurement, not two. This removes the only apparent replication of the unbound state.

### 2.4 The proposed switch and hinge

| Claim | Repo | Recomputed | Verdict |
|---|---|---|---|
| Tyr176:OH–Gln794:OE1 in bound 7TRE | 2.338 Å | 2.338 Å | VERIFIED |
| Same distance in apo 7TRD / 7TRF | 4.372 Å | 4.372 / 4.372 Å | VERIFIED |
| Lys757 B-factor 7TRE / 7TRD | 25.8 / 121.5 Å² | 25.78 / 121.53 Å² | VERIFIED |
| "78.8 % quench" | 78.8 % | 78.79 % | VERIFIED numerically |
| **The quench is site-specific** | implied | Lys757 window ratio 0.24 vs TRAP region 0.26 vs global 0.19 | **CONTRADICTED** |
| Lys757→TPP1 | 25.43 Å | 25.43 Å (whole residue) / 27.94 Å (NZ) | VERIFIED (convention) |
| Lys757→DT12:OP1 | 5.75 Å | 5.75 Å (backbone O) / 6.31 Å (NZ) | VERIFIED (convention) |
| Lys757→hTR A55:N6 | 8.33 Å | 8.33 Å (whole residue) / 11.49 Å (NZ) | VERIFIED (convention) |
| Arg756→hTR A55:C2 | 3.12 Å | 3.12 Å **via NH2** (NH1–C2 = 4.13 Å) | VERIFIED (atom label wrong) |

### 2.5 MD and ligand

| Claim | Repo | Recomputed | Verdict |
|---|---|---|---|
| Cleft distance 3.820 ± 0.633 vs 5.620 ± 0.196 Å | as stated | means reproduce exactly; **± is replica-level (n=3)** while the mean is frame-level (n=360) | PARTIALLY_VERIFIED |
| Sampling = 3.0 ps total / 2.4 ps production | as stated | 3.00 / 2.36 ps | VERIFIED (numerically) |
| Replica-level p = 0.00928, t = −4.705 | as stated | p = 0.0093, t = −4.705 (n = 3 vs 3) | VERIFIED (statistic) |
| **MD system = TPP1-bound holoenzyme** | implied | 97 residues, 4 chains, **no hTR, no DNA**, no active site | **CONTRADICTED** |
| **A K757A mutant system exists** | implied | `junction_bound_k757a_fixed.pdb` is byte-identical to WT (Lys757 intact) | **CONTRADICTED** |
| TELIN-1 = 36 heavy atoms, C26N2O6SF | as stated | 36 atoms, C26 N2 O6 S1 F1 | VERIFIED |
| TELIN-1 = two disconnected fragments | as stated | 2 fragments | VERIFIED |
| Dock: 160 / 21 / 8 clashes, min 0.17 Å | as stated | 160 / 21 / 8, min 0.17 Å (TERT Ala130 backbone N) | VERIFIED |
| TPP1 absent from docked complex | absent | absent | VERIFIED |
| Shape complementarity quantified | claimed (Aims) | no Sc value exists anywhere in `results/` | UNTESTED |

---

## 3. What the structural work genuinely establishes

These statements survive independent re-derivation and should be kept, in this form:

1. The 7TRE TERT–TPP1 interface buries **857 Å² per partner** (ΔSASA 1 715 Å²; Shrake–Rupley,
   Bondi radii, 1.4 Å probe, protein heavy atoms only, `/2` convention). Report the convention
   with the number: the same interface is 683 Å² at a 1.2 Å probe and 1 036 Å² at 1.6 Å.
2. The interface core (≤3.5 Å) is **5 residue pairs / 9 heavy-atom contacts**: TERT Ala45–TPP1
   Glu171, Ala46–Glu169, Tyr772–Leu211, **Arg774–Glu168**, Asn799–Arg92. The 4.0 Å shell adds
   16 pairs (including Lys78–Glu215 = 3.584 Å); the 4.5 Å shell adds 11 threshold-sensitive pairs.
3. **Arg774–Glu168 = 2.320 Å** is the tightest contact in the interface and is a positive control
   (Liu 2022; Sekne 2022), not a discovery. At 3.5 Å map resolution this distance is shorter than
   a physical H-bond and should be reported as refinement-compressed.
4. The TPP1 OB-fold is **pre-organised at Leu183 (0.512 Å)** and **disordered at the DE-loop
   (Lys170 4.293 Å, Glu169 3.809 Å, Glu171 3.320 Å, Glu168 1.253 Å)** between the bound (7TRE)
   and unbound (2I46) states — a real, quantified disorder-to-order difference of a loop, which
   is a legitimate (if conventional) observation.
5. The TERT catalytic core is **coordinate-invariant** between 7TRE and 7TRD (0.3377 Å over 770
   Cα, i.e. below the 0.35–0.5 Å coordinate uncertainty of a 3.3–3.5 Å cryo-EM map).

---

## 4. Where the reasoning breaks

### 4.1 Atom-pair inflation

Three of the four "key hydrophobic/electrostatic contacts" in `report.md` Table 6.1 quote a
distance that does not exist between the named atoms (K78:NZ–E215:OE2 2.89 vs 3.58; L769:CD2–L183:CD1
3.52 vs 5.53; P771:CG–L212:CD1 3.86 vs 5.63). In each case a *different* atom pair at the *right*
distance exists in the same residue pair, and the evidence ledger already contains the correct
value. The effect is to make a 3.5–4 Å van der Waals contact look like a 2.9–3.5 Å near-H-bond.

### 4.2 B-factors cannot carry the mechanistic weight placed on them

* Per-residue ADP collapse: 100 % of residues in 7TRE, 7TRD and 7TRC carry a single B-value
  copied across all atoms (group ADP refinement). B-factors are therefore per-residue, not
  per-atom, and cannot be read as local flexibility at the atomic level.
* Cross-entry scale: mean B = **16.5 Å² (7TRE)** vs **85.9 Å² (7TRD)**. 7TRE is *lower*-resolution
  (3.5 Å) yet has 5× lower ADPs, so the two entries are not on a common scale. Ratios, not
  differences, are the only meaningful comparison — and by ratio the Lys757 window (0.24) is
  indistinguishable from its own domain (TRAP 0.26) and not far from the global value (0.19).
* Consequence: the "67–82 % B-factor collapse at the DNA exit pore" is a restatement of the global
  ADP-scale difference between the two depositions, not evidence of TPP1-induced quenching at
  Lys757.

### 4.3 Mixed levels of replication

The headline cleft-distance table reports a **frame-pooled mean** (± SD computed at replica level,
n = 3) with a **replica-level p-value**, and then describes the result as "p < 0.01" in the README.
Frame-pooled SDs are 0.978 / 0.758 Å, not 0.633 / 0.196 Å. Two files in the repository
(`data/md_replicate_summary.txt` and `results/experimental_validation/RECALCULATED_REPLICA_STATISTICS.csv`)
give different summary values for the same six CSVs (3.496/5.296 vs the recomputable 3.820/5.620).

### 4.4 The MD is not an MD of anything relevant

| Check | Result |
|---|---|
| System composition (`junction_bound_fixed.pdb`) | 97 residues, 4 chains: TERT 165–178, TERT 750–795, TPP1 165–184, TPP1 210–225 |
| Nucleic acid present | **No** (no hTR, no DNA) |
| Active site / catalytic triad present | No |
| Simulated time | 2.36 ps production (60 frames × 20 steps × 2 fs) |
| "K757A" system | byte-identical to WT (MD5 `1b284b2f59095bb44edbe50fe86d1b04`) |
| Force field / method | AMBER ff14SB + OBC2, Langevin, 300 K, no cut-off — disclosed correctly |
| Script portability | hard-coded `C:/Users/hp/...` path in `scripts/run_md_trajectories.py` |

A 2.4 ps trajectory samples ~10⁹-fold less time than one repeat-addition cycle (0.5–2 s).
No conclusion about processivity, allostery or TPP1 stimulation can be drawn from it, and the
K757A predictions were never simulated at all.

### 4.5 Circular "validation"

`scripts/generate_independent_benchmarks.py` fabricates the SPR sensorgrams (k_on = 4.3×10⁴ M⁻¹s⁻¹,
k_off = 6.2×10⁻³ s⁻¹ → Kd = 144 nM, `np.random.seed(42)`) and the RAP densitometries, and
`results/experimental_validation/spr_benchmark_results.csv` then reports "CRITERION 2A PASSED
(Nanomolar Binding Preserved)" — on data that were generated from the answer. The "WT reference
Kd = 150 ± 25 nM" that anchors the pre-registered criteria has **no cited experimental source**;
it is the mean of the three synthetic WT replicates. Pipeline validation is legitimate; presenting
it as biological evidence is not.

---

## 5. Data-provenance and repository defects (all independently confirmed)

| Defect | Detail | Status |
|---|---|---|
| Report hashes wrong | all four SHA-256 values in `FINAL_RESEARCH_REPORT.md` §3.1 differ from the shipped files; `data_manifest.csv` values match the files | **ERRATA added** |
| Missing manifest copy | `data/manifest/data_manifest.csv` was absent, breaking 7 tests and 9 fixtures on a clean clone | **FIXED** (file added; suite now 338 passed / 1 failed / 3 errors in a minimal env) |
| Case-sensitive filename | `tests/test_adversarial_b2_challenger2.py` opens `data/raw/7tre.cif`; the repository ships `7TRE.cif` — passes only on Windows/macOS | **FIXED** |
| Label-chain path broken | `compute_bsa(..., use_auth=False)` returns 331.65 Å² instead of 856.82 Å²: the freeSASA branch honours `use_auth`, but the Biopython fallback indexes by *author* chain id | open (test `test_compute_bsa_label_chains` fails without freeSASA) |
| Non-deterministic table | `results/tables/tel_patch_contacts.csv` row order changes between runs (unstable sort on tied distances) — regenerating the file dirties the git tree | open |
| "Optional" dependency is required | `freesasa` is listed in `requirements.txt` but fails to build in this environment; without it the BSA tests silently exercise a different code path | note |

**Not demonstrable either way (blocked offline):** whether `data/raw/*.cif` are pristine RCSB
downloads. The per-residue ADP grouping is consistent with group-B refinement and is not itself
evidence of tampering, but the report's mismatched hashes mean authenticity must be re-established
by a fresh download + hash check before publication.

---

## 6. TELIN-1: a molecule that does not exist

| Check | Result |
|---|---|
| Composition | C₂₆N₂O₆SF, 36 heavy atoms — matches the label |
| Connectivity | **2 disconnected covalent fragments** (a 29-atom polycyclic core + an unbonded 7-atom acid tail); longest nominal "bond" 1.85 Å |
| Provenance | 0 hits in PubChem/ChEMBL/PubMed/patent space (as reported); no synthesis route, no certificate |
| Docking | 160 clashes < 2.8 Å, 21 < 1.5 Å, 8 < 1.0 Å, minimum 0.17 Å into the TERT Ala130 backbone N; **no TPP1 in the model at all** |
| Affinity / uptake / selectivity | **none measured, none computed with a scoring function that is described** |

**Verdict:** TELIN-1 is not a "computational lead" in any usable sense. It is an unparameterised
sketch whose docked pose is physically impossible (a 0.17 Å interatomic distance is a coordinate
collision, not a contact). Any downstream document that repeats "TELIN-1" as a candidate should be
quarantined, which `report.md` §10 already does — keep that, and delete the ligand files from the
"models" section of the README.

---

## 7. The prior art that decides the question

The repository's central proposal is: *"TERT Lys757 is a hinge; K757A selectively uncouples
TPP1-stimulated RAP while preserving recruitment and basal catalysis."* Three published results
bear directly on it:

1. **Sekne et al., Science 375, 1173–1176 (2022)** — the paper that deposited the TPP1–POT1
   structures — reconstituted telomerase with TERT **Y176A/Q177A, K757A/F759A and Q794A** and
   reported that all **"significantly reduced telomerase activity and RAP"** (Fig. 4E–G). The
   residues the repository calls a "novel allosteric lock" (Tyr176, Gln794) and its "premier
   separation-of-function candidate" (Lys757) are the paper's **DNA anchor site**.
2. **Tesmer et al., Life Science Alliance 6, e202201727 (2023)** — "Interaction hub critical for
   telomerase recruitment and primer-template handling for catalysis". TERT **K757E and F759E**
   are defective in activity *and* processivity, yet **localise to telomeres like WT**
   (recruitment intact). That is precisely the separation-of-function phenotype the repository
   predicts for K757A — already published, with the opposite sign for basal activity.
3. **Kocak et al., Genes Dev. 28, 2090–2102 (2014)** and **Bisht et al., PNAS 113, 13021–13026
   (2016)** — TPP1 K170Δ (TEL patch) causes Hoyeraal–Hreidarsson syndrome and reduces both
   recruitment and processivity; the TEL patch is a *separation-of-function* module on the TPP1
   side.

Additional citation errors to correct: Latrick & Cech 2010 is **EMBO J. 29, 924–933** (not NSMB);
Wang et al. 2007 (2I46) is **Nature 445, 506–510** (not NSMB 14, 1011); Liu et al. 2022 is
**604, 578–583** (not 578–585); Sekne 2022 is **1173–1176** (not 1173–1179); the 2I46 construct is
TPP1 90–301 at **2.70 Å**.

---

## 8. The corrected claim

Replacing "TERT Lys757 is an allosteric hinge that can be mutated to uncouple processivity from
recruitment" with a claim that is (i) not already published, (ii) consistent with the published
anchor-site data, and (iii) falsifiable with the assays this repository already specifies:

> **The TERT IFD-TRAP DNA-anchor surface (Lys757/Phe759/Gln794) and the TPP1 TEL-patch docking
> surface (Arg774, Lys78) are functionally coupled: the fold-stimulation of repeat-addition
> processivity by POT1–TPP1 (χ_RAP) is not a simple multiplicative consequence of recruitment,
> but depends on an intact anchor surface, and the coupling can be measured as a non-additivity
> (epistasis) term between interface and anchor mutations.**

The full statement, with decision criteria, controls, and the experiments that would kill it, is in
[`corrected_breakthrough_claim.md`](corrected_breakthrough_claim.md).

---

## 9. Recommended repo actions

**Retract / relabel**
1. Delete TELIN-1 from the "Structural Models" section; keep only the quarantine note.
2. Remove the MD statistics table from the README headline position; keep it as a methods
   benchmark, relabelled "2.4 ps relaxation of a 97-residue fragment".
3. Remove the words "breakthrough", "discovery" and "validated" from anything resting on
   `data/benchmarks/*` (synthetic) or on `data/trajectories/*` (fragment MD).
4. Replace the K757A "separation-of-function" claim with §8 and cite Sekne 2022 / Tesmer 2023
   as the prior tests of K757.

**Fix**
5. `data/manifest/data_manifest.csv` (done), the lowercase CIF path in the challenger test (done),
   the SHA-256 block in `FINAL_RESEARCH_REPORT.md` (errata added), the 2I46 resolution and range,
   the three atom-pair distances in `report.md` Table 6.1, the Glu168 displacement, the citation
   list, and the BSA label-chain code path.

**Add**
6. `results/verification/claim_recheck.csv` as the single source of truth for headline numbers,
   regenerated by `scripts/verify_all_claims.py` in CI.
7. A `PROVENANCE.md` recording the download date, URL and hash for each CIF, produced by a script
   that re-downloads and re-hashes (blocked in this sandbox; run where RCSB is reachable).
