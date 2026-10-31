# Telomerase Allosteric Switch (`telomerase-allosteric-switch`)

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Test Suite](https://img.shields.io/badge/pytest-376%20passed-success.svg)](tests/)
[![Code Hygiene](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Simulation Engine](https://img.shields.io/badge/OpenMM-8.6.1-orange.svg)](https://openmm.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Mechanism%20audit%20complete%20--%20two%20experiments%20proposed-blue.svg)](results/mechanism/MECHANISM_ASSESSMENT.md)

**A structural and experimental-evidence audit of how TPP1 stimulates human telomerase repeat-addition processivity — including the record of which claims in this repository failed and why.**

> ## ⚠️ Status (2026-10-31): mechanism audit complete, no mechanism proven
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
> The **pytest badge (376 passed)** is verified on a clean clone with `pip install -e .`; the
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
├── results/mechanism/                   # CURRENT deliverable: the mechanism assessment (WP1-WP5)
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
│   └── generate_code_freeze_manifest.py # Computes SHA-256 manifest and creates freeze document
│
├── tests/                               # Comprehensive test suite (376 passing tests)
│   ├── e2e/                             # End-to-end pipeline execution tests
│   ├── test_adversarial_b2_challenger.py# Adversarial superimposition and SVD stress tests
│   ├── test_adversarial_m1.py ... m5.py # Data integrity, corrupted CIF, and figure generation tests
│   ├── test_experimental_scripts.py     # 10 unit tests for SPR, RAP, and MD analysis pipelines
│   ├── test_ligand_telin1.py            # 2 unit tests verifying Telin-1 clashes and composition
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

### Reproduce Key Analyses

```bash
# 1. Run full project test suite (357 passing tests)
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
- **Test Suite**: `357 passed` across 22 test modules.
- **Linting & Formatting**: Zero errors via `ruff check src tests scripts`.
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
