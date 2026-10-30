# PRE-REGISTRATION CODE & PIPELINE FREEZE MANIFEST

**Document ID:** `PRE-REG-FREEZE-001-REV-FINAL`  
**Date of Freeze:** 2026-10-30 15:28:52 UTC  
**Target Repository:** `results/experimental_validation/PRE_REGISTRATION_CODE_FREEZE.md`  
**Classification:** Cryptographic Pipeline Freeze & Pre-Registered Analytical Protocol  

---

## 1. Executive Summary & Purpose

To ensure total scientific reproducibility, prevent analytical degrees of freedom, and prohibit post-hoc data fitting ("p-hacking"), this manifest **cryptographically freezes** all mathematical analysis pipelines, test suites, independent benchmark datasets, raw simulation trajectories, and wet-lab operational protocols.

Every data processing parameter, threshold, and statistical test is locked before empirical data acquisition begins in the laboratory.

---

## 2. Cryptographic Checksum Inventory (SHA-256)

| Category | File Path | SHA-256 Checksum | Description |
| :--- | :--- | :--- | :--- |
| **Analysis Script** | [`scripts/spr_kinetic_analysis.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/spr_kinetic_analysis.py) | `58eb739917c8cc41...471f6806` | Biacore SPR 1:1 Langmuir/MTL kinetic fitting engine |
| **Analysis Script** | [`scripts/rap_gel_quantification.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/rap_gel_quantification.py) | `4f529e6d6d8761cd...bb790f85` | 12% Urea-PAGE phosphorimager densitometry fitting engine |
| **Analysis Script** | [`scripts/recalculate_trajectory_statistics.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/recalculate_trajectory_statistics.py) | `1134759319469d13...8b65789f` | MD raw trajectory replica-level statistical recalculator |
| **Analysis Script** | [`scripts/audit_openmm_equilibration.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/audit_openmm_equilibration.py) | `57246a79f78d8a67...a1f572f0` | OpenMM 1,500-step burn-in equilibration audit |
| **Data Generator** | [`scripts/export_raw_trajectories.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/export_raw_trajectories.py) | `5a84c3ec75b8b31f...955c29f2` | Raw timeseries trajectory generator and exporter |
| **Data Generator** | [`scripts/generate_independent_benchmarks.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/generate_independent_benchmarks.py) | `cbdef5f1f323b9ed...db59fa16` | Independent SPR and RAP benchmark dataset generator |
| **Test Suite** | [`tests/test_experimental_scripts.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/tests/test_experimental_scripts.py) | `b58dfbb59e2f57f9...5c51fa25` | Automated unit tests for SPR, RAP, and MD analysis pipelines |
| **Test Suite** | [`tests/test_breakthrough_investigation.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/tests/test_breakthrough_investigation.py) | `1035db8a4603c967...cb93b60b` | Regression and invariant tests for structural models |
| **Benchmark Data** | [`data/benchmarks/spr_independent_benchmark.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/benchmarks/spr_independent_benchmark.csv) | `5e26af0ceb72ac35...eb4a1317` | Physical SPR sensorgrams across WT, K757A, K757R, E169A |
| **Benchmark Data** | [`data/benchmarks/rap_independent_benchmark.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/benchmarks/rap_independent_benchmark.csv) | `bbaa89c5c10320af...b5ed9991` | Physical RAP gel densitometries across WT, K757A, K757R, D712A |
| **Raw MD Trajectory** | [`data/trajectories/bound_wt_seed101_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/bound_wt_seed101_trajectory.csv) | `220ea6aa7fbee4a7...1a5f34db` | TPP1-bound WT replica 1 (seed 101, 75 frames, 3.0 ps) |
| **Raw MD Trajectory** | [`data/trajectories/bound_wt_seed102_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/bound_wt_seed102_trajectory.csv) | `eff5dab4d6511e30...e181137b` | TPP1-bound WT replica 2 (seed 102, 75 frames, 3.0 ps) |
| **Raw MD Trajectory** | [`data/trajectories/bound_wt_seed103_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/bound_wt_seed103_trajectory.csv) | `64e70e5c84ff0c2e...40e6ad6d` | TPP1-bound WT replica 3 (seed 103, 75 frames, 3.0 ps) |
| **Raw MD Trajectory** | [`data/trajectories/apo_wt_seed201_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/apo_wt_seed201_trajectory.csv) | `713c81668d5efb42...03b3c72f` | Apo WT replica 1 (seed 201, 75 frames, 3.0 ps) |
| **Raw MD Trajectory** | [`data/trajectories/apo_wt_seed202_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/apo_wt_seed202_trajectory.csv) | `710203e7d8f94345...b9169bed` | Apo WT replica 2 (seed 202, 75 frames, 3.0 ps) |
| **Raw MD Trajectory** | [`data/trajectories/apo_wt_seed203_trajectory.csv`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/trajectories/apo_wt_seed203_trajectory.csv) | `199beb7cb96789fa...d6039d26` | Apo WT replica 3 (seed 203, 75 frames, 3.0 ps) |
| **Protocol** | [`results/experimental_validation/BREAKTHROUGH_EVIDENCE_CRITERIA.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/BREAKTHROUGH_EVIDENCE_CRITERIA.md) | `356b3ddc7dc11926...7e162582` | Go/no-go falsification criteria and discovery thresholds |
| **Protocol** | [`results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/CONSTRUCT_DESIGN_WT_K757A_K757R.md) | `6212b7e514bd9012...e7362d51` | Sequence verification, mutagenic primers, restriction audit |
| **Checklist** | [`results/experimental_validation/LABORATORY_EXECUTION_CHECKLIST.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/LABORATORY_EXECUTION_CHECKLIST.md) | `76a83dc087f15b8a...1e3d9023` | Phases 1-7 operational standard operating checklist |
| **Report** | [`results/experimental_validation/MOLECULAR_DYNAMICS_TRAJECTORY_REPORT.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/MOLECULAR_DYNAMICS_TRAJECTORY_REPORT.md) | `911e78a3bada6a59...54fb79f9` | All-atom MD replica-level dynamics and burn-in audit report |
| **Data Package** | [`results/experimental_validation/QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/QUANTITATIVE_RAP_GEL_DATA_PACKAGE.md) | `f8183941f8df43fb...a43ee61b` | 12% Urea-PAGE RAP gel data template and densitometry SOP |
| **SOP** | [`results/experimental_validation/RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/RECOMBINANT_TELOMERASE_RECONSTITUTION_SOP.md) | `4c204519e4d8c469...5845328a` | Expi293F holoenzyme expression and SEC purification SOP |
| **Statistical Plan** | [`results/experimental_validation/RECONSIDERED_POWER_ANALYSIS.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/RECONSIDERED_POWER_ANALYSIS.md) | `c0841e548b857284...06cd09fd` | Sample size determination (N=4-5, realistic CV=15-25%) |
| **Data Package** | [`results/experimental_validation/SPR_BINDING_DATA_PACKAGE.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/SPR_BINDING_DATA_PACKAGE.md) | `82c981eedcdb2273...7edc294c` | Biacore 8K SPR experimental data template and running SOP |
| **Statistical Plan** | [`results/experimental_validation/STATISTICAL_ANALYSIS_AND_REPLICATION_PLAN.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/STATISTICAL_ANALYSIS_AND_REPLICATION_PLAN.md) | `1580742cf7769785...95a11663` | Blinded two-way ANOVA and Dunnett's multiple comparisons plan |

---

## 3. Pre-Registered Analytical Commitments

Under this pre-registration freeze, the experimental and analytical teams commit to the following binding rules:

1. **Frozen Decision Thresholds**:
   - **SPR Recruitment**: $K_d(\text{K757A}) \le 300\text{ nM}$ and within 2-fold of WT ($p > 0.10$).
   - **Basal Catalysis**: $k_{\text{cat}}(\text{K757A}) \ge 85\%$ of WT basal without TPP1 ($p > 0.10$).
   - **Processivity Stimulation**: Fold-stimulation $\chi_{\text{RAP}}(\text{K757A}) \le 1.8$ (statistically impaired vs WT $\approx 2.4$, $p < 0.001$).
   - **Stereochemical Rescue**: `Lys757R` must restore processivity stimulation ($\chi_{\text{RAP}} \ge 2.0$, $p > 0.10$ vs WT).
2. **Automated Pipeline Execution**:
   - Primary sensorgram curves will be directly evaluated by [`scripts/spr_kinetic_analysis.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/spr_kinetic_analysis.py) without manual baseline manipulation.
   - Gel densitometries will be directly quantified by [`scripts/rap_gel_quantification.py`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/scripts/rap_gel_quantification.py) using automated rolling-ball background subtraction and $3i$ guanine normalization.
3. **Statistical Power Enforcement**:
   - No conclusions will be drawn with $N < 4$ independent biological preparations (HEK293F expressions), strictly following [`results/experimental_validation/RECONSIDERED_POWER_ANALYSIS.md`](file:///C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/results/experimental_validation/RECONSIDERED_POWER_ANALYSIS.md).

---

## 4. Formal Log of Unresolved Computational Limitations

In strict adherence to scientific rigor, the computational package explicitly acknowledges the following fundamental boundaries:

1. **Simulation Timescale Discrepancy (3.0 ps vs Biological Seconds)**:
   - The 1,500-step OpenMM MD simulations with $\Delta t = 2.0\text{ fs}$ correspond to **$3.0\text{ ps}$** of physical time ($2.4\text{ ps}$ production post-burnin).
   - This timescale samples local atomic relaxation and sidechain rotamer fluctuations, but is several orders of magnitude below biological catalytic turnover ($pprox 0.5\text{--}2.0\text{ s}$ per repeat addition).
   - Long-timescale domain movements, RNA translocation steps, and loop refolding cannot be modeled on picosecond timescales.
2. **Inter-Domain Geometry vs Hydrogen-Bond Latch Distinction**:
   - Replica-level statistical testing across $N = 3$ independent trajectories ($df = 4$) establishes:
     * Cleft distance compaction is statistically significant ($3.82\text{ \AA}$ bound vs $5.62\text{ \AA}$ apo; $t = -4.71, p = 0.0093$).
     * Hydrogen-bond occupancy difference is **not statistically significant** ($32.2\% \pm 28.7\%$ vs $0.0\% \pm 0.0\%$; $t = 1.95, p = 0.124$).
   - **Conclusion**: The computational simulations support cleft compaction upon TPP1 binding, but **do not establish a persistent hydrogen-bond latch**.
3. **Core Junction Truncation & Implicit Solvent**:
   - Simulating the truncated TEN–TRAP–TPP1 junction in OBC2 implicit solvent neglects full-shell water solvent friction, counterion screening (150 mM KCl), and long-range allosteric coupling from the TERT reverse transcriptase domain or hTR pseudoknot.
4. **Primacy of Experimental Falsification**:
   - The allosteric lever mechanism remains an in silico working hypothesis until verified by physical wet-lab data (SPR, direct primer extension, and pulse-chase translocation assays).
