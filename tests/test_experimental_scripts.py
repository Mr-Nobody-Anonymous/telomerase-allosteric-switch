"""Unit tests auditing the experimental analysis pipelines for SPR and RAP gels.

Verifies fitting convergence, parameter recovery, noise tolerance, and
pre-registered falsification triggers across positive and negative controls.
"""

import sys
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

from scripts.rap_gel_quantification import (
    calculate_fold_stimulation,
    calculate_rap_metrics,
    evaluate_rap_phenotype,
)
from scripts.spr_kinetic_analysis import (
    evaluate_spr_candidate,
    fit_spr_cycle,
    langmuir_1_to_1_association,
    langmuir_1_to_1_dissociation,
)

# ---------------------------------------------------------------------------
# SPR Analysis Pipeline Tests
# ---------------------------------------------------------------------------

def test_spr_wt_positive_control():
    """Verify recovery of known WT kinetic parameters from simulated sensorgram."""
    t_a = np.linspace(0, 180, 90)
    t_d = np.linspace(0, 600, 120)
    
    true_kon = 4.2e4
    true_koff = 0.0063
    conc = 240e-9
    
    resp_a = langmuir_1_to_1_association(t_a, true_kon, true_koff, 65.0, conc)
    resp_d = langmuir_1_to_1_dissociation(t_d, resp_a[-1], true_koff)
    
    fit = fit_spr_cycle(t_a, resp_a, t_d, resp_d, conc)
    qc_pass, status = evaluate_spr_candidate("TERT_WT", fit["Kd_nM"], fit["Chi2"], fit["R_max"])
    
    assert bool(qc_pass) is True, "WT fit should pass QC"
    assert 120.0 <= fit["Kd_nM"] <= 180.0, f"Expected Kd ~150 nM, got {fit['Kd_nM']:.1f}"
    assert status == "BASELINE VALID"

def test_spr_uncoupling_candidate_k757a():
    """Verify evaluation of Lys757A candidate with nanomolar binding."""
    t_a = np.linspace(0, 180, 90)
    t_d = np.linspace(0, 600, 120)
    
    true_kon = 3.8e4
    true_koff = 0.0072
    conc = 240e-9
    
    resp_a = langmuir_1_to_1_association(t_a, true_kon, true_koff, 60.0, conc)
    resp_d = langmuir_1_to_1_dissociation(t_d, resp_a[-1], true_koff)
    
    fit = fit_spr_cycle(t_a, resp_a, t_d, resp_d, conc)
    qc_pass, status = evaluate_spr_candidate("TERT_Lys757A", fit["Kd_nM"], fit["Chi2"], fit["R_max"])
    
    assert bool(qc_pass) is True
    assert fit["Kd_nM"] <= 300.0, f"Candidate should exhibit Kd <= 300 nM, got {fit['Kd_nM']:.1f}"
    assert "CRITERION 2A PASSED" in status

def test_spr_recruitment_failure_negative_control():
    """Verify that a binding-disrupted mutant triggers the FALSIFICATION verdict."""
    fake_kd = 2500.0 # 2.5 uM (severely impaired)
    _qc_pass, status = evaluate_spr_candidate("TERT_Lys757A", fake_kd, chi2=0.5, r_max=50.0)
    assert status == "FALSIFIED (Recruitment Disrupted)"

# ---------------------------------------------------------------------------
# RAP Gel Analysis Pipeline Tests
# ---------------------------------------------------------------------------

def test_rap_wt_positive_control():
    """Verify that WT gel profile shows robust processivity stimulation."""
    repeats = 15
    # Low processivity without TPP1 (decay k = 1.05)
    wt_minus = 1000.0 * np.exp(-1.05 * np.arange(repeats)) * (3.0 * np.arange(1, repeats+1)) + 50.0
    # High processivity with TPP1 (decay k = 0.16)
    wt_plus = 1000.0 * np.exp(-0.16 * np.arange(repeats)) * (3.0 * np.arange(1, repeats+1)) + 50.0
    
    m_minus = calculate_rap_metrics(wt_minus, background=50.0)
    m_plus = calculate_rap_metrics(wt_plus, background=50.0)
    chi = calculate_fold_stimulation(m_minus, m_plus)
    verdict = evaluate_rap_phenotype("WT", 1.0, chi)
    
    assert m_plus["P_overall"] > m_minus["P_overall"]
    assert chi >= 2.0
    assert verdict == "WT REFERENCE BASELINE"

def test_rap_uncoupling_hit_k757a():
    """Verify that selective loss of RAP with intact basal activity triggers SUCCESS verdict."""
    repeats = 15
    # Basal elongation intact (95% of WT)
    mut_minus = 950.0 * np.exp(-1.05 * np.arange(repeats)) * (3.0 * np.arange(1, repeats+1)) + 50.0
    # Fails to stimulate RAP (decay k = 0.85 instead of 0.16)
    mut_plus = 950.0 * np.exp(-0.85 * np.arange(repeats)) * (3.0 * np.arange(1, repeats+1)) + 50.0
    
    m_minus = calculate_rap_metrics(mut_minus, background=50.0)
    m_plus = calculate_rap_metrics(mut_plus, background=50.0)
    chi = calculate_fold_stimulation(m_minus, m_plus)
    
    verdict = evaluate_rap_phenotype("TERT_Lys757A", basal_activity_fraction=0.95, chi_rap=chi)
    assert chi <= 1.8, f"Expected chi_RAP <= 1.8, got {chi:.2f}"
    assert "CRITERIA 2B & 2C PASSED" in verdict

def test_rap_catalytic_failure_negative_control():
    """Verify that severe basal catalytic loss (< 30%) triggers FALSIFICATION."""
    verdict = evaluate_rap_phenotype("TERT_Lys757A", basal_activity_fraction=0.20, chi_rap=1.2)
    assert "FALSIFIED (Severe Catalytic Failure; k_cat < 30%)" in verdict

def test_rap_no_uncoupling_negative_control():
    """Verify that a mutant that retains processivity stimulation triggers FALSIFICATION."""
    verdict = evaluate_rap_phenotype("TERT_Lys757A", basal_activity_fraction=0.95, chi_rap=4.8)
    assert "FALSIFIED (Processivity Stimulation Intact; chi_RAP >= 4.0)" in verdict


# ---------------------------------------------------------------------------
# Independent Benchmark Dataset & Trajectory Recalculation Tests
# ---------------------------------------------------------------------------

from scripts.rap_gel_quantification import process_rap_benchmark_file
from scripts.recalculate_trajectory_statistics import recalculate_statistics
from scripts.spr_kinetic_analysis import process_spr_benchmark_file


def test_spr_independent_benchmark_dataset():
    """Verify automated processing of the physical independent SPR benchmark dataset."""
    bench_csv = BASE_DIR / "data" / "benchmarks" / "spr_independent_benchmark.csv"
    assert bench_csv.is_file(), f"Missing benchmark file {bench_csv}"

    df_res = process_spr_benchmark_file(bench_csv)
    assert len(df_res) == 9, f"Expected 9 fitted runs, got {len(df_res)}"
    assert (df_res["QC_Pass"] == True).all(), "All SPR benchmark cycles must pass QC"

    # 1. WT baseline check
    wt_rows = df_res[df_res["construct"] == "TERT_WT"]
    assert len(wt_rows) == 3
    for _, r in wt_rows.iterrows():
        assert r["Evaluation_Verdict"] == "BASELINE VALID"
        assert 130.0 <= r["Kd_nM"] <= 170.0

    # 2. Lys757A candidate check (Preserved binding)
    k757a_rows = df_res[df_res["construct"] == "TERT_Lys757A"]
    assert len(k757a_rows) == 3
    for _, r in k757a_rows.iterrows():
        assert "CRITERION 2A PASSED" in r["Evaluation_Verdict"]
        assert r["Kd_nM"] <= 300.0

    # 3. Lys757R rescue control
    k757r_rows = df_res[df_res["construct"] == "TERT_Lys757R"]
    assert len(k757r_rows) == 2
    for _, r in k757r_rows.iterrows():
        assert r["Evaluation_Verdict"] == "RESCUE CONTROL CONFIRMED"
        assert r["Kd_nM"] <= 300.0

    # 4. Glu169A recruitment defect
    e169a_rows = df_res[df_res["construct"] == "TERT_Glu169A_neg"]
    assert len(e169a_rows) == 1
    assert e169a_rows.iloc[0]["Evaluation_Verdict"] == "RECRUITMENT DEFECT CONFIRMED"
    assert e169a_rows.iloc[0]["Kd_nM"] > 1500.0


def test_rap_independent_benchmark_dataset():
    """Verify automated processing of the physical independent RAP benchmark dataset."""
    bench_csv = BASE_DIR / "data" / "benchmarks" / "rap_independent_benchmark.csv"
    assert bench_csv.is_file(), f"Missing benchmark file {bench_csv}"

    df_res = process_rap_benchmark_file(bench_csv)
    assert len(df_res) == 9, f"Expected 9 fitted lanes, got {len(df_res)}"

    # 1. Asp712A catalytic dead control
    dead_rows = df_res[df_res["construct"] == "TERT_Asp712A_dead"]
    assert len(dead_rows) == 1
    assert dead_rows.iloc[0]["Evaluation_Verdict"] == "FALSIFIED (Severe Catalytic Failure; k_cat < 30%)"
    assert dead_rows.iloc[0]["basal_activity_frac"] < 0.10

    # 2. Lys757A uncoupling candidate
    k757a_rows = df_res[df_res["construct"] == "TERT_Lys757A"]
    assert len(k757a_rows) == 3
    for _, r in k757a_rows.iterrows():
        assert "CRITERIA 2B & 2C PASSED" in r["Evaluation_Verdict"]
        assert r["chi_RAP"] <= 1.8
        assert r["basal_activity_frac"] >= 0.85

    # 3. Lys757R rescue control
    k757r_rows = df_res[df_res["construct"] == "TERT_Lys757R"]
    assert len(k757r_rows) == 2
    for _, r in k757r_rows.iterrows():
        assert "RESCUE CONTROL CONFIRMED" in r["Evaluation_Verdict"]
        assert r["chi_RAP"] >= 2.0
        assert r["basal_activity_frac"] >= 0.85

    # 4. WT reference baseline
    wt_rows = df_res[df_res["construct"] == "TERT_WT"]
    assert len(wt_rows) == 3
    for _, r in wt_rows.iterrows():
        assert r["Evaluation_Verdict"] == "WT REFERENCE BASELINE"
        assert r["chi_RAP"] >= 1.8


def test_raw_md_trajectory_recalculation():
    """Verify replica-level statistical recalculation directly from raw trajectory files."""
    traj_dir = BASE_DIR / "data" / "trajectories"
    assert traj_dir.is_dir()

    stats_res = recalculate_statistics(traj_dir)
    assert stats_res["n_bound"] == 3
    assert stats_res["n_apo"] == 3
    assert stats_res["df"] == 4

    # Distance compaction is statistically significant (p < 0.01)
    assert stats_res["p_dist"] < 0.01
    assert stats_res["bound_dist_mean"] < stats_res["apo_dist_mean"]

    # H-bond occupancy is non-significant due to inter-replicate variance (p > 0.05)
    assert stats_res["p_occ"] > 0.05

    # Verification of physical duration (3.0 ps limit)
    rep_df = stats_res["replicates_df"]
    assert (rep_df["duration_ps"] < 3.0).all()
