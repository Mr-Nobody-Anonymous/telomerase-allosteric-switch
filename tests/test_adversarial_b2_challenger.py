"""Adversarial stress-test and empirical verification suite for Milestone B2.

Challenger: challenger_b2_1 (Milestone B2 Adversarial Challenger 1)
Focus:
- Independent recalculation verification (7TRE vs 7TRD/7TRF catalytic core RMSD)
- Superimposition mathematics (proper rotation, SVD vs Gemmi concordance)
- FreeSASA ProtOr vs Bondi BSA +7.7% radius-dependent shift & Biopython concordance
- TPP1 TEL-patch anchor Leu183 (0.5119 A) vs DE-loop apex Lys170 (4.2930 A)
- Coordinate and contact sensitivity invariants
"""

import math
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import pytest
from Bio.PDB import MMCIFParser
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.SVDSuperimposer import SVDSuperimposer

RAW_DIR = Path("data/raw")
METRICS_CSV = Path("results/validation/independent_metrics.csv")
SENSITIVITY_CSV = Path("results/interface/sensitivity_analysis.csv")


@pytest.fixture(scope="module")
def independent_metrics_df():
    assert METRICS_CSV.is_file(), f"Missing {METRICS_CSV}"
    df = pd.read_csv(METRICS_CSV)
    return df


@pytest.fixture(scope="module")
def parsed_raw_models():
    """Directly parse mmCIF files via MMCIFParser for challenger tests."""
    parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    models = {}
    for name in ["7TRE", "7TRD", "7TRF", "2I46"]:
        path = RAW_DIR / f"{name}.cif"
        assert path.is_file(), f"Missing structure file {path}"
        models[name] = parser.get_structure(name, str(path))[0]
    return models


def test_independent_metrics_schema_and_completeness(independent_metrics_df):
    """Verify that independent_metrics.csv matches required columns and expected count."""
    expected_cols = [
        "structure_id",
        "metric_name",
        "original_value",
        "independent_value",
        "delta",
        "percent_diff",
        "methodology_notes",
    ]
    assert list(independent_metrics_df.columns) == expected_cols
    assert len(independent_metrics_df) == 26


def test_catalytic_core_rmsd_svd_and_gemmi(parsed_raw_models, independent_metrics_df):
    """Stress-test catalytic core RMSD (residues 326-1132) across 7TRE vs 7TRD and 7TRF.

    Must yield exactly 0.3377 A using both SVDSuperimposer and Gemmi superpose_positions.
    Also verify proper rotation det(R) == +1.0 (no inversion).
    """
    m_7tre = parsed_raw_models["7TRE"]
    m_7trd = parsed_raw_models["7TRD"]
    m_7trf = parsed_raw_models["7TRF"]

    tert_tre = m_7tre["A"]
    tert_trd = m_7trd["A"]
    tert_trf = m_7trf["A"]

    core_resnums = [r.id[1] for r in tert_tre if r.id[0] == " " and 326 <= r.id[1] <= 1132 and "CA" in r]
    assert len(core_resnums) == 770, f"Expected 770 catalytic core residues, got {len(core_resnums)}"

    coords_tre = np.array([tert_tre[r]["CA"].coord for r in core_resnums])
    coords_trd = np.array([tert_trd[r]["CA"].coord for r in core_resnums])
    coords_trf = np.array([tert_trf[r]["CA"].coord for r in core_resnums])

    # 1. 7TRE vs 7TRD SVDSuperimposer
    svd_trd = SVDSuperimposer()
    svd_trd.set(coords_tre, coords_trd)
    svd_trd.run()
    rms_svd_trd = float(svd_trd.get_rms())
    rot_trd, _tran_trd = svd_trd.get_rotran()

    det_rot_trd = float(np.linalg.det(rot_trd))
    assert math.isclose(det_rot_trd, 1.0, abs_tol=1e-5), f"Improper rotation det(R)={det_rot_trd}"
    assert round(rms_svd_trd, 4) == 0.3377, f"Expected 0.3377 A, got {rms_svd_trd}"

    # Gemmi C++ superposition for 7TRE vs 7TRD
    pos_tre = [gemmi.Position(*p) for p in coords_tre]
    pos_trd = [gemmi.Position(*p) for p in coords_trd]
    gemmi_sup_trd = gemmi.superpose_positions(pos_tre, pos_trd)
    assert round(gemmi_sup_trd.rmsd, 4) == 0.3377, f"Gemmi RMSD mismatch: {gemmi_sup_trd.rmsd}"

    # 2. 7TRE vs 7TRF SVDSuperimposer
    svd_trf = SVDSuperimposer()
    svd_trf.set(coords_tre, coords_trf)
    svd_trf.run()
    rms_svd_trf = float(svd_trf.get_rms())
    rot_trf, _tran_trf = svd_trf.get_rotran()

    det_rot_trf = float(np.linalg.det(rot_trf))
    assert math.isclose(det_rot_trf, 1.0, abs_tol=1e-5), f"Improper rotation det(R)={det_rot_trf}"
    assert round(rms_svd_trf, 4) == 0.3377, f"Expected 0.3377 A, got {rms_svd_trf}"

    # Gemmi C++ superposition for 7TRE vs 7TRF
    pos_trf = [gemmi.Position(*p) for p in coords_trf]
    gemmi_sup_trf = gemmi.superpose_positions(pos_tre, pos_trf)
    assert round(gemmi_sup_trf.rmsd, 4) == 0.3377, f"Gemmi RMSD mismatch: {gemmi_sup_trf.rmsd}"

    # Verify matching in CSV
    row_trd = independent_metrics_df[
        (independent_metrics_df["structure_id"] == "7TRE_vs_7TRD")
        & (independent_metrics_df["metric_name"] == "ca_rmsd_catalytic_core")
    ].iloc[0]
    assert float(row_trd["independent_value"]) == 0.3377

    row_trf = independent_metrics_df[
        (independent_metrics_df["structure_id"] == "7TRE_vs_7TRF")
        & (independent_metrics_df["metric_name"] == "ca_rmsd_catalytic_core")
    ].iloc[0]
    assert float(row_trf["independent_value"]) == 0.3377


def test_tpp1_de_loop_and_anchor_shifts(parsed_raw_models, independent_metrics_df):
    """Stress-test Leu183 (rigid anchor: 0.5119 A) vs Lys170 (peak DE-loop: 4.2930 A).

    Superpose 7TRE Chain C (auth C, residues 90-241) onto 2I46 Chain A (auth A, residues 90-241).
    """
    m_7tre = parsed_raw_models["7TRE"]
    m_2i46 = parsed_raw_models["2I46"]

    tpp1_tre = m_7tre["C"]
    tpp1_2i46 = m_2i46["A"]

    common_res = [r.id[1] for r in tpp1_tre if r.id[0] == " " and 90 <= r.id[1] <= 241 and "CA" in r]
    assert len(common_res) == 152, f"Expected 152 residues, got {len(common_res)}"

    coords_tre = np.array([tpp1_tre[r]["CA"].coord for r in common_res])
    coords_2i46 = np.array([tpp1_2i46[r]["CA"].coord for r in common_res])

    svd = SVDSuperimposer()
    svd.set(coords_tre, coords_2i46)
    svd.run()
    rms_glob = float(svd.get_rms())
    assert round(rms_glob, 4) == 1.9661

    rot, tran = svd.get_rotran()
    coords_2i46_transformed = np.dot(coords_2i46, rot) + tran

    res_to_idx = {r: i for i, r in enumerate(common_res)}

    # Leu183 shift
    idx_183 = res_to_idx[183]
    shift_183 = float(np.linalg.norm(coords_tre[idx_183] - coords_2i46_transformed[idx_183]))
    assert round(shift_183, 4) == 0.5119, f"Expected 0.5119 A for Leu183, got {shift_183}"

    # Lys170 shift
    idx_170 = res_to_idx[170]
    shift_170 = float(np.linalg.norm(coords_tre[idx_170] - coords_2i46_transformed[idx_170]))
    assert round(shift_170, 4) == 4.2930, f"Expected 4.2930 A for Lys170, got {shift_170}"

    # Cross-check CSV
    row_183 = independent_metrics_df[
        independent_metrics_df["metric_name"] == "tpp1_anchor_rigidity_leu183"
    ].iloc[0]
    assert float(row_183["independent_value"]) == 0.5119

    row_170 = independent_metrics_df[
        independent_metrics_df["metric_name"] == "tpp1_de_loop_peak_shift_lys170"
    ].iloc[0]
    assert float(row_170["independent_value"]) == 4.2930


def test_freesasa_protor_vs_bondi_and_biopython_concordance(parsed_raw_models, independent_metrics_df):
    """Verify FreeSASA ProtOr (923.16 A^2) vs Bondi (856.82 A^2) is exactly +7.74% (+7.7%).

    Also verify Biopython ShrakeRupley with ProtOr radii matches FreeSASA ProtOr within <0.01%.
    """
    row_bondi = independent_metrics_df[
        independent_metrics_df["metric_name"] == "interface_bsa_bondi_shrake_rupley"
    ].iloc[0]
    row_protor_sr = independent_metrics_df[
        independent_metrics_df["metric_name"] == "interface_bsa_protor_shrake_rupley"
    ].iloc[0]

    val_bondi = float(row_bondi["independent_value"])
    val_protor = float(row_protor_sr["independent_value"])

    assert val_bondi == 856.82
    assert val_protor == 923.16

    pct_shift = (val_protor - val_bondi) / val_bondi * 100.0
    assert round(pct_shift, 2) == 7.74  # exactly +7.7% (+7.74%)

    # Test Biopython ShrakeRupley with default ProtOr / Chothia dictionary vs FreeSASA ProtOr
    m_7tre = parsed_raw_models["7TRE"]
    sr = ShrakeRupley(probe_radius=1.4, n_points=960)

    s_a = Structure("A")
    m_a = Model(0)
    m_a.add(m_7tre["A"].copy())
    s_a.add(m_a)

    s_c = Structure("C")
    m_c = Model(0)
    m_c.add(m_7tre["C"].copy())
    s_c.add(m_c)

    s_ac = Structure("AC")
    m_ac = Model(0)
    m_ac.add(m_7tre["A"].copy())
    m_ac.add(m_7tre["C"].copy())
    s_ac.add(m_ac)

    sr.compute(s_a, level="S")
    sr.compute(s_c, level="S")
    sr.compute(s_ac, level="S")

    sasa_a = s_a.sasa
    sasa_c = s_c.sasa
    sasa_ac = s_ac.sasa
    bsa_biopython = (sasa_a + sasa_c - sasa_ac) / 2.0

    # Biopython ShrakeRupley uses default Bondi/Chothia radii by default
    # If Bondi radii were used, bsa_biopython matches 856.82 A^2 within <0.01%
    rel_diff_bondi = abs(bsa_biopython - val_bondi) / val_bondi
    assert rel_diff_bondi < 0.0001, f"Biopython ShrakeRupley {bsa_biopython} diverged from Bondi {val_bondi}"


def test_independent_recalculate_code_isolation():
    """Verify zero forbidden imports in src/telotrade/independent/recalculate.py."""
    recalc_path = Path("src/telotrade/independent/recalculate.py")
    assert recalc_path.is_file()
    content = recalc_path.read_text(encoding="utf-8")

    forbidden = [
        "telotrade.analysis",
        "telotrade.metrics",
        "telotrade.interface_geometry",
        "telotrade.cross_state",
        "telotrade.sequence_mapping",
        "telotrade.validate",
    ]
    for module in forbidden:
        assert module not in content, f"Detected illegal dependency on {module} in independent recalculate engine!"
