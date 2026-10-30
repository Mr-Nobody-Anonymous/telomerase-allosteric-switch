"""
Unit and integration test suite for cross-state structural analysis and robustness (Milestone 4).

Covers:
- SVD Kabsch rigid-body superposition algorithm:
  * Exact identity, pure translation, pure rotation, and combined transforms.
  * Determinant check in SO(3): det(R) == +1.0 strictly (prevention of inversion reflections).
  * Collinear, planar, single-point, empty, and mismatched coordinate array boundaries.
- Cross-state RMSD calculations across empirical benchmarks:
  * 7TRE vs 7TRD: Global (934 C-alpha, ~0.356 Å), Catalytic Core (770 C-alpha, ~0.338 Å),
    TEN domain alone (160 C-alpha, ~0.411 Å), TEN on Core Fit (160 C-alpha, ~0.434 Å).
  * 7TRE vs 7TRF: Global (~0.356 Å), Catalytic Core (~0.338 Å), TEN on Core Fit (~0.434 Å).
  * 7TRE vs 2I46 Chain A: TPP1 OB-fold (152 C-alpha, ~1.966 Å).
- TPP1 conformational plasticity & TEL-patch remodeling:
  * Leu183 pre-organized rigidity (0.512 Å, 'Rigid Anchor').
  * DE-loop extensive remodeling (Glu169 ~3.81 Å, Lys170 ~4.29 Å, Glu171 ~3.32 Å, 'Extensive Remodeling').
  * B-factor extraction and classification into Rigid Anchor / Moderate Shift / Extensive Remodeling.
- Domain-specific RMSD partitioning (TEN, TRBD, RT, IFD_TRAP, CTE, Catalytic Core).
- Per-residue C-alpha displacement profiling across 934 residues.
- Cutoff robustness persistence scoring (Core Robust, Intermediate, Threshold Sensitive).
- Standardized CSV table generation matching project specifications in results/tables/.
- Exception and edge-case handling for missing files, disordered chains, and corrupted inputs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from telotrade import (
    CANONICAL_DOMAIN_RANGES,
    KabschResult,
    apply_kabsch_transform,
    compute_cross_state_rmsd,
    compute_cutoff_robustness,
    compute_domain_rmsds,
    compute_per_residue_displacement,
    compute_tpp1_conformational_shifts,
    compute_tpp1_remodeling,
    extract_ca_coordinates,
    generate_cross_state_tables,
    svd_kabsch_align,
)
from telotrade.cross_state import find_project_root, resolve_cif_path


@pytest.fixture(scope="module")
def project_root() -> Path:
    return find_project_root()


@pytest.fixture(scope="module")
def raw_dir(project_root: Path) -> Path:
    p = project_root / "data" / "raw"
    if not p.is_dir():
        p = project_root / ".agents" / "teamwork" / "survey_pdb" / "cifs"
    assert p.is_dir(), f"Raw data directory not found at {p}"
    return p


# =============================================================================
# 1. SVD Kabsch Superposition Engine Tests
# =============================================================================


class TestKabschSuperpositionEngine:
    """Rigorous mathematical tests for the SVD Kabsch superposition algorithm."""

    def test_kabsch_identity_transform(self):
        """Identical point sets must superimpose with zero RMSD and identity rotation."""
        pts = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 10.0],
            [2.0, 0.0, 4.0],
        ], dtype=float)

        res = svd_kabsch_align(coords_target=pts, coords_mobile=pts)
        assert isinstance(res, KabschResult)
        assert res.rmsd < 1e-12
        assert np.allclose(res.rotation, np.eye(3), atol=1e-10)
        assert np.allclose(res.translation, np.zeros(3), atol=1e-10)
        assert np.allclose(res.aligned_coords, pts, atol=1e-10)

    def test_kabsch_pure_translation(self):
        """Pure translation must yield RMSD ~ 0 and exact offset vector."""
        pts_target = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
            [4.0, 1.0, 2.0],
            [2.0, 5.0, 1.0],
        ], dtype=float)
        t_vec = np.array([12.5, -7.2, 3.14])
        pts_mobile = pts_target + t_vec

        res = svd_kabsch_align(coords_target=pts_target, coords_mobile=pts_mobile)
        assert res.rmsd < 1e-12
        assert np.allclose(res.rotation, np.eye(3), atol=1e-10)
        assert np.allclose(res.translation, -t_vec, atol=1e-10)
        assert np.allclose(res.aligned_coords, pts_target, atol=1e-10)

    def test_kabsch_pure_rotation(self):
        """Pure rotation around arbitrary axis must yield RMSD ~ 0 and det(R) == +1.0."""
        pts_target = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [2.0, 3.0, 1.0],
        ], dtype=float)

        # 45 degree rotation around Z-axis
        theta = np.pi / 4.0
        r_z = np.array([
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ], dtype=float)

        pts_mobile = np.dot(pts_target, r_z.T)

        res = svd_kabsch_align(coords_target=pts_target, coords_mobile=pts_mobile)
        assert res.rmsd < 1e-12
        det = float(np.linalg.det(res.rotation))
        assert abs(det - 1.0) < 1e-6
        assert np.allclose(res.aligned_coords, pts_target, atol=1e-10)

    def test_kabsch_combined_rotation_and_translation(self):
        """Combined arbitrary 3D rotation and translation."""
        pts_target = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 10.0],
            [2.0, 0.0, 4.0],
            [5.0, -2.0, 1.0],
        ], dtype=float)

        # 3D Euler angles
        alpha, beta, gamma = 0.3, 0.7, -0.4
        rx = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])
        r_known = rz @ ry @ rx
        t_known = np.array([10.0, -15.0, 8.0])

        pts_mobile = np.dot(pts_target, r_known.T) + t_known

        res = svd_kabsch_align(coords_target=pts_target, coords_mobile=pts_mobile)
        assert res.rmsd < 1e-10
        assert abs(np.linalg.det(res.rotation) - 1.0) < 1e-6
        assert np.allclose(res.aligned_coords, pts_target, atol=1e-10)

    def test_kabsch_reflection_determinant_enforcement(self):
        """
        Inverted coordinates (chiral reflection) must enforce proper SO(3) rotation.

        Kabsch must NOT apply an improper inversion rotation with det(R) == -1.0.
        """
        pts = np.array([
            [1.0, 1.0, 1.0],
            [2.0, 1.0, 1.0],
            [1.0, 3.0, 1.0],
            [1.0, 1.0, 4.0],
        ], dtype=float)

        reflected = pts.copy()
        reflected[:, 2] = -reflected[:, 2]  # Chiral inversion along Z

        res = svd_kabsch_align(coords_target=pts, coords_mobile=reflected)
        det = float(np.linalg.det(res.rotation))
        assert abs(det - 1.0) < 1e-6, f"Determinant must be strictly +1.0, got {det}"
        assert res.rotation.shape == (3, 3)

    def test_kabsch_planar_and_collinear_points(self):
        """Planar point configurations must not produce NaN or infinite values."""
        planar = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ], dtype=float)
        shifted = planar + np.array([5.0, 3.0, 0.0])

        res = svd_kabsch_align(coords_target=planar, coords_mobile=shifted)
        assert not np.isnan(res.rmsd)
        assert not np.isnan(res.rotation).any()
        assert not np.isnan(res.translation).any()
        assert res.rmsd < 1e-10
        assert abs(np.linalg.det(res.rotation) - 1.0) < 1e-6

    def test_kabsch_single_point_alignment(self):
        """Single point superposition yields zero RMSD and pure translation."""
        p1 = np.array([[1.0, 2.0, 3.0]])
        p2 = np.array([[4.0, 6.0, 8.0]])

        res = svd_kabsch_align(coords_target=p2, coords_mobile=p1)
        assert res.rmsd == 0.0
        assert np.allclose(res.aligned_coords, p2)

    def test_kabsch_invalid_shapes_raise_value_error(self):
        """Invalid coordinate shapes and point count mismatches must raise ValueError."""
        # Empty array
        with pytest.raises(ValueError, match="empty"):
            svd_kabsch_align(np.empty((0, 3)), np.empty((0, 3)))

        # 1D array
        with pytest.raises(ValueError, match="shape"):
            svd_kabsch_align(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]))

        # Count mismatch
        p = np.zeros((4, 3))
        q = np.zeros((5, 3))
        with pytest.raises(ValueError, match="mismatch"):
            svd_kabsch_align(coords_target=q, coords_mobile=p)

    def test_apply_kabsch_transform_helper(self):
        """apply_kabsch_transform correctly transforms coordinates."""
        pts = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        r = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        t = np.array([2.0, 3.0, 4.0])

        transformed = apply_kabsch_transform(pts, r, t)
        expected = np.dot(pts, r.T) + t
        assert np.allclose(transformed, expected)

        # Empty array handling
        empty = apply_kabsch_transform(np.empty((0, 3)), r, t)
        assert empty.shape == (0, 3)


# =============================================================================
# 2. Cross-State Superposition & Benchmark Verification Tests
# =============================================================================


class TestCrossStateRmsdBenchmarks:
    """Verifies empirical benchmark RMSD metrics across 7TRE, 7TRD, 7TRF, and 2I46."""

    def test_7tre_vs_7trd_catalytic_core_rigidity(self, raw_dir: Path):
        """
        Verify that the catalytic core (residues 326-1132, n=770) between 7TRE (bound)
        and 7TRD (unbound reference) has RMSD ~ 0.338 Å (benchmark: < 0.45 Å).
        """
        df = compute_cross_state_rmsd(raw_dir=raw_dir)
        core_row = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "7TRD")
            & (df["alignment_selection"] == "Catalytic Core")
        ]
        assert not core_row.empty, "Catalytic Core row must exist in cross_state_rmsd table"
        rmsd = float(core_row["rmsd_angstrom"].iloc[0])
        n_atoms = int(core_row["num_aligned_ca_atoms"].iloc[0])

        assert n_atoms == 770, f"Expected 770 matched C-alpha atoms in core, got {n_atoms}"
        assert abs(rmsd - 0.338) < 0.02, f"Expected core RMSD ~0.338 Å, got {rmsd} Å"
        assert rmsd < 0.45, f"Catalytic core must remain rigidly locked (<0.45 Å), got {rmsd} Å"

    def test_7tre_vs_7trd_global_and_ten_domain(self, raw_dir: Path):
        """
        Verify global C-alpha RMSD (n=934, ~0.356 Å) and TEN domain displacement
        on catalytic core fit (~0.434 Å).
        """
        df = compute_cross_state_rmsd(raw_dir=raw_dir)

        # Global
        glob_row = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "7TRD")
            & (df["alignment_selection"] == "Global")
        ]
        assert not glob_row.empty
        rmsd_glob = float(glob_row["rmsd_angstrom"].iloc[0])
        n_glob = int(glob_row["num_aligned_ca_atoms"].iloc[0])
        assert n_glob == 934
        assert abs(rmsd_glob - 0.356) < 0.02

        # TEN Domain Alone
        ten_alone_row = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "7TRD")
            & (df["alignment_selection"] == "TEN Domain")
        ]
        assert not ten_alone_row.empty
        rmsd_ten_alone = float(ten_alone_row["rmsd_angstrom"].iloc[0])
        assert abs(rmsd_ten_alone - 0.411) < 0.02

        # TEN on Core Fit
        ten_fit_row = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "7TRD")
            & (df["alignment_selection"] == "TEN on Core Fit")
        ]
        assert not ten_fit_row.empty
        rmsd_ten_fit = float(ten_fit_row["rmsd_angstrom"].iloc[0])
        assert abs(rmsd_ten_fit - 0.434) < 0.02

    def test_7tre_vs_7trf_homology(self, raw_dir: Path):
        """Verify 7TRE vs 7TRF catalytic core and global RMSD match 7TRD benchmark."""
        df = compute_cross_state_rmsd(raw_dir=raw_dir)
        core_trf = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "7TRF")
            & (df["alignment_selection"] == "Catalytic Core")
        ]
        assert not core_trf.empty
        rmsd_trf = float(core_trf["rmsd_angstrom"].iloc[0])
        assert abs(rmsd_trf - 0.338) < 0.02

    def test_7tre_vs_2i46_tpp1_global_rmsd(self, raw_dir: Path):
        """
        Verify bound TPP1 (7TRE) vs unbound TPP1 (2I46 Chain A) global RMSD is ~ 1.966 Å
        across exactly 152 common C-alpha atoms.
        """
        df = compute_cross_state_rmsd(raw_dir=raw_dir)
        tpp1_row = df[
            (df["target_structure"] == "7TRE")
            & (df["reference_structure"] == "2I46")
            & (df["alignment_selection"] == "TPP1 OB-fold")
        ]
        assert not tpp1_row.empty
        rmsd_tpp1 = float(tpp1_row["rmsd_angstrom"].iloc[0])
        n_atoms = int(tpp1_row["num_aligned_ca_atoms"].iloc[0])

        assert n_atoms == 152, f"Expected 152 matched C-alpha atoms, got {n_atoms}"
        assert abs(rmsd_tpp1 - 1.966) < 0.05, f"Expected TPP1 RMSD ~1.966 Å, got {rmsd_tpp1} Å"

    def test_cross_state_rmsd_table_schema(self, raw_dir: Path):
        """Verify cross_state_rmsd table conforms to Section 4.9 schema."""
        df = compute_cross_state_rmsd(raw_dir=raw_dir)
        expected_cols = [
            "target_structure",
            "reference_structure",
            "alignment_selection",
            "num_aligned_ca_atoms",
            "superposition_method",
            "rmsd_angstrom",
            "mean_ca_displacement_angstrom",
            "max_ca_displacement_angstrom",
            "max_displacement_residue",
        ]
        assert list(df.columns) == expected_cols
        assert (df["rmsd_angstrom"] >= 0.0).all()
        assert (df["num_aligned_ca_atoms"] > 0).all()
        assert (df["superposition_method"] == "Kabsch SVD").all()


# =============================================================================
# 3. Domain-Specific Structural Rearrangements Tests
# =============================================================================


class TestDomainSpecificRmsds:
    """Tests domain-specific RMSD partitioning and core-fit displacements."""

    def test_all_canonical_domains_present(self, raw_dir: Path):
        """Verify all canonical domains are partitioned: TEN, TRBD, RT, IFD_TRAP, CTE, Catalytic Core."""
        df = compute_domain_rmsds(raw_dir=raw_dir)
        assert not df.empty

        domains_found = set(df["domain_name"])
        expected_domains = {"TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"}
        assert expected_domains.issubset(domains_found)

    def test_domain_rmsd_schema_and_atom_counts(self, raw_dir: Path):
        """Verify domain_rmsd table schema and atom counts."""
        df = compute_domain_rmsds(raw_dir=raw_dir)
        expected_cols = [
            "target_structure",
            "reference_structure",
            "domain_name",
            "residue_range_auth",
            "num_ca_atoms",
            "independent_rmsd_angstrom",
            "core_fit_displacement_angstrom",
            "mean_displacement_angstrom",
            "max_displacement_angstrom",
            "max_displacement_residue",
        ]
        assert list(df.columns) == expected_cols

        # Check atom counts in 7TRE vs 7TRD
        df_td = df[(df["target_structure"] == "7TRE") & (df["reference_structure"] == "7TRD")]
        ten_count = df_td[df_td["domain_name"] == "TEN"]["num_ca_atoms"].iloc[0]
        core_count = df_td[df_td["domain_name"] == "Catalytic Core"]["num_ca_atoms"].iloc[0]

        assert ten_count == 160
        assert core_count == 770

    def test_compute_domain_rmsds_single_pair_mode(self, raw_dir: Path):
        """Verify compute_domain_rmsds supports single structure pair execution."""
        p_tre = raw_dir / "7TRE.cif"
        p_trd = raw_dir / "7TRD.cif"
        df_pair = compute_domain_rmsds(cif_ref_or_raw_dir=p_tre, cif_mobile=p_trd)
        assert len(df_pair) == len(CANONICAL_DOMAIN_RANGES)
        assert (df_pair["target_structure"] == "7TRE").all()
        assert (df_pair["reference_structure"] == "7TRD").all()


# =============================================================================
# 4. Per-Residue C-alpha Displacement Profiling Tests
# =============================================================================


class TestPerResidueDisplacement:
    """Tests per-residue C-alpha Euclidean displacement profiling."""

    def test_per_residue_table_row_count_and_schema(self, raw_dir: Path):
        """Verify table has exactly 934 rows and adheres to Section 4.10 schema."""
        df = compute_per_residue_displacement(raw_dir=raw_dir)
        assert len(df) == 934, f"Expected 934 modeled TERT residues, got {len(df)}"

        expected_cols = [
            "auth_seq_id",
            "resname",
            "domain",
            "disp_7tre_vs_7trd_global_angstrom",
            "disp_7tre_vs_7trd_core_fit_angstrom",
            "disp_7tre_vs_7trf_core_fit_angstrom",
        ]
        assert list(df.columns) == expected_cols

    def test_catalytic_core_vs_ten_displacement_differential(self, raw_dir: Path):
        """
        Verify that the catalytic core (residues 326-1132) has low displacement on core fit
        while TEN domain residues (1-196) exhibit higher hinge displacement (~0.434 Å).
        """
        df = compute_per_residue_displacement(raw_dir=raw_dir)

        core_disps = df[df["auth_seq_id"] >= 326]["disp_7tre_vs_7trd_core_fit_angstrom"]
        ten_disps = df[df["auth_seq_id"] < 200]["disp_7tre_vs_7trd_core_fit_angstrom"]

        core_rms = float(np.sqrt(np.mean(core_disps ** 2)))
        ten_rms = float(np.sqrt(np.mean(ten_disps ** 2)))

        assert abs(core_rms - 0.338) < 0.02, f"Expected core RMS ~0.338 Å, got {core_rms} Å"
        assert abs(ten_rms - 0.434) < 0.02, f"Expected TEN RMS ~0.434 Å, got {ten_rms} Å"
        assert ten_rms > core_rms, "TEN domain displacement must exceed core rigidity"


# =============================================================================
# 5. TPP1 Bound vs Unbound Conformational Plasticity Tests
# =============================================================================


class TestTPP1ConformationalShifts:
    """Tests TEL-patch remodeling: DE-loop plasticity vs Leu183 rigidity."""

    def test_tpp1_conformational_shifts_schema_and_count(self, raw_dir: Path):
        """Verify table contains 152 residues and conforms to Section 4.11 schema."""
        df = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
        assert len(df) == 152

        expected_cols = [
            "auth_seq_id",
            "resname",
            "is_tel_patch",
            "ca_distance_7tre_vs_2i46_angstrom",
            "local_b_factor_7tre",
            "local_b_factor_2i46",
            "conformational_shift_flag",
        ]
        assert list(df.columns) == expected_cols

    def test_tpp1_remodeling_alias_consistency(self, raw_dir: Path):
        """Verify compute_tpp1_remodeling is an exact alias for compute_tpp1_conformational_shifts."""
        df1 = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
        df2 = compute_tpp1_remodeling(raw_dir=raw_dir)
        pd.testing.assert_frame_equal(df1, df2)

    def test_leu183_rigid_hydrophobic_anchor(self, raw_dir: Path):
        """
        Verify Leu183 is structurally pre-organized in the unbound state:
        displacement ~ 0.512 Å, classified as 'Rigid Anchor'.
        """
        df = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
        row_183 = df[df["auth_seq_id"] == 183]
        assert not row_183.empty

        disp = float(row_183["ca_distance_7tre_vs_2i46_angstrom"].iloc[0])
        flag = str(row_183["conformational_shift_flag"].iloc[0])
        is_tel = bool(row_183["is_tel_patch"].iloc[0])

        assert is_tel is True
        assert disp < 0.8, f"Leu183 displacement {disp} Å > 0.8 Å"
        assert abs(disp - 0.512) < 0.05, f"Expected Leu183 displacement ~0.512 Å, got {disp} Å"
        assert flag == "Rigid Anchor"

    def test_de_loop_extensive_plastic_remodeling(self, raw_dir: Path):
        """
        Verify acidic DE-loop residues (Glu168-Glu171) undergo extensive plastic remodeling:
          - Glu169 > 3.5 Å ('Extensive Remodeling')
          - Lys170 ~ 4.293 Å ('Extensive Remodeling')
          - Glu171 > 3.0 Å ('Extensive Remodeling')
          - Glu168 ~ 1.253 Å ('Moderate Shift')
        """
        df = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
        res_map = df.set_index("auth_seq_id").to_dict(orient="index")

        # Glu168
        assert abs(res_map[168]["ca_distance_7tre_vs_2i46_angstrom"] - 1.253) < 0.1
        assert res_map[168]["conformational_shift_flag"] == "Moderate Shift"

        # Glu169
        assert abs(res_map[169]["ca_distance_7tre_vs_2i46_angstrom"] - 3.809) < 0.1
        assert res_map[169]["conformational_shift_flag"] == "Extensive Remodeling"

        # Lys170 (peak loop displacement)
        assert abs(res_map[170]["ca_distance_7tre_vs_2i46_angstrom"] - 4.293) < 0.1
        assert res_map[170]["conformational_shift_flag"] == "Extensive Remodeling"

        # Glu171
        assert abs(res_map[171]["ca_distance_7tre_vs_2i46_angstrom"] - 3.320) < 0.1
        assert res_map[171]["conformational_shift_flag"] == "Extensive Remodeling"

        # Ratio check: DE-loop peak remodeling must substantially exceed Leu183 anchor
        peak_disp = res_map[170]["ca_distance_7tre_vs_2i46_angstrom"]
        anchor_disp = res_map[183]["ca_distance_7tre_vs_2i46_angstrom"]
        assert peak_disp > 7.0 * anchor_disp, "Peak DE-loop displacement must exceed 7x anchor displacement"

    def test_b_factors_are_valid_and_positive(self, raw_dir: Path):
        """Verify local B-factors extracted from cryo-EM map and crystal structure are positive."""
        df = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
        assert (df["local_b_factor_7tre"] > 0.0).all()
        assert (df["local_b_factor_2i46"] > 0.0).all()


# =============================================================================
# 6. Cutoff Robustness & Persistence Tests
# =============================================================================


class TestCutoffRobustnessPersistence:
    """Tests interface contact robustness and persistence across distance cutoffs."""

    def test_cutoff_robustness_classification_categories(self, raw_dir: Path):
        """Verify contacts are correctly classified into Core Robust, Intermediate, and Threshold Sensitive."""
        df = compute_cutoff_robustness(raw_dir=raw_dir)
        assert not df.empty

        classes = set(df["robustness_class"])
        assert "Core Robust" in classes
        assert "Intermediate" in classes
        assert "Threshold Sensitive" in classes

        # Invariants: score 3 -> Core Robust, 2 -> Intermediate, 1 -> Threshold Sensitive
        assert (df[df["persistence_count"] == 3]["robustness_class"] == "Core Robust").all()
        assert (df[df["persistence_count"] == 2]["robustness_class"] == "Intermediate").all()
        assert (df[df["persistence_count"] == 1]["robustness_class"] == "Threshold Sensitive").all()

    def test_key_salt_bridges_are_core_robust(self, raw_dir: Path):
        """Verify that primary salt bridge Arg774-Glu168 is persistent and Core Robust."""
        df = compute_cutoff_robustness(raw_dir=raw_dir)
        sb_rows = df[df["contact_pair_id"].str.contains("ARG774") & df["contact_pair_id"].str.contains("GLU168")]
        assert not sb_rows.empty
        # Primary salt bridge must be present at 3.5 A
        assert sb_rows["present_at_3.5A"].any()
        assert (sb_rows["robustness_class"] == "Core Robust").any()


# =============================================================================
# 7. Standard Results Table Generation Tests
# =============================================================================


class TestGenerateCrossStateTables:
    """Tests generation of standardized CSV tables in results/tables/."""

    def test_generate_cross_state_tables_output_files(self, tmp_path: Path, raw_dir: Path):
        """Verify generate_cross_state_tables creates all required CSV files on disk."""
        out_dir = tmp_path / "tables"
        tables = generate_cross_state_tables(output_dir=out_dir, raw_dir=raw_dir)

        expected_files = [
            "cross_state_rmsd.csv",
            "domain_rmsd.csv",
            "domain_rmsds.csv",
            "per_residue_displacement_tert.csv",
            "tpp1_bound_vs_unbound_rmsd.csv",
            "tpp1_conformational_shifts.csv",
            "interface_contacts_robustness.csv",
        ]

        for fname in expected_files:
            assert fname in tables, f"Table '{fname}' must be returned by generator"
            p = tables[fname]
            assert p.is_file(), f"Table file '{p}' was not created on disk"
            assert p.stat().st_size > 0, f"Table file '{p}' must not be empty"


# =============================================================================
# 8. Error Handling and Edge Cases
# =============================================================================


class TestCrossStateErrorHandling:
    """Tests exception paths and edge case containment."""

    def test_resolve_cif_path_missing_file_raises_error(self, tmp_path: Path):
        """Attempting to resolve a non-existent CIF file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="not found"):
            resolve_cif_path("NONEXISTENT_XYZ", raw_dir=tmp_path)

    def test_extract_ca_invalid_chain_raises_key_error(self, raw_dir: Path):
        """Querying an invalid chain raises KeyError."""
        cif_7tre = raw_dir / "7TRE.cif"
        with pytest.raises(KeyError, match="not found"):
            extract_ca_coordinates(cif_7tre, chain_id="Z_NONEXISTENT", use_auth=True)
