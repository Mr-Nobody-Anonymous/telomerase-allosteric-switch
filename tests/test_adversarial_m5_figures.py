"""
Adversarial stress-testing suite for Milestone 5 Remediation:
Publication Figures Generator Dynamic Sensitivity & Resource Integrity.
Authored by Challenger 1 (challenger_m5_re1).

Empirically verifies:
1. Dynamic Runtime Sensitivity of Figures to Table Perturbations:
   - Modifying values in cross_state_rmsd.csv alters figure2_cross_state_rmsd.png.
   - Modifying domain_rmsd.csv alters figure2_cross_state_rmsd.png.
   - Modifying contact row counts in interface_contacts_7tre_3.5A.csv alters figure4_cutoff_robustness.png.
   - Modifying robustness classes in interface_contacts_robustness.csv alters figure4_cutoff_robustness.png.
   - Modifying polar interaction distances alters figure1_interface_contacts.png.
   - Modifying TPP1 displacements alters figure3_tpp1_conformational_shifts.png.
2. Resource Management & Integrity:
   - Repeated figure generations leave zero unclosed Matplotlib figures (plt.get_fignums() == []).
   - Repeated runs produce valid non-zero images (> 1024 bytes).
   - Atomic replacement leaves zero residual temp files (.tmp).
"""

import hashlib
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from telotrade.figures import (
    generate_all_figures,
    generate_figure1_interface_contacts,
    generate_figure2_cross_state_rmsd,
    generate_figure3_tpp1_conformational_shifts,
    generate_figure4_cutoff_robustness,
    save_figure_safely,
)


def _file_hash(path: Path) -> str:
    """Compute SHA-256 checksum of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def tables_copy(tmp_path: Path) -> Path:
    """Create a temporary isolated copy of results/tables/ for perturbation testing."""
    repo_root = Path(__file__).resolve().parent.parent
    src_tables = repo_root / "results" / "tables"
    dest_tables = tmp_path / "tables"
    shutil.copytree(src_tables, dest_tables)
    return dest_tables


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    """Create an isolated output directory for figure generation."""
    out = tmp_path / "figures"
    out.mkdir(parents=True, exist_ok=True)
    return out


class TestFigureDynamicSensitivity:
    """Proves figures dynamically ingest and reflect changes from results/tables CSVs."""

    def test_figure2_sensitive_to_cross_state_rmsd_perturbation(
        self, tables_copy: Path, output_dir: Path
    ):
        """Altering values in cross_state_rmsd.csv must alter figure2_cross_state_rmsd.png."""
        fig_baseline = output_dir / "fig2_baseline.png"
        fig_perturbed = output_dir / "fig2_perturbed.png"

        # 1. Baseline generation
        generate_figure2_cross_state_rmsd(tables_copy, fig_baseline)
        assert fig_baseline.exists()
        assert fig_baseline.stat().st_size > 1024
        baseline_hash = _file_hash(fig_baseline)

        # 2. Perturb cross_state_rmsd.csv (Panel B data source)
        csv_path = tables_copy / "cross_state_rmsd.csv"
        df = pd.read_csv(csv_path)

        # Alter Global TERT RMSD and TPP1 OB-fold RMSD
        mask_global = (df["reference_structure"] == "7TRD") & (df["alignment_selection"] == "Global")
        assert mask_global.any(), "Expected 7TRD Global row in cross_state_rmsd.csv"
        df.loc[mask_global, "rmsd_angstrom"] = 8.888

        mask_tpp1 = (df["reference_structure"] == "2I46") & (df["alignment_selection"] == "TPP1 OB-fold")
        assert mask_tpp1.any(), "Expected 2I46 TPP1 OB-fold row in cross_state_rmsd.csv"
        df.loc[mask_tpp1, "rmsd_angstrom"] = 0.042

        df.to_csv(csv_path, index=False)

        # 3. Regenerate figure
        generate_figure2_cross_state_rmsd(tables_copy, fig_perturbed)
        assert fig_perturbed.exists()
        assert fig_perturbed.stat().st_size > 1024
        perturbed_hash = _file_hash(fig_perturbed)

        # 4. Prove file content changed
        assert baseline_hash != perturbed_hash, (
            "figure2_cross_state_rmsd.png did not change after perturbing cross_state_rmsd.csv! "
            "Figure is unresponsive to data changes."
        )

        # Verify pixel-level array difference
        img_baseline = plt.imread(fig_baseline)
        img_perturbed = plt.imread(fig_perturbed)
        assert not np.array_equal(img_baseline, img_perturbed), "Image pixel matrices are identical!"

    def test_figure2_sensitive_to_domain_rmsd_perturbation(
        self, tables_copy: Path, output_dir: Path
    ):
        """Altering values in domain_rmsd.csv must alter Panel A of Figure 2."""
        fig_baseline = output_dir / "fig2_baseline_d.png"
        fig_perturbed = output_dir / "fig2_perturbed_d.png"

        generate_figure2_cross_state_rmsd(tables_copy, fig_baseline)
        baseline_hash = _file_hash(fig_baseline)

        # Perturb domain_rmsd.csv
        csv_path = tables_copy / "domain_rmsd.csv"
        df = pd.read_csv(csv_path)
        mask = (df["reference_structure"] == "7TRD") & (df["domain_name"] == "TEN")
        assert mask.any()
        df.loc[mask, "independent_rmsd_angstrom"] = 4.750
        df.to_csv(csv_path, index=False)

        generate_figure2_cross_state_rmsd(tables_copy, fig_perturbed)
        perturbed_hash = _file_hash(fig_perturbed)

        assert baseline_hash != perturbed_hash
        img1 = plt.imread(fig_baseline)
        img2 = plt.imread(fig_perturbed)
        assert not np.array_equal(img1, img2)

    def test_figure4_sensitive_to_contact_count_perturbation(
        self, tables_copy: Path, output_dir: Path
    ):
        """Modifying row counts in interface_contacts_7tre_3.5A.csv alters figure4_cutoff_robustness.png."""
        fig_baseline = output_dir / "fig4_baseline.png"
        fig_perturbed = output_dir / "fig4_perturbed.png"

        generate_figure4_cutoff_robustness(tables_copy, fig_baseline)
        assert fig_baseline.exists()
        assert fig_baseline.stat().st_size > 1024
        baseline_hash = _file_hash(fig_baseline)

        # Perturb 3.5 A contacts: duplicate rows 10x to simulate massive contact increase
        p35 = tables_copy / "interface_contacts_7tre_3.5A.csv"
        df35 = pd.read_csv(p35)
        orig_count = len(df35)
        df35_inflated = pd.concat([df35] * 8, ignore_index=True)
        # Give distinct resnums so residue pairs count also scales
        df35_inflated["tert_resnum"] = [1000 + i for i in range(len(df35_inflated))]
        df35_inflated.to_csv(p35, index=False)

        generate_figure4_cutoff_robustness(tables_copy, fig_perturbed)
        assert fig_perturbed.exists()
        assert fig_perturbed.stat().st_size > 1024
        perturbed_hash = _file_hash(fig_perturbed)

        assert baseline_hash != perturbed_hash, (
            "figure4_cutoff_robustness.png did not change after perturbing interface_contacts_7tre_3.5A.csv! "
            f"Original row count: {orig_count}, perturbed: {len(df35_inflated)}"
        )
        img1 = plt.imread(fig_baseline)
        img2 = plt.imread(fig_perturbed)
        assert not np.array_equal(img1, img2)

    def test_figure4_sensitive_to_robustness_class_distribution(
        self, tables_copy: Path, output_dir: Path
    ):
        """Modifying robustness_class in interface_contacts_robustness.csv alters Panel B pie chart."""
        fig_baseline = output_dir / "fig4_pie_base.png"
        fig_perturbed = output_dir / "fig4_pie_pert.png"

        generate_figure4_cutoff_robustness(tables_copy, fig_baseline)
        baseline_hash = _file_hash(fig_baseline)

        rob_path = tables_copy / "interface_contacts_robustness.csv"
        df_rob = pd.read_csv(rob_path)
        # Shift all classes to 'Threshold Sensitive'
        df_rob["robustness_class"] = "Threshold Sensitive"
        df_rob.to_csv(rob_path, index=False)

        generate_figure4_cutoff_robustness(tables_copy, fig_perturbed)
        perturbed_hash = _file_hash(fig_perturbed)

        assert baseline_hash != perturbed_hash
        img1 = plt.imread(fig_baseline)
        img2 = plt.imread(fig_perturbed)
        assert not np.array_equal(img1, img2)

    def test_figure1_sensitive_to_polar_interactions_perturbation(
        self, tables_copy: Path, output_dir: Path
    ):
        """Modifying polar_interactions_7tre.csv alters Figure 1."""
        fig_baseline = output_dir / "fig1_baseline.png"
        fig_perturbed = output_dir / "fig1_perturbed.png"

        generate_figure1_interface_contacts(tables_copy, fig_baseline)
        baseline_hash = _file_hash(fig_baseline)

        polar_path = tables_copy / "polar_interactions_7tre.csv"
        df_pol = pd.read_csv(polar_path)
        # Invert all distances
        df_pol["distance_angstrom"] = df_pol["distance_angstrom"] + 1.5
        df_pol.to_csv(polar_path, index=False)

        generate_figure1_interface_contacts(tables_copy, fig_perturbed)
        perturbed_hash = _file_hash(fig_perturbed)

        assert baseline_hash != perturbed_hash
        img1 = plt.imread(fig_baseline)
        img2 = plt.imread(fig_perturbed)
        assert not np.array_equal(img1, img2)

    def test_figure3_sensitive_to_tpp1_displacements_perturbation(
        self, tables_copy: Path, output_dir: Path
    ):
        """Modifying tpp1_bound_vs_unbound_rmsd.csv alters Figure 3."""
        fig_baseline = output_dir / "fig3_baseline.png"
        fig_perturbed = output_dir / "fig3_perturbed.png"

        generate_figure3_tpp1_conformational_shifts(tables_copy, fig_baseline)
        baseline_hash = _file_hash(fig_baseline)

        tpp1_path = tables_copy / "tpp1_bound_vs_unbound_rmsd.csv"
        df_t = pd.read_csv(tpp1_path)
        # Set all displacements to 6.8 A
        df_t["ca_distance_7tre_vs_2i46_angstrom"] = 6.8
        df_t.to_csv(tpp1_path, index=False)

        generate_figure3_tpp1_conformational_shifts(tables_copy, fig_perturbed)
        perturbed_hash = _file_hash(fig_perturbed)

        assert baseline_hash != perturbed_hash
        img1 = plt.imread(fig_baseline)
        img2 = plt.imread(fig_perturbed)
        assert not np.array_equal(img1, img2)


class TestFigureResourceIntegrityAndRepetition:
    """Stress-tests repeated execution, resource leak prevention, and atomic saving."""

    def test_repeated_generation_leaves_zero_open_figures(
        self, tables_copy: Path, output_dir: Path
    ):
        """Repeatedly calling generate_all_figures must close all figures and create valid files."""
        # Ensure clean initial state
        plt.close("all")
        assert len(plt.get_fignums()) == 0

        num_iterations = 5
        for i in range(num_iterations):
            run_dir = output_dir / f"run_{i}"
            figs = generate_all_figures(tables_copy, run_dir)

            # Assert all 4 figures exist and are non-empty (> 1024 bytes)
            assert len(figs) == 4
            for name, path in figs.items():
                assert path.exists(), f"Figure {name} missing in iteration {i}"
                size = path.stat().st_size
                assert size > 1024, f"Figure {name} has size {size} <= 1024 bytes"

            # Assert no open figures in matplotlib manager (prevents memory leak)
            open_figs = plt.get_fignums()
            assert len(open_figs) == 0, f"Iteration {i} leaked open figures: {open_figs}"

            # Assert no dangling temporary files (.tmp) left behind
            tmp_files = list(run_dir.glob(".*.tmp"))
            assert len(tmp_files) == 0, f"Residual temporary files found: {tmp_files}"

    def test_save_figure_safely_prevents_zero_byte_corruption(self, output_dir: Path):
        """save_figure_safely guarantees atomic replacement and rejects empty figures."""
        target = output_dir / "test_atomic.png"

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot([0, 1], [0, 1])
        try:
            saved = save_figure_safely(fig, target)
            assert saved == target
            assert target.exists()
            assert target.stat().st_size > 1024
        finally:
            plt.close(fig)

        # Residual check
        assert len(list(output_dir.glob(".*.tmp"))) == 0
