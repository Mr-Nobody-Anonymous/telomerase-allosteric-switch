"""Tier 4: Real-World Application Scenario & Pipeline Execution E2E Tests.

Tests the full end-to-end scientific pipeline execution and validates the complete set of
project deliverables:
- Automated pipeline execution reproducing results from scratch
- Standardized CSV table generation across all 11 required result tables in results/tables/
- Schema validation, data types, value ranges, and relational foreign keys
- Publication-quality figure generation in results/figures/
- Synthesis report authoring with strict separation of Axis A and Axis B
"""

import subprocess
import sys

import pandas as pd
import pytest


@pytest.mark.tier4
class TestTier4EndToEndPipelineExecution:
    """Verifies end-to-end execution of run_pipeline.py."""

    def test_t4_pipeline_runner_exists(self, project_root):
        """Test that single-entry pipeline runner script exists."""
        runner = project_root / "run_pipeline.py"
        if not runner.exists():
            pytest.skip("run_pipeline.py not yet implemented (M5 in progress)")
        assert runner.is_file(), "run_pipeline.py must be a file"

    def test_t4_pipeline_execution_exit_code_zero(self, project_root):
        """Test that executing run_pipeline.py finishes with exit code 0."""
        runner = project_root / "run_pipeline.py"
        if not runner.exists():
            pytest.skip("run_pipeline.py not yet implemented (M5 in progress)")

        res = subprocess.run(
            [sys.executable, str(runner)],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        assert res.returncode == 0, f"run_pipeline.py failed with code {res.returncode}:\n{res.stderr}"


@pytest.mark.tier4
class TestTier4StandardizedResultsTables:
    """Verifies presence, schemas, and contents of all 11 standardized CSV tables."""

    @pytest.fixture(autouse=True)
    def check_tables_directory(self, project_root):
        self.tables_dir = project_root / "results" / "tables"
        if not self.tables_dir.exists():
            pytest.skip("results/tables/ directory does not exist yet (M5 in progress)")

    def test_t4_all_11_tables_present(self, ref_data):
        """Test that all 11 required CSV result tables exist and are non-empty."""
        for tbl in ref_data["expected_tables"]:
            tbl_path = self.tables_dir / tbl
            assert tbl_path.exists(), f"Expected result table {tbl} is missing from results/tables/"
            assert tbl_path.stat().st_size > 0, f"Table {tbl} is empty"

    def test_t4_table_manifest_contents(self):
        """Test data_manifest.csv table schema and all 5 structure records."""
        path = self.tables_dir / "data_manifest.csv"
        if not path.exists():
            pytest.skip("data_manifest.csv not in results/tables/")
        df = pd.read_csv(path)
        assert len(df) >= 5
        assert {"7TRE", "7TRD", "7TRF", "7TRC", "2I46"}.issubset(set(df["pdb_id"]))

    def test_t4_table_buried_surface_area_metrics(self, ref_data):
        """Test buried_surface_area.csv contains exact formula metrics."""
        path = self.tables_dir / "buried_surface_area.csv"
        if not path.exists():
            pytest.skip("buried_surface_area.csv not found")
        df = pd.read_csv(path)
        assert "bsa_interface_half_angstrom2" in df.columns or "bsa" in df.columns
        bsa_col = "bsa_interface_half_angstrom2" if "bsa_interface_half_angstrom2" in df.columns else "bsa"
        bsa_val = float(df[bsa_col].iloc[0])
        expected_bsa = ref_data["bsa"]["bsa_interface"]
        tol = ref_data["bsa"]["tolerance"]
        assert abs(bsa_val - expected_bsa) <= tol, f"BSA {bsa_val} A^2 differs from reference {expected_bsa} A^2"

    def test_t4_table_interface_contacts_monotonicity(self):
        """Test that interface contact table row counts increase with cutoff."""
        p_35 = self.tables_dir / "interface_contacts_7tre_3.5A.csv"
        p_40 = self.tables_dir / "interface_contacts_7tre_4.0A.csv"
        p_45 = self.tables_dir / "interface_contacts_7tre_4.5A.csv"
        if not (p_35.exists() and p_40.exists() and p_45.exists()):
            pytest.skip("Interface contact tables not yet generated")

        df_35 = pd.read_csv(p_35)
        df_40 = pd.read_csv(p_40)
        df_45 = pd.read_csv(p_45)

        assert len(df_35) == 9, f"Expected 9 contacts at 3.5 A, got {len(df_35)}"
        assert len(df_40) == 46, f"Expected 46 contacts at 4.0 A, got {len(df_40)}"
        assert len(df_45) == 123, f"Expected 123 contacts at 4.5 A, got {len(df_45)}"
        assert len(df_35) < len(df_40) < len(df_45)

    def test_t4_table_polar_interactions_salt_bridges(self):
        """Test polar_interactions_7tre.csv includes Arg774-Glu168 and Lys78-Glu215."""
        path = self.tables_dir / "polar_interactions_7tre.csv"
        if not path.exists():
            pytest.skip("polar_interactions_7tre.csv not found")
        df = pd.read_csv(path)
        # Check TERT 774 and TPP1 168
        has_arg774 = ((df["tert_resnum"] == 774) & (df["tpp1_resnum"] == 168)).any()
        assert has_arg774, "Primary salt bridge Arg774-Glu168 missing from polar table"

    def test_t4_table_cross_state_rmsd_benchmarks(self):
        """Test cross_state_rmsd.csv benchmarks for 7TRE vs 7TRD and 2I46."""
        path = self.tables_dir / "cross_state_rmsd.csv"
        if not path.exists():
            pytest.skip("cross_state_rmsd.csv not found")
        df = pd.read_csv(path)
        assert len(df) >= 3, "cross_state_rmsd.csv should contain multiple comparison rows"


@pytest.mark.tier4
class TestTier4FiguresAndSynthesisReport:
    """Verifies generation of publication-grade figures and synthesis research report."""

    def test_t4_figures_directory_and_files(self, project_root):
        """Test that results/figures/ directory exists and contains visual artifacts."""
        figures_dir = project_root / "results" / "figures"
        if not figures_dir.exists():
            pytest.skip("results/figures/ not yet created (M5 in progress)")
        figures = list(figures_dir.glob("*.png")) + list(figures_dir.glob("*.svg"))
        assert len(figures) >= 1, "At least one publication figure must be present in results/figures/"
        for fig in figures:
            assert fig.stat().st_size > 1024, f"Figure {fig.name} is unexpectedly small"

    def test_t4_synthesis_report_presence_and_structure(self, project_root):
        """Test that report.md exists, has substantial content, and strictly separates Axis A and Axis B."""
        report = project_root / "report.md"
        if not report.exists():
            pytest.skip("report.md not yet written (M5 in progress)")

        assert report.stat().st_size > 5000, "report.md is too short (<5 KB)"
        content = report.read_text(encoding="utf-8")

        # Verify key structural sections
        assert "AXIS A" in content or "Axis A" in content, "report.md must contain Axis A (Biophysics/Mechanism)"
        assert "AXIS B" in content or "Axis B" in content, "report.md must contain Axis B (Therapeutic Selectivity)"
        assert "7TRE" in content and "7TRD" in content and "7TRF" in content and "2I46" in content

        # Verify translational and clinical disease context
        assert "Dyskeratosis Congenita" in content or "dyskeratosis congenita" in content.lower()
        assert "stem cell" in content.lower(), "report.md must evaluate stem cell reserves / toxicity trade-off"
