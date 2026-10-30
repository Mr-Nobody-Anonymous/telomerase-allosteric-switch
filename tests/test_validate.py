"""
Unit and integration tests for telotrade baseline validation engine (Milestone 1).

Covers:
- SHA-256 calculation and authoritative checksum verification across all 5 structures.
- Manifest schema, integrity, and immutability checks.
- Zero-crash mmCIF parsing via Biopython and Gemmi on missing loops, disordered atoms, and altlocs.
- Asymmetric disorder handling in 2I46 (Chain A vs Chain B).
- Validation quality metrics extraction (cryo-EM Q-score, clashscore, X-ray R-free, Wilson B).
- CLI runner execution and exit code verification.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from telotrade.validate import (
    AUTHORITATIVE_HASHES,
    TARGET_PDB_IDS,
    compute_sha256,
    find_project_root,
    parse_structure_summary,
    parse_validation_metrics,
    validate_all,
    verify_checksums,
)


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Fixture providing verified project root directory."""
    return find_project_root()


@pytest.fixture(scope="session")
def raw_dir(project_root: Path) -> Path:
    """Fixture providing verified data/raw directory."""
    path = project_root / "data" / "raw"
    assert path.is_dir(), f"Raw data directory does not exist: {path}"
    return path


@pytest.fixture(scope="session")
def manifest_path(project_root: Path) -> Path:
    """Fixture providing verified data_manifest.csv path."""
    path = project_root / "data" / "manifest" / "data_manifest.csv"
    assert path.is_file(), f"Manifest file does not exist: {path}"
    return path


class TestChecksumVerification:
    """Tests for SHA-256 hashing and checksum verification."""

    def test_compute_sha256_known_content(self, tmp_path: Path):
        """Test SHA-256 calculation on known string content."""
        test_file = tmp_path / "hello.txt"
        test_file.write_bytes(b"hello world\n")
        # SHA256('hello world\n') = a948904f2f0f479b8f8197694b30184b0d2ed1c1cd2a1ec0fb85d299a192a447
        h = compute_sha256(test_file)
        assert h == "a948904f2f0f479b8f8197694b30184b0d2ed1c1cd2a1ec0fb85d299a192a447"

    def test_compute_sha256_missing_file_raises(self, tmp_path: Path):
        """Test that compute_sha256 raises FileNotFoundError for missing paths."""
        with pytest.raises(FileNotFoundError):
            compute_sha256(tmp_path / "nonexistent.cif")

    @pytest.mark.parametrize("pid", TARGET_PDB_IDS)
    def test_raw_cif_checksums_match_authoritative_values(self, raw_dir: Path, pid: str):
        """Verify that every raw mmCIF file in data/raw matches its authoritative SHA-256."""
        cif_file = raw_dir / f"{pid}.cif"
        assert cif_file.is_file(), f"Missing mmCIF file: {cif_file}"
        actual_hash = compute_sha256(cif_file)
        expected_hash = AUTHORITATIVE_HASHES[pid]
        assert actual_hash == expected_hash, (
            f"Checksum mismatch for {pid}: actual={actual_hash}, expected={expected_hash}"
        )

    def test_verify_checksums_passes_on_real_dataset(self, manifest_path: Path, raw_dir: Path):
        """Test verify_checksums returns True on real project dataset."""
        assert verify_checksums(manifest_path=manifest_path, raw_dir=raw_dir) is True

    def test_verify_checksums_fails_on_corrupted_file(self, tmp_path: Path, manifest_path: Path):
        """Test verify_checksums detects and flags a corrupted mmCIF file."""
        fake_raw = tmp_path / "raw"
        fake_raw.mkdir()
        for pid in TARGET_PDB_IDS:
            (fake_raw / f"{pid}.cif").write_text("corrupted content", encoding="utf-8")

        assert verify_checksums(manifest_path=manifest_path, raw_dir=fake_raw) is False

    def test_verify_checksums_raises_on_missing_manifest(self, tmp_path: Path, raw_dir: Path):
        """Test verify_checksums raises FileNotFoundError when given invalid manifest path."""
        with pytest.raises(FileNotFoundError):
            verify_checksums(manifest_path=tmp_path / "missing_manifest.csv", raw_dir=raw_dir)


class TestDataManifest:
    """Tests for data manifest schema, completeness, and values."""

    def test_manifest_exists_in_both_locations(self, project_root: Path):
        """Verify data_manifest.csv exists both at root and in data/manifest/."""
        root_manifest = project_root / "data_manifest.csv"
        sub_manifest = project_root / "data" / "manifest" / "data_manifest.csv"

        assert root_manifest.is_file()
        assert sub_manifest.is_file()
        assert root_manifest.read_text(encoding="utf-8") == sub_manifest.read_text(encoding="utf-8")

    def test_manifest_required_columns_present(self, manifest_path: Path):
        """Verify all required columns are present in data_manifest.csv."""
        df = pd.read_csv(manifest_path)
        required_cols = [
            "pdb_id",
            "experimental_method",
            "resolution",
            "sha256",
            "source_url",
            "chain_mappings",
            "unmodeled_segments",
            "uncertainty_notes",
        ]
        for col in required_cols:
            assert col in df.columns, f"Missing required column: {col}"

    def test_manifest_covers_all_five_targets(self, manifest_path: Path):
        """Verify manifest contains entries for all 5 target PDB IDs."""
        df = pd.read_csv(manifest_path)
        pids = set(df["pdb_id"].str.upper())
        assert pids == set(TARGET_PDB_IDS)
        assert len(df) == 5

    def test_manifest_hashes_match_authoritative_dict(self, manifest_path: Path):
        """Verify each hash in the manifest matches authoritative constants."""
        df = pd.read_csv(manifest_path)
        for _, row in df.iterrows():
            pid = row["pdb_id"].upper()
            assert row["sha256"].lower() == AUTHORITATIVE_HASHES[pid]


class TestStructureParsingZeroCrash:
    """Tests for Biopython and Gemmi parsing with zero-crash guarantee."""

    @pytest.mark.parametrize("pid", TARGET_PDB_IDS)
    def test_parse_structure_summary_succeeds(self, raw_dir: Path, pid: str):
        """Verify parse_structure_summary parses cleanly via both Biopython and Gemmi."""
        cif_path = raw_dir / f"{pid}.cif"
        summary = parse_structure_summary(cif_path)

        assert summary["status"] == "SUCCESS"
        assert summary["biopython_status"] == "SUCCESS"
        assert summary["gemmi_status"] == "SUCCESS"
        assert summary["models_count"] == 1
        assert summary["chains_count"] > 0
        assert summary["modeled_residues_count"] > 0
        assert summary["modeled_atoms_count"] > 0

    def test_2i46_asymmetry_and_disorder_handling(self, raw_dir: Path):
        """
        Verify that 2I46 parses without error despite Chain B having unobserved TEL-patch
        residues (168-171 unmodeled).
        """
        cif_path = raw_dir / "2I46.cif"
        summary = parse_structure_summary(cif_path)

        assert summary["pdb_id"] == "2I46"
        assert summary["chains_count"] == 2
        assert set(summary["chain_ids"]) == {"A", "B"}
        # Modeled residues: 374 (Chain A: 187, Chain B: 187 including waters/ligands)
        assert summary["modeled_residues_count"] == 374
        assert summary["modeled_atoms_count"] == 2304
        assert summary["unobserved_residues_count"] == 36

    def test_7tre_unmodeled_loops_and_tag_handling(self, raw_dir: Path):
        """
        Verify that 7TRE parses cleanly despite large disordered loops (TEN-TRBD linker 179-321
        and Twin-Strep tag 1-35).
        """
        cif_path = raw_dir / "7TRE.cif"
        summary = parse_structure_summary(cif_path)

        assert summary["pdb_id"] == "7TRE"
        assert summary["chains_count"] == 4
        # Author chains: B (RNA), A (TERT), D (DNA), C (TPP1)
        assert set(summary["chain_ids"]) == {"A", "B", "C", "D"}
        assert summary["modeled_residues_count"] == 1310
        assert summary["modeled_atoms_count"] == 13411
        assert summary["unobserved_residues_count"] == 489


class TestValidationMetrics:
    """Tests for extracting validation metrics from XML and mmCIF."""

    def test_7tre_validation_metrics(self, raw_dir: Path):
        """Verify 7TRE metrics extracted from validation XML and mmCIF."""
        metrics = parse_validation_metrics("7TRE", raw_dir=raw_dir)
        assert metrics["pdb_id"] == "7TRE"
        assert metrics["resolution"] == 3.50
        assert metrics["clashscore"] == 7.45
        assert metrics["emdb_id"] == "EMD-26087"
        assert metrics["atom_inclusion_all_atoms"] == 0.7940
        assert metrics["atom_inclusion_backbone"] == 0.8660
        assert metrics["q_score"] == 0.392

    def test_2i46_validation_metrics(self, raw_dir: Path):
        """Verify 2I46 X-ray crystallographic metrics extracted from validation XML."""
        metrics = parse_validation_metrics("2I46", raw_dir=raw_dir)
        assert metrics["pdb_id"] == "2I46"
        assert metrics["resolution"] == 2.70
        assert metrics["clashscore"] == 12.97
        assert metrics["r_free"] == 0.2509
        assert pytest.approx(metrics["wilson_b"], 0.01) == 43.362


class TestEndToEndValidation:
    """Tests for complete validation execution and CLI."""

    def test_validate_all_end_to_end(self, manifest_path: Path, raw_dir: Path):
        """Verify validate_all returns 100% pass on project workspace."""
        results = validate_all(raw_dir=raw_dir, manifest_path=manifest_path)
        assert results["passed"] is True
        assert results["checksums_valid"] is True
        assert results["all_structures_parsed"] is True
        assert len(results["structures"]) == 5

    def test_cli_execution_direct_script(self, project_root: Path):
        """Test running python src/telotrade/validate.py returns exit code 0."""
        script = project_root / "src" / "telotrade" / "validate.py"
        res = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0
        assert "ALL INTEGRITY CRITERIA SATISFIED" in res.stdout
        assert "Status: PASS (All SHA-256 hashes verified)" in res.stdout

    def test_cli_execution_module(self, project_root: Path):
        """Test running python -m telotrade.validate returns exit code 0."""
        res = subprocess.run(
            [sys.executable, "-m", "telotrade.validate"],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0
        assert "ALL INTEGRITY CRITERIA SATISFIED" in res.stdout

    def test_cli_execution_quiet_mode(self, project_root: Path):
        """Test running validate.py in quiet mode."""
        script = project_root / "src" / "telotrade" / "validate.py"
        res = subprocess.run(
            [sys.executable, str(script), "--quiet"],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0
        assert res.stdout.strip() == ""
