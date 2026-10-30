"""
Adversarial stress-testing suite for Milestone 1 baseline integrity and validation engine.
Authored by Challenger 1 (challenger_m1_1).

Empirically probes:
- Single-byte, multi-byte, and truncated hash corruption across all 5 targets.
- Missing CIF files and missing manifest file exception paths.
- Empty (zero-byte) files and corrupted binary files.
- Graceful error containment in Gemmi and Biopython parsers (no unhandled tracebacks).
- Corrupted validation XML handling.
- Manifest schema corruption.
- CLI exit code contracts (non-zero on integrity failure).
"""

import shutil
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


@pytest.fixture(scope="module")
def project_root() -> Path:
    return find_project_root()


@pytest.fixture(scope="module")
def clean_manifest(project_root: Path) -> Path:
    p = project_root / "data_manifest.csv"
    assert p.is_file()
    return p


@pytest.fixture(scope="module")
def clean_raw_dir(project_root: Path) -> Path:
    p = project_root / "data" / "raw"
    assert p.is_dir()
    return p


class TestAdversarialHashCorruption:
    """Stress-test hash corruption detection across all targets."""

    @pytest.mark.parametrize("target_pid", TARGET_PDB_IDS)
    def test_single_byte_flip_detected(
        self, tmp_path: Path, clean_manifest: Path, clean_raw_dir: Path, target_pid: str
    ):
        """Verify that flipping a single byte in any target CIF is detected by verify_checksums."""
        raw_copy = tmp_path / "raw"
        raw_copy.mkdir()
        for pid in TARGET_PDB_IDS:
            shutil.copy(clean_raw_dir / f"{pid}.cif", raw_copy / f"{pid}.cif")

        # Corrupt 1 byte in the selected target
        target_file = raw_copy / f"{target_pid}.cif"
        data = bytearray(target_file.read_bytes())
        data[250] = (data[250] + 1) % 256
        target_file.write_bytes(bytes(data))

        # Check that hash changed
        assert compute_sha256(target_file) != AUTHORITATIVE_HASHES[target_pid]

        # Check that verify_checksums returns False cleanly
        result = verify_checksums(manifest_path=clean_manifest, raw_dir=raw_copy)
        assert result is False

        # Check that validate_all reflects failure
        v_res = validate_all(raw_dir=raw_copy, manifest_path=clean_manifest)
        assert v_res["passed"] is False
        assert v_res["checksums_valid"] is False

    def test_all_files_truncated_detected(
        self, tmp_path: Path, clean_manifest: Path, clean_raw_dir: Path
    ):
        """Verify that truncated files are detected and rejected."""
        raw_copy = tmp_path / "raw"
        raw_copy.mkdir()
        for pid in TARGET_PDB_IDS:
            full = (clean_raw_dir / f"{pid}.cif").read_bytes()
            (raw_copy / f"{pid}.cif").write_bytes(full[:500])

        assert verify_checksums(manifest_path=clean_manifest, raw_dir=raw_copy) is False


class TestAdversarialMissingFiles:
    """Stress-test behavior when required files or directories are missing."""

    @pytest.mark.parametrize("missing_pid", TARGET_PDB_IDS)
    def test_missing_single_target_detected(
        self, tmp_path: Path, clean_manifest: Path, clean_raw_dir: Path, missing_pid: str
    ):
        """Verify that missing any single target file causes validation to fail cleanly."""
        raw_copy = tmp_path / "raw"
        raw_copy.mkdir()
        for pid in TARGET_PDB_IDS:
            if pid != missing_pid:
                shutil.copy(clean_raw_dir / f"{pid}.cif", raw_copy / f"{pid}.cif")

        result = verify_checksums(manifest_path=clean_manifest, raw_dir=raw_copy)
        assert result is False

        v_res = validate_all(raw_dir=raw_copy, manifest_path=clean_manifest)
        assert v_res["passed"] is False
        assert v_res["all_structures_parsed"] is False

    def test_missing_manifest_raises_cleanly(self, tmp_path: Path, clean_raw_dir: Path):
        """Verify that a nonexistent manifest path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError) as exc_info:
            verify_checksums(manifest_path=tmp_path / "nonexistent.csv", raw_dir=clean_raw_dir)
        assert "not found" in str(exc_info.value).lower()

    def test_missing_raw_dir_raises_cleanly(self, tmp_path: Path, clean_manifest: Path):
        """Verify that a nonexistent raw directory raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError) as exc_info:
            verify_checksums(manifest_path=clean_manifest, raw_dir=tmp_path / "nonexistent_dir")
        assert "not found" in str(exc_info.value).lower()


class TestAdversarialMalformedInputs:
    """Stress-test parser containment on invalid, empty, or garbage files."""

    def test_zero_byte_file_fails_cleanly(self, tmp_path: Path):
        """Verify that a 0-byte file does not crash parse_structure_summary."""
        empty_file = tmp_path / "empty.cif"
        empty_file.write_bytes(b"")

        summary = parse_structure_summary(empty_file)
        assert summary["status"] == "FAILED"
        assert "ERROR" in summary["gemmi_status"]

    def test_binary_garbage_fails_cleanly(self, tmp_path: Path):
        """Verify that random binary bytes do not crash parse_structure_summary."""
        junk_file = tmp_path / "junk.cif"
        junk_file.write_bytes(b"\x00\xff\xfe\x01\x80\x99" * 300)

        summary = parse_structure_summary(junk_file)
        assert summary["status"] == "FAILED"
        assert "ERROR" in summary["gemmi_status"]

    def test_valid_syntax_empty_datablock_handled(self, tmp_path: Path):
        """Verify that valid CIF syntax with zero coordinates/models does not crash."""
        cif_file = tmp_path / "empty_block.cif"
        cif_file.write_text("data_EMPTY\n# no coordinates\n", encoding="utf-8")

        summary = parse_structure_summary(cif_file)
        # Should parse without exception, reporting 0 models and 0 residues
        assert summary["models_count"] == 0
        assert summary["modeled_residues_count"] == 0
        assert summary["modeled_atoms_count"] == 0

    def test_corrupted_validation_xml_handled_gracefully(self, tmp_path: Path):
        """Verify that a corrupted (non-gzip / non-XML) validation file does not raise an unhandled exception."""
        bad_xml = tmp_path / "7TRE_validation.xml.gz"
        bad_xml.write_bytes(b"NOT A GZIP ARCHIVE")

        # parse_validation_metrics should catch the error and record it in xml_error
        metrics = parse_validation_metrics("7TRE", raw_dir=tmp_path)
        assert "xml_error" in metrics
        assert metrics["pdb_id"] == "7TRE"

    def test_corrupted_manifest_schema_raises_valueerror(self, tmp_path: Path, clean_raw_dir: Path):
        """Verify that a manifest missing required columns raises ValueError."""
        bad_manifest = tmp_path / "bad_manifest.csv"
        bad_manifest.write_text("wrong_col1,wrong_col2\n1,2\n", encoding="utf-8")

        with pytest.raises(ValueError) as exc_info:
            verify_checksums(manifest_path=bad_manifest, raw_dir=clean_raw_dir)
        assert "must contain 'pdb_id' and 'sha256'" in str(exc_info.value)


class TestAdversarialCLIContracts:
    """Stress-test CLI exit code contracts and failure propagation."""

    def test_cli_fails_on_corrupt_data(self, tmp_path: Path, clean_manifest: Path, project_root: Path):
        """Verify CLI returns exit code 1 when checksum verification fails."""
        raw_copy = tmp_path / "raw"
        raw_copy.mkdir()
        # Put 0-byte files
        for pid in TARGET_PDB_IDS:
            (raw_copy / f"{pid}.cif").write_bytes(b"corrupted")

        res = subprocess.run(
            [
                sys.executable,
                "-m",
                "telotrade.validate",
                "--raw-dir",
                str(raw_copy),
                "--manifest",
                str(clean_manifest),
            ],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 1
        assert "OVERALL RESULT: VALIDATION FAILED" in res.stdout

    def test_cli_quiet_fails_silently_with_code_1(self, tmp_path: Path, clean_manifest: Path, project_root: Path):
        """Verify CLI in quiet mode with corrupted data produces no stdout and exit code 1."""
        raw_copy = tmp_path / "raw"
        raw_copy.mkdir()
        for pid in TARGET_PDB_IDS:
            (raw_copy / f"{pid}.cif").write_bytes(b"corrupted")

        res = subprocess.run(
            [
                sys.executable,
                "-m",
                "telotrade.validate",
                "--raw-dir",
                str(raw_copy),
                "--manifest",
                str(clean_manifest),
                "--quiet",
            ],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 1
        assert res.stdout.strip() == ""


class TestIndependentManifestIntegrity:
    """Independent verification of the manifest file itself."""

    def test_root_and_sub_manifests_identical(self, project_root: Path):
        """Ensure bit-for-bit identity between root and subfolder manifests."""
        root_m = project_root / "data_manifest.csv"
        sub_m = project_root / "data" / "manifest" / "data_manifest.csv"
        assert root_m.read_bytes() == sub_m.read_bytes()

    def test_all_manifest_urls_are_https_rcsb(self, clean_manifest: Path):
        """Verify all download source URLs in manifest point to HTTPS RCSB."""
        df = pd.read_csv(clean_manifest)
        for _, row in df.iterrows():
            url = row["source_url"]
            assert url.startswith("https://files.rcsb.org/download/")
            assert url.endswith(f"{row['pdb_id']}.cif")

    def test_manifest_resolutions_match_xml_reports(self, clean_manifest: Path, clean_raw_dir: Path):
        """Cross-check manifest resolutions against parse_validation_metrics."""
        df = pd.read_csv(clean_manifest)
        for _, row in df.iterrows():
            pid = row["pdb_id"]
            manifest_res = float(row["resolution"])
            metrics = parse_validation_metrics(pid, raw_dir=clean_raw_dir)
            assert metrics["resolution"] == pytest.approx(manifest_res, rel=1e-3)
