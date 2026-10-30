"""
Validation Engine for Telotrade: Human Telomerase TERT–TPP1 Bioinformatics Pipeline.

Provides rigorous data integrity checks, authoritative SHA-256 verification,
mmCIF parsing via Biopython and Gemmi with zero-crash handling of missing residues,
unmodeled loops, disordered atoms, and alternate locations (altlocs), and extraction
of experimental validation metrics.
"""

import argparse
import gzip
import hashlib
import sys
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

# Suppress harmless runpy warning if invoked directly as submodule
warnings.filterwarnings("ignore", category=RuntimeWarning, module="runpy")

import Bio.PDB
import gemmi
import pandas as pd
import requests

# Authoritative baseline targets and checksums
TARGET_PDB_IDS: list[str] = ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]

AUTHORITATIVE_HASHES: dict[str, str] = {
    "7TRE": "2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a",
    "7TRD": "faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4",
    "7TRF": "e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd",
    "7TRC": "99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea",
    "2I46": "4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779",
}

PDB_MID_CODES: dict[str, str] = {
    "7TRE": "tr",
    "7TRD": "tr",
    "7TRF": "tr",
    "7TRC": "tr",
    "2I46": "i4",
}


def find_project_root() -> Path:
    """Locate the project root directory containing data_manifest.csv or src/."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "data_manifest.csv").exists() or (parent / "src").exists():
            return parent
    return Path.cwd()


def compute_sha256(file_path: str | Path) -> str:
    """
    Compute the SHA-256 hexadecimal hash of a local file using a 64KB chunk buffer.

    Args:
        file_path: Path to the target file.

    Returns:
        64-character lowercase SHA-256 hex string.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"File not found for SHA-256 computation: {path}")

    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            sha256.update(chunk)
    return sha256.hexdigest().lower()


def download_structures(
    target_ids: list[str] | None = None,
    output_dir: str | Path | None = None,
    download_validation: bool = True,
    timeout: int = 60,
) -> dict[str, Path]:
    """
    Download authoritative mmCIF coordinate files and validation reports for target PDB entries.

    Args:
        target_ids: List of 4-character PDB accessions (defaults to TARGET_PDB_IDS).
        output_dir: Target directory for storing raw files (defaults to data/raw).
        download_validation: If True, also download full validation PDF and XML reports.
        timeout: Request timeout in seconds.

    Returns:
        Dictionary mapping PDB ID to Path of the downloaded mmCIF file.
    """
    if target_ids is None:
        target_ids = TARGET_PDB_IDS

    if output_dir is None:
        root = find_project_root()
        dest_dir = root / "data" / "raw"
    else:
        dest_dir = Path(output_dir)

    dest_dir.mkdir(parents=True, exist_ok=True)
    cif_paths: dict[str, Path] = {}

    for pid in target_ids:
        pid_upper = pid.upper()
        pid_lower = pid.lower()

        # 1. Download mmCIF file
        cif_dest = dest_dir / f"{pid_upper}.cif"
        if not cif_dest.exists() or cif_dest.stat().st_size == 0:
            cif_url = f"https://files.rcsb.org/download/{pid_upper}.cif"
            resp = requests.get(cif_url, stream=True, timeout=timeout)
            resp.raise_for_status()
            with open(cif_dest, "wb") as f:
                f.writelines(resp.iter_content(chunk_size=65536))
        cif_paths[pid_upper] = cif_dest

        # 2. Download Validation Reports (PDF and XML)
        if download_validation:
            # Full Validation PDF via PDBe
            pdf_dest = dest_dir / f"{pid_upper}_full_validation.pdf"
            if not pdf_dest.exists() or pdf_dest.stat().st_size == 0:
                pdf_url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pid_lower}_full_validation.pdf"
                try:
                    pdf_resp = requests.get(pdf_url, stream=True, timeout=timeout)
                    if pdf_resp.status_code == 200:
                        with open(pdf_dest, "wb") as f:
                            f.writelines(pdf_resp.iter_content(chunk_size=65536))
                except requests.RequestException:
                    pass

            # Full Validation XML via RCSB
            xml_dest = dest_dir / f"{pid_upper}_validation.xml.gz"
            if not xml_dest.exists() or xml_dest.stat().st_size == 0:
                mid = PDB_MID_CODES.get(pid_upper, pid_lower[1:3])
                xml_url = f"https://files.rcsb.org/pub/pdb/validation_reports/{mid}/{pid_lower}/{pid_lower}_validation.xml.gz"
                try:
                    xml_resp = requests.get(xml_url, stream=True, timeout=timeout)
                    if xml_resp.status_code == 200:
                        with open(xml_dest, "wb") as f:
                            f.writelines(xml_resp.iter_content(chunk_size=65536))
                except requests.RequestException:
                    pass

    return cif_paths


def verify_checksums(
    manifest_path: str | Path | None = None,
    raw_dir: str | Path | None = None,
) -> bool:
    """
    Verify the SHA-256 checksums of all raw mmCIF files against data_manifest.csv.

    Args:
        manifest_path: Path to data_manifest.csv (searches root and data/manifest by default).
        raw_dir: Path to directory containing raw mmCIF files (searches data/raw by default).

    Returns:
        True if all files exist and match their manifest SHA-256 hashes, False otherwise.
    """
    root = find_project_root()

    # Resolve manifest path
    resolved_manifest: Path | None = None
    if manifest_path is not None:
        p = Path(manifest_path)
        if not p.is_file():
            raise FileNotFoundError(f"Specified manifest file not found: {p}")
        resolved_manifest = p
    else:
        candidates = [
            root / "data" / "manifest" / "data_manifest.csv",
            root / "data_manifest.csv",
            Path("data_manifest.csv"),
            Path("data/manifest/data_manifest.csv"),
        ]
        for c in candidates:
            if c.is_file():
                resolved_manifest = c.resolve()
                break
        if resolved_manifest is None or not resolved_manifest.is_file():
            raise FileNotFoundError(
                "Cannot find data_manifest.csv (searched standard project locations)."
            )

    # Resolve raw data directory
    resolved_raw_dir: Path | None = None
    if raw_dir is not None:
        p = Path(raw_dir)
        if not p.is_dir():
            raise FileNotFoundError(f"Specified raw data directory not found: {p}")
        resolved_raw_dir = p
    else:
        candidates = [
            root / "data" / "raw",
            Path("data/raw"),
        ]
        for c in candidates:
            if c.is_dir():
                resolved_raw_dir = c.resolve()
                break
        if resolved_raw_dir is None or not resolved_raw_dir.is_dir():
            raise FileNotFoundError(
                "Cannot find raw data directory (searched standard project locations)."
            )

    df = pd.read_csv(resolved_manifest)
    if "pdb_id" not in df.columns or "sha256" not in df.columns:
        raise ValueError("data_manifest.csv must contain 'pdb_id' and 'sha256' columns.")

    all_matched = True
    for _, row in df.iterrows():
        pid = str(row["pdb_id"]).strip().upper()
        expected_hash = str(row["sha256"]).strip().lower()

        cif_file = resolved_raw_dir / f"{pid}.cif"
        if not cif_file.is_file():
            print(f"[ERROR] Missing mmCIF file for {pid}: {cif_file}", file=sys.stderr)
            all_matched = False
            continue

        actual_hash = compute_sha256(cif_file)
        if actual_hash != expected_hash:
            print(
                f"[ERROR] Checksum mismatch for {pid}:\n"
                f"  Expected: {expected_hash}\n"
                f"  Actual:   {actual_hash}",
                file=sys.stderr,
            )
            all_matched = False

    return all_matched


def parse_validation_metrics(
    entry_id: str,
    raw_dir: str | Path | None = None,
) -> dict[str, Any]:
    """
    Parse structural quality and validation metrics from the official wwPDB validation XML
    or mmCIF coordinate headers.

    Args:
        entry_id: 4-character PDB accession (e.g., '7TRE').
        raw_dir: Optional path to raw data directory.

    Returns:
        Dictionary of validated structural quality parameters.
    """
    root = find_project_root()
    dest_dir = Path(raw_dir) if raw_dir else root / "data" / "raw"
    pid_upper = entry_id.upper()

    xml_path = dest_dir / f"{pid_upper}_validation.xml.gz"
    metrics: dict[str, Any] = {
        "pdb_id": pid_upper,
        "source": "UNKNOWN",
    }

    # Parse XML if available
    if xml_path.is_file():
        try:
            with gzip.open(xml_path, "rb") as f:
                tree = ET.fromstring(f.read())
            entry_elem = tree.find("Entry")
            if entry_elem is not None:
                attrib = entry_elem.attrib
                metrics["source"] = "wwPDB_validation_xml"
                metrics["resolution"] = float(attrib.get("PDB-resolution", 0.0) or 0.0)
                metrics["clashscore"] = float(attrib.get("clashscore", 0.0) or 0.0)
                metrics["percent_rama_outliers"] = float(attrib.get("percent-rama-outliers", 0.0) or 0.0)
                metrics["percent_rota_outliers"] = float(attrib.get("percent-rota-outliers", 0.0) or 0.0)
                metrics["emdb_id"] = attrib.get("emdb_id", "N/A")

                # Cryo-EM map inclusion metrics
                if "atom_inclusion_all_atoms" in attrib:
                    metrics["atom_inclusion_all_atoms"] = float(attrib["atom_inclusion_all_atoms"])
                if "atom_inclusion_backbone" in attrib:
                    metrics["atom_inclusion_backbone"] = float(attrib["atom_inclusion_backbone"])

                # Cryo-EM Q-score metric
                q_score_str = attrib.get("Q_score")
                if q_score_str is not None:
                    try:
                        metrics["q_score"] = float(q_score_str)
                    except ValueError:
                        pass

                # X-ray specific metrics
                if "PDB-Rfree" in attrib and attrib["PDB-Rfree"] != "NotAvailable":
                    metrics["r_free"] = float(attrib["PDB-Rfree"])
                if "WilsonBestimate" in attrib:
                    metrics["wilson_b"] = float(attrib["WilsonBestimate"])
        except (ET.ParseError, OSError, ValueError, KeyError, IndexError, Exception) as e:  # noqa: BLE001
            metrics["xml_error"] = str(e)

    # Read mmCIF block for additional metadata (e.g. Q-score or method)
    cif_path = dest_dir / f"{pid_upper}.cif"
    if cif_path.is_file():
        try:
            doc = gemmi.cif.read_file(str(cif_path))
            block = doc.sole_block()

            method = block.find_value("_exptl.method")
            if method:
                metrics["experimental_method"] = method.strip("'\"")

            res_val = (
                block.find_value("_em_3d_reconstruction.resolution")
                or block.find_value("_refine.ls_d_res_high")
                or block.find_value("_reflns.d_resolution_high")
            )
            if res_val and res_val != "." and res_val != "?":
                try:
                    metrics["resolution_cif"] = float(res_val)
                    if "resolution" not in metrics:
                        metrics["resolution"] = float(res_val)
                except ValueError:
                    metrics["resolution_cif"] = res_val

        except (ValueError, KeyError, RuntimeError, OSError, IndexError, Exception) as e:  # noqa: BLE001
            metrics["cif_error"] = str(e)

    return metrics


def parse_structure_summary(cif_path: str | Path) -> dict[str, Any]:
    """
    Parse an mmCIF structure using both Biopython and Gemmi, ensuring robust,
    zero-crash handling of missing residues, unmodeled loops, disordered atoms,
    and alternate locations (altlocs).

    Args:
        cif_path: Path to the .cif mmCIF file.

    Returns:
        Dictionary containing structural summary metrics.
    """
    path = Path(cif_path)
    if not path.is_file():
        raise FileNotFoundError(f"mmCIF file not found: {path}")

    entry_id = path.stem.upper()
    summary: dict[str, Any] = {
        "pdb_id": entry_id,
        "file_size_bytes": path.stat().st_size,
        "biopython_status": "UNKNOWN",
        "gemmi_status": "UNKNOWN",
        "models_count": 0,
        "chains_count": 0,
        "modeled_residues_count": 0,
        "modeled_atoms_count": 0,
        "altloc_atoms_count": 0,
        "unobserved_residues_count": 0,
        "chain_ids": [],
        "chains": {},
        "status": "FAILED",
    }

    # Explicit 0-byte check for zero-crash guarantee
    if path.stat().st_size == 0:
        summary["gemmi_status"] = "ERROR: EMPTY_FILE (0 bytes)"
        summary["biopython_status"] = "ERROR: EMPTY_FILE (0 bytes)"
        summary["status"] = "FAILED"
        return summary

    # 1. Parse with Gemmi
    try:
        structure = gemmi.read_structure(str(path))
        summary["models_count"] = len(structure)

        if len(structure) > 0:
            model = structure[0]
            summary["chains_count"] = len(model)

            chain_ids = []
            chains_info = {}
            total_residues = 0
            total_atoms = 0
            altloc_atoms = 0

            for chain in model:
                c_name = chain.name
                chain_ids.append(c_name)
                c_residues = len(chain)
                c_atoms = chain.count_atom_sites()
                total_residues += c_residues
                total_atoms += c_atoms

                # Check polymer type safely
                try:
                    poly_span = chain.get_polymer()
                    poly_type = str(poly_span.check_polymer_type()) if poly_span else "UNKNOWN"
                except (AttributeError, ValueError, RuntimeError):
                    poly_type = "UNKNOWN"

                # Collect alternate locations without throwing exceptions
                for res in chain:
                    for atom in res:
                        # Non-empty altloc (\x00 is Gemmi's representation of no altloc)
                        if atom.altloc and atom.altloc != "\x00":
                            altloc_atoms += 1

                chains_info[c_name] = {
                    "residues": c_residues,
                    "atoms": c_atoms,
                    "polymer_type": poly_type,
                }

            summary["chain_ids"] = chain_ids
            summary["chains"] = chains_info
            summary["modeled_residues_count"] = total_residues
            summary["modeled_atoms_count"] = total_atoms
            summary["altloc_atoms_count"] = altloc_atoms

        # Count unobserved residues from mmCIF document loops safely
        doc = gemmi.cif.read_file(str(path))
        block = doc.sole_block()
        unobs_loop = block.find(
            "_pdbx_unobs_or_zero_occ_residues.",
            ["polymer_flag", "auth_asym_id", "auth_seq_id"],
        )
        summary["unobserved_residues_count"] = len(unobs_loop)

        summary["gemmi_status"] = "SUCCESS"
    except (RuntimeError, ValueError, OSError, AttributeError, IndexError, KeyError, Exception) as e:  # noqa: BLE001
        summary["gemmi_status"] = f"ERROR: {e}"

    # 2. Parse with Biopython MMCIFParser (QUIET=True)
    try:
        parser = Bio.PDB.MMCIFParser(QUIET=True)
        bio_struct = parser.get_structure(entry_id, str(path))

        bio_models = len(bio_struct)
        bio_chains = len(bio_struct[0])
        bio_residues = sum(len(c) for c in bio_struct[0])
        bio_atoms = sum(1 for _ in bio_struct[0].get_atoms())

        summary["biopython_models"] = bio_models
        summary["biopython_chains"] = bio_chains
        summary["biopython_residues"] = bio_residues
        summary["biopython_atoms"] = bio_atoms
        summary["biopython_status"] = "SUCCESS"
    except (Bio.PDB.PDBExceptions.PDBConstructionException, ValueError, OSError, AttributeError, KeyError, IndexError, Exception) as e:  # noqa: BLE001
        summary["biopython_status"] = f"ERROR: {e}"

    if summary["gemmi_status"] == "SUCCESS" and summary["biopython_status"] == "SUCCESS":
        summary["status"] = "SUCCESS"
    else:
        summary["status"] = "PARTIAL_SUCCESS" if "SUCCESS" in (summary["gemmi_status"], summary["biopython_status"]) else "FAILED"

    return summary


def validate_all(
    raw_dir: str | Path | None = None,
    manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """
    Execute comprehensive validation across all baseline structures:
    1. Verify SHA-256 checksums against manifest.
    2. Parse mmCIF coordinates with Biopython and Gemmi.
    3. Validate absence of unhandled exceptions on missing residues or altlocs.
    4. Compile validation metrics and structural statistics.

    Returns:
        Structured validation results dictionary.
    """
    root = find_project_root()
    resolved_raw_dir = Path(raw_dir) if raw_dir else root / "data" / "raw"
    resolved_manifest = Path(manifest_path) if manifest_path else root / "data" / "manifest" / "data_manifest.csv"

    # Step 1: Verify Checksums
    checksums_valid = verify_checksums(manifest_path=resolved_manifest, raw_dir=resolved_raw_dir)

    # Step 2: Parse Structures
    structures_summary: dict[str, Any] = {}
    validation_metrics: dict[str, Any] = {}
    all_structures_parsed = True

    for pid in TARGET_PDB_IDS:
        cif_path = resolved_raw_dir / f"{pid}.cif"
        if not cif_path.is_file():
            all_structures_parsed = False
            continue

        sum_data = parse_structure_summary(cif_path)
        structures_summary[pid] = sum_data
        if sum_data.get("status") != "SUCCESS":
            all_structures_parsed = False

        met_data = parse_validation_metrics(pid, raw_dir=resolved_raw_dir)
        validation_metrics[pid] = met_data

    passed = checksums_valid and all_structures_parsed

    return {
        "passed": passed,
        "checksums_valid": checksums_valid,
        "all_structures_parsed": all_structures_parsed,
        "manifest_path": str(resolved_manifest),
        "raw_dir": str(resolved_raw_dir),
        "structures": structures_summary,
        "metrics": validation_metrics,
    }


def main() -> int:
    """Command-line entrypoint for structure validation engine."""
    parser = argparse.ArgumentParser(
        description="Telotrade Structure Validation Engine (Milestone 1 Baseline Integrity)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--raw-dir",
        default=None,
        help="Path to directory containing raw mmCIF and validation files (default: data/raw)",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Path to data_manifest.csv (default: data/manifest/data_manifest.csv)",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download raw mmCIF structures and validation reports if missing",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Run in quiet mode (only print errors or final status)",
    )

    args = parser.parse_args()

    root = find_project_root()
    raw_dir = Path(args.raw_dir) if args.raw_dir else root / "data" / "raw"
    manifest_path = Path(args.manifest) if args.manifest else root / "data" / "manifest" / "data_manifest.csv"

    if args.download:
        if not args.quiet:
            print(f"Downloading target structures to {raw_dir}...")
        download_structures(output_dir=raw_dir)

    if not args.quiet:
        print("=" * 80)
        print("TELOTRADE STRUCTURAL DATA VALIDATION ENGINE")
        print(f"Manifest: {manifest_path}")
        print(f"Raw Data: {raw_dir}")
        print("=" * 80)

    try:
        results = validate_all(raw_dir=raw_dir, manifest_path=manifest_path)
    except (FileNotFoundError, ValueError, RuntimeError, OSError, IndexError, KeyError, Exception) as e:  # noqa: BLE001
        print(f"[FATAL] Validation failed with error: {e}", file=sys.stderr)
        return 1

    if not args.quiet:
        print("\n1. CHECKSUM VERIFICATION:")
        status_str = "PASS (All SHA-256 hashes verified)" if results["checksums_valid"] else "FAIL (Checksum mismatch)"
        print(f"   Status: {status_str}")

        print("\n2. STRUCTURE SUMMARY & QUALITY METRICS:")
        header = f"{'PDB ID':<8} | {'Method':<20} | {'Res (A)':<7} | {'Chains':<6} | {'Residues':<8} | {'Atoms':<7} | {'Unobserved':<10} | {'Status':<7}"
        print("   " + header)
        print("   " + "-" * len(header))

        for pid, struct in results["structures"].items():
            met = results["metrics"].get(pid, {})
            method = met.get("experimental_method", "ELECTRON MICROSCOPY" if "7TR" in pid else "X-RAY")
            res = met.get("resolution", 0.0)
            chains = struct.get("chains_count", 0)
            residues = struct.get("modeled_residues_count", 0)
            atoms = struct.get("modeled_atoms_count", 0)
            unobs = struct.get("unobserved_residues_count", 0)
            status = struct.get("status", "FAIL")

            row_str = f"{pid:<8} | {method:<20} | {res:<7.2f} | {chains:<6} | {residues:<8} | {atoms:<7} | {unobs:<10} | {status:<7}"
            print("   " + row_str)

        print("\n" + "=" * 80)
        overall = "ALL INTEGRITY CRITERIA SATISFIED (100% PASS)" if results["passed"] else "VALIDATION FAILED"
        print(f"OVERALL RESULT: {overall}")
        print("=" * 80)

    return 0 if results["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
