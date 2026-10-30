"""
End-to-End Automated Pipeline Engine for Telotrade.

Executes all structural bioinformatics workflows:
1. Data Integrity & Checksum Verification (R1)
2. Sequence-Aware Chain & Residue Correspondence Mapping (R2)
3. Quantitative Interface Geometry & Interaction Profiling (R3)
4. Cross-State Structural Superposition & Robustness Testing (R4)
5. Standardized CSV Table Compilation in results/tables/ (R5)
6. Publication-Grade Figure Generation in results/figures/ (R5)
"""

import argparse
import shutil
import sys
import time
from pathlib import Path
from typing import Any

from .cross_state import generate_cross_state_tables
from .figures import generate_all_figures
from .interface_geometry import generate_interface_tables
from .sequence_mapping import (
    generate_chain_entity_mapping,
    generate_tert_correspondence_table,
    generate_tpp1_correspondence_table,
)
from .validate import find_project_root, verify_checksums

REQUIRED_TABLES: list[str] = [
    "data_manifest.csv",
    "chain_entity_mapping.csv",
    "residue_correspondence_tert.csv",
    "residue_correspondence_tpp1.csv",
    "interface_contacts_7tre_3.5A.csv",
    "interface_contacts_7tre_4.0A.csv",
    "interface_contacts_7tre_4.5A.csv",
    "interface_contacts_robustness.csv",
    "buried_surface_area.csv",
    "polar_interactions_7tre.csv",
    "cross_state_rmsd.csv",
]

REQUIRED_FIGURES: list[str] = [
    "figure1_interface_contacts.png",
    "figure2_cross_state_rmsd.png",
    "figure3_tpp1_conformational_shifts.png",
    "figure4_cutoff_robustness.png",
]


def resolve_directories(
    raw_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    tables_dir: str | Path | None = None,
    figures_dir: str | Path | None = None,
) -> tuple[Path, Path, Path, Path]:
    """Resolve and create standardized project directories."""
    root = find_project_root()

    # Raw mmCIF directory
    if raw_dir is not None:
        p_raw = Path(raw_dir)
    else:
        p_raw = root / "data" / "raw"
        if not p_raw.exists():
            survey_cifs = root / ".agents" / "teamwork" / "survey_pdb" / "cifs"
            if survey_cifs.exists():
                p_raw = survey_cifs

    # Output base directory
    if output_dir is not None:
        p_out = Path(output_dir)
    else:
        p_out = root / "results"

    # Tables directory
    if tables_dir is not None:
        p_tables = Path(tables_dir)
    else:
        p_tables = p_out / "tables"

    # Figures directory
    if figures_dir is not None:
        p_figures = Path(figures_dir)
    else:
        p_figures = p_out / "figures"

    p_tables.mkdir(parents=True, exist_ok=True)
    p_figures.mkdir(parents=True, exist_ok=True)

    return p_raw, p_out, p_tables, p_figures


def run_pipeline(
    raw_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
    tables_dir: str | Path | None = None,
    figures_dir: str | Path | None = None,
    manifest_path: str | Path | None = None,
    skip_validation: bool = False,
    generate_figs: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """
    Execute the full end-to-end telomerase bioinformatics pipeline.

    Args:
        raw_dir: Path to directory containing raw mmCIF files.
        output_dir: Base output directory (results/).
        tables_dir: Standardized CSV output directory (results/tables/).
        figures_dir: Publication figures directory (results/figures/).
        manifest_path: Path to data_manifest.csv.
        skip_validation: If True, bypass SHA-256 validation checks.
        generate_figs: If True, generate all 4 publication figures.
        verbose: If True, print progress logs to stdout.

    Returns:
        dict[str, Any]: Execution status, runtime metrics, generated artifact paths.
    """
    start_time = time.time()
    root = find_project_root()
    p_raw, _p_out, p_tables, p_figures = resolve_directories(raw_dir, output_dir, tables_dir, figures_dir)

    if verbose:
        print("=" * 78)
        print("Telotrade: Human Telomerase TERT-TPP1 Bioinformatics Pipeline")
        print(f"Project Root: {root}")
        print(f"Raw mmCIF Data: {p_raw}")
        print(f"Tables Target: {p_tables}")
        print(f"Figures Target: {p_figures}")
        print("=" * 78)

    # -------------------------------------------------------------------------
    # Stage 1: Data Integrity & Checksum Verification
    # -------------------------------------------------------------------------
    if not skip_validation:
        if verbose:
            print("\n[Stage 1/5] Verifying Raw Data Integrity & SHA-256 Checksums...")
        try:
            checksum_ok = verify_checksums(manifest_path=manifest_path, raw_dir=p_raw)
            if not checksum_ok:
                print("[WARNING] One or more checksums failed verification.", file=sys.stderr)
        except (FileNotFoundError, ValueError, OSError) as exc:
            if verbose:
                print(f"[NOTE] Checksum verification note: {exc}")
    else:
        if verbose:
            print("\n[Stage 1/5] Checksum verification skipped by user request.")

    # Copy data_manifest.csv to tables_dir
    src_manifest = manifest_path or root / "data_manifest.csv"
    if not Path(src_manifest).exists():
        src_manifest = root / "data" / "manifest" / "data_manifest.csv"
    dest_manifest = p_tables / "data_manifest.csv"
    if Path(src_manifest).exists() and not dest_manifest.exists():
        shutil.copy2(src_manifest, dest_manifest)

    # -------------------------------------------------------------------------
    # Stage 2: Sequence-Aware Chain & Residue Correspondence Mapping
    # -------------------------------------------------------------------------
    if verbose:
        print("\n[Stage 2/5] Establishing Sequence-Aware Residue Correspondence...")

    # 1. Chain & Entity Mapping
    df_chain = generate_chain_entity_mapping(raw_dir=p_raw)
    chain_csv = p_tables / "chain_entity_mapping.csv"
    df_chain.to_csv(chain_csv, index=False)
    if verbose:
        print(f"  [OK] {chain_csv.name} ({len(df_chain)} records)")

    # 2. TERT Correspondence
    df_tert = generate_tert_correspondence_table(raw_dir=p_raw)
    tert_csv = p_tables / "residue_correspondence_tert.csv"
    df_tert.to_csv(tert_csv, index=False)
    if verbose:
        print(f"  [OK] {tert_csv.name} ({len(df_tert)} residues)")

    # 3. TPP1 Correspondence
    df_tpp1 = generate_tpp1_correspondence_table(raw_dir=p_raw)
    tpp1_csv = p_tables / "residue_correspondence_tpp1.csv"
    df_tpp1.to_csv(tpp1_csv, index=False)
    if verbose:
        print(f"  [OK] {tpp1_csv.name} ({len(df_tpp1)} residues)")

    # -------------------------------------------------------------------------
    # Stage 3: Quantitative Interface Geometry & Interaction Profiling
    # -------------------------------------------------------------------------
    if verbose:
        print("\n[Stage 3/5] Profiling 7TRE TERT-TPP1 Interface Geometry...")

    cif_7tre = p_raw / "7TRE.cif"
    if not cif_7tre.exists():
        raise FileNotFoundError(f"7TRE.cif not found in {p_raw}")

    gen_interface = generate_interface_tables(cif_path=cif_7tre, output_dir=p_tables)
    if verbose:
        for fname, fpath in gen_interface.items():
            print(f"  [OK] {fname} ({fpath.stat().st_size} bytes)")

    # -------------------------------------------------------------------------
    # Stage 4: Cross-State Structural Superposition & Robustness Testing
    # -------------------------------------------------------------------------
    if verbose:
        print("\n[Stage 4/5] Computing Cross-State Superpositions (Kabsch SVD) & Robustness...")

    gen_cross = generate_cross_state_tables(output_dir=p_tables, raw_dir=p_raw)
    if verbose:
        for fname, fpath in gen_cross.items():
            print(f"  [OK] {fname} ({fpath.stat().st_size} bytes)")

    # -------------------------------------------------------------------------
    # Stage 5: Publication Figures Generation
    # -------------------------------------------------------------------------
    gen_figs: dict[str, Path] = {}
    if generate_figs:
        if verbose:
            print("\n[Stage 5/5] Generating Publication-Grade Structural Figures...")
        gen_figs = generate_all_figures(tables_dir=p_tables, figures_dir=p_figures)
        if verbose:
            for fname, fpath in gen_figs.items():
                print(f"  [OK] {fname} ({fpath.stat().st_size} bytes, 300 DPI)")

    elapsed = time.time() - start_time

    # Validate deliverables
    missing_tables = [t for t in REQUIRED_TABLES if not (p_tables / t).exists()]
    missing_figures = [f for f in REQUIRED_FIGURES if not (p_figures / f).exists()]

    if missing_tables:
        print(f"[ERROR] Missing required tables: {missing_tables}", file=sys.stderr)
    if missing_figures and generate_figs:
        print(f"[ERROR] Missing required figures: {missing_figures}", file=sys.stderr)

    success = len(missing_tables) == 0 and (not generate_figs or len(missing_figures) == 0)

    if verbose:
        print("\n" + "=" * 78)
        status_str = "SUCCESS" if success else "COMPLETED WITH WARNINGS"
        print(f"Pipeline Finished: {status_str} in {elapsed:.2f} seconds")
        print(f"Total Standardized Tables in {p_tables}: {len(list(p_tables.glob('*.csv')))}")
        print(f"Total Publication Figures in {p_figures}: {len(list(p_figures.glob('*.png')))}")
        print("=" * 78 + "\n")

    return {
        "success": success,
        "elapsed_seconds": elapsed,
        "tables_dir": p_tables,
        "figures_dir": p_figures,
        "missing_tables": missing_tables,
        "missing_figures": missing_figures,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments for the pipeline runner."""
    parser = argparse.ArgumentParser(
        description="Telotrade: Accelerated Human Telomerase Structural Bioinformatics Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--raw-dir",
        type=str,
        default=None,
        help="Path to directory containing raw mmCIF coordinate files (default: data/raw/)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Base output directory for pipeline deliverables (default: results/)",
    )
    parser.add_argument(
        "--tables-dir",
        type=str,
        default=None,
        help="Target directory for standardized CSV tables (default: results/tables/)",
    )
    parser.add_argument(
        "--figures-dir",
        type=str,
        default=None,
        help="Target directory for publication figures (default: results/figures/)",
    )
    parser.add_argument(
        "--manifest",
        type=str,
        default=None,
        help="Path to data_manifest.csv (default: data_manifest.csv)",
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip cryptographic SHA-256 integrity validation checks",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip publication figure rendering",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress informational stdout logging",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for telotrade.pipeline."""
    args = parse_args(argv)
    res = run_pipeline(
        raw_dir=args.raw_dir,
        output_dir=args.output_dir,
        tables_dir=args.tables_dir,
        figures_dir=args.figures_dir,
        manifest_path=args.manifest,
        skip_validation=args.skip_validation,
        generate_figs=not args.no_figures,
        verbose=not args.quiet,
    )
    return 0 if res["success"] else 1


if __name__ == "__main__":
    sys.exit(main())
