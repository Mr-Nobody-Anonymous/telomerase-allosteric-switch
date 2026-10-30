"""Test suite verifying deliverables for Benchmark Milestones B4 and B5.

Covers:
- results/translation/cancer_aging_evidence_matrix.csv
- results/translation/therapeutic_window_limits.md
- results/validation_experiments/falsifiable_predictions.md
- results/final_adversarial_review.md
- report.md (13 mandatory sections, Axis A/B separation, de-hyped catalog)
"""

from pathlib import Path

import pandas as pd

RESULTS_DIR = Path("results")
REPORT_PATH = Path("report.md")


def test_cancer_aging_evidence_matrix_schema_and_contents():
    """Verify cancer_aging_evidence_matrix.csv schema and targets."""
    matrix_path = RESULTS_DIR / "translation" / "cancer_aging_evidence_matrix.csv"
    assert matrix_path.is_file(), f"Missing {matrix_path}"
    df = pd.read_csv(matrix_path)

    expected_cols = [
        "target_or_interface",
        "cancer_dependency",
        "stem_cell_dependency",
        "therapeutic_window_evidence",
        "experimental_status",
        "citation",
    ]
    for col in expected_cols:
        assert col in df.columns, f"Missing column {col} in {matrix_path}"

    assert len(df) >= 5, f"Expected at least 5 target rows, got {len(df)}"

    targets = set(df["target_or_interface"])
    assert any("catalytic" in t.lower() for t in targets), "Missing catalytic active site target"
    assert any("recruitment" in t.lower() or "patch" in t.lower() for t in targets), "Missing recruitment interface target"
    assert any("alt" in t.lower() for t in targets), "Missing ALT target"


def test_therapeutic_window_limits_report():
    """Verify therapeutic_window_limits.md existence, size, and critical content."""
    doc_path = RESULTS_DIR / "translation" / "therapeutic_window_limits.md"
    assert doc_path.is_file(), f"Missing {doc_path}"
    assert doc_path.stat().st_size > 3000, "therapeutic_window_limits.md is too short"

    content = doc_path.read_text(encoding="utf-8")
    assert "Dyskeratosis Congenita" in content or "dyskeratosis congenita" in content.lower()
    assert "Hoyeraal" in content, "Missing Hoyeraal-Hreidarsson syndrome reference"
    assert "lag" in content.lower(), "Missing kinetic lag paradox discussion"
    assert "stem cell" in content.lower(), "Missing stem cell analysis"


def test_falsifiable_predictions_report():
    """Verify falsifiable_predictions.md specifies required experimental protocols."""
    doc_path = RESULTS_DIR / "validation_experiments" / "falsifiable_predictions.md"
    assert doc_path.is_file(), f"Missing {doc_path}"
    assert doc_path.stat().st_size > 3000, "falsifiable_predictions.md is too short"

    content = doc_path.read_text(encoding="utf-8")
    assert "stopped-flow" in content.lower() or "spr" in content.lower(), "Missing stopped-flow / SPR protocol"
    assert "hdx-ms" in content.lower() or "hydrogen-deuterium" in content.lower(), "Missing HDX-MS protocol"
    assert "fret" in content.lower(), "Missing single-molecule FRET protocol"
    assert "trap" in content.lower() or "rap" in content.lower(), "Missing processivity / TRAP protocol"
    assert "falsification" in content.lower(), "Missing quantitative falsification criteria"


def test_final_adversarial_review():
    """Verify final_adversarial_review.md red-team audit conclusions."""
    doc_path = RESULTS_DIR / "final_adversarial_review.md"
    assert doc_path.is_file(), f"Missing {doc_path}"
    assert doc_path.stat().st_size > 3000, "final_adversarial_review.md is too short"

    content = doc_path.read_text(encoding="utf-8")
    assert "CLM-021" in content or "two-step" in content.lower()
    assert "CLM-022" in content or "allosteric" in content.lower()
    assert "TELIN-1" in content or "telin1" in content.lower()
    assert "quarantine" in content.lower(), "TELIN-1 must be marked as quarantined"


def test_rewritten_report_13_sections_and_standards():
    """Verify report.md has all 13 mandatory sections and passes tier 4 standards."""
    assert REPORT_PATH.is_file(), f"Missing {REPORT_PATH}"
    assert REPORT_PATH.stat().st_size > 10000, f"report.md unexpectedly small ({REPORT_PATH.stat().st_size} bytes)"

    content = REPORT_PATH.read_text(encoding="utf-8")

    # Verify 13 mandatory numbered sections
    mandatory_sections = [
        "1. Executive Summary",
        "2. Repository Inventory",
        "3. Prior-Art & Novelty Analysis",
        "4. Falsifiable Hypotheses",
        "5. Primary Structural Data",
        "6. Quantitative Interface Geometry",
        "7. Parameter Sensitivity",
        "8. Cross-State Comparative Analysis",
        "9. Mechanistic Evaluation of Allostery",
        "10. Candidate Ligand Reality Check",
        "11. Translation: The Cancer vs. Healthy Aging Dilemma",
        "12. Falsifiable Wet-Lab Experimental Protocols",
        "13. Explicit Catalog of Overstated / Rejected Claims",
    ]
    for section in mandatory_sections:
        assert section in content, f"Missing mandatory section header: {section}"

    # Verify Tier 4 required keywords
    assert "AXIS A" in content or "Axis A" in content, "Missing Axis A"
    assert "AXIS B" in content or "Axis B" in content, "Missing Axis B"
    assert "7TRE" in content and "7TRD" in content and "7TRF" in content and "2I46" in content
    assert "Dyskeratosis Congenita" in content or "dyskeratosis congenita" in content.lower()
    assert "stem cell" in content.lower()
