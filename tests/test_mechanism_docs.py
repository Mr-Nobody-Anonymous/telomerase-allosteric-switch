"""Guards for the mechanism-assessment deliverables (results/mechanism, WP1-WP5).

Documentation drifts faster than analysis code. These tests pin the things the assessment
must keep saying, and the two things it must keep *not* saying: that allostery is absent, and
that 756-775 is one continuous helix.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MECH = ROOT / "results" / "mechanism"
README = ROOT / "README.md"

DOCS = {
    "structural_validation.md": ["## 1. Chain identity", "## 3. What does TPP1 actually touch",
                                 "## 5. What secondary structure", "## 7. Can TPP1 reach"],
    "experimental_evidence.md": ["## 2. The measurements", "## 3. What is now established",
                                 "## 5. The experiment that does not exist yet"],
    "model_comparison.md": ["## 3. Predictions", "## 4. What the repository's own mechanical analysis"],
    "experimental_decision_matrix.md": ["## 1. Ranking", "H762E is now a required negative control"],
    "MECHANISM_ASSESSMENT.md": ["## 1. Findings ranked by confidence", "## 5. Claims rejected",
                                "## 6. Standing hypothesis"],
}


@pytest.mark.parametrize("name,sections", sorted(DOCS.items()))
def test_deliverable_exists_with_required_sections(name, sections):
    path = MECH / name
    assert path.exists(), f"missing deliverable {name}"
    text = path.read_text()
    for section in sections:
        assert section in text, f"{name} is missing section {section!r}"


def test_models_are_scored_not_asserted():
    text = (MECH / "model_comparison.md").read_text()
    # each model must be stated up front so that it can fail
    for model in ("Model A", "Model B", "Model C"):
        assert f"**{model}" in text
    # allostery must not be declared absent
    assert "neither supported nor excluded" in text


def test_assessment_forbids_the_two_overstatements():
    text = (MECH / "MECHANISM_ASSESSMENT.md").read_text()
    assert "inconclusive" in text
    assert "not evidence against allostery" in text
    # the retraction must be explicit, not implied
    assert "strand" in text and "turn" in text and "loop" in text


def test_experimental_evidence_is_distinguished_from_simulation():
    text = (MECH / "experimental_evidence.md").read_text()
    assert "no" in text.lower() and "digitised" in text
    assert "H762E" in text and "86%" in text
    assert "rescued" in text.lower() or "rescue" in text.lower()
    assert "independent transfection" in text


def test_older_documents_carry_the_correction_banner():
    for rel in ("results/allostery/ALLOSTERIC_COUPLING_ANALYSIS.md",
                "results/verification/corrected_breakthrough_claim.md"):
        text = (ROOT / rel).read_text()
        assert "Correction (2026-10-31)" in text, rel
        assert "results/mechanism/structural_validation.md" in text, rel


def test_readme_points_at_the_assessment_and_separates_the_four_framings():
    text = README.read_text()
    for name in DOCS:
        assert name in text, name
    for phrase in ("bipartite clamp", "structural audit", "conduit hypothesis",
                   "inconclusive, not evidence against allostery"):
        assert phrase in text.lower() or phrase in text, phrase


def test_readme_test_badge_matches_the_suite():
    proc = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q",
                           "-p", "no:cacheprovider"],
                          cwd=ROOT, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        pytest.skip("collection failed in this environment; badge check not applicable")
    m = re.search(r"(\d+) tests? collected", proc.stdout)
    assert m, proc.stdout[-500:]
    collected = int(m.group(1))
    badge = re.search(r"pytest-(\d+)%20passed", README.read_text())
    assert badge, "README test badge is missing or has an unexpected format"
    assert int(badge.group(1)) == collected, (
        f"README advertises {badge.group(1)} passing tests but {collected} are collected; "
        "update the badge when the suite changes")
