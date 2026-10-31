"""Guards for the affinity-processivity preregistration and its data provenance.

The point of these tests is that a preregistration is worthless if it can drift silently. They
pin (i) that the transcribed source data really are the published numbers, (ii) that nobody has
edited the provenance files since the hashes were recorded, (iii) that the preregistration still
contains the inferential rules that constrain what the investigation is allowed to conclude, and
(iv) that the quantities quoted in the preregistration are reproducible from the data.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external" / "tesmer2023"
CAUSALITY = ROOT / "results" / "causality"
PREREG = CAUSALITY / "PREREGISTRATION.md"


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, comment="#")


@pytest.fixture(scope="module")
def processivity() -> pd.DataFrame:
    return read(EXT / "fig3b_processivity_plusPT.csv")


@pytest.fixture(scope="module")
def activity() -> pd.DataFrame:
    return read(EXT / "fig3c_activity_plusPT.csv")


# ---------------------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------------------
def test_source_files_exist_and_carry_provenance_headers():
    for name in ("fig3b_processivity_plusPT.csv", "fig3c_activity_plusPT.csv",
                 "fig3e_motif3n_activity.csv"):
        text = (EXT / name).read_text()
        assert text.startswith("#"), f"{name} lost its provenance header"
        assert "LSA-2022-01727" in text, name
        assert "NOT technical replicates" in text, name


def test_manifest_hashes_match_the_files_on_disk():
    manifest = read(EXT / "SOURCE_MANIFEST.csv")
    assert len(manifest) == 3
    for _, row in manifest.iterrows():
        path = EXT / row["file"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == row["sha256_local"], (
            f"{row['file']} was edited after its hash was recorded; either revert it or update "
            "SOURCE_MANIFEST.csv and the corresponding entry in PREREGISTRATION.md section 3")


def test_upstream_hash_is_marked_unavailable_rather_than_blank():
    manifest = read(EXT / "SOURCE_MANIFEST.csv")
    assert set(manifest["sha256_upstream"]) == {"TESMER_FIG3B_PT", "TESMER_FIG3C_PT", "TESMER_FIG3E"}
    assert manifest["upstream_retrievable"].str.startswith("NO").all()


def test_replicate_semantics_are_three_independent_experiments(processivity, activity):
    for df in (processivity, activity):
        assert list(df["experiment"]) == [1, 2, 3]
        assert df.shape[1] == 11  # experiment + 10 genotypes, no technical-replicate axis


# ---------------------------------------------------------------------------------------
# the transcription is the published data
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("column,expected_pct", [("H762E", 86), ("D788R", 64)])
def test_transcription_reproduces_published_processivity(processivity, column, expected_pct):
    """Tesmer 2023 text states 86% and 64% of WT processivity for these two mutants."""
    fraction = processivity[column].mean() / processivity["WT"].mean()
    assert round(fraction * 100) == expected_pct, (
        f"{column} processivity = {fraction:.3f}; the paper states {expected_pct}% - the "
        "transcription is wrong or the source file changed")


@pytest.mark.parametrize("column,expected_pct", [("H762E", 68), ("D788R", 24)])
def test_transcription_reproduces_published_activity(activity, column, expected_pct):
    """One percentage point of tolerance.

    The mean of the D788R rows is 23.5% while the text says 24%, so exact equality would fail on
    a value that is correct. The discrepancy is recorded in PREREGISTRATION.md section 3 rather
    than tuned away in the test.
    """
    fraction = activity[column].mean() * 100
    assert abs(fraction - expected_pct) <= 1.0, (
        f"{column} activity = {fraction:.1f}%; the paper states {expected_pct}%")


def test_normalisation_is_within_experiment(activity):
    """Fig 3C is normalised to WT inside each experiment; WT must be exactly 1 in every row."""
    assert np.allclose(activity["WT"].to_numpy(), 1.0)


def test_negative_impossible_value_is_retained_verbatim_not_clipped(processivity):
    """F776E replicate 3 is negative; clipping it would fabricate a measurement."""
    assert processivity.loc[2, "F776E"] < 0
    assert "F776E" in (EXT / "README.md").read_text()


# ---------------------------------------------------------------------------------------
# the precision claims in the preregistration are reproducible
# ---------------------------------------------------------------------------------------
def _detectable_ratio(cv: float, n: int = 3) -> float:
    """Two-sided alpha=0.05, 80% power, per-group n."""
    return math.exp((1.96 + 0.8416) * math.sqrt(2.0 / n) * cv)


def test_quoted_precision_limits_match_the_data(processivity, activity):
    # WT is excluded (activity is normalised to it, so its CV is 0 by construction) and so is
    # F776E (near-zero signal; see PREREGISTRATION.md section 3).
    drop = ["WT", "F776E"]
    cv_p = (processivity.drop(columns=drop).std(ddof=1) / processivity.drop(columns=drop).mean())
    cv_a = (activity.drop(columns=drop).std(ddof=1) / activity.drop(columns=drop).mean())
    assert round(float(np.median(cv_p)) * 100) == 35   # prereg section 4.3 says 34.5%
    assert round(float(np.median(cv_a)) * 100) == 90   # prereg section 4.3 says 90.5%
    assert 2.0 < _detectable_ratio(float(np.median(cv_p))) < 2.4   # ~2.2x
    assert 7.4 < _detectable_ratio(float(np.median(cv_a))) < 8.4   # ~7.9x


# ---------------------------------------------------------------------------------------
# the preregistration still says what it must say
# ---------------------------------------------------------------------------------------
REQUIRED = [
    "## 1. The single question",
    "## 2. Frozen directional predictions",
    "## 4. Inferential rules",
    "**4.4 Non-identifiability rule (binding).**",
    "**4.6 A geometric distribution is not evidence for affinity.**",
    "**4.2 Two response variables, both required.**",
    "## 5. Feasibility verdict (gate G0)",
    "## 8. What would change the plan",
]


def _flat(text: str) -> str:
    """Collapse whitespace so that markdown line-wrapping cannot break a content check."""
    return " ".join(text.split())


@pytest.mark.parametrize("section", REQUIRED)
def test_preregistration_contains_required_rule(section):
    assert _flat(section) in _flat(PREREG.read_text()), f"preregistration lost: {section}"


def test_preregistration_keeps_the_three_acceptable_outcomes():
    text = _flat(PREREG.read_text())
    for outcome in ("affinity sufficient", "affinity insufficient",
                    "not identifiable with available data"):
        assert outcome in text, outcome


def test_preregistration_states_that_delta_rap_is_not_in_the_source_data():
    text = _flat(PREREG.read_text())
    assert "cannot be computed from them at all" in text
    assert "ΔRAP **cannot be computed from them at all**" in text


def test_preregistration_does_not_declare_a_mechanism():
    text = PREREG.read_text().lower()
    for forbidden in ("proves affinity", "affinity is the mechanism",
                      "allostery is excluded", "allostery is absent"):
        assert forbidden not in text, f"preregistration asserts {forbidden!r}"


def test_frozen_hash_or_amendment():
    """PREREGISTRATION.md is frozen by hash; changes require an amendments file."""
    frozen = (CAUSALITY / "PREREGISTRATION.sha256")
    assert frozen.exists(), "the frozen hash file is missing"
    recorded = frozen.read_text().split()[0]
    actual = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    if actual != recorded:
        amendments = CAUSALITY / "PREREGISTRATION_AMENDMENTS.md"
        assert amendments.exists(), (
            "PREREGISTRATION.md was edited after it was frozen and no "
            "PREREGISTRATION_AMENDMENTS.md exists. Revert the edit, or record why it changed.")


# ---------------------------------------------------------------------------------------
# housekeeping that keeps the deliverable directory coherent
# ---------------------------------------------------------------------------------------
def test_legacy_documents_moved_but_preserved():
    legacy = ROOT / "results" / "mechanism"
    for name in ("allostery_evidence_assessment.md", "de_loop_mechanism_report.md"):
        assert not (legacy / name).exists(), f"{name} should live under legacy/"
        assert (legacy / "legacy" / name).exists(), f"{name} was deleted, not relocated"
        assert (legacy / "legacy" / name).stat().st_size > 1000


def test_current_assessment_is_not_mixed_with_legacy():
    current = {p.name for p in (ROOT / "results" / "mechanism").glob("*.md")}
    assert {"MECHANISM_ASSESSMENT.md", "structural_validation.md"} <= current
    assert "de_loop_mechanism_report.md" not in current
