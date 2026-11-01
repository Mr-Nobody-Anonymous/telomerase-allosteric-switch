"""Guards for paper/MANUSCRIPT.md.

The manuscript quotes numbers that are computed by
`scripts/mechanism_classification_checks.py`. This module regenerates them and fails if the text
and the computation disagree, and it fails if any retracted claim is reintroduced.

The paper's own standard, stated in its §8: an exact result within a model is not an
experimentally established biological fact. These tests cannot verify biology — only that the
paper does not overstate what it has, and that its arithmetic still matches its code.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "MANUSCRIPT.md"
FIGDIR = ROOT / "paper" / "figures"
CHECKS = ROOT / "scripts" / "mechanism_classification_checks.py"
JSON_OUT = ROOT / "results" / "causality" / "classification_checks.json"


@pytest.fixture(scope="module")
def text() -> str:
    assert PAPER.exists(), "paper/MANUSCRIPT.md is missing"
    return PAPER.read_text()


@pytest.fixture(scope="module")
def flat(text) -> str:
    return " ".join(text.split())


@pytest.fixture(scope="module")
def results() -> dict:
    proc = subprocess.run([sys.executable, str(CHECKS)], cwd=ROOT, check=False,
                          capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return json.loads(JSON_OUT.read_text())


# ---------------------------------------------------------------------------------------
# authorship and structure
# ---------------------------------------------------------------------------------------
def test_author_is_named(text):
    assert "Beamlak Sintayehu" in text


def test_has_abstract_and_the_required_sections(text):
    assert "## Abstract" in text
    for section in ("## 1.", "## 2.", "## 3.", "## 4.", "## 5.", "## 6.", "## 7.", "## 8."):
        assert re.search(rf"^{re.escape(section)} ", text, re.MULTILINE), section
    assert "## Figure legends" in text


def test_abstract_states_no_mechanism_is_demonstrated(text):
    abstract = text.split("## Abstract")[1].split("## 1.")[0]
    assert "No mechanism is demonstrated" in abstract


def test_references_are_complete(text):
    refs = text.split("## References")[1].split("## Figure legends")[0]
    numbered = re.findall(r"^\d+\. ", refs, re.MULTILINE)
    assert len(numbered) == 7, f"expected 7 references, found {len(numbered)}"
    assert refs.count("doi:") == 7, "every reference must carry a verifiable DOI"


# ---------------------------------------------------------------------------------------
# every quoted number must match the computation
# ---------------------------------------------------------------------------------------
def test_ratio_standard_errors_match(results, flat):
    se = results["identifiability"]["ratio_relative_SE"]
    quoted = {"k_loss/k_ext": 44, "k_off/k_ext": 12, "k_on/k_ext": 60}
    for key, expected in quoted.items():
        assert abs(100 * se[key] - expected) < 0.6, (key, 100 * se[key], expected)
        assert f"{expected}%" in flat, f"{expected}% not quoted for {key}"


def test_condition_numbers_match(results, flat):
    ident = results["identifiability"]
    assert abs(ident["cond_number_ratios"] - 60.3) < 0.1
    assert "60.3" in flat
    assert "2.4×10⁸" in flat or "2.4x10⁸" in flat


def test_absolute_fisher_table_matches(results, flat):
    af = results["absolute_fisher"]
    withbind = af["plus_binding_10pct"]["relative_SE_pct"]
    assert round(withbind["k_loss"], 1) == 49.0
    assert round(withbind["k_off"], 1) == 9.9
    only = af["processivity_only"]["relative_SE_pct"]
    assert round(only["k_loss"], 1) == 34.3
    assert "49.0%" in flat and "34.3%" in flat and "9.9%" in flat


def test_ranks_match(results, flat):
    ident = results["identifiability"]
    assert ident["rank_processivity_only"] == 3
    assert ident["rank_with_turnover_assumed"] == 4
    assert ident["rank_with_binding"] == 4
    assert ident["rank_with_free_coupling"] == [4, 5]
    assert "rank **3 of 4**" in flat or "rank 3 of 4" in flat


def test_power_table_matches(results, flat):
    rows = {(r["delta_log"], r["cv"], r["n_tests"]): r["n"] for r in results["power"]["tost_table"]}
    assert rows[(0.1, 0.345, 3)] == 105
    assert rows[(0.1, 0.345, 1)] == 74
    assert rows[(0.1, 0.15, 3)] == 20
    assert rows[(0.2, 0.345, 3)] == 27
    assert "**105**" in flat and "74 (1 test)" in flat


def test_closed_form_agreement_claim(results, flat):
    err = results["tail_slope"]["closed_form_max_abs_error"]
    assert err < 1e-15
    assert "4.4×10⁻¹⁶" in flat


def test_fabricated_c2_range_matches(results, flat):
    implied = [f["implied_k_loss_over_k_ext"] for f in results["tail_slope"]["fabricated_C2"]]
    assert 0.015 < min(implied) < 0.020
    assert 0.030 < max(implied) < 0.040
    assert "0.017–0.036" in flat


def test_onset_length_matches(results, flat):
    row = results["mixture_onset"]["rows"][0]
    assert row["s_majority"] == 0.60 and row["phi"] == 0.01 and row["s2"] == 0.95
    assert abs(row["n_star"] - 15.52) < 0.05
    assert row["midpoint_crossed_at"] == 15
    assert "15.5" in flat


def test_window_sweep_numbers_match(results, flat):
    row = results["tail_slope"]["window_sweep"][0]
    by = row["by_window"]
    assert abs(by["3-10"] - 0.6098) < 1e-3
    assert abs(by["40-80"] - 0.95) < 1e-6
    for token in ("0.610", "0.909", "0.950"):
        assert token in flat
    near, far = (1 - by["3-10"]) / by["3-10"], (1 - by["40-80"]) / by["40-80"]
    assert 11.0 < near / far < 13.0
    assert "12-fold" in flat


def test_interpolation_bound_is_reported(results, flat):
    mb = results["mixture_bound"]
    assert mb["monotone"] is True
    assert mb["bound_worst_relative_violation"] < 1e-10
    assert mb["r_as_u_to_0"] == "s"
    assert "min(s, s₂) ≤ r(n) ≤ max(s, s₂)" in flat


# ---------------------------------------------------------------------------------------
# the paper must not resurrect anything that was retracted
# ---------------------------------------------------------------------------------------
def test_delta_s_approximation_is_retired(flat):
    assert "φ(1−*s*)" in flat
    assert "retired" in flat.lower() or "withdrawn" in flat.lower()
    assert "50–90×" in flat


def test_c1a_c1b_equivalence_is_never_asserted(text):
    """Every mention of 'identical at matched' must be negated or inside the retraction table."""
    for line in text.splitlines():
        plain = line.replace("*", "")
        if "identical at matched" not in plain:
            continue
        if "claimed" in plain or "error" in plain.lower():
            continue
        assert "not identical" in plain, f"unguarded equivalence claim: {line[:160]}"
    flat = " ".join(text.split())
    assert "They are not, and the reason matters" in flat


def test_saturation_minus_0_040_is_withdrawn(flat):
    assert "not reproducible" in flat
    assert "withdrawn here" in flat


def test_b_vs_c1_claim_is_conditional_not_global(text, flat):
    """B and C1 agree only at matched s for a given [P] — with one fixed parameter pair each they
    agree at no [P] at all. An earlier phrasing ('globally identical') overstated this."""
    assert "globally identical" not in flat
    assert "at any given primer concentration" in flat
    # the wording that shows the limitation is present
    assert "agree at *no* primer concentration" in flat
    assert "0.889" in flat and "0.800" in flat


def test_scope_disclaimer_present(flat):
    assert "mathematically exact result within a model is not an experimentally established" in flat
    assert "at any precision" in flat


def test_does_not_claim_a_mechanism(flat):
    for overclaim in ("TPP1 acts by", "demonstrates that TPP1 changes", "we prove that TPP1"):
        assert overclaim not in flat
    assert "We do not claim TPP1 changes primer retention" in flat


def test_next_step_is_experimental(flat):
    assert "The next experiment is not a simulation" in flat
    assert "densitometry" in flat
    assert "falsif" in flat.lower()


# ---------------------------------------------------------------------------------------
# artefacts
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", [
    "fig1_collapse", "fig2_affinity_variants", "fig3_crossing",
    "fig4_mixture_ratio", "fig5_window_dependence", "fig6_identifiability",
])
def test_figures_exist_in_both_formats(name):
    assert (FIGDIR / f"{name}.png").exists(), f"{name}.png missing"
    assert (FIGDIR / f"{name}.pdf").exists(), f"{name}.pdf missing"
    assert (FIGDIR / f"{name}.png").stat().st_size > 20_000, "suspiciously small figure"


def test_every_figure_is_cited_in_the_legends(text):
    legends = text.split("## Figure legends")[1]
    for n in range(1, 7):
        assert f"**Figure {n}." in legends, f"Figure {n} has no legend"


def test_repository_url_is_present(flat):
    assert "https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch" in flat
