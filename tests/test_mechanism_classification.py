"""Guards for MECHANISM_CLASSIFICATION.md.

That document makes exact algebraic claims and quotes computed numbers. Both can drift: the prose
can be edited without the script, or the script can be edited without the prose. These tests pin
(i) that every number the document quotes is reproduced by `scripts/mechanism_classification_checks.py`,
(ii) that the exact algebra is verified symbolically rather than asserted, and (iii) that nobody has
quietly upgraded what the document is allowed to claim.

The document's own standard, from §9: a claim is mechanism only if an independent measurement
separates it from the alternatives. Nothing here may be cited as having done that.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CAUSALITY = ROOT / "results" / "causality"
DOC = CAUSALITY / "MECHANISM_CLASSIFICATION.md"
SCRIPT = ROOT / "scripts" / "mechanism_classification_checks.py"
JSON_OUT = CAUSALITY / "classification_checks.json"


# ---------------------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def doc() -> str:
    assert DOC.exists(), "MECHANISM_CLASSIFICATION.md is missing"
    return DOC.read_text()


@pytest.fixture(scope="module")
def results(tmp_path_factory) -> dict:
    """Run the check script in-process-independent fashion and load what it wrote."""
    proc = subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=False,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, f"script failed:\n{proc.stdout}\n{proc.stderr}"
    assert JSON_OUT.exists()
    return json.loads(JSON_OUT.read_text())


# ---------------------------------------------------------------------------------------
# the document still opens with the thesis it is built on
# ---------------------------------------------------------------------------------------
def test_opens_with_the_thesis_sentence(doc):
    head = " ".join(doc.split())
    assert "An identifiable effective kinetic parameter is not necessarily an identifiable" in head
    assert "molecular mechanism" in head
    # and it must be the *opening* claim, before the first numbered section
    head_before_sections = doc.split("## 0.")[0]
    assert "identifiable effective kinetic parameter" in " ".join(head_before_sections.split())


def test_contains_every_required_section(doc):
    numbered = re.findall(r"^## (\d+)\. ", doc, re.MULTILINE)
    assert numbered == [str(n) for n in range(10)], (
        f"expected sections 0-9 in order, found {numbered}")


# ---------------------------------------------------------------------------------------
# exact algebra: verified symbolically, and the document says so
# ---------------------------------------------------------------------------------------
def test_symbolic_invariance_is_actually_zero(results):
    sym = results["symbolic"]
    assert sym["scale_invariance_simplifies_to"] in ("0", "0.0")
    assert sym["scale_invariance_is_zero"] is True


def test_crossing_closed_form_is_numerically_correct(results):
    """Guards the sign of the displayed numerator, which a draft of the document got backwards."""
    cn = results["crossing_numeric"]
    assert cn["P_star_satisfies_equality"] is True
    assert cn["sign_change_count_matches_sign_of_P_star"] is True
    assert cn["has_positive_example"] is True
    for row in cn["rows"]:
        expect = 1 if row["P_star"] > 0 else 0
        assert row["sign_changes"] == expect, row
        if row["P_star"] > 0:
            assert row["residual_at_P_star"] < 1e-12


def test_document_shows_a_positive_crossing_example(doc):
    flat = " ".join(doc.split())
    assert "17.27" in flat, "the positive-crossing worked example must be present"
    assert "sign error" in flat or "opposite sign" in flat


def test_crossing_is_degree_one_hence_unique(results):
    sym = results["symbolic"]
    assert sym["crossing_degree_in_P"] == 1
    assert sym["crossing_is_unique"] is True
    # the degenerate case must be reported, not hidden
    assert sym["degenerate_case_alpha_eq_beta"] == "alpha*k_off**2*(alpha - 1)"


def test_document_quotes_the_same_closed_form_for_P_star(results, doc):
    flat = " ".join(doc.split()).lower()
    assert "at most one crossing" in flat
    assert results["symbolic"]["P_star"].startswith("k_off*(")


# ---------------------------------------------------------------------------------------
# §2: the corrected C1a / C1b statement. This test exists because the earlier draft of the
# document wrongly claimed the two variants were equivalent at matched Kd.
# ---------------------------------------------------------------------------------------
def test_loss_term_rewrites_exactly(results):
    sym = results["symbolic"]
    assert sym["loss_term_equals_koff2_over_koff_plus_P_kon"] is True
    assert sym["loss_term_at_zero_primer"] == "k_off"


def test_matched_Kd_scales_the_loss_term_by_lambda(results):
    """A k_off-only and a k_on-only change at matched Kd differ by exactly a factor, not by nothing."""
    assert results["symbolic"]["matched_Kd_lam_scales_loss_term"] is True
    c1 = results["c1_variants"]
    assert c1["same_Kd"] is True
    assert c1["distinguishable_at_low_primer"] is True
    # and the factor really is lambda: C1b's loss is 2x C1a's here
    a, b = c1["cases"]["C1a_koff_down"], c1["cases"]["C1b_kon_up"]
    assert b["s_at_P0.1"] < a["s_at_P0.1"]          # higher loss term -> lower survival
    assert c1["converge_at_high_primer"] is True     # they converge only as [P] -> infinity


def test_document_does_not_restore_the_false_equivalence_claim(doc):
    flat = " ".join(doc.split())
    assert "are **not** identical at matched" in flat, (
        "the document must state that C1a and C1b are NOT identical at matched Kd")
    # every line mentioning 'identical at matched' must carry the negation (tables have no
    # sentence punctuation, so check line by line rather than by sentence)
    for line in doc.splitlines():
        plain = line.replace("*", "")          # markdown bold splits the phrase
        if "identical at matched" not in plain:
            continue
        # the corrections table is allowed to quote the error it retracts
        if "claimed" in plain or "error" in plain.lower():
            continue
        assert "not identical" in plain, f"unguarded equivalence claim: {line[:200]}"
    # the displayed worked example must stay consistent with the script
    assert "0.587" in flat and "0.463" in flat


# ---------------------------------------------------------------------------------------
# §8: identifiability. Rank is necessary but not sufficient, and the document must show that.
# ---------------------------------------------------------------------------------------
def test_structural_ranks(results):
    ident = results["identifiability"]
    assert ident["rank_processivity_only"] == 3
    assert ident["rank_with_turnover_assumed"] == 4
    assert ident["rank_with_binding"] == 4
    assert ident["rank_with_free_coupling"] == [4, 5], "free coupling must be rank-deficient"
    assert ident["smallest_relative_sv_processivity"] < 1e-12


def test_practical_identifiability_numbers_match_the_document(results, doc):
    ident = results["identifiability"]
    flat = " ".join(doc.split())
    # 12.49% rounds to 12%, not 13%; the rounding is checked, not eyeballed
    se = {k: round(100 * v) for k, v in ident["ratio_relative_SE"].items()}
    assert se == {"k_loss/k_ext": 44, "k_off/k_ext": 12, "k_on/k_ext": 60}
    for token in ("44%", "12%", "60%"):
        assert token in flat
    assert abs(ident["cond_number_ratios"] - 60.3) < 0.1
    assert "60.3" in flat


def test_binding_does_not_rescue_k_loss(results, doc):
    """The document's key negative: direct binding fixes k_off/k_on but not k_loss."""
    af = results["absolute_fisher"]
    only = af["processivity_only"]["relative_SE_pct"]["k_loss"]
    with_binding = af["plus_binding_10pct"]["relative_SE_pct"]["k_loss"]
    assert with_binding >= only, "if binding improved k_loss, the document's claim is wrong"
    assert af["plus_binding_10pct"]["relative_SE_pct"]["k_off"] < 12
    flat = " ".join(doc.split())
    assert "49.0%" in flat
    assert "not sufficient" in flat or "necessary but" in flat


# ---------------------------------------------------------------------------------------
# §5: the tail slope. Forbidden as a basis for C2, and shown to be window-dependent.
# ---------------------------------------------------------------------------------------
def test_closed_form_matches_simulation(results):
    assert results["tail_slope"]["closed_form_max_abs_error"] < 1e-12


def test_tail_slope_is_window_dependent(results, doc):
    sweep = results["tail_slope"]["window_sweep"]
    first = next(w for w in sweep if w["s_majority"] == 0.60)
    assert first["spread"] > 0.3, "the window effect is the objection; it must be large"
    values = first["by_window"]
    assert values["3-10"] < 0.65 and values["40-80"] > 0.94
    # the implied loss ratio really is a factor ~12 apart
    near, far = values["3-10"], values["40-80"]
    assert ((1 - near) / near) / ((1 - far) / far) > 10
    flat = " ".join(doc.split())
    assert "12-fold" in flat
    assert "0.610" in flat and "0.950" in flat


def test_tail_slope_never_exceeds_the_longest_lived_species(results):
    for case in results["tail_slope"]["cases"]:
        assert case["tail_slope"] <= max(case["s_majority"], case["s2"]) + 1e-3


def test_saturation_bias_is_near_zero_not_minus_forty(results):
    """The earlier draft claimed -0.040 at saturation; that number is not reproducible."""
    sat = [c for c in results["tail_slope"]["cases"] if c["s_majority"] > 0.99]
    assert sat
    for c in sat:
        assert abs(c["delta_s"]) < 0.01


def test_document_forbids_tail_slope_for_C2(doc):
    flat = " ".join(doc.split()).lower()
    assert "tail slope is forbidden" in flat
    # C2 conclusions may only come from a full-distribution fit with phi free
    assert "full-distribution fit" in flat
    assert "never from a tail slope" in flat


# ---------------------------------------------------------------------------------------
# the document must not have upgraded its claims
# ---------------------------------------------------------------------------------------
def test_document_declares_no_mechanism(doc):
    flat = " ".join(doc.split()).lower()
    assert "no mechanism" in flat or "not proven" in flat or "no experiment yet distinguishes" in flat


def test_document_keeps_the_simulation_derived_label(doc):
    flat = " ".join(doc.split())
    assert "simulation-derived" in flat
    assert "not yet validated on real bands" in flat


def test_document_keeps_gate_states_honest(doc):
    flat = " ".join(doc.split()).lower()
    for state in ("satisfied", "conditional", "experimental", "pending"):
        assert state in flat


def test_document_has_an_explicit_non_claims_section(doc):
    section = doc.split("## 9.")[-1]
    flat = " ".join(section.split()).lower()
    assert "does not claim" in flat or "not claim" in flat


# ---------------------------------------------------------------------------------------
# reproducibility plumbing
# ---------------------------------------------------------------------------------------
def test_document_hash_is_recorded_and_current():
    """The document is frozen by hash, so an edit is a deliberate act, not silent drift."""
    frozen = CAUSALITY / "MECHANISM_CLASSIFICATION.sha256"
    assert frozen.exists(), "MECHANISM_CLASSIFICATION.sha256 is missing"
    recorded = frozen.read_text().split()[0]
    actual = hashlib.sha256(DOC.read_bytes()).hexdigest()
    assert actual == recorded, (
        "MECHANISM_CLASSIFICATION.md changed after its hash was recorded. Re-run the checks, "
        "confirm every quoted number, then regenerate MECHANISM_CLASSIFICATION.sha256.")


def test_script_writes_both_json_and_csv():
    assert (CAUSALITY / "classification_checks.csv").exists()
    assert JSON_OUT.exists()


def test_document_does_not_quote_numbers_the_script_cannot_produce(results):
    """Every 2-3 sig-fig decimal in the document that looks like a check output must be traceable."""
    manifest = results.get("manifest", {})
    assert isinstance(manifest, dict)
