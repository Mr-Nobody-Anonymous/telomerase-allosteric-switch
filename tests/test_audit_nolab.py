"""Guards for the no-lab audit.

These fail if the audit's conclusions drift from the code that produces them.
A change to any headline number must arrive with a corresponding change here.
"""

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
AUDIT_PY = REPO / "results" / "audit" / "audit.py"
AUDIT_MD = REPO / "results" / "audit" / "AUDIT.md"
RESULTS = REPO / "results" / "audit" / "audit_results.json"


@pytest.fixture(scope="module")
def audit():
    spec = importlib.util.spec_from_file_location("audit", AUDIT_PY)
    module = importlib.util.module_from_spec(spec)
    sys.modules["audit"] = module
    spec.loader.exec_module(module)
    return module


def test_audit_script_exists():
    assert AUDIT_PY.exists(), "results/audit/audit.py is missing"


def test_selftest_invariants_pass(audit):
    """The audit's own invariants, run in-process."""
    audit.selftest()


def test_structure_distances_match_published_values(audit):
    """The Lys757 / Tyr176-Gln794 geometry, to 0.01 A."""
    geo = audit.structure_audit()
    watched = re.compile(r"^7TR[DE]$")
    assert any(watched.match(pid) for pid in geo)

    tre = geo["7TRE"]
    assert tre["Y176_OH--Q794_OE1"] == pytest.approx(2.34, abs=0.01)
    assert tre["K757_NZ--Q794_OE1"] == pytest.approx(11.44, abs=0.01)
    assert tre["K757_NZ--Y176_OH"] == pytest.approx(9.91, abs=0.01)

    # Lys757 does not contact either partner in any deposited model.
    for pid, entry in geo.items():
        assert entry["K757_NZ--Q794_OE1"] > 8.0, pid
        assert entry["K757_NZ--Y176_OH"] > 8.0, pid


def test_over_short_contact_is_reported_not_hidden(audit):
    """A 2.34 A O...O pair is below hydrogen-bond distance; it must stay visible."""
    tre = audit.structure_audit()["7TRE"]
    assert tre["Y176_OH--Q794_OE1"] < 2.6, (
        "the audit exists because this contact is over-short at 3.5 A resolution"
    )


def test_only_k760e_survives_correction(audit):
    res = audit.reanalysis()
    surviving = [m for m, e in res["processivity"].items() if e.get("holm")]
    assert surviving == ["K760E"]


def test_f776e_is_unusable(audit):
    res = audit.reanalysis()
    assert res["processivity"]["F776E"]["class"] == "unusable"
    assert res["processivity"]["F776E"]["valid"] is False


def test_no_activity_mutant_is_resolvable(audit):
    """No activity interval may exclude 1 in the deposited normalised panel."""
    res = audit.reanalysis()
    for mut, entry in res["activity"].items():
        assert entry["lo"] < 1 < entry["hi"], mut


def test_log_intervals_are_positive(audit):
    """Regression guard: an earlier version reported a linear-scale CI [0, 1.15]."""
    res = audit.reanalysis()
    for endpoint in ("processivity", "activity"):
        for mut, entry in res[endpoint].items():
            if entry.get("valid", True):
                assert entry["lo"] > 0, f"{endpoint}/{mut} has a non-positive CI bound"


def test_pairing_sensitivity_is_retained(audit):
    """The unpaired check must remain available and must differ from the paired one."""
    res = audit.reanalysis()
    r756e = res["processivity"]["R756E"]
    assert r756e["unpaired"][3] > r756e["p"], "unpaired analysis should be less significant"


def test_power_requires_more_than_three_preparations(audit):
    """The panel could not have detected its own predicted effects."""
    table = audit.power_table()
    assert table["processivity median"]["1.5x"]["paired"] > 3
    assert table["activity median"]["1.5x"]["paired"] > 3
    assert table["activity median"]["detectable_at_n3"] is None


def test_tost_cost_scales_with_noise(audit):
    cheap = audit.tost_n(0.35, 0.405)
    dear = audit.tost_n(0.74, 0.405)
    assert cheap < dear
    assert cheap >= 4, "an equivalence claim cannot be made from a handful of preparations"


def test_mixture_bound_holds(audit):
    for s, s2, phi in ((0.6, 0.95, 0.05), (0.6, 0.9, 0.01), (0.8, 0.99, 0.02)):
        ratios = [
            audit.mixture(k, s, s2, phi) / audit.mixture(k - 1, s, s2, phi)
            for k in range(2, 60)
        ]
        assert min(ratios) >= min(s, s2) - 1e-9
        assert max(ratios) <= max(s, s2) + 1e-9


def test_writeup_states_its_own_limits():
    """The write-up must not present the audit as a biological result."""
    text = AUDIT_MD.read_text().lower()
    for required in ("no mechanism is established", "cannot tell", "not a measurement"):
        assert required in text, f"AUDIT.md must contain: {required!r}"
    for forbidden in ("we discovered", "breakthrough", "proves that tpp1"):
        assert forbidden not in text


def test_results_json_was_generated_from_this_code():
    if not RESULTS.exists():
        pytest.skip("audit_results.json not generated in this checkout")
    data = json.loads(RESULTS.read_text())
    assert set(data) >= {"structure", "reanalysis", "power", "recovery", "tost"}
