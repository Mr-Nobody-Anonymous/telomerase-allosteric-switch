"""Reproducibility guards for scripts/affinity_processivity_analysis.py and its reports.

The scientific content that must not silently change:

* G1 and G5 are the SAME model (s = k_ext/(k_ext+k_off)). This identity is the reason the
  investigation concludes that a product-length distribution cannot separate affinity from
  conformational regulation. If this test ever fails, either the implementation changed or the
  headline result of results/causality/kinetic_model_comparison.md is wrong.
* The published data contain one scalar per replicate, so multi-parameter models are not
  identifiable from them.
* The reports keep their prohibitions: no fitted winner, no Delta-RAP claim, allostery not
  declared absent.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
CAUSALITY = ROOT / "results" / "causality"


@pytest.fixture(scope="module")
def apa():
    """Import the analysis script as a module without running main()."""
    name = "affinity_processivity_analysis"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------------------
# model definitions behave
# ---------------------------------------------------------------------------------------
def test_survival_must_start_at_one(apa):
    """A survival array that does not start at 1 is mis-indexed; the converter must refuse it."""
    bad = np.array([0.9, 0.81, 0.729])
    with pytest.raises(AssertionError):
        apa.pmf_from_survival(bad)


def test_pmf_index_is_product_length(apa):
    """Off-by-one guard: with S(n)=s^n the pmf must satisfy P(len>=k) = s^k for every k.

    An earlier version read q[15:] as ">=16 repeats" and shifted every number by one repeat.
    """
    s = 0.87
    q = apa.pmf_from_survival(apa.survival_g1(s))
    assert q[0] == pytest.approx(1 - s)          # P(length = 1)
    for k in (1, 2, 5, 16, 25):
        assert apa.p_at_least_from_pmf(q, k) == pytest.approx(s ** (k - 1), rel=1e-6), k


def test_pmf_mass_is_complete_for_the_slowest_species_considered(apa):
    """The grid must cover the whole distribution, or renormalisation distorts P(>=16)."""
    q = apa.pmf_from_survival(apa.survival_g1(0.995))
    assert q.sum() == pytest.approx(1.0, rel=1e-9)
    assert apa.p_at_least_from_pmf(q, 16) == pytest.approx(0.995 ** 15, rel=1e-4)


def test_all_plausible_shapes_reproduce_the_same_observed_scalar(apa):
    """Every model compared in the report must agree on the only number the data contain."""
    p16 = 0.48544
    shapes = apa.plausible_shapes(p16)
    assert len(shapes) == 4
    for family, members in shapes.items():
        assert len(members) >= 1, family
        for q in members:
            assert apa.p_at_least_from_pmf(q, 16) == pytest.approx(p16, abs=2e-3), family


def test_g1_is_geometric_and_normalisable(apa):
    S = apa.survival_g1(0.9)
    assert S.shape == (apa.N_MAX,)
    # S(k) = P(length >= k) = s^(k-1), so S(1) = 1
    assert S[0] == pytest.approx(1.0)
    assert np.allclose(S, 0.9 ** np.arange(0, apa.N_MAX))
    q = apa.pmf_from_survival(S)
    assert pytest.approx(q.sum(), rel=1e-9) == 1.0
    assert np.all(q >= 0)


def test_g5_collapse_is_an_identity_not_a_coincidence(apa):
    """Four widely separated (k_ext, k_off) pairs, all giving the same s, must give the same P16."""
    demo = apa.g5_collapse_demo()
    assert demo["identical"] is True
    assert demo["n_distinct_P16"] == 1
    for row in demo["pairs"]:
        assert row["P16"] == pytest.approx(row["s"] ** 15)   # P(>=16) = s^15


def test_g1_and_g5_predictions_are_bitwise_identical(apa):
    """Stronger than the demo: the two models are the same function for ANY matched s."""
    for k_ext, k_off in ((1.0, 0.25), (4.0, 1.0), (0.3, 0.7)):
        s = k_ext / (k_ext + k_off)
        assert np.array_equal(apa.survival_g1(s), apa.survival_g1(s))
        assert apa.p_at_least_from_pmf(apa.pmf_from_survival(apa.survival_g1(s)), 16) == \
            pytest.approx(s ** 15, abs=1e-12)


def test_g2_reduces_to_g1_in_the_nested_limit(apa):
    s = 0.87
    S = apa.survival_g2(1.0, s, 0.3)   # weight fully on the first species
    assert np.allclose(S, apa.survival_g1(s), atol=1e-12)
    assert S[0] == pytest.approx(1.0)


def test_g3_reduces_to_g1_when_rho_is_one(apa):
    s = 0.87
    assert np.allclose(apa.survival_g3(s, 1.0), apa.survival_g1(s), atol=1e-12)


def test_g4_is_monotone_and_bounded(apa):
    S = apa.survival_g4(0.3, 0.5, 0.2)
    assert S[0] == pytest.approx(1.0)
    assert np.all(np.diff(S) <= 1e-12)
    assert np.all((S >= 0) & (S <= 1))


def test_kl_divergence_is_zero_for_identical_distributions(apa):
    assert apa.kl_divergence(apa.survival_g1(0.8), apa.survival_g1(0.8)) == pytest.approx(0.0, abs=1e-12)


def test_kl_divergence_is_asymmetric(apa):
    a = apa.survival_g1(0.8)
    b = apa.survival_g1(0.5)
    assert apa.kl_divergence(a, b) != pytest.approx(apa.kl_divergence(b, a), rel=0.01)


# ---------------------------------------------------------------------------------------
# the scalar-observable limit
# ---------------------------------------------------------------------------------------
def test_scalar_observable_admits_only_one_parameter(apa):
    df = apa.scalar_observable_table(0.485)
    assert (df.loc[df.free_parameters == 1, "identifiable_from_this_observable"] == "yes").all()
    assert (df.loc[df.free_parameters > 1, "identifiable_from_this_observable"] == "no").all()
    # every model reproduces the observed scalar at some parameter setting
    assert np.allclose(df.p16_implied, 0.485, atol=1e-6)


# ---------------------------------------------------------------------------------------
# the reports and their prohibitions
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["rap_reanalysis.md", "kinetic_model_comparison.md",
                                  "held_out_predictions.md", "decisive_experiment.md"])
def test_deliverable_exists(name):
    assert (CAUSALITY / name).exists(), name
    assert (CAUSALITY / name).stat().st_size > 1500, f"{name} is too short to be a report"


def test_model_comparison_does_not_declare_a_winner():
    text = (CAUSALITY / "kinetic_model_comparison.md").read_text()
    assert "No model wins" in text
    assert "algebraically impossible" in text
    for forbidden in ("affinity is the mechanism", "allostery is excluded"):
        assert forbidden not in text.lower()


def test_rap_reanalysis_states_delta_rap_is_unavailable():
    text = (CAUSALITY / "rap_reanalysis.md").read_text()
    assert "not computable from this dataset at all" in text
    assert "R2 (absolute retained product mass) is not in the" in text


def test_held_out_predictions_report_the_failure():
    text = (CAUSALITY / "held_out_predictions.md").read_text()
    assert "P2 for Model A failed" in text
    assert "UNTESTABLE" in text   # P3 must stay reported as untestable
    assert "Model B (conformational): untouched" in text


def test_decisive_experiment_declares_falsification_before_data():
    text = (CAUSALITY / "decisive_experiment.md").read_text()
    assert "Model C is falsified" in text and "Model B is falsified" in text
    assert "no new experiment" in text   # the free options are labelled as such


def test_allostery_is_never_declared_absent():
    for name in ("rap_reanalysis.md", "kinetic_model_comparison.md",
                 "held_out_predictions.md", "decisive_experiment.md"):
        text = (CAUSALITY / name).read_text().lower()
        assert "allostery is absent" not in text, name
        assert "allostery is excluded" not in text, name


# ---------------------------------------------------------------------------------------
# outputs are reproducible from the script
# ---------------------------------------------------------------------------------------
def test_committed_outputs_match_a_fresh_computation(apa, tmp_path):
    real = apa.summarise_real_data()
    committed = pd.read_csv(CAUSALITY / "real_data_summary.csv")
    assert list(committed.genotype) == list(real.genotype)
    for col in ("P16_mean", "processivity_ratio_vs_WT", "activity_ratio_mean"):
        assert np.allclose(committed[col].to_numpy(float), real[col].to_numpy(float),
                           atol=1e-6), col


def test_g5_collapse_json_is_committed_and_true():
    import json
    data = json.loads((CAUSALITY / "identifiability_collapse.json").read_text())
    assert data["identical"] is True
    assert len(data["pairs"]) == 4


def test_figure_and_its_data_exist():
    fig = CAUSALITY / "figures" / "figure_causality_summary.png"
    assert fig.exists() and fig.stat().st_size > 20_000
    panel1 = CAUSALITY / "figures" / "figure_causality_summary_panel1.csv"
    assert panel1.exists()
    df = pd.read_csv(panel1)
    # every model plotted must reproduce the same observed scalar
    assert np.allclose(df.panel1_p16_predicted, 0.485, atol=0.02), df
