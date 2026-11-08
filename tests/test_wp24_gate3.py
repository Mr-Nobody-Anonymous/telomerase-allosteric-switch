"""Guards for WP24: the threshold fraction cannot satisfy Gate 3."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
WP24 = REPO / "scripts" / "solve" / "wp24_gate3_harness.py"


@pytest.fixture(scope="module")
def m():
    spec = importlib.util.spec_from_file_location("wp24", WP24)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_collision_is_exact_at_the_recorded_threshold(m):
    """The two distributions must be indistinguishable to the recorded scalar."""
    a, b, T = m.make_collision()
    assert m.thresh_fraction(a) == pytest.approx(T, abs=1e-12)
    assert m.thresh_fraction(b) == pytest.approx(T, abs=1e-12)
    assert m.thresh_fraction(a) == pytest.approx(m.thresh_fraction(b), abs=1e-12)


def test_collision_actually_differs_in_shape(m):
    """...while differing substantially below the threshold - else it proves nothing."""
    a, b, _ = m.make_collision()
    early_a, early_b = sum(a[:4]), sum(b[:4])
    assert early_b / early_a > 3.0, "the two shapes must be clearly different"


def test_both_distributions_are_valid_probability_vectors(m):
    a, b, _ = m.make_collision()
    for p in (a, b):
        assert all(x >= 0 for x in p)
        assert sum(p) == pytest.approx(1.0, abs=1e-12)


def test_no_per_band_ladder_exists_in_the_repo(m):
    """If someone adds a ladder, this test should be updated deliberately."""
    ladders = []
    for f in sorted(m.TESMER.glob("*.csv")):
        with f.open() as fh:
            rows = [r for r in fh if not r.startswith("#")]
        # a ladder would have one row per repeat number, i.e. many rows
        if len(rows) > 25:
            ladders.append(f.name)
    assert ladders == [], f"per-band data appeared: {ladders} - update Gate 3 status"
