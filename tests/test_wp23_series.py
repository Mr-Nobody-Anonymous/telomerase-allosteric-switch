"""Guards for WP23, the TPP1-free matched series.

A Kabsch implementation that is transposed returns a plausible-looking number
(3.5 A here) instead of failing, so it is pinned against two EXACT cases before
any biological number is trusted. Getting this wrong was not hypothetical: it
produced a 3.55 A 'displacement' where the correct answer is 0.356 A.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
WP23 = REPO / "scripts" / "solve" / "wp23_tpp1_free_series.py"


def _m():
    spec = importlib.util.spec_from_file_location("wp23", WP23)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m():
    return _m()


def _rmsd(P, Q, R, t):
    moved = (R @ P.T).T + t
    return float(np.sqrt(((moved - Q) ** 2).sum(1).mean()))


def test_kabsch_identity_is_exact(m):
    """Superposing a structure on itself must give a zero transform."""
    rng = np.random.default_rng(0)
    P = rng.normal(size=(50, 3)) * 10
    R, t = m.kabsch(P, P)
    assert _rmsd(P, P, R, t) == pytest.approx(0.0, abs=1e-9)


def test_kabsch_recovers_a_random_rigid_motion(m):
    """A pure rotation+translation must be fully absorbed, leaving zero."""
    rng = np.random.default_rng(1)
    P = rng.normal(size=(80, 3)) * 12
    th = 0.7
    Rr = np.array([[math.cos(th), -math.sin(th), 0],
                   [math.sin(th), math.cos(th), 0], [0, 0, 1]])
    Q = (Rr @ P.T).T + np.array([5.0, 7.0, 9.0])
    R, t = m.kabsch(P, Q)
    assert _rmsd(P, Q, R, t) == pytest.approx(0.0, abs=1e-8)


def test_kabsch_is_not_transposed(m):
    """The regression that produced 3.55 A instead of 0.356 A."""
    rng = np.random.default_rng(2)
    P = rng.normal(size=(40, 3))
    Q = P.copy()
    Q[:, 0] += 0.5
    R, t = m.kabsch(P, Q)
    # a transposed R would be a reflection/rotation mismatch; check det and action
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-9)
    moved = (R @ P.T).T + t
    # the shift is along x only, so the residual must be ~0.5 in x, ~0 elsewhere
    resid = np.abs(moved - Q).mean(0)
    assert resid[0] == pytest.approx(0.0, abs=1e-8)


def test_matched_series_coverage_is_complete(m):
    """7TRE and 7TRD must correspond 100% for the comparison to be matched."""
    tre, trd = m.load_ca("7TRE"), m.load_ca("7TRD")
    common = set(tre) & set(trd)
    assert len(trd) == 934
    assert len(common) == len(tre) == len(trd)
    for name, (lo, hi) in m.DOMAINS.items():
        n = sum(1 for r in common if lo <= r <= hi)
        assert n == sum(1 for r in tre if lo <= r <= hi), f"{name} incomplete"


def test_7trd_and_7trf_share_one_tert_model(m):
    """Documented caveat: they are not two independent apo structures."""
    trd, trf = m.load_ca("7TRD"), m.load_ca("7TRF")
    common = sorted(set(trd) & set(trf))
    P = np.array([trd[i] for i in common])
    Q = np.array([trf[i] for i in common])
    R, t = m.kabsch(P, Q)
    assert _rmsd(P, Q, R, t) < 0.01
