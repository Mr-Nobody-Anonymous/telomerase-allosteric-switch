"""Regression tests for the four open findings of the RAP-state / mechanism audit.

These encode results that were previously only in prose (or only in a Colab
notebook that was never committed) so they cannot silently rot.

Three are theorems about the kinetic model. One documents *why* the RAP-state
superposition was circular and locks in the corrected frames. One guards the
domain-range agreement between the two Kabsch implementations. One guards the
anchor-site report against contradicting itself.

Provenance of the mechanism result: the model is the one executed in the
Colab notebook (cell 26, `evaluate_hypotheses`); its saved output read
gap = 0.000478 at A0 = 0.99 and gap = 0.02178 at A0 = 0.70, both reproduced
here to <5%.
"""

from __future__ import annotations

import importlib.util
import math
import pathlib
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
for _p in (REPO, REPO / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from colab.tpp1_colab import DOMAINS, rigid_transform  # noqa: E402
from telotrade.cross_state import CANONICAL_DOMAIN_RANGES  # noqa: E402


# --------------------------------------------------------------------------
# 1. Mechanism: R1 (affinity) vs R2 (extension) under eq. 1
# --------------------------------------------------------------------------
def folds(a0: float, p: float, n: float = 5.0, a: float = 4.0) -> tuple[float, float]:
    """Fold-stimulation for the affinity arm (R1) and extension arm (R2).

        theta = P/(1+P)                      primer occupancy
        r     = (1-A0)/A0                    resting k_off / k_pol
        s     = 1/(1 + r*(1-theta))          per-repeat survival
        R1    : theta -> P*n/(1+P*n)         TPP1 lowers K_d by n
        R2    : r     -> r/a                 TPP1 raises k_pol by a
    """
    r = (1.0 - a0) / a0
    theta = p / (1.0 + p)
    theta_r1 = (p * n) / (1.0 + p * n)
    s_base = 1.0 / (1.0 + r * (1.0 - theta))
    r1 = 1.0 / (1.0 + r * (1.0 - theta_r1)) / s_base
    r2 = 1.0 / (1.0 + (r / a) * (1.0 - theta)) / s_base
    return r1, r2


def test_reproduces_colab_notebook_gaps() -> None:
    """The two numbers in the Colab notebook's saved output."""
    r1, r2 = folds(0.99, 0.9326)
    assert abs(r2 - r1) == pytest.approx(0.000478, abs=0.0001)
    r1, r2 = folds(0.70, 0.9326)
    assert abs(r2 - r1) == pytest.approx(0.02178, rel=0.05)


def test_r1_r2_degenerate_at_high_processivity() -> None:
    """At the assay's A0 = 0.99 the two hypotheses are indistinguishable."""
    r1, r2 = folds(0.99, 1.0)
    assert abs(r1 - r2) < 0.005


def test_gap_grows_as_processivity_falls() -> None:
    """The R1/R2 separation is set by A0, not by primer concentration."""
    gaps = [abs(folds(a0, 1.0)[1] - folds(a0, 1.0)[0]) for a0 in (0.99, 0.90, 0.70)]
    assert gaps[0] < gaps[1] < gaps[2]


def test_low_processivity_cannot_host_the_observed_effect() -> None:
    """The decisive negative result.

    The regime that separates R1 from R2 caps the TPP1 effect near 1.29x,
    whereas the effect to be explained is 2.7 +/- 0.6x.  A primer titration
    into that regime therefore does not test the mechanism of the effect that
    was actually measured -- it tests a much smaller effect in a crippled
    enzyme.  This is why no combination of primer titration and processivity
    lowering cleanly separates R1 from R2.
    """
    best = max(max(folds(0.70, p)) for p in np.logspace(-2, 2, 400))
    assert best < 1.5, f"A0=0.70 caps the effect at {best:.3f}x, far below 2.7x"


# --------------------------------------------------------------------------
# 2. Superposition: the frame must not contain the domain being measured
# --------------------------------------------------------------------------
def test_ifd_trap_is_inside_the_rt_frame() -> None:
    """Documents the hazard: IFD-TRAP sits entirely inside RT."""
    rt = set(DOMAINS["RT 541-935"])
    ifd = set(DOMAINS["IFD-TRAP 755-785"])
    assert ifd and ifd <= rt


def _recover(lo: int, hi: int, shift: float = 2.0) -> dict[str, float]:
    """Inject a known displacement; recover it using the same rule as the code."""
    rng = np.random.default_rng(0)
    ca: dict[int, np.ndarray] = {}
    for i in range(1, 197):
        ca[i] = rng.normal(0, 12, 3)
    for i in range(326, 541):
        ca[i] = rng.normal(0, 20, 3)
    for i in range(541, 1000):
        ca[i] = rng.normal(0, 25, 3)

    tgt = {k: v.copy() for k, v in ca.items()}
    for k in tgt:
        if lo <= k <= hi:
            tgt[k] = tgt[k] + np.array([shift, 0.0, 0.0])

    rt = list(DOMAINS["RT 541-935"])
    out: dict[str, float] = {}
    for dname, rng_ in DOMAINS.items():
        common = [i for i in rng_ if i in ca and i in tgt]
        frame = [i for i in rt if i in ca and i in tgt and i not in set(rng_)]
        if len(common) < 10 or len(frame) < 20:
            continue
        rot, tr = rigid_transform(
            np.array([ca[i] for i in frame]), np.array([tgt[i] for i in frame])
        )
        moved = (rot @ np.array([ca[i] for i in common]).T).T + tr
        out[dname] = float(
            np.linalg.norm(moved - np.array([tgt[i] for i in common]), axis=1).mean()
        )
    return out


def test_clean_frame_recovers_ground_truth() -> None:
    """Excluding the measured domain from the frame recovers 2.000 A."""
    assert _recover(755, 785)["IFD-TRAP 755-785"] == pytest.approx(2.0, abs=0.05)


def test_clean_frame_reports_no_motion_for_rigid_domains() -> None:
    assert _recover(755, 785)["TRBD 326-540"] < 0.5


def test_circular_frame_would_understate_the_truth() -> None:
    """Pin the size of the bug that was fixed: the old frame loses ~16%."""
    rng = np.random.default_rng(0)
    ca = {}
    for i in range(1, 197):
        ca[i] = rng.normal(0, 12, 3)
    for i in range(326, 541):
        ca[i] = rng.normal(0, 20, 3)
    for i in range(541, 1000):
        ca[i] = rng.normal(0, 25, 3)
    tgt = {k: v.copy() for k, v in ca.items()}
    for k in tgt:
        if 755 <= k <= 785:
            tgt[k] = tgt[k] + np.array([2.0, 0.0, 0.0])

    rt = list(DOMAINS["RT 541-935"])
    circ = [i for i in rt if i in ca]                      # the old, circular frame
    rot, tr = rigid_transform(np.array([ca[i] for i in circ]),
                              np.array([tgt[i] for i in circ]))
    common = [i for i in range(755, 786)]
    moved = (rot @ np.array([ca[i] for i in common]).T).T + tr
    got = float(np.linalg.norm(moved - np.array([tgt[i] for i in common]), axis=1).mean())
    assert 1.5 < got < 1.9, f"expected ~1.68 A attenuation, got {got:.3f}"


# --------------------------------------------------------------------------
# 3. Domain definitions must agree across the two implementations
# --------------------------------------------------------------------------
def test_ifd_trap_range_agrees_between_implementations() -> None:
    colab = DOMAINS["IFD-TRAP 755-785"]
    src_lo, src_hi = CANONICAL_DOMAIN_RANGES["IFD_TRAP"]
    assert (min(colab), max(colab)) == (src_lo, src_hi), (
        f"colab {min(colab)}-{max(colab)} != src {src_lo}-{src_hi}"
    )


# --------------------------------------------------------------------------
# 4. The anchor-site report must not contradict itself
# --------------------------------------------------------------------------
def test_wp16_9shy_row_says_contact() -> None:
    md = (REPO / "results/solve/WP16_ANCHOR_SITE_RESULT.md").read_text(encoding="utf-8")
    row = next(l for l in md.splitlines() if l.startswith("| 9SHY"))
    assert "NOT RESOLVED" not in row, f"stale verdict row: {row.strip()[:90]}"
    assert "CONTACT" in row


# ---------------------------------------------------------------- WP18 ----
WP18 = Path(__file__).resolve().parents[1] / "scripts" / "solve" / "wp18_candidate_model.py"


def _wp18():
    spec = importlib.util.spec_from_file_location("wp18", WP18)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_wp18_exists_and_is_importable():
    assert WP18.exists(), "WP18 script missing"
    _wp18()


def test_wp18_translocation_rates_match_wp3():
    m = _wp18()
    assert m.KOFF_TR_NO == pytest.approx(0.05281, abs=5e-5)
    assert m.KOFF_TR_PT == pytest.approx(0.00101, abs=5e-5)
    assert m.KOFF_TR_NO / m.KOFF_TR_PT == pytest.approx(52.3, abs=0.5)


def test_wp18_degeneracy_is_real():
    """The whole point: k_off_tr drop and k_tr rise are one observable."""
    m = _wp18()
    eps = m.KOFF_TR_NO / m.KOFF_TR_PT
    # r = COUPLING/eps reproduces the measured pair ...
    assert m.ktr_ratio_needed(eps) == pytest.approx(m.KTR_PT / m.KTR_NO_PT, rel=1e-9)
    # ... and r = 1 with a larger eps reproduces the SAME phi.
    eps_r1 = m.COUPLING
    assert m.ktr_ratio_needed(eps_r1) == pytest.approx(1.0, rel=1e-9)
    assert eps_r1 != pytest.approx(eps, rel=0.2)


def test_wp18_identity_holds():
    """log(s) - log(phi) = -6 log(1+a), independent of k_tr."""
    m = _wp18()
    a = 0.03
    for phi in (0.57, 0.99):
        s = m.s_of(a, phi)
        assert math.log(s) - math.log(phi) == pytest.approx(-6 * math.log(1 + a))


def test_wp18_only_translocation_candidates_survive():
    m = _wp18()
    eps = m.KOFF_TR_NO / m.KOFF_TR_PT
    s0 = m.s_of(m.A0, m.PHI_NO_PT)
    # translocation family: resting k_off untouched, phi reaches 0.99
    assert m.s_of(m.A0, m.PHI_PT) / s0 == pytest.approx(1.7368, abs=1e-3)
    # retention: global affinity change necessarily moves the resting k_off
    assert (m.KOFF_RESTING / eps) != pytest.approx(m.KOFF_RESTING)


def test_wp18_wp3_bracket_claim_is_retracted():
    """RESULTS.md / wp3 must not claim the range brackets the observation."""
    m = _wp18()
    lo = m.PHI_PT / m.PHI_NO_PT
    hi = lo * (m.KTR_PT / m.KTR_NO_PT)
    assert not (lo <= m.RAP_FOLD_OBS <= hi), "range would bracket - claim invalid"
    root = Path(__file__).resolve().parents[1]
    for f in ("results/solve/RESULTS.md", "scripts/solve/wp3_step_resolved_model.py"):
        lines = (root / f).read_text(encoding="utf-8").splitlines()
        for i, line in enumerate(lines):
            if "always brackets" not in line:
                continue
            # allowed only inside an explicit retraction (+/- 2 lines)
            window = " ".join(lines[max(0, i - 2): i + 3]).lower()
            assert "correct" in window or "earlier version" in window, (
                f"{f}:{i+1} still overclaims: {line.strip()}")


# --------------------------------------------------- WP19 / circular frame ----
def test_src_fit_frame_excludes_every_domain_it_reports():
    """The src/ frame must not be the catalytic core.

    TRBD 326-540, RT 541-935, IFD-TRAP 755-785 and CTE 936-1132 are each
    100% inside the catalytic core 326-1132.  Using 326-1132 as the
    superposition frame therefore lets every reported domain vote on the
    frame that reports it.
    """
    from telotrade.cross_state import CANONICAL_DOMAIN_RANGES as C

    rt_lo, rt_hi = C["RT"]
    ifd_lo, ifd_hi = C["IFD_TRAP"]
    frame = {r for r in range(rt_lo, rt_hi + 1) if not (ifd_lo <= r <= ifd_hi)}

    for name in ("TEN", "TRBD", "IFD_TRAP", "CTE"):
        lo, hi = C[name]
        inside = sum(1 for r in range(lo, hi + 1) if r in frame)
        assert inside == 0, (
            f"{name} {lo}-{hi} overlaps the fit frame in {inside} residues; "
            "the frame is still circular")
    assert len(frame) > 300, "frame unexpectedly small"


def test_wp19_published_measurements_are_encoded():
    """The values WP19 relies on must match the primary text."""
    import importlib.util as ilu

    spec = ilu.spec_from_file_location(
        "wp19", REPO / "scripts" / "solve" / "wp19_literature_closure.py")
    m = ilu.module_from_spec(spec)
    spec.loader.exec_module(m)

    assert (m.PHI0, m.PHI1) == (0.57, 0.99)          # L&C Fig 2F, A
    assert (m.KTR0, m.KTR1) == (0.07, 0.10)          # L&C Fig 2F, k
    assert m.KOFF_REST == 0.012                       # resting, unchanged
    assert m.T12_CYCLING_FOLD == 4.0                  # "at least four-fold"
    assert (m.R_HALF0, m.R_HALF1) == (0.78, 3.3)      # Wang 2007

    # s = 2^(-1/R_1/2)
    assert m.s_from_rhalf(m.R_HALF0) == pytest.approx(0.4112, abs=1e-3)
    assert m.s_from_rhalf(m.R_HALF1) == pytest.approx(0.8105, abs=1e-3)

    # the rebinding gap that WP19 turns into a prediction
    a0 = m.a_from_identity(m.s_from_rhalf(m.R_HALF0), m.PHI0)
    a1 = m.a_from_identity(m.s_from_rhalf(m.R_HALF1), m.PHI1)
    assert a0 / a1 == pytest.approx(1.65, abs=0.02)
    assert m.T12_CYCLING_FOLD / (a0 / a1) == pytest.approx(2.42, abs=0.05)


def test_wp19_retracts_the_never_measured_claim():
    """Docs must no longer assert that cycling k_off was never measured."""
    root = Path(__file__).resolve().parents[1]
    body = (root / "README.md").read_text(encoding="utf-8")
    assert "has never been measured\" was wrong" in body or \
           "Correction (WP19)" in body, "README still lacks the WP19 correction"
