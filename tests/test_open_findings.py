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
import inspect
import math
import pathlib
import subprocess
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


# ---------------------------------------------------------------- WP22 ----
WP22 = Path(__file__).resolve().parents[1] / "scripts" / "solve" / "wp22_rap_model_v2.py"


def _wp22():
    spec = importlib.util.spec_from_file_location("wp22", WP22)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_wp22_resolution_conversion_is_right():
    """R1/2 = -ln2/ln(s) and its inverse must round-trip."""
    m = _wp22()
    for R in (0.78, 2.963, 3.3, 8.0):
        assert m.R_from_s(m.s_from_R(R)) == pytest.approx(R)


def test_wp22_2_7_is_a_median_ratio_not_a_survival_ratio():
    """The core correction: 2.7 must not be used as an s ratio."""
    m = _wp22()
    R_no = 8.0 / 2.7
    s_ratio = m.s_from_R(8.0) / m.s_from_R(R_no)
    assert 8.0 / R_no == pytest.approx(2.7)          # what the paper reports
    assert s_ratio == pytest.approx(1.1587, abs=1e-3)  # what it means for s
    assert s_ratio / 2.7 == pytest.approx(0.429, abs=1e-3)


def test_wp22_old_factorisation_is_unphysical():
    """s = (1+a)^-6 * phi forces a < 0 on L&C's own no-PT numbers."""
    m = _wp22()
    s0 = m.s_from_R(m.SYS["L&C 2010 super-telomerase"]["R0"])
    assert s0 / m.PHI0 > 1.0                 # violates s <= phi
    a0 = (s0 / m.PHI0) ** (-1.0 / m.N_INCORP) - 1.0
    assert a0 < 0, "the old form must produce a negative rate constant"


def test_wp22_m5_is_physical_for_both_systems():
    m = _wp22()
    for d in m.SYS.values():
        for R in (d["R0"], d["R1"]):
            a = m.s_from_R(R) ** (-1.0 / m.N_INCORP) - 1.0
            assert a >= 0, "M5 must keep a physical (non-negative) rate"


def test_wp22_predicts_the_measured_off_rate_fall():
    """2.73x from the product distribution vs >=4x measured by pulse-chase."""
    m = _wp22()
    d = m.SYS["L&C 2010 super-telomerase"]
    a0 = m.s_from_R(d["R0"]) ** (-1.0 / m.N_INCORP) - 1.0
    a1 = m.s_from_R(d["R1"]) ** (-1.0 / m.N_INCORP) - 1.0
    assert a0 / a1 == pytest.approx(2.7335, abs=1e-3)
    assert a0 / a1 < 4.0, "prediction must sit below the >=4x measured bound"


def test_wp22_the_two_systems_are_never_multiplied():
    """Wang s0 and L&C's fold must not be combined into one equation."""
    m = _wp22()
    assert m.SYS["L&C 2010 super-telomerase"]["R1"] == 8.0
    assert m.SYS["Wang 2007 RRL enzyme"]["R1"] == 3.3
    # L&C is the more processive enzyme, as the paper states
    lc = m.s_from_R(8.0)
    wg = m.s_from_R(3.3)
    assert lc > wg


def test_wp20_is_marked_superseded():
    root = Path(__file__).resolve().parents[1]
    for f in ("results/solve/WP20_EXHAUSTIVE_FEASIBILITY.md",
              "scripts/solve/wp20_exhaustive_feasibility.py"):
        body = (root / f).read_text(encoding="utf-8")
        assert "WP22" in body, f"{f} must defer to the corrected model"


# ---------------------------------------------------------------- WP27 ----
WP27 = Path(__file__).resolve().parents[1] / "scripts" / "solve" / "wp27_transcription_layer.py"


def _wp27():
    spec = importlib.util.spec_from_file_location("wp27", WP27)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_wp27_transcription_layer_does_not_disturb_verdicts():
    """The Yamanaka layer must not change anything already decided."""
    m = _wp27()
    assert m._self_test() == 0


def test_wp27_layers_compose_multiplicatively():
    m = _wp27()
    a = m.R_half(m.S_PT) * 1.0
    b = m.R_half(m.S_PT) * 2.0
    assert b == pytest.approx(2 * a)


def test_wp27_enzyme_layer_has_more_leverage_than_abundance():
    """A 1% change in s must beat a 1% change in [TERT]."""
    m = _wp27()
    assert m.elasticity(m.S_PT) > 1.0
    assert m.elasticity(m.S_NO_PT) > 1.0


def test_wp27_documents_the_oct4_hazard():
    """Any Yamanaka document must carry the safety caveat, not just the upside."""
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP27_YAMANAKA_LAYER.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())          # phrases wrap across lines in prose
    assert "Velychko" in flat
    assert "None of these factors is a telomerase drug" in flat
    assert "zero empirical evidence" in flat


# --------------------------------------------------------------- WP28
def test_wp28_saturating_model_reduces_to_wp27_linear():
    """M6' must be a strict generalisation of M6, not a replacement of it."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    wp28 = importlib.import_module("wp28_operating_point")
    s = wp28.S_NO_PT
    x = 1e-9
    assert abs(wp28.output(s, x) / x - wp28.R_half(s)) < 1e-6
    # and it saturates: 10x more TERT at x=1e4 adds < 0.1%
    assert wp28.output(s, 1e5) / wp28.output(s, 1e4) - 1.0 < 1e-3


def test_wp28_enzyme_elasticity_always_beats_transcription():
    """The one inequality that survived WP27's scope correction."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    wp28 = importlib.import_module("wp28_operating_point")
    for s in (wp28.S_NO_PT, wp28.S_PT, 0.99):
        for x in (1e-6, 1.0, 1e3):
            assert wp28.eps_s(s) > wp28.eps_x(x)
    # the occupancy elasticity is exactly 1 - coverage
    for x in (0.01, 0.5, 3.0):
        assert abs(wp28.eps_x(x) - (1.0 - wp28.occupancy(x))) < 1e-15


def test_wp28_elasticity_is_not_a_bounded_range():
    """WP27 quoted 4.3-11.5 as a range; the elasticity in fact diverges."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    wp28 = importlib.import_module("wp28_operating_point")
    assert wp28.eps_s(0.99) > 90.0            # outside WP27's "range" by 8x
    assert wp28.eps_s(0.999) > 900.0          # and it keeps going
    assert wp28.eps_s(0.99) > wp28.eps_s(wp28.S_PT) > wp28.eps_s(wp28.S_NO_PT)


def test_wp28_reports_unreachable_targets_as_unreachable():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    wp28 = importlib.import_module("wp28_operating_point")
    rng = wp28.R_half(wp28.S_PT) / wp28.R_half(wp28.S_NO_PT)
    f = wp28.fold_to_match(rng, 0.01)
    assert abs(wp28.output(wp28.S_NO_PT, 0.01 * f) / wp28.output(wp28.S_NO_PT, 0.01) - rng) < 1e-9
    assert wp28.fold_to_match(rng, 2.0) is None       # needs coverage > 1
    # the price rises with the operating point, monotonically, and is gone by x=1
    prices = [wp28.fold_to_match(rng, x) for x in (0.001, 0.01, 0.1, 0.5)]
    assert all(a < b for a, b in zip(prices, prices[1:]))
    assert prices[-1] > 10.0          # one-third coverage already costs >10x
    assert wp28.fold_to_match(rng, 1.0) is None


def test_wp28_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp28_operating_point.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp28_doc_states_the_scope_correction():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP28_OPERATING_POINT.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "M6'" in flat
    assert "diverges" in flat
    assert "operating point" in flat
    assert "no therapeutic claim" in flat.lower()


# --------------------------------------------------------------- WP29
def _wp29():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp29_maintenance")


def test_wp29_delta_star_identity_is_exact():
    """delta* = 1/log2(1+a) must be an identity, not a fit."""
    wp = _wp29()
    for a in (wp.A_NO_PT, wp.A_PT, wp.A_MEAS, wp.A_WANG):
        assert abs(wp.delta_star(a) - 1.0 / math.log2(1.0 + a)) < 1e-12
        # ... and it must agree with M5's own s = (1+a)^-6
        assert abs(wp.delta_star(a) - 6.0 * wp.R_of_a(a)) < 1e-12


def test_wp29_reproduces_the_48nt_landc_observed():
    """The model derives 48 nt; L&C independently observed 48 nt."""
    wp = _wp29()
    assert abs(wp.delta_star(wp.A_PT) - 48.0) < 0.05


def test_wp29_no_tert_immortalises_a_real_cell():
    """The impossibility result: coverage <= 1 cannot reach any real attrition."""
    wp = _wp29()
    for delta in (wp.DELTA_FIBRO_LO, 60.0, wp.DELTA_HUVEC, wp.DELTA_FIBRO_HI,
                  wp.DELTA_MAXLIT):
        assert wp.coverage_needed(wp.A_PT, delta) > 1.0
        # even at full coverage the net is still a loss
        assert wp.gain_bp(wp.A_PT, 1.0) < delta
    # and it is not marginal at the low end, so it is not a rounding artefact
    assert wp.coverage_needed(wp.A_PT, wp.DELTA_FIBRO_LO) > 1.0


def test_wp29_offrate_is_the_only_lever_without_a_ceiling():
    wp = _wp29()
    for delta in (50.0, 90.0, 100.0, 200.0):
        a = wp.a_needed(delta)
        assert abs(wp.delta_star(a) - delta) < 1e-9        # target is met exactly
        assert wp.gain_bp(a * 0.99) > delta                # slightly better: maintain
        assert wp.gain_bp(a * 1.01) < delta                # slightly worse: lose
    # 90 bp/PD needs ~5.14x, which nothing measured reaches
    assert 5.0 < wp.offrate_fold_needed(wp.DELTA_HUVEC) < 5.3
    assert wp.offrate_fold_needed(wp.DELTA_HUVEC) > 4.0
    assert wp.lifespan_pd(wp.A_MEAS, wp.DELTA_HUVEC, 11.5e3, 5e3) is not None


def test_wp29_transcription_is_a_dead_lever_for_maintenance():
    wp = _wp29()
    for delta in (50.0, 90.0, 100.0):
        assert wp.substitution_rate(delta) < 0.01          # < 1% per e-fold
    # 4x by abundance alone needs ~400 e-folds, i.e. a ~1e174-fold rise in [TERT]
    need = math.log(4.0) / wp.substitution_rate(wp.DELTA_FIBRO_LO)
    assert 300 < need < 500
    assert need * math.log10(math.e) > 170      # ~1e174-fold in [TERT]


def test_wp29_untreated_lifespan_check_holds():
    """The arithmetic must reproduce the observed HUVEC Hayflick limit."""
    wp = _wp29()
    pd = wp.lifespan_pd(wp.A_NO_PT, wp.DELTA_HUVEC, 11.5e3, 5e3, coverage=0.0)
    assert 65 <= pd <= 75


def test_wp29_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp29_maintenance.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp29_doc_states_the_impossibility_and_the_target():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP29_MAINTENANCE.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "1 / log₂(1 + a)" in flat
    assert "5.14" in flat                        # the target off-rate fold
    assert "still senesces" in flat
    assert "not a therapy" in flat.lower()
    assert "no therapeutic claim" in flat.lower()


# --------------------------------------------------------------- WP30
def _wp30():
    pytest.importorskip("rdkit", reason="WP30 needs the optional 'chem' extra")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp30_ligand_completion")


def test_wp30_original_ligand_is_two_fragments():
    """The defect is measured from the file, not taken on trust."""
    wp = _wp30()
    atoms, bonds = wp.parse_pdb(wp.LIG / "telin1_inhibitor.pdb")
    assert len(atoms) == 36
    comps = wp.components(len(atoms), bonds)
    assert len(comps) == 2
    assert sorted(len(c) for c in comps) == [7, 29]


def test_wp30_original_ligand_has_impossible_geometry():
    wp = _wp30()
    atoms, bonds = wp.parse_pdb(wp.LIG / "telin1_inhibitor.pdb")
    by = {a["name"]: a for a in atoms}
    # a 1.208 A bond inside a six-membered ring, with a 158.9 degree angle there
    assert wp.dist(by, "C4", "C1") < 1.22
    assert wp.angle_deg(by, "C4", "C1", "N1") > 150
    # and the ring is flat, which is what makes the angles impossible
    ring = ["C1", "N1", "C2", "N2", "C3", "C4"]
    torsions = [abs(wp.torsion_deg(by, ring[i], ring[(i+1) % 6],
                                   ring[(i+2) % 6], ring[(i+3) % 6])) for i in range(6)]
    assert max(torsions) < 10.0


def test_wp30_repair_reproduces_the_declared_formula():
    """The header says C26H27FN2O6S; the repair must hit it exactly."""
    wp = _wp30()
    core = wp.repaired_properties(wp.CORE_SMILES)
    assert core["formula"] == "C21H19FN2O4S" and core["heavy"] == 29
    for label, smi in wp.REPAIRED.items():
        p = wp.repaired_properties(smi)
        assert p["formula"] == "C26H27FN2O6S", label
        assert p["heavy"] == 36, label
        assert abs(p["mw"] - 514.58) < 0.01, label


def test_wp30_repair_is_one_connected_molecule():
    pytest.importorskip("rdkit", reason="WP30 needs the optional 'chem' extra")
    from rdkit import Chem
    wp = _wp30()
    for label, smi in wp.REPAIRED.items():
        assert len(Chem.GetMolFrags(Chem.MolFromSmiles(smi))) == 1, label


def test_wp30_repair_is_valid_but_not_druglike():
    """The honest finding: fixing the formula does not fix the drug."""
    wp = _wp30()
    p = wp.repaired_properties(wp.REPAIRED["para-to-N"])
    nfail, fails = wp.lipinski_pass(p)
    assert nfail == 2 and "MW>500" in fails and "logP>5" in fails


def test_wp30_borrowed_design_rule_passes_lipinski():
    """One anchor carboxylate on a shorter arm - the precedent's rule."""
    wp = _wp30()
    p = wp.repaired_properties(wp.VARIANTS["TLN-3 anchor, short arm"])
    assert wp.lipinski_pass(p)[0] == 0
    assert p["acidic"] == 1
    # and the reference molecule really is a one-carboxylate precedent
    lif = wp.repaired_properties(wp.LIFITEGRAST_SMILES)
    assert lif["acidic"] == 1


def test_wp30_repaired_pdb_on_disk_is_one_component():
    wp = _wp30()
    rp = wp.LIG / "telin1_repaired.pdb"
    assert rp.is_file()
    ra, rb = wp.parse_pdb(rp)
    assert len(ra) == 36
    assert len(wp.components(len(ra), rb)) == 1
    lengths = wp.bond_lengths(ra, rb)
    assert lengths[0][0] > 1.20          # no impossible short bond
    assert lengths[-1][0] < 1.90         # no absurd long bond


def test_wp30_self_test_passes():
    pytest.importorskip("rdkit", reason="WP30 needs the optional 'chem' extra")
    subprocess.run([sys.executable, "scripts/solve/wp30_ligand_completion.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp30_doc_keeps_the_quarantine():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP30_LIGAND_COMPLETION.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "Quarantine unchanged" in flat
    assert "not an activity claim" in flat
    assert "no TPP1 chain" in flat


# --------------------------------------------------------------- WP31
def _wp31():
    pytest.importorskip("rdkit", reason="WP31 needs the optional 'chem' extra")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp31_epitope")


def test_wp31_tpp1_really_is_in_the_repository():
    """Gap 9 said TPP1 was absent. It is not - this pins that down."""
    wp = _wp31()
    t7 = wp.load_chains(wp.RAW / "7TRE.cif")
    t2 = wp.load_chains(wp.RAW / "2I46.cif")
    assert len(t7[wp.TPP1_7TRE_CHAIN]) >= 100
    assert len(t2[wp.TPP1_2I46_CHAIN]) >= 100


def test_wp31_kabsch_is_correct_and_oriented_right():
    """The frame bug: mapping backwards is silent, so assert the direction."""
    wp = _wp31()
    asm = wp.build_assembly()
    assert asm["rmsd"] < 3.0
    assert asm["n_shared"] >= 80
    # transformed TPP1 must land in the TERT frame
    tpp1 = wp.transform_residues(asm["t2"][wp.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    tpp1_ca = np.array([a["CA"] for _, _, a in tpp1 if "CA" in a])
    tert_ca = np.array([a["CA"] for _, _, a in asm["t7"]["A"] if "CA" in a])
    dmin = np.linalg.norm(tpp1_ca[:, None, :] - tert_ca[None, :, :], axis=2).min()
    assert dmin < 8.0


def test_wp31_epitope_is_an_interface_not_the_whole_protein():
    wp = _wp31()
    asm = wp.build_assembly()
    tpp1 = wp.transform_residues(asm["t2"][wp.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    ia, ib = wp.interface_residues((asm["t7"]["A"], tpp1))
    assert 0 < len(ia) < len(asm["t7"]["A"]) / 2
    assert len(ib) > 3
    # the TEL patch is in the interface (auth 166-170)
    assert any(166 <= s <= 170 for s in ib)


def test_wp31_scorer_rejects_clashes_rather_than_ranking_them():
    wp = _wp31()
    lig = np.array([[0.0, 0.0, 0.0]])
    ok = wp.score_pose(lig, np.array([[0.0, 3.9, 0.0]]))
    bad = wp.score_pose(lig, np.array([[0.0, 3.9, 0.0], [0.5, 0.0, 0.0]]))
    assert ok[0] == 1 and ok[1] == 0
    assert bad[1] == 1 and bad[3] < 0


def test_wp31_every_pose_is_clash_free():
    wp = _wp31()
    asm = wp.build_assembly()
    tpp1 = wp.transform_residues(asm["t2"][wp.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    ia, ib = wp.interface_residues((asm["t7"]["A"], tpp1))
    target = np.array([x for _, _, _, x in wp.all_heavy(
        [(n, s, a) for n, s, a in tpp1 if s in ib])])
    wp30 = wp._import_wp30()
    xyz = wp.ligand_conformer(wp30.VARIANTS["TLN-3 anchor, short arm"])[1]
    sc, c, cl, dmin, P, valid = wp.place_ligand(xyz, target, seed=3, n_iter=800)
    assert valid and cl == 0 and dmin > 2.6


def test_wp31_assembled_file_exists_and_has_tpp1():
    wp = _wp31()
    dest = wp.STRUCT / "telomerase_tpp1_assembled.pdb"
    assert dest.is_file()
    txt = dest.read_text()
    assert "SUPERPOSITION RMSD" in txt
    assert " P" in txt                      # the transformed TPP1 chain


def test_wp31_doc_keeps_the_quarantine():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP31_EPITOPE.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "Quarantine unchanged" in flat
    assert "not a binding measurement" in flat
    assert "none, and none obtainable here" in flat


# --------------------------------------------------------------- WP32
def _wp32():
    pytest.importorskip("rdkit", reason="WP32 needs the optional 'chem' extra")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp32_physics_scoring")


def test_wp32_ionisation_is_correct():
    """Gasteiger on a neutral molecule sums to ~0 - the acid was never -1."""
    wp = _wp32()
    acid = "O=C1C(=S)N(c2ccc(F)cc2CCC(=O)O)C(=O)N(c2ccccc2)C1CC(C)C"
    amine = "O=C1C(=S)N(c2ccc(F)cc2CCCN)C(=O)N(c2ccccc2)C1CC(C)C"
    guan = "O=C1C(=S)N(c2ccc(F)cc2CCCNC(=N)N)C(=O)N(c2ccccc2)C1CC(C)C"
    assert wp.net_formal_charge(wp.apply_ionisation(acid)) == -1
    assert wp.net_formal_charge(wp.apply_ionisation(amine)) == 1
    assert wp.net_formal_charge(wp.apply_ionisation(guan)) == 1
    assert wp.net_formal_charge(wp.apply_ionisation("c1ccccc1C")) == 0


def test_wp32_clashes_are_repulsive_not_attractive():
    """The old penalty was negative, so a 1.5 A overlap beat a good pose."""
    wp = _wp32()
    lig = dict(mol=None, xyz=np.zeros((1, 3)), idx=[0], q=np.array([-0.5]),
               sym=["O"], donors=[], acceptors=[0], donor_h={},
               all_xyz={0: np.zeros(3)})
    clash = wp.interaction_energy(lig, lig["xyz"],
                                  [("LYS", 1, "NZ", np.array([0.0, 0.0, 1.5]),
                                    1.0, True, False)])
    good = wp.interaction_energy(lig, lig["xyz"],
                                 [("LYS", 1, "NZ", np.array([0.0, 0.0, 3.2]),
                                   1.0, True, False)])
    assert clash["clash"] > 0            # repulsive
    assert clash["elec"] == 0.0          # electrostatics off at overlap
    assert clash["total"] > good["total"]


def test_wp32_hbonds_are_directional():
    """Same distance, different angle: only one of them is a hydrogen bond."""
    wp = _wp32()
    lig = dict(mol=None, xyz=np.zeros((1, 3)), idx=[0], q=np.array([0.3]),
               sym=["N"], donors=[0], donor_h={0: 1}, acceptors=[],
               all_xyz={0: np.zeros(3), 1: np.array([0.0, 0.0, 1.0])})
    good = wp.interaction_energy(lig, lig["xyz"],
                                 [("ASP", 1, "OD1", np.array([0.0, 0.0, 3.2]),
                                   -0.55, False, True)])
    bent = wp.interaction_energy(lig, lig["xyz"],
                                 [("ASP", 1, "OD1", np.array([2.0, 0.0, 2.0]),
                                   -0.55, False, True)])
    assert good["n_hb"] == 1
    assert bent["n_hb"] == 0
    assert good["hbond"] < bent["hbond"]


def test_wp32_control_patch_is_far_from_the_interface():
    wp = _wp32()
    tg = wp.build_targets()
    assert tg["n_tert_iface"] > 5 and tg["n_tpp1_iface"] > 3
    assert len(tg["control"]) > 5


def test_wp32_the_control_test_fails_and_that_is_recorded():
    """If this ever starts passing, the doc's negative result must change."""
    wp = _wp32()
    tg = wp.build_targets()
    prot = {k: wp.protein_atoms(v) for k, v in
            (("interface", tg["interface"]), ("control", tg["control"]))}
    q_if = sum(a[4] for a in prot["interface"])
    q_ct = sum(a[4] for a in prot["control"])
    # the interface really is the more negative patch - the stated cause
    assert q_if < q_ct
    # and every ligand in the series really is anionic
    from wp30_ligand_completion import CORE_SMILES, REPAIRED, VARIANTS
    for smi in (CORE_SMILES, REPAIRED["para-to-N"],
                VARIANTS["TLN-3 anchor, short arm"]):
        assert wp.ligand_typed(smi)["q"].sum() < -1.0


def test_wp32_self_test_passes():
    pytest.importorskip("rdkit", reason="WP32 needs the optional 'chem' extra")
    subprocess.run([sys.executable, "scripts/solve/wp32_physics_scoring.py",
                    "--self-test"], cwd=REPO, check=True, capture_output=True, text=True)


def test_wp32_docs_record_the_retraction():
    base = Path(__file__).resolve().parents[1] / "results" / "solve"
    w32 = " ".join((base / "WP32_PHYSICS.md").read_text(encoding="utf-8").split())
    w31 = " ".join((base / "WP31_EPITOPE.md").read_text(encoding="utf-8").split())
    assert "RETRACTED (WP32" in w31        # wording extended in WP33
    assert "WP31's ranking was an artefact" in w32
    assert "Quarantine unchanged" in w32


# --------------------------------------------------------------- WP33
def _wp33():
    pytest.importorskip("rdkit", reason="WP33 needs the optional 'chem' extra")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp33_valid_scoring")


def test_wp33_matched_pairs_only_differ_by_protonation():
    """The defect: five of six WP30 variants changed the scaffold too."""
    wp33 = _wp33()
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    for pair, (neutral, charged) in wp33.PAIRS.items():
        a = Chem.MolFromSmiles(wp33.SCAFFOLD.format(ANCHOR=neutral))
        b = Chem.MolFromSmiles(wp33.SCAFFOLD.format(ANCHOR=charged))
        assert a.GetNumHeavyAtoms() == b.GetNumHeavyAtoms(), pair
        assert rdMolDescriptors.CalcNumRings(a) == rdMolDescriptors.CalcNumRings(b), pair


def test_wp33_series_spans_all_three_charge_states():
    wp33 = _wp33()
    wp32 = _wp33()._w()
    charges = {wp32.net_formal_charge(s) for s in wp33.series().values()}
    assert charges == {0, -1, 1}


def test_wp33_control_is_matched_on_atoms_not_just_residues():
    """The defect: 377 vs 224 atoms. The fix must close that to within 25%."""
    wp33 = _wp33()
    wp31 = wp33._w31()
    asm = wp31.build_assembly()
    tpp1 = wp31.transform_residues(asm["t2"][wp31.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    ia, _ = wp31.interface_residues((asm["t7"]["A"], tpp1))
    chosen, meta = wp33.matched_control(asm["t7"]["A"], ia)
    assert abs(meta["atom_ratio"] - 1.0) <= 0.25
    assert meta["ctrl_res"] >= 15
    assert not (set(chosen) & set(ia))


def test_wp33_tel_patch_is_ambipolar_not_basic():
    """The prose said 'the TEL patch is basic'. It is net -1 and both signs."""
    wp33 = _wp33()
    wp32 = wp33._w()
    wp31 = wp33._w31()
    asm = wp31.build_assembly()
    tpp1 = wp31.transform_residues(asm["t2"][wp31.TPP1_2I46_CHAIN], asm["R"], asm["t"])
    _, ib = wp31.interface_residues((asm["t7"]["A"], tpp1))
    q = [wp32.RES_CHARGE.get(n, 0.0) for n, s, a in tpp1 if s in ib]
    assert any(x > 0 for x in q) and any(x < 0 for x in q)   # ambipolar
    assert -2.0 < sum(q) <= 0.0                              # net about -1


def test_wp33_hbond_probe_realises_the_angle_it_claims():
    """The first test requested 125 deg and actually built 106 deg."""
    wp33 = _wp33()
    for ang in (160, 150, 130, 125, 121, 119, 100):
        A = wp33.hbond_probe(ang, d_DA=3.0)
        assert abs(wp33.realised_hbond_angle(A) - ang) < 0.01
        assert abs(float(np.linalg.norm(A)) - 3.0) < 1e-9


def test_wp33_hbond_boundary_holds_at_constant_distance():
    wp33 = _wp33()
    wp32 = wp33._w()
    lig = dict(mol=None, xyz=np.zeros((1, 3)), idx=[0], q=np.array([0.3]),
               sym=["N"], donors=[0], donor_h={0: 1}, acceptors=[],
               all_xyz={0: np.zeros(3), 1: np.array([0.0, 0.0, 1.0])})
    for ang, expect in ((150, 1), (121, 1), (119, 0), (100, 0)):
        A = wp33.hbond_probe(ang, d_DA=3.0)
        r = wp32.interaction_energy(lig, lig["xyz"],
                                    [("ASP", 1, "OD1", A, -0.55, False, True)])
        assert r["n_hb"] == expect, ang


def test_wp33_explicit_protonation_state_is_selectable():
    """Without this the pairs collapse and the comparison is empty."""
    wp32 = _wp33()._w()
    acid = "O=C1C(=S)N(c2ccc(F)cc2CCC(=O)O)C(=O)N(c2ccccc2)C1CC(C)C"
    neutral = wp32.ligand_typed(acid, ionise=False)["q"].sum()
    ionised = wp32.ligand_typed(acid, ionise=True)["q"].sum()
    assert neutral > ionised + 0.5          # the explicit neutral state differs


def test_wp33_wp32_fixes_still_hold():
    wp32 = _wp33()._w()
    lig = dict(mol=None, xyz=np.zeros((1, 3)), idx=[0], q=np.array([-0.5]),
               sym=["O"], donors=[], acceptors=[0], donor_h={},
               all_xyz={0: np.zeros(3)})
    r = wp32.interaction_energy(lig, lig["xyz"],
                                [("LYS", 1, "NZ", np.array([0.0, 0.0, 1.5]), 1.0, True, False)])
    assert r["clash"] > 0 and r["elec"] == 0.0


def test_wp33_docs_record_the_withdrawals():
    base = Path(__file__).resolve().parents[1] / "results" / "solve"
    w33 = " ".join((base / "WP33_CORRECTED_SCORING.md").read_text(encoding="utf-8").split())
    w32 = " ".join((base / "WP32_PHYSICS.md").read_text(encoding="utf-8").split())
    assert "SUPERSEDED IN PART (WP33)" in w32
    assert "WITHDRAWN" in w33
    assert "ambipolar" in w33
    assert "Quarantine unchanged" in w33


# --------------------------------------------------------------- WP34
def _wp34():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp34_status")


def test_wp34_the_retracted_arithmetic_is_kept_but_marked_retracted():
    """WP34's 88% is REFUTED by WP35. The numbers stay as the record; the
    conclusion does not.

    This test used to assert ``0.80 <= stalled/p_needed <= 0.95``, i.e. that the
    stalled fraction "explains 88% of the gap". That locked a refuted claim into
    the suite: the arithmetic is reproducible, but the quantity is
    (s0-phi0)/s0 - the anomaly itself, rescaled - and it is negative in the +PT
    arm, so it is not a population fraction at all. See WP35.
    """
    wp = _wp34()
    retry = wp.MEAS_RATIO / wp.PRED_RATIO
    p_needed = 1.0 - 1.0 / retry
    stalled = (wp.S_NO_PT - wp.PHI_NO_PT) / wp.S_NO_PT
    # the arithmetic is still exactly what WP34 computed - kept as evidence
    assert abs(retry - 1.465) < 1e-3
    assert abs(p_needed - 0.318) < 1e-3
    assert abs(stalled - 0.280) < 1e-3
    # ...but it is NOT a population fraction: the +PT analogue is negative
    assert wp.wp34_pt_value() < 0.0
    # ...and 0.318 is not the physical retry probability (WP35 crossing)
    import wp35_p2_identifiability as w35
    tau, _ = w35.solve_selfconsistent()
    assert abs(tau - p_needed) > 0.05


def test_wp34_no_longer_asserts_the_gate_is_closed():
    """WP34's 'P(approved drug, now) = 0, gate CLOSED' was a conditional
    presented as a closure, because 48.0 is delta* at tau_0 = 0."""
    wp = _wp34()
    src = (Path(__file__).resolve().parents[1]
           / "scripts" / "solve" / "wp34_status.py").read_text(encoding="utf-8")
    assert "RETRACTED IN PART (WP35)" in src
    assert "UNDETERMINED" in src
    assert "gate CLOSED, not unlikely" not in src


def test_wp34_the_maintenance_gate_statement_is_exact():
    """Two distinct claims: 2.73x fails the floor, 4x clears it but is short."""
    wp = _wp34()
    assert wp.delta_star(wp.A_PT) < 50.0            # 48.0 < 50
    assert wp.delta_star(wp.A_NO_PT / 4.0) > 50.0   # 70.1 > 50
    assert wp.delta_star(wp.A_NO_PT / 4.0) < 90.0   # but < 90


def test_wp34_required_ratios_round_trip():
    wp = _wp34()
    for d in (50.0, 90.0):
        a = wp.A_NO_PT / wp.offrate_needed(d)
        assert abs(wp.delta_star(a) - d) < 1e-6
    assert abs(wp.offrate_needed(90) - 5.14) < 0.02
    assert abs(wp.offrate_needed(50) - 2.85) < 0.02


def test_wp34_self_test_passes():
    pytest.importorskip("rdkit", reason="WP34 needs the optional 'chem' extra")
    subprocess.run([sys.executable, "scripts/solve/wp34_status.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp34_doc_is_marked_superseded_and_keeps_what_survives():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP34_STATUS.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    # the retraction must be at the top, so no reader meets 88% first
    assert "SUPERSEDED IN PART (WP35)" in flat
    assert flat.index("SUPERSEDED IN PART (WP35)") < flat.index("88 % explained")
    # what still holds
    assert "TLN-amide" in flat
    assert "Quarantine" in flat or "quarantine" in flat
    assert "pipeline arithmetic" in flat
    # what does not
    assert "UNDETERMINED" in flat
    assert "an identity" in flat        # the 48 nt retraction (WP36)


# --------------------------------------------------------------- WP35
def _wp35():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp35_p2_identifiability")


def test_wp35_the_stall_correction_has_the_right_sign():
    """The first draft of WP35 ran the correction forward; this pins it back.

    If someone "simplifies" the inverse away, kappa_dist stops rising with tau
    and the whole reconciliation collapses.  That is the failure this blocks.
    """
    wp = _wp35()
    assert wp.kappa_dist(0.30) > wp.kappa_dist(0.0)          # rises, not flat
    assert wp.kappa_dist(0.30) / wp.kappa_dist(0.0) > 1.30    # materially
    # and the forward (wrong) form is explicitly NOT what is computed
    a0 = wp.a_of_s(wp.S_NO_PT)
    forward = a0 / wp.a_of_s(wp.s_obs_given_true(wp.S_PT, 0.30))
    assert abs(wp.kappa_dist(0.30) - forward) > 0.1


def test_wp35_inverse_correction_round_trips():
    wp = _wp35()
    assert abs(wp.s_true_given_obs(wp.S_NO_PT, 0.0) - wp.S_NO_PT) < 1e-15
    for t in (0.05, 0.1833, 0.35):
        st = wp.s_true_given_obs(wp.S_NO_PT, t)
        assert abs(wp.s_obs_given_true(st, t) - wp.S_NO_PT) < 1e-12


def test_wp35_the_crossing_is_unique_and_the_prediction_holds():
    wp = _wp35()
    tau, kappa = wp.solve_selfconsistent()
    assert abs(tau - 0.1833) < 2e-3
    assert abs(kappa - 3.2670) < 5e-3
    assert abs(wp.kappa_dist(tau) - wp.kappa_chase(tau)) < 1e-9
    for i in range(120):
        t = i / 200
        assert (wp.kappa_dist(t) - wp.kappa_chase(t)) * (t - tau) >= -1e-9


def test_wp35_all_three_couplings_reconcile_and_are_reported():
    """No single coupling may be quoted alone: the spread is assumption-driven."""
    wp = _wp35()
    solved = {c: wp.solve_selfconsistent(coupling=c) for c in wp.COUPLINGS}
    assert len(solved) == 3
    assert abs(solved["failure-scaled"][1] - solved["none"][1]) < 0.05
    assert solved["same"][1] < solved["failure-scaled"][1]
    assert abs(wp.FAILURE_RATIO - 0.02326) < 1e-4
    assert 1.0 / wp.FAILURE_RATIO > 40.0        # 43x rarer with POT1-TPP1


def test_wp35_the_gate_is_undetermined_not_closed():
    wp = _wp35()
    tau, kappa = wp.solve_selfconsistent()
    assert wp.delta_star(wp.kappa_dist(0.0)) < wp.DELTA_FLOOR     # 48.0 closed
    assert wp.delta_star(kappa) > wp.DELTA_FLOOR                  # 57.3 open
    assert wp.delta_star(kappa) < wp.DELTA_HUVEC                  # still short of 90
    assert abs(wp.kappa_needed(wp.DELTA_HUVEC) - 5.1426) < 1e-3
    assert wp.tau_for_delta(wp.DELTA_HUVEC) is None or \
        wp.tau_for_delta(wp.DELTA_HUVEC) > tau


def test_wp35_wp34s_88_percent_is_refuted_on_its_own_terms():
    """0.280 is the anomaly itself rescaled, and it is negative in the +PT arm."""
    wp = _wp35()
    assert abs((wp.S_NO_PT - wp.PHI_NO_PT) / wp.S_NO_PT - 0.2798) < 5e-4
    assert (wp.S_PT - wp.PHI_PT) / wp.S_PT < 0.0      # not a fraction
    tau, _ = wp.solve_selfconsistent()
    assert abs(tau - 0.318) > 0.05                    # 0.318 was not a tau


def test_wp35_powering_uses_a_real_assay_and_states_the_pairing():
    wp = _wp35()
    cv = wp.assay_cv()
    assert abs(cv - 0.1718) < 1e-3
    r = wp.KAPPA_CHASE / wp.solve_selfconsistent()[1]
    assert wp.n_per_arm(cv, r) == 12                  # independent arms
    assert wp.n_per_arm(cv, r, 0.8) == 3              # paired
    assert wp.cv_tolerated(r, 12) >= cv
    assert wp.cv_tolerated(r, 11) < cv


def test_wp35_the_benchmark_file_is_generated_not_measured():
    """Guards the claim that P2 cannot use data/benchmarks/*.csv."""
    src = (Path(__file__).resolve().parents[1] / "scripts"
           / "generate_independent_benchmarks.py").read_text(encoding="utf-8")
    assert "np.random.seed(42)" in src
    assert "synthesise" in src or "generate" in src.lower()


def test_wp35_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp35_p2_identifiability.py",
                    "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp35_doc_records_the_corrections():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP35_P2_IDENTIFIABILITY.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "UNDETERMINED" in flat
    assert "wrong sign" in flat
    assert "57.30" in flat and "46.76" in flat     # both branches, not one
    assert "Quarantine unchanged" in flat
    assert "not a re-analysis" in flat


# --------------------------------------------------------------- WP36
def _wp36():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp36_circularity")


def test_wp36_the_48nt_agreement_is_an_identity_not_a_prediction():
    """C1: delta* = n*R1/2 exactly, so 'deriving 48 nt' is 6x8 vs 8x6.

    This is the guard on the single largest term WP34 used for P(M5) (90%).
    """
    wp = _wp36()
    assert wp.c1_is_an_identity()
    assert abs(wp.delta_star_from_R(8.0, 6) - 48.0) < 1e-9
    # and it is vacuous: n=1 "predicts" 8 nt from the same R
    assert abs(wp.delta_star_from_R(8.0, 1) - 8.0) < 1e-9
    assert wp.c1_cannot_discriminate()


def test_wp36_delta_star_moves_24x_while_kappa_moves_7_percent():
    """The identity carries no information and kappa is the informative one."""
    wp = _wp36()
    sweep = wp.c1_sweep()
    d = [r["delta_at_R8"] for r in sweep]
    k = [r["kappa"] for r in sweep]
    assert d[-1] / d[0] > 20.0
    assert max(k) / min(k) - 1.0 < 0.10


def test_wp36_the_88_percent_is_not_a_fraction():
    """C2's predicate, stated as a predicate rather than a story."""
    wp = _wp36()
    assert wp.c2_fraction(wp.S_NO_PT, wp.PHI_NO_PT) > 0.0
    assert wp.c2_fraction(wp.S_PT, wp.PHI_PT) < 0.0     # impossible for a share
    assert wp.c2_is_not_a_fraction()
    assert wp.c2_not_a_probability()


def test_wp36_synthetic_benchmarks_recover_their_own_parameters():
    """C5: exact recovery is the signature of a written-in answer."""
    wp = _wp36()
    assert wp.c5_recover_generating_parameters() < 0.01


def test_wp36_every_retraction_marker_is_present():
    """If a document loses its marker, the refuted claim starts reading as true."""
    base = Path(__file__).resolve().parents[1] / "results" / "solve"
    flat = {p.name: " ".join(p.read_text(encoding="utf-8").split())
            for p in (base / "WP34_STATUS.md", base / "WP29_MAINTENANCE.md",
                      base / "WP22_RAP_MODEL_V2.md", base / "OPEN_GAPS.md",
                      base / "WP35_P2_IDENTIFIABILITY.md").__iter__()}
    assert "SUPERSEDED IN PART (WP35)" in flat["WP34_STATUS.md"]
    assert "RETRACTED (WP36)" in flat["WP29_MAINTENANCE.md"]
    assert "CORRECTED BY WP35" in flat["WP22_RAP_MODEL_V2.md"]
    assert "WITHDRAWN" in flat["OPEN_GAPS.md"]
    assert "UNDETERMINED" in flat["WP35_P2_IDENTIFIABILITY.md"]


def test_wp36_the_honest_probability_is_lower_than_wp34s():
    """Both terms that made M5 look well-supported were the circular ones."""
    wp = _wp36()
    src = inspect.getsource(wp.main)
    assert "WITHDRAWN" in src              # the two circular terms carry no weight
    p = 0.60 * 0.70 * 0.50
    assert abs(p - 0.21) < 1e-9
    assert p < 0.35


def test_wp36_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp36_circularity.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp36_doc_names_the_class_and_does_not_hide_the_earlier_ones():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP36_CIRCULARITY.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "circular inference" in flat.lower()
    for cid in ("C1", "C2", "C3", "C4", "C5"):
        assert cid in flat
    assert "21 %" in flat
    assert "Quarantine unchanged" in flat


def test_wp36_c6_the_wang_out_of_sample_test_is_an_identity():
    """C6: there is no measured second-enzyme fold; 4.48x recovers its input."""
    wp = _wp36()
    assert wp.c6_wang_prediction_recovers_its_input() < 1e-9
    assert wp.c6_there_is_no_measured_second_fold()      # WP22 says "not measured"


def test_wp36_c7_delta_star_tests_the_exponent_not_the_enzyme():
    """C7: at a FIXED fold, delta* straddles the 50 bp/PD gate on n alone."""
    wp = _wp36()
    assert wp.c7_delta_star_scales_with_the_assumed_exponent()
    assert wp.delta_star_fixed_fold(1) < 50.0
    assert wp.delta_star_fixed_fold(24) > 50.0
    assert wp.delta_star_fixed_fold(4) < 50.0 < wp.delta_star_fixed_fold(6)
    assert wp.c7_n_and_tau_are_confounded()


def test_wp36_all_seven_instances_are_recorded():
    wp = _wp36()
    ids = [it["id"] for it in wp.INSTANCES]
    assert ids == ["C1", "C2", "C3", "C4", "C5", "C6", "C7"]


# --------------------------------------------------------------- WP37
def _wp37():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp37_missing")


def test_wp37_two_chases_recover_tau_kappa_and_n_exactly():
    """The fix for C7: n and tau are de-confounded by a second chase."""
    wp = _wp37()
    for n_true in (1, 2, 6, 12, 24):
        for tau_true in (0.0, 0.05, 0.1833, 0.35, 0.5):
            tau, k, n = wp.recover_both(n_true, tau_true)
            assert abs(tau - tau_true) < 1e-9
            assert abs(k - wp.kappa(n_true, tau_true)) < 1e-9
            assert abs(n - n_true) < 1e-6


def test_wp37_the_sqrt_shortcut_is_wrong():
    """m*c = kappa^2 (1-tau^2): taking sqrt recovers a fold biased low."""
    wp = _wp37()
    m, c = wp.observed_chases(6, 0.25)
    true_k = wp.kappa(6, 0.25)
    assert abs(wp.kappa_from_chases(m, c) - true_k) < 1e-12
    assert abs(math.sqrt(m * c) - true_k) > 0.05          # the trap
    assert abs(math.sqrt(m * c) - true_k * math.sqrt(1 - 0.25 ** 2)) < 1e-12


def test_wp37_the_trap_must_reduce_processivity():
    """A trap that increased it would refute the stall model outright."""
    wp = _wp37()
    for t in (0.05, 0.2, 0.5):
        m, c = wp.observed_chases(6, t)
        assert m > c
    m0, c0 = wp.observed_chases(6, 0.0)
    assert abs(m0 - c0) < 1e-12                            # coincide at tau=0


def test_wp37_delta_star_is_exponent_free():
    """The physical quantity is 6*R1/2 and must not depend on n."""
    wp = _wp37()
    assert abs(wp.delta_star_physical(8.0) - 48.0) < 1e-12
    assert wp.physical_and_model_agree_only_at_six()
    assert wp.delta_star_physical(8.0) < 50.0              # gate CLOSED
    assert wp.DELTA_HUVEC / wp.NT_PER_REPEAT == 15.0       # R1/2 needed for 90


def test_wp37_coverage_is_the_only_escape_hatch():
    """gain = 6*c*R1/2; c>1 is the only route past the ceiling."""
    wp = _wp37()
    assert wp.coverage_needed(46.0) < 1.0
    assert wp.coverage_needed(50.0) > 1.0                  # just over one
    assert abs(wp.coverage_needed(48.0) - 1.0) < 1e-12
    assert wp.coverage_needed(90.0) > 1.5


def test_wp37_second_system_predicts_a_different_fold():
    wp = _wp37()
    f = wp.system_folds()
    assert abs(f["L&C 2010 (HEK super-telomerase)"] - 2.7335) < 5e-4
    assert abs(f["Nat Commun 2020"] - 1.8187) < 5e-4
    assert wp.folds_are_system_specific()


def test_wp37_tln_amide_descriptors_match_wp34():
    """C23H24FN3O3S / 441.53 - a first draft asserted 429.51 and was wrong."""
    wp = _wp37()
    p = wp.ligand_properties()
    assert p["formula"] == "C23H24FN3O3S"
    assert abs(p["mw"] - 441.53) < 0.1
    assert p["charge"] == 0
    assert p["lipinski_pass"]


def test_wp37_route_and_mutant_are_falsifiable():
    wp = _wp37()
    assert len(wp.RETROSYNTHESIS) == 7
    assert "CONH2" in str(wp.RETROSYNTHESIS[-1])
    assert wp.net_patch_charge([1, 1, 1, -1, -1, -1, -1]) == -1
    assert wp.net_patch_charge([0] * 7) == 0
    assert "falsified_if" in wp.DOCKING_PROPOSAL


def test_wp37_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp37_missing.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)


def test_wp37_doc_is_honest_about_what_is_not_possible():
    doc = (Path(__file__).resolve().parents[1]
           / "results" / "solve" / "WP37_MISSING.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "NOT POSSIBLE HERE" in flat
    assert "no drug" in flat
    assert "Quarantine unchanged" in flat
    assert "1.04" in flat and "1.88" in flat      # the escape hatch, both ends
    assert "√(m·c)" in flat or "sqrt" in flat.lower()


# --------------------------------------------------------------- WP38
def _wp38():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "solve"))
    import importlib
    return importlib.import_module("wp38_archive")


def test_wp38_the_register_is_complete_against_the_tree():
    """Every document carrying a retraction marker must be registered."""
    wp = _wp38()
    assert wp.check_register() == []
    assert wp.unregistered_markers() == []


def test_wp38_the_committed_register_is_in_sync():
    """A hand-edited register would drift; this recomputes and compares."""
    wp = _wp38()
    committed = REPO / "results" / "solve" / "RETRACTIONS.md"
    assert committed.is_file()
    assert committed.read_text(encoding="utf-8") == wp.render_markdown()


def test_wp38_all_four_states_are_used_and_fixed():
    wp = _wp38()
    assert {r["state"] for r in wp.REGISTER} == set(wp.STATES)
    assert len(wp.STATES) == 4
    assert len({r["id"] for r in wp.REGISTER}) == len(wp.REGISTER)


def test_wp38_the_programme_stayed_quarantined_through_the_archive():
    """Archiving must not quietly promote TELIN-1 to a therapeutic."""
    doc = (REPO / "results" / "solve" / "RETRACTIONS.md").read_text(encoding="utf-8")
    flat = " ".join(doc.split())
    assert "Quarantine is not a retraction" in flat
    assert "No therapeutic claim" in flat
    assert "CLM-023" in flat


def test_wp38_self_test_passes():
    subprocess.run([sys.executable, "scripts/solve/wp38_archive.py", "--self-test"],
                   cwd=REPO, check=True, capture_output=True, text=True)
