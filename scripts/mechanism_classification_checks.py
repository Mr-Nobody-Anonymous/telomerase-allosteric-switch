#!/usr/bin/env python3
"""Reproducible checks for results/causality/MECHANISM_CLASSIFICATION.md.

Every claim in the classification document is regenerated here, so the document does not depend
on numbers typed by hand. Five independent checks:

  1. scale invariance      -- symbolic, sympy returns exactly 0
  2. crossing formula      -- symbolic, degree in [P] must be 1 (hence unique)
  3. Jacobian rank         -- complex-step derivatives (finite differences are not precise enough)
  4. tail-slope estimator  -- closed-form r(n) for a two-species mixture
  5. TOST power            -- paired design, Bonferroni over the number of doublings

Outputs
-------
results/causality/classification_checks.json
results/causality/classification_checks.csv

Usage
-----
    python scripts/mechanism_classification_checks.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "causality"

# Reference parameter set used throughout (arbitrary units, not fitted to anything).
REF = dict(k_ext=2.0, k_loss=0.5, k_off=2.0, k_on=1.0)
RESULTS: dict[str, object] = {}
ROWS: list[dict] = []


def record(section: str, item: str, value, note: str = "") -> None:
    ROWS.append({"section": section, "item": item, "value": value, "note": note})


# =======================================================================================
# 1-2. Symbolic checks
# =======================================================================================
def symbolic_checks() -> dict:
    import sympy as sp

    k, kL, koff, kon, P, lam, al, be = sp.symbols(
        "k k_L k_off k_on P lambda alpha beta", positive=True)
    Kd = koff / kon
    R = koff * Kd / (Kd + P)
    s = k / (k + kL + R)

    # --- scale invariance: multiplies every rate by lambda, s must not change
    sub = {k: lam * k, kL: lam * kL, koff: lam * koff, kon: lam * kon}
    invariance = sp.simplify(s.subs(sub) - s)

    # --- B vs C2 crossing: linear in P (hence at most one crossing)
    sB = al * k / (al * k + kL + R)
    sC2 = k / (k + kL / be + R)
    num, _ = sp.fraction(sp.cancel(sp.together(sB - sC2)))
    bracket = sp.simplify(num / (k * (P * kon + koff)))
    degree_P = sp.degree(sp.Poly(sp.expand(bracket), P), P)
    P_star = sp.solve(sp.Eq(bracket, 0), P)
    P_star = sp.simplify(P_star[0]) if P_star else None

    # --- exact form of the primer-dependent loss term
    R_sym = koff * Kd / (Kd + P)
    R_rewritten = sp.simplify(R_sym - koff**2 / (koff + P * kon))
    loss_term_rewritten = bool(R_rewritten == 0)
    R_at_zero_primer = sp.simplify(R_sym.subs(P, 0))          # = k_off, independent of Kd
    # --- at matched Kd, scaling BOTH rates by lam scales the entire loss term by lam
    lam = sp.Symbol("lam", positive=True)
    scaled = R_sym.subs({koff: lam * koff, kon: lam * kon})
    lam_scaling_law = bool(sp.simplify(scaled - lam * R_sym) == 0)

    # --- at fixed Kd (koff = Kd*kon) the survival must not depend on the split between them
    Kd_sym = sp.Symbol("Kd", positive=True)
    s_fixed_Kd = s.subs(koff, Kd_sym * kon)
    # at fixed Kd, koff = Kd*kon, so changing kon also changes koff. The survival DOES change:
    # this is why C1a (k_off down) and C1b (k_on up) are not equivalent at matched Kd.
    fixed_kd_derivative = sp.simplify(sp.diff(s_fixed_Kd, kon))
    fixed_kd_depends_on_kon = bool(fixed_kd_derivative != 0)

    # --- three cases
    degenerate = sp.simplify(bracket.subs(be, al))          # alpha == beta
    no_extension = sp.simplify(bracket.subs(al, 1))         # alpha == 1

    results = {
        "scale_invariance_simplifies_to": str(invariance),
        "scale_invariance_is_zero": invariance == 0,
        "crossing_degree_in_P": int(degree_P),
        "crossing_is_unique": int(degree_P) == 1,
        "P_star": str(P_star),
        "degenerate_case_alpha_eq_beta": str(degenerate),
        "case_alpha_eq_1": str(no_extension),
        "loss_term_equals_koff2_over_koff_plus_P_kon": loss_term_rewritten,
        "loss_term_at_zero_primer": str(R_at_zero_primer),
        "matched_Kd_lam_scales_loss_term": lam_scaling_law,
        "at_fixed_Kd_s_depends_on_kon": fixed_kd_depends_on_kon,
        "at_fixed_Kd_ds_dkon": str(fixed_kd_derivative),
    }
    record("symbolic", "scale_invariance", results["scale_invariance_simplifies_to"])
    record("symbolic", "crossing_degree_in_P", results["crossing_degree_in_P"])
    record("symbolic", "P_star", str(P_star))
    record("symbolic", "alpha_eq_beta_numerator", str(degenerate))
    return results


# =======================================================================================
# 3. Jacobian rank with complex-step derivatives
# =======================================================================================
def s_of(P: float, k: float, kL: float, koff: float, kon: float) -> float:
    Kd = koff / kon
    return k / (k + kL + koff * Kd / (Kd + P))


def _cs_div(f, x: list[complex], i: int, h: float = 1e-20) -> float:
    xc = [complex(v) for v in x]
    xc[i] += 1j * h
    return f(*xc).imag / h


def jacobian(P_grid, theta, columns: int = 4) -> np.ndarray:
    rows = []
    for P in P_grid:
        rows.append([_cs_div(lambda *a: s_of(P, *a), theta, i) for i in range(columns)])
    return np.array(rows)


def rank_of(J: np.ndarray) -> tuple[int, np.ndarray]:
    sv = np.linalg.svd(J, compute_uv=False)
    rel = sv / sv[0]
    return int((rel > 1e-8).sum()), rel


def identifiability_checks() -> dict:
    theta = [REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"]]
    P = np.logspace(-2, 3, 12)

    r_proc, rel_proc = rank_of(jacobian(P, theta))
    J_turn = np.vstack([jacobian(P, theta), [[1, 0, 0, 0]]])
    r_turn, _ = rank_of(J_turn)
    J_free = np.vstack([np.hstack([jacobian(P, theta), np.zeros((len(P), 1))]),
                        [0, 0, 0, 0, 1]])
    r_free, _ = rank_of(J_free)
    J_bind = np.vstack([jacobian(P, theta), [[0, 0, 1, 0]], [[0, 0, 0, 1]]])
    r_bind, _ = rank_of(J_bind)

    # practical SEs on the identifiable ratios. Scale removed by construction: k_ext = 1, and
    # the ratio reference point (0.25, 2.0, 0.5) is the one quoted in the classification document.
    RATIO_REF = {"k_loss_over_k_ext": 0.25, "k_off_over_k_ext": 2.0, "k_on_over_k_ext": 0.5}
    ratio_theta = [1.0, RATIO_REF["k_loss_over_k_ext"],
                   RATIO_REF["k_off_over_k_ext"], RATIO_REF["k_on_over_k_ext"]]

    def s_ratio(P, kL, koff, kon):
        return s_of(P, 1.0, kL, koff, kon)

    s0 = np.array([s_ratio(p, *ratio_theta[1:]) for p in P])
    W = 1.0 / (0.15 * s0) ** 2                      # 15% assay CV
    # differentiate in LOG parameters, so the SEs are relative (log-scale) and the
    # condition number is scale-free
    log_theta = list(np.log(ratio_theta[1:]))
    Jr = np.array([[_cs_div(lambda *a: s_ratio(p, *np.exp(np.array(a))), log_theta, i)
                    for i in range(3)] for p in P])
    F = Jr.T @ (Jr * W[:, None])
    C = np.linalg.inv(F)
    se = np.sqrt(np.diag(C))
    cond_ratio = float(np.linalg.cond(F))

    # conditioning of the raw four-parameter form
    s_raw = np.array([s_of(p, *theta) for p in P])
    Wr = 1.0 / (0.15 * s_raw) ** 2
    Jraw = jacobian(P, theta)
    cond_raw = float(np.linalg.cond(Jraw.T @ (Jraw * Wr[:, None])))

    results = {
        "rank_processivity_only": r_proc,
        "rank_with_turnover_assumed": r_turn,
        "rank_with_free_coupling": [r_free, 5],
        "rank_with_binding": r_bind,
        "smallest_relative_sv_processivity": float(rel_proc[-1]),
        "cond_number_raw_four_rate": cond_raw,
        "cond_number_ratios": cond_ratio,
        "ratio_relative_SE": {
            "k_loss/k_ext": float(se[0]),
            "k_off/k_ext": float(se[1]),
            "k_on/k_ext": float(se[2]),
        },
        "ratio_estimates": {k: v for k, v in RATIO_REF.items()},
    }
    for key in ("rank_processivity_only", "rank_with_turnover_assumed",
                "rank_with_binding", "cond_number_raw_four_rate", "cond_number_ratios"):
        record("identifiability", key, results[key])
    return results


# =======================================================================================
# 3a. Practical identifiability in the absolute (unscaled) four-rate parametrisation.
#     Condition numbers here are NOT comparable with the log-ratio numbers below unless the
#     parametrisation is stated: everything in this block is a log-space Fisher matrix.
# =======================================================================================
def absolute_fisher(bind_cv: float | None = None) -> tuple[float, np.ndarray]:
    theta = [REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"]]
    log_theta = list(np.log(theta))
    P = np.logspace(-2, 3, 12)
    rows = []
    for p in P:
        m = s_of(p, *theta)
        f = lambda *a: s_of(p, *np.exp(np.array(a)))
        rows.append([_cs_div(f, log_theta, i) / (0.15 * m) for i in range(4)])
    if bind_cv is not None:
        for j in (2, 3):          # direct binding measures k_off and k_on, each at bind_cv
            r = [0.0, 0.0, 0.0, 0.0]
            r[j] = 1.0 / bind_cv
            rows.append(r)
    J = np.array(rows)
    F = J.T @ J
    ev = np.linalg.eigvalsh(F)
    ev = ev[ev > 0]
    cond = float(np.sqrt(ev.max() / ev.min()))
    se = 100.0 * np.sqrt(np.diag(np.linalg.pinv(F)))
    return cond, se


# =======================================================================================
# 3b. C1a vs C1b: are they really equivalent at matched Kd?
# =======================================================================================
def c1_variant_check() -> dict:
    """R = k_off*Kd/(Kd+P) is LINEAR in k_off but only weakly dependent on k_on.

    So a k_off-only change and a k_on-only change that produce the same Kd do NOT produce the
    same processivity. An earlier draft of the classification document claimed they did; this
    check is what caught it.
    """
    k, kL = REF["k_ext"], REF["k_loss"]
    cases = {}
    for label, (koff, kon) in {"wild_type": (2.0, 1.0), "C1a_koff_down": (1.0, 1.0),
                               "C1b_kon_up": (2.0, 2.0)}.items():
        cases[label] = {
            "Kd": koff / kon,
            "s_at_P0.1": s_of(0.1, k, kL, koff, kon),
            "s_at_P1": s_of(1.0, k, kL, koff, kon),
            "s_at_P100": s_of(100.0, k, kL, koff, kon),
        }
    same_kd = abs(cases["C1a_koff_down"]["Kd"] - cases["C1b_kon_up"]["Kd"]) < 1e-12
    distinguishable = abs(cases["C1a_koff_down"]["s_at_P0.1"] -
                          cases["C1b_kon_up"]["s_at_P0.1"]) > 0.01
    converge_at_high_P = abs(cases["C1a_koff_down"]["s_at_P100"] -
                             cases["C1b_kon_up"]["s_at_P100"]) < 0.01
    record("c1_variants", "same_Kd", same_kd)
    record("c1_variants", "distinguishable_at_low_P", distinguishable)
    record("c1_variants", "converge_at_high_P", converge_at_high_P)
    return {"cases": cases, "same_Kd": same_kd, "distinguishable_at_low_primer": distinguishable,
            "converge_at_high_primer": converge_at_high_P}


# =======================================================================================
# 3c. The crossing theorem, checked numerically: P* must satisfy s_B = s_C2, and the number of
#     sign changes in [P] must be 0 when P* < 0 and exactly 1 when P* > 0. A sign error in the
#     displayed closed form (which this script once contained) fails both tests.
# =======================================================================================
def P_star_numeric(a, b, k, kL, koff, kon) -> float:
    return koff * (-a * b * koff - a * kL + b * kL + b * koff) / (kL * kon * (a - b))


def _sB_minus_sC2(P, a, b, k, kL, koff, kon):
    Kd = koff / kon
    R = koff * Kd / (Kd + P)
    return a * k / (a * k + kL + R) - k / (k + kL / b + R)


def crossing_numeric_check() -> dict:
    grid = np.logspace(-3, 4, 20001)
    rows = []
    for a in (0.5, 2.0, 5.0, 30.0):
        for b in (0.5, 1.0, 3.42, 10.0):
            if abs(a - b) < 1e-12:
                continue
            ps = P_star_numeric(a, b, REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"])
            d = _sB_minus_sC2(grid, a, b, REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"])
            n_sign = int((np.diff(np.sign(d)) != 0).sum())
            residual = abs(_sB_minus_sC2(ps, a, b, REF["k_ext"], REF["k_loss"],
                                        REF["k_off"], REF["k_on"])) if ps > 0 else None
            rows.append({"alpha": a, "beta": b, "P_star": ps, "sign_changes": n_sign,
                         "residual_at_P_star": residual})
    ok_residual = all(r["residual_at_P_star"] < 1e-10 for r in rows if r["P_star"] > 0)
    ok_count = all((r["sign_changes"] == (1 if r["P_star"] > 0 else 0)) for r in rows)
    record("crossing_numeric", "P_star_satisfies_equality", ok_residual)
    record("crossing_numeric", "sign_change_count_matches_sign_of_P_star", ok_count)
    return {"rows": rows, "P_star_satisfies_equality": ok_residual,
            "sign_change_count_matches_sign_of_P_star": ok_count,
            "has_positive_example": any(r["P_star"] > 0 for r in rows)}


# =======================================================================================
# 4. Tail-slope estimator
# =======================================================================================
def mixture_pmf(n_max: int, s: float, phi: float, s2: float) -> np.ndarray:
    n = np.arange(1, n_max + 1)
    q = (1 - phi) * (1 - s) * s ** (n - 1) + phi * (1 - s2) * s2 ** (n - 1)
    return q / q.sum()


def tail_slope_checks() -> dict:
    """Show that the tail slope converges to the longest-lived species, not the majority."""
    cases = []
    for s, phi, s2 in ((0.60, 0.01, 0.95), (0.70, 0.01, 0.90), (0.85, 0.05, 0.95),
                       (0.999, 0.01, 0.95), (0.999, 0.005, 0.95)):
        q = mixture_pmf(400, s, phi, s2)
        r = q[1:] / q[:-1]
        tail = float(np.mean(r[40:80]))
        cases.append({"s_majority": s, "phi": phi, "s2": s2, "tail_slope": tail,
                      "delta_s": tail - s, "predicted_by_max": max(s, s2),
                      "naive_phi_1_minus_s": phi * (1 - s)})
        record("tail_slope", f"s={s:.3g},phi={phi},s2={s2}", round(tail, 4))

    # the closed form r(n) = s*[1+kappa*rho**n]/[1+kappa*rho**(n-1)] must match the simulation
    s, phi, s2 = 0.60, 0.01, 0.95
    q = mixture_pmf(400, s, phi, s2)
    r_sim = q[1:] / q[:-1]
    kappa = phi * (1 - s2) / ((1 - phi) * (1 - s))
    rho = s2 / s
    n = np.arange(1, 400)
    r_closed = s * (1 + kappa * rho ** n) / (1 + kappa * rho ** (n - 1))
    closed_form_max_abs_error = float(np.max(np.abs(r_sim - r_closed[: len(r_sim)])))

    # The tail slope is a WINDOW-DEPENDENT function, not a property of the distribution.
    # Same truth, different read window => wildly different apparent s.
    windows = [(3, 10), (10, 16), (16, 25), (25, 40), (40, 80)]
    window_sweep = []
    for s_w, phi_w, s2_w in ((0.60, 0.01, 0.95), (0.70, 0.01, 0.90), (0.85, 0.05, 0.95),
                             (0.999, 0.01, 0.95)):
        q_w = mixture_pmf(400, s_w, phi_w, s2_w)
        r_w = q_w[1:] / q_w[:-1]
        row = {"s_majority": s_w, "phi": phi_w, "s2": s2_w, "by_window": {}}
        for n0, n1 in windows:
            row["by_window"][f"{n0}-{n1}"] = float(np.mean(r_w[n0:n1]))
        row["spread"] = max(row["by_window"].values()) - min(row["by_window"].values())
        window_sweep.append(row)
        record("tail_slope_window", f"s={s_w},phi={phi_w},s2={s2_w}",
               f"spread {row['spread']:.4f} over windows {min(row['by_window'].values()):.4f}-"
               f"{max(row['by_window'].values()):.4f}")

    # fabricated C2 from zero true loss
    fabricated = []
    for phi, s2 in ((0.005, 0.95), (0.01, 0.90), (0.02, 0.90)):
        q = mixture_pmf(400, 1 - 1e-6, phi, s2)
        r = q[1:] / q[:-1]
        s_obs = float(np.mean(r[60:120]))
        fabricated.append({"phi": phi, "s2": s2, "s_obs": s_obs,
                           "implied_k_loss_over_k_ext": (1 - s_obs) / s_obs})
        record("fabricated_C2", f"phi={phi},s2={s2}", round((1 - s_obs) / s_obs, 4))

    return {"cases": cases, "closed_form_max_abs_error": closed_form_max_abs_error,
            "fabricated_C2": fabricated, "window_sweep": window_sweep}


# =======================================================================================
# 5. TOST power
# =======================================================================================
def tost_replicates(delta: float, cv: float, power: float = 0.80,
                    alpha: float = 0.05, paired: bool = True,
                    n_tests: int = 1) -> int:
    from scipy.stats import norm

    z = norm.ppf(1 - alpha / n_tests) + norm.ppf(power)
    factor = 1.0 if paired else np.sqrt(2.0)
    return int(np.ceil((z * cv * factor / delta) ** 2))


def power_checks() -> dict:
    table = []
    for delta in (0.10, 0.20):
        for cv in (0.345, 0.15):
            for n_tests in (1, 3):
                n = tost_replicates(delta, cv, n_tests=n_tests)
                table.append({"delta_log": delta, "cv": cv, "n_tests": n_tests, "n": n})
                record("power", f"delta={delta},cv={cv},tests={n_tests}", n)
    return {"tost_table": table}


# =======================================================================================
def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    RESULTS["symbolic"] = symbolic_checks()
    RESULTS["identifiability"] = identifiability_checks()
    for label, cv in (("processivity_only", None), ("plus_binding_10pct", 0.10)):
        cond_a, se_a = absolute_fisher(cv)
        record("absolute_fisher", f"{label}_cond", cond_a)
        RESULTS.setdefault("absolute_fisher", {})[label] = {
            "condition_number": cond_a,
            "relative_SE_pct": dict(zip(["k_ext", "k_loss", "k_off", "k_on"], se_a.tolist())),
        }
    RESULTS["c1_variants"] = c1_variant_check()
    RESULTS["crossing_numeric"] = crossing_numeric_check()
    RESULTS["tail_slope"] = tail_slope_checks()
    RESULTS["power"] = power_checks()

    (OUT / "classification_checks.json").write_text(json.dumps(RESULTS, indent=2))
    with (OUT / "classification_checks.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["section", "item", "value", "note"])
        w.writeheader()
        w.writerows(ROWS)

    print(json.dumps({k: (v if k != "tail_slope" else
                          {"closed_form_max_abs_error": v["closed_form_max_abs_error"],
                           "cases": v["cases"][:3]})
                      for k, v in RESULTS.items()}, indent=2)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
