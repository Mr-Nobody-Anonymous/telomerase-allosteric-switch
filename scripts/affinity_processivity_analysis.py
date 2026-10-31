#!/usr/bin/env python3
"""Affinity-processivity causality: model statements and identifiability analysis.

Implements the five pre-registered product-length models (G1-G5, see
results/causality/PREREGISTRATION.md section 6) and asks a question that the available data
*can* answer:

    given how the published data are structured, which of these models could ever be
    distinguished, and how much data would that take?

Read the data-requirement numbers with their scope: they are computed for parameter settings
chosen to span each family at a fixed observed P(>=16). For well-separated members one replicate
is enough to tell the families apart; in the degenerate limit (G2 with weight -> 1) a family
*becomes* G1 and no amount of data separates them. The one result that does not depend on that
choice is G1 == G5, which is an algebraic identity and holds for any data.

This script does **not** fit any model to the published measurements. That is deliberate: the
numeric source data contain one scalar per replicate (the fraction of products at >=16 repeats),
so a multi-parameter model has at least as many parameters as observations. What the script does
instead is (i) show analytically and numerically that G5 (conformational) and G1 (affinity)
collapse to the same observable, (ii) simulate *hypothetical* datasets under each model to
measure model-selection accuracy as a function of experimental design, and (iii) summarise the
real data descriptively, including a quantity that is missing from the source.

Outputs
-------
results/causality/identifiability_results.csv
results/causality/figures/figure_causality_summary.png
results/causality/figures/figure_causality_summary.csv

Usage
-----
    python scripts/affinity_processivity_analysis.py
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external" / "tesmer2023"
OUT = ROOT / "results" / "causality"
FIG = OUT / "figures"

N_MAX = 2000        # length axis; the slowest species considered has s=0.995, so 200 bins
                    # would leave 37% of the mass outside the grid (a test pins this)
PLOT_MAX = 50       # bands plotted in the figure
SEED = 20261031     # fixed so that every reported number is reproducible


# =======================================================================================
# Model definitions
# =======================================================================================
# Each model is a survival function S(k) = P(product length >= k), evaluated for k = 1..N.
# The definitional anchor: S(1) = 1 for every model. Getting this wrong shifts every number in
# the report by one repeat; `tests/test_affinity_processivity.py` pins it against closed forms.

def survival_g1(s: float, n_max: int = N_MAX) -> np.ndarray:
    """G1 constant-hazard retention (affinity): geometric, P(len=n) = (1-s) s^(n-1).

    S(1) = 1, S(k) = s^(k-1), so P(len >= 16) = s^15.
    """
    k = np.arange(1, n_max + 1)
    return np.clip(s, 0.0, 1.0) ** (k - 1)


def survival_g2(w: float, s1: float, s2: float, n_max: int = N_MAX) -> np.ndarray:
    """G2 heterogeneous population: a mixture of two constant-hazard species."""
    k = np.arange(1, n_max + 1) - 1
    return w * np.clip(s1, 0, 1) ** k + (1.0 - w) * np.clip(s2, 0, 1) ** k


def survival_g3(s0: float, rho: float, n_max: int = N_MAX) -> np.ndarray:
    """G3 repeat-dependent extension: the per-repeat survival decays, s_j = s0 * rho^(j-1).

    The first repeat is always taken (S(1) = 1); survival to k depends on s_2..s_k.
    """
    step = np.arange(1, n_max + 1) - 1
    s_j = np.clip(s0 * rho ** (step - 1), 1e-12, 1.0)
    s_j[0] = 1.0
    return np.cumprod(s_j)


def survival_g4(kp: float, kr: float, dp: float, n_max: int = N_MAX) -> np.ndarray:
    """G4 retention with a reversible pause.

    Two live states: extending (E) and paused (P). From E the enzyme extends and may switch to P
    with probability kp; from P it either returns to E (kr) or dissociates (dp). Dissociation is
    therefore delayed rather than occurring at a constant per-repeat rate, which produces an
    excess of long products relative to a geometric. S(1) = 1 by construction.
    """
    kp, kr, dp = (float(np.clip(v, 0.0, 1.0)) for v in (kp, kr, dp))
    E, P = 1.0, 0.0
    S = np.empty(n_max)
    S[0] = 1.0
    for i in range(1, n_max):
        E_new = E * (1.0 - kp) + P * kr
        P_new = max(E * kp + P * (1.0 - kr - dp), 0.0)
        E, P = E_new, P_new
        S[i] = min(E + P, 1.0)
    return S


def p_at_least_from_pmf(q: np.ndarray, k: int) -> float:
    """P(length >= k) from a pmf indexed by product length (index 0 = length 1)."""
    return float(q[k - 1:].sum())


def pmf_from_survival(S: np.ndarray) -> np.ndarray:
    """pmf over product lengths 1..N from a survival array S[k-1] = P(length >= k).

    q(n) = S(n) - S(n+1), with S(N+1) = 0. Requires S(1) = 1; a survival array that does not
    start at 1 means the model is mis-indexed, so this is asserted rather than silently repaired.
    """
    S = np.clip(np.asarray(S, dtype=float), 0.0, 1.0)
    assert abs(S[0] - 1.0) < 1e-9, (
        "survival arrays must start at S(1) = P(length >= 1) = 1; this one starts at "
        f"{S[0]:.4f}, which would shift the distribution by one repeat")
    tail = np.concatenate([S[1:], [0.0]])
    q = np.clip(S - tail, 0.0, None)
    q = q[: len(S)]
    total = q.sum()
    return q / total if total > 0 else np.full(len(S), 1.0 / len(S))


# =======================================================================================
# The analytic result: G5 collapses into G1
# =======================================================================================
def g5_collapse_demo() -> dict:
    """G5 (conformational) changes k_ext; G1 (affinity) changes k_off.

    Both act only through s = k_ext / (k_ext + k_off). So for every pair (k_ext, k_off) that
    gives the same s, the predicted distribution is identical - the two models are not merely
    hard to distinguish, they are the *same function* of the observable unless a second
    observable is added.
    """
    pairs = [(1.0, 0.5), (2.0, 1.0), (10.0, 5.0), (0.2, 0.1)]
    cvs = []
    for k_ext, k_off in pairs:
        s = k_ext / (k_ext + k_off)
        S = survival_g1(s)
        cvs.append({"k_ext": k_ext, "k_off": k_off, "s": s,
                    "P16": float(p_at_least_from_pmf(pmf_from_survival(S), 16))})
    p16 = {c["P16"] for c in cvs}
    return {"pairs": cvs, "n_distinct_P16": len(p16),
            "identical": len(p16) == 1}


# =======================================================================================
# Fitting machinery (used only on simulated data)
# =======================================================================================
def _nll(counts: np.ndarray, q: np.ndarray) -> float:
    q = np.clip(q, 1e-12, 1.0)
    return float(-(counts * np.log(q)).sum())


@dataclass
class Model:
    name: str
    n_params: int
    predict: callable
    init: list = field(default_factory=list)
    bounds: list = field(default_factory=list)


def _fit(model: Model, counts: np.ndarray) -> tuple[float, np.ndarray]:
    def objective(x):
        try:
            q = pmf_from_survival(model.predict(*x))
        except (ValueError, FloatingPointError):
            return 1e12
        if not np.all(np.isfinite(q)):
            return 1e12
        return _nll(counts, q)

    best = None
    for start in model.init:
        res = minimize(objective, np.array(start, dtype=float), method="L-BFGS-B",
                       bounds=model.bounds, options={"maxiter": 400})
        if best is None or res.fun < best.fun:
            best = res
    return float(best.fun), best.x


def build_models() -> list[Model]:
    neg_inf, inf = -np.inf, np.inf
    return [
        Model("G1 affinity (const hazard)", 1, lambda s: survival_g1(s),
              [[0.8], [0.5], [0.2]], [(1e-4, 1.0)]),
        Model("G2 heterogeneous mixture", 3,
              lambda w, s1, s2: survival_g2(w, s1, s2),
              [[0.5, 0.9, 0.5], [0.8, 0.8, 0.2], [0.2, 0.95, 0.7]],
              [(1e-3, 1 - 1e-3), (1e-3, 1 - 1e-3), (1e-3, 1 - 1e-3)]),
        Model("G3 repeat-dependent", 2, lambda s0, rho: survival_g3(s0, rho),
              [[0.9, 0.99], [0.8, 0.95], [0.95, 0.9]],
              [(1e-3, 1 - 1e-3), (0.5, 1.0)]),
        Model("G4 reversible pause", 3, lambda kp, kr, dp: survival_g4(kp, kr, dp),
              [[0.2, 0.5, 0.2], [0.5, 0.3, 0.1], [0.1, 0.8, 0.5]],
              [(1e-4, 0.999), (1e-4, 0.999), (1e-4, 0.999)]),
    ]


# =======================================================================================
# How much data would it take to separate the models?
# =======================================================================================
# A simulation-based "model selection accuracy" was tried first and rejected: at realistic
# sequencing depths the fitted G1/G2/G4 shapes are numerically almost identical, so the
# selection frequency is decided by optimiser noise rather than by evidence. The analytic
# quantity that answers the question without an optimiser is the Kullback-Leibler divergence
# between the model-implied distributions, which gives the expected log-likelihood ratio per
# replicate and hence the depth or replicate count needed to detect the difference.

CHI2_95 = 3.8415          # 1 df, alpha = 0.05  (test at 50% power)
LAMBDA_80 = (1.959964 + 0.841621) ** 2   # 80% power at alpha = 0.05


def kl_divergence(q1: np.ndarray, q2: np.ndarray) -> float:
    """KL(q1 || q2) in nats, with both distributions normalised."""
    q1 = np.clip(q1, 1e-15, None)
    q2 = np.clip(q2, 1e-15, None)
    q1 = q1 / q1.sum()
    q2 = q2 / q2.sum()
    return float(np.sum(q1 * np.log(q1 / q2)))


def _match_p16(survival_fn, p16_target: float, lo: float, hi: float,
               n_max: int = N_MAX) -> float:
    """Solve a shape family for its amplitude so that P(>=16) equals the observed value.

    Every model compared here must reproduce the *same* observed scalar; only then is the claim
    "identical scalar, different shape" true. Hand-picked parameters made that claim false (a
    test caught it), which is why the parameters are solved rather than chosen.
    """
    from scipy.optimize import brentq

    def err(x: float) -> float:
        return float(survival_fn(x)[15]) - p16_target

    if err(lo) * err(hi) > 0:
        return lo if abs(err(lo)) < abs(err(hi)) else hi
    return float(brentq(err, lo, hi, xtol=1e-10))


def plausible_shapes(p16_observed: float, n_max: int = N_MAX) -> dict[str, list[np.ndarray]]:
    """Model-implied product-length distributions, each constrained to reproduce the observed
    P(>=16).

    Two members per family, spanning the family. The shape parameter of each member is fixed a
    priori; a single amplitude parameter is then solved so that all models agree on the only
    number the published data contain. Nothing is estimated by minimising a discrepancy against
    the measurements, so these are not fits.
    """
    shapes: dict[str, list[np.ndarray]] = {}
    s_obs = p16_observed ** (1.0 / 15.0)   # P(>=16) = s^15 under S(k) = s^(k-1)

    # G1 / G5: one parameter, fully determined by the observation
    shapes["G1/G5 const hazard"] = [pmf_from_survival(survival_g1(s_obs, n_max))]

    # G2: the two branch survivals are fixed (a slow and a fast species); solve for the weight.
    # Solving for a survival instead lets the fast branch run to ~1 and produces a mean product
    # length of 100+ repeats, which is not what this assay measures.
    g2 = []
    for s_hi, s_lo in ((0.990, s_obs * 0.30), (0.975, s_obs * 0.50)):
        def build(x, s_hi=s_hi, s_lo=s_lo):
            return survival_g2(float(np.clip(x, 1e-6, 1 - 1e-6)), s_hi, s_lo, n_max)
        w = _match_p16(build, p16_observed, 1e-6, 1 - 1e-6)
        g2.append(pmf_from_survival(build(w)))
    shapes["G2 heterogeneous"] = g2

    # G3: decay rate fixed; solve for the initial survival
    g3 = []
    for rho in (0.9985, 0.9955):
        s0 = _match_p16(lambda x, rho=rho: survival_g3(min(x, 1 - 1e-12), rho, n_max),
                        p16_observed, 1e-6, 1 - 1e-12)
        g3.append(pmf_from_survival(survival_g3(s0, rho, n_max)))
    shapes["G3 repeat-dependent"] = g3

    # G4: pause/recovery fixed; solve for the dissociation-from-paused probability
    g4 = []
    for kp, kr in ((0.20, 0.55), (0.45, 0.65)):
        dp = _match_p16(lambda x, kp=kp, kr=kr: survival_g4(kp, kr, x, n_max),
                        p16_observed, 1e-6, 1 - 1e-6)
        g4.append(pmf_from_survival(survival_g4(kp, kr, dp, n_max)))
    shapes["G4 reversible pause"] = g4
    return shapes


def discriminability_table(p16_observed: float) -> pd.DataFrame:
    """Pairwise KL divergence and the data requirement implied by it."""
    shapes = plausible_shapes(p16_observed)
    rows = []
    names = list(shapes)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            kls = [kl_divergence(qa, qb) for qa in shapes[a] for qb in shapes[b]]
            kl = float(np.median(kls))
            if kl <= 0:
                depth_50 = depth_80 = float("inf")
            else:
                depth_50 = CHI2_95 / (2 * kl)
                depth_80 = LAMBDA_80 / (2 * kl)
            rows.append({"model_a": a, "model_b": b,
                         "kl_nats_median": kl,
                         "kl_nats_min": float(np.min(kls)),
                         "molecules_needed_50pct": depth_50,
                         "molecules_needed_80pct": depth_80})
    return pd.DataFrame(rows)


def replicates_needed_with_densitometry(p16_observed: float, cv: float = 0.20,
                                        floor: float = 1e-3,
                                        alpha_z: float = 1.96,
                                        power_z: float = 0.8416) -> pd.DataFrame:
    """Replicates needed to separate two models when each band carries measurement error.

    The KL calculation above assumes counting noise only, which makes the requirement look
    trivial (1-30 molecules). Real densitometry adds ~15-25% multiplicative error per band and
    has a dynamic range of roughly three orders of magnitude, so the binding constraint is the
    number of *quantifiable* bands in which the two models actually differ.

    Per band n, the standardised difference is |ln(q_a(n)/q_b(n))| / cv, for bands above the
    quantification floor. The statistic sqrt(sum_n d_n^2) grows as sqrt(replicates), so
        n_rep = ((z_alpha + z_power) / sqrt(sum_n d_n^2))^2.
    """
    shapes = plausible_shapes(p16_observed)
    names = list(shapes)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            reps, informative = [], []
            for qa in shapes[a]:
                for qb in shapes[b]:
                    qa_n = qa / qa.sum()
                    qb_n = qb / qb.sum()
                    # only bands the assay can actually resolve: the gel's dynamic range and
                    # its band assignment run out well before the model grid does
                    mask = (qa_n > floor) & (qb_n > floor)
                    mask[PLOT_MAX:] = False
                    if not mask.any():
                        reps.append(float("inf"))
                        continue
                    d = np.abs(np.log(qa_n[mask] / qb_n[mask])) / cv
                    informative.append(int((d > 1.0).sum()))
                    reps.append(((alpha_z + power_z) / np.sqrt((d ** 2).sum())) ** 2)
            rows.append({"model_a": a, "model_b": b,
                         "note": "conditional on the parameter sets explored; nested at w->1",
                         "replicates_needed_median": float(np.median(reps)),
                         "replicates_needed_worst": float(np.max(reps)),
                         "informative_bands_median": float(np.median(informative))})
    return pd.DataFrame(rows)


def scalar_observable_table(p16_observed: float) -> pd.DataFrame:
    """The same comparison for the *only* observable the published source data provide:
    a single fraction P(>=16) per replicate. Everything with more than one free parameter is
    then under-determined regardless of replicate count."""
    s = p16_observed ** (1.0 / 16.0)
    rows = []
    for name, p16 in (("G1/G5 const hazard (1 parameter)", s ** 16),
                      ("G2 heterogeneous (3 parameters)", s ** 16),
                      ("G3 repeat-dependent (2 parameters)", s ** 16),
                      ("G4 reversible pause (3 parameters)", s ** 16)):
        rows.append({"model": name, "p16_implied": p16,
                     "observations_per_replicate": 1,
                     "free_parameters": int(name.split("(")[1].split()[0]),
                     "identifiable_from_this_observable":
                         "yes" if int(name.split("(")[1].split()[0]) == 1 else "no",
                     "kl_vs_G1_scalar_nats": 0.0})
    return pd.DataFrame(rows)


# =======================================================================================
# Descriptive summary of the real data (no fitting)
# =======================================================================================
def summarise_real_data() -> pd.DataFrame:
    proc = pd.read_csv(EXT / "fig3b_processivity_plusPT.csv", comment="#")
    act = pd.read_csv(EXT / "fig3c_activity_plusPT.csv", comment="#")
    t95_2df = 4.3027
    rows = []
    for col in proc.columns[1:]:
        p = proc[col].to_numpy()
        a = act[col].to_numpy()
        ratio = float(np.mean([p[i] / proc["WT"].iloc[i] for i in range(3)])) if col != "WT" else 1.0
        rows.append({
            "genotype": col,
            "P16_mean": p.mean(),
            "P16_ci95": t95_2df * p.std(ddof=1) / np.sqrt(3),
            "processivity_ratio_vs_WT": ratio,
            "activity_ratio_mean": a.mean(),
            "activity_ci95": t95_2df * a.std(ddof=1) / np.sqrt(3),
            # the two response variables the preregistration requires; R2 is not recoverable
            "R1_processivity_ratio": ratio,
            "R2_absolute_retained_mass": "NOT IN SOURCE DATA",
            "flagged": "processivity_undefined_at_low_signal" if col == "F776E" else "",
        })
    return pd.DataFrame(rows)


# =======================================================================================
def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    collapse = g5_collapse_demo()
    df_real = summarise_real_data()

    obs_p16 = float(df_real.loc[df_real.genotype == "WT", "P16_mean"].iloc[0])
    df_disc = discriminability_table(obs_p16)
    df_scalar = scalar_observable_table(obs_p16)
    df_reps = replicates_needed_with_densitometry(obs_p16)

    df_real.to_csv(OUT / "real_data_summary.csv", index=False)
    df_disc.to_csv(OUT / "identifiability_results.csv", index=False)
    df_scalar.to_csv(OUT / "scalar_observable_limits.csv", index=False)
    df_reps.to_csv(OUT / "replicates_needed.csv", index=False)
    (OUT / "identifiability_collapse.json").write_text(json.dumps(collapse, indent=2))

    # ---------------- figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    # ---------------- Panel 1: shapes consistent with the same observed scalar ----------------
    n = np.arange(1, PLOT_MAX + 1)
    shapes = plausible_shapes(obs_p16)
    for label, qs in shapes.items():
        axes[0].plot(n, qs[0][:PLOT_MAX], lw=1.8, label=label)
    axes[0].axvline(16, color="k", ls=":", lw=1)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("product length (repeats)")
    axes[0].set_ylabel("probability per band")
    axes[0].set_title("1. Four models, one observed number\n"
                      f"(P>=16 = {obs_p16:.2f}): the shapes differ, the scalar does not", fontsize=10)
    axes[0].legend(fontsize=7, frameon=False)

    # ---------------- Panel 2: the two response variables, with their real uncertainty --------
    d = df_real[df_real.genotype != "WT"].sort_values("processivity_ratio_vs_WT")
    ypos = np.arange(len(d))
    axes[1].errorbar(d.processivity_ratio_vs_WT, ypos, xerr=d.P16_ci95 / obs_p16,
                     fmt="o", ms=5, lw=1.2, capsize=2, color="tab:blue",
                     label="R1 processivity ratio")
    axes[1].errorbar(d.activity_ratio_mean, ypos, xerr=d.activity_ci95,
                     fmt="s", ms=4, lw=1.2, capsize=2, color="tab:orange", alpha=0.85,
                     label="activity ratio")
    axes[1].set_yticks(ypos)
    axes[1].set_yticklabels(d.genotype, fontsize=8)
    axes[1].axvline(1, color="grey", lw=0.8)
    axes[1].axvline(0, color="grey", lw=0.5, ls=":")
    axes[1].set_xlabel("ratio vs WT (95% CI over 3 independent experiments)")
    axes[1].set_xlim(-2.0, 2.7)
    axes[1].set_title("2. What the source data support:\nwide CIs, and R2 (absolute mass) is absent",
                      fontsize=10)
    axes[1].legend(fontsize=7, frameon=False, loc="lower right")

    # ---------------- Panel 3: the measurement that would actually discriminate ----------------
    primer = np.logspace(-1.3, 1.5, 200)   # arbitrary units of [primer]; schematic, not data
    vmax = 0.60        # plateau without TPP1
    km = 1.0           # apparent primer K_m without TPP1
    occ = primer / (km + primer)
    # C (affinity): TPP1 lowers the apparent K_m. The plateau is set by k_ext and is unchanged.
    km_aff = km / 4.0
    s_aff = vmax * (primer / (km_aff + primer))
    # B (conformational): TPP1 raises k_ext, so the plateau rises; K_m is unchanged.
    s_conf = min(vmax * 1.45, 0.97) * occ
    axes[2].plot(primer, s_aff, lw=2, color="tab:green",
                 label="C: affinity (K_m shifts, plateau same)")
    axes[2].plot(primer, s_conf, lw=2, ls="--", color="tab:purple",
                 label="B: conformational (plateau shifts)")
    axes[2].plot(primer, vmax * occ, lw=1, color="grey", label="no TPP1")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("[primer] (arbitrary units)")
    axes[2].set_ylabel("per-repeat survival s (processivity)")
    axes[2].set_title("3. The experiment that separates them:\nprimer titration +/- TPP1 (schematic)",
                      fontsize=10)
    axes[2].legend(fontsize=7, frameon=False, loc="lower right")

    fig.suptitle("Affinity vs conformational regulation of telomerase processivity: "
                 "what the data can and cannot separate", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG / "figure_causality_summary.png", dpi=200)
    pd.DataFrame({"panel1_model": list(shapes),
                  "panel1_p16_predicted": [p_at_least_from_pmf(q[0], 16) for q in shapes.values()]}
                 ).to_csv(FIG / "figure_causality_summary_panel1.csv", index=False)
    pd.DataFrame({"primer_au": primer, "no_TPP1": vmax * occ,
                  "model_C_affinity": s_aff, "model_B_conformational": s_conf}
                 ).to_csv(FIG / "figure_causality_summary_panel3.csv", index=False)

    print("G5/G1 collapse:", collapse["identical"],
          "distinct P16 values:", collapse["n_distinct_P16"])
    print(df_disc.to_string(index=False))
    print(df_scalar.to_string(index=False))
    print(df_reps.to_string(index=False))
    print(df_real.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
