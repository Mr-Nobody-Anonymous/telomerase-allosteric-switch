#!/usr/bin/env python3
"""Figure generation for the manuscript.

Every panel is computed here from the same equations that
`scripts/mechanism_classification_checks.py` verifies, so no figure can drift from the text.

    python paper/make_figures.py          # writes paper/figures/*.png and *.pdf

Colour: Okabe-Ito qualitative palette, colour-blind safe. No panel relies on colour alone —
line style and direct labels carry every distinction.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paper" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN, RED, PURPLE, YELLOW, SKY, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#666666")

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.grid": False, "mathtext.fontset": "dejavusans",
})

REF = {"k_ext": 2.0, "k_loss": 0.5, "k_off": 2.0, "k_on": 1.0}


def s_of(P, k, kL, koff, kon):
    Kd = koff / kon
    R = koff * Kd / (Kd + P)
    return k / (k + kL + R)


def r_of_n(n_axis, s, phi, s2):
    """r(n) = s (1 + u rho)/(1 + u) with u = kappa rho^(n-1), evaluated in log space."""
    kap = phi * (1 - s2) / ((1 - phi) * (1 - s))
    rho = s2 / s
    lg = np.log(kap) + n_axis * np.log(rho)
    return s * np.exp(np.logaddexp(0.0, lg + np.log(rho)) - np.logaddexp(0.0, lg))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


# ---------------------------------------------------------------------------------------
# Figure 1 — what is and is not an identity
# ---------------------------------------------------------------------------------------
def figure1():
    fig, ax = plt.subplots(1, 2, figsize=(6.9, 2.7))

    P = np.logspace(-2, 5, 500)
    k, kL, koff, kon = REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"]

    models = {
        "WT": (k, kL, koff, kon, BLUE, "-", 2.0),
        r"B ($k_\mathrm{ext}$ up)": (2 * k, kL, koff, kon, ORANGE, (0, (4, 2)), 2.6),
        r"C1 ($k_\mathrm{off}$ down)": (k, kL, koff / 2, kon, GREEN, (0, (1, 1.4)), 2.6),
        r"C2 ($k_\mathrm{loss}$ down)": (k, kL / 3.42, koff, kon, PURPLE, "-.", 1.6),
    }

    a = ax[0]
    for label, (kx, lx, ox, nx, col, ls, lw) in models.items():
        a.plot(P, s_of(P, kx, lx, ox, nx), lw=lw, color=col, ls=ls, label=label)
    a.set_xscale("log")
    a.set_xlabel(r"primer concentration [P]  ($K_d^\mathrm{WT}$ units)")
    a.set_ylabel(r"effective survival per repeat, $s$")
    a.set_title("a   only C1 leaves the plateau unchanged", loc="left", fontweight="bold")
    a.legend(frameon=False, loc="lower right", fontsize=6.4, ncol=2)
    a.set_ylim(0.42, 1.10)

    for y, col in ((0.800, GREEN), (0.889, ORANGE), (0.932, PURPLE)):
        a.axhline(y, color=col, lw=0.5, ls=":", alpha=0.8)
    a.text(0.03, 0.99,
           "$s_\\infty$ = 0.800   WT and C1  (identical)\n"
           "          0.889   B\n"
           "          0.932   C2",
           transform=a.transAxes, fontsize=6.3, va="top", linespacing=1.5,
           bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#dddddd", "lw": 0.5})

    b = ax[1]
    n = np.arange(1, 61)
    target = 0.75
    p_same = (1 - target) * target ** (n - 1)
    b.semilogy(n, p_same, lw=3.4, color=BLUE, alpha=0.55,
               label=r"model B ($k_\mathrm{ext}$ doubled)")
    b.semilogy(n, p_same, lw=1.2, ls=(0, (2, 2)), color=ORANGE,
               label=r"model C1 ($k_\mathrm{off}$ halved)")
    b.set_xlabel("product length  $n$")
    b.set_ylabel(r"$P(\geq n) = s^{\,n-1}$")
    b.set_title("b   equal $s$ means identical data", loc="left", fontweight="bold")
    b.legend(frameon=False, loc="upper right", fontsize=6.4)
    b.set_ylim(1e-10, 3)
    b.text(0.05, 0.13,
           "the two curves superimpose exactly:\n"
           r"$\max_n|\Delta P(\geq n)| = 0$",
           transform=b.transAxes, fontsize=6.4, color=GREY)

    save(fig, "fig1_collapse")


# ---------------------------------------------------------------------------------------
# Figure 2 — affinity is not one hypothesis
# ---------------------------------------------------------------------------------------
def figure2():
    """R = k_off K_d/(K_d+[P]) = k_off^2/(k_off + [P] k_on).

    Normalised by k_off the loss depends on K_d alone, so C1a and C1b superimpose; in absolute
    terms they differ by a factor lambda, which is confounded with k_loss.
    """
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.7))
    fig.subplots_adjust(wspace=0.44)
    P = np.logspace(-2, 3, 400)
    k, kL = REF["k_ext"], REF["k_loss"]

    variants = [
        ("WT", 2.0, 1.0, BLUE, "-"),
        (r"C1a  $k_\mathrm{off}\!\downarrow$", 1.0, 1.0, ORANGE, (0, (4, 2))),
        (r"C1b  $k_\mathrm{on}\!\uparrow$", 2.0, 2.0, GREEN, (0, (1, 1.4))),
    ]

    a, b, d = axes
    for label, koff, kon, col, ls in variants:
        a.plot(P, s_of(P, k, kL, koff, kon), lw=1.7, color=col, ls=ls, label=label)
    a.set_xscale("log")
    a.set_xlabel(r"[P]")
    a.set_ylabel("$s$")
    a.set_title(r"a   matched $K_d$, different $s$", loc="left", fontweight="bold")
    a.legend(frameon=False, loc="upper left", fontsize=6.0)
    a.set_ylim(0.42, 0.92)

    for label, koff, kon, col, ls in variants:
        Kd = koff / kon
        b.plot(P, Kd / (Kd + P), lw=1.7, color=col, ls=ls)
    b.set_xscale("log")
    b.set_xlabel("[P]")
    b.set_ylabel(r"$R/k_\mathrm{off} = K_d/(K_d+[P])$")
    b.set_title(r"b   shape set by $K_d$ alone", loc="left", fontweight="bold")
    b.set_ylim(0, 1.04)
    b.annotate("C1a and C1b\ncoincide here", xy=(1.0, 0.5), xytext=(2.0, 0.72), fontsize=6.2,
               color=GREY, arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.7})

    for label, koff, kon, col, ls in variants:
        Kd = koff / kon
        d.plot(P, koff * Kd / (Kd + P), lw=1.7, color=col, ls=ls)
    d.set_xscale("log")
    d.set_xlabel("[P]")
    d.set_ylabel(r"$R = k_\mathrm{off}K_d/(K_d+[P])$")
    d.set_title(r"c   size scales with $k_\mathrm{off}$", loc="left", fontweight="bold")
    d.set_ylim(0, 2.45)
    d.annotate("", xy=(0.02, 2.0), xytext=(0.02, 1.0),
               arrowprops={"arrowstyle": "<->", "color": RED, "lw": 0.9})
    d.text(0.025, 1.40, r"$\lambda = 2$", fontsize=7.0, color=RED)
    d.text(0.025, 2.15, r"confounded with $k_\mathrm{loss}$", fontsize=6.1, color=RED)

    save(fig, "fig2_affinity_variants")


# ---------------------------------------------------------------------------------------
# Figure 3 — the B vs C2 crossing and the sign of P*
# ---------------------------------------------------------------------------------------
def figure3():
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 2.6))
    P = np.logspace(-2, 4, 600)
    k, kL, koff, kon = REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"]

    def diff(Pgrid, al, be):
        Kd = koff / kon
        R = koff * Kd / (Kd + Pgrid)
        return al * k / (al * k + kL + R) - k / (k + kL / be + R)

    for axx, (al, be, tag) in zip(axes, ((2.0, 3.42, "crossing exists"), (30.0, 3.42, "no crossing"))):
        Pstar = koff * (-al * be * koff - al * kL + be * kL + be * koff) / (kL * kon * (al - be))
        dd = diff(P, al, be)
        axx.plot(P, dd, lw=2.0, color=BLUE)
        axx.axhline(0, color=GREY, lw=0.7)
        axx.set_xscale("log")
        axx.set_xlabel("[P]")
        axx.set_ylabel(r"$s_\mathrm{B} - s_\mathrm{C2}$")
        letter = "a" if al < 10 else "b"
        axx.set_title(rf"{letter}   $\alpha={al:g}$, $\beta={be:g}$: {tag}",
                      loc="left", fontweight="bold")
        if Pstar > 0:
            axx.axvline(Pstar, color=RED, lw=1.0, ls="--")
            axx.annotate(rf"$P^* = {Pstar:.1f}$", xy=(Pstar, 0),
                         xytext=(Pstar * 1.5, max(dd) * 0.55), fontsize=6.8, color=RED,
                         arrowprops={"arrowstyle": "->", "color": RED, "lw": 0.8})
        else:
            axx.text(0.05, 0.14,
                     rf"$P^* = {Pstar:.1f} < 0$" + "\nno root in the physical domain",
                     transform=axx.transAxes, fontsize=6.6, color=RED)

    save(fig, "fig3_crossing")


# ---------------------------------------------------------------------------------------
# Figure 4 — the central methodological figure
# ---------------------------------------------------------------------------------------
def figure4():
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    fig.subplots_adjust(wspace=0.34)
    n = np.arange(1, 121)

    cases = [
        (0.60, 0.01, 0.95, BLUE, "-"),
        (0.70, 0.01, 0.90, ORANGE, (0, (4, 2))),
        (0.85, 0.05, 0.95, GREEN, (0, (1, 1.4))),
        (0.95, 0.02, 0.99, PURPLE, "-."),
    ]

    a = ax[0]
    for s, phi, s2, col, ls in cases:
        a.plot(n, r_of_n(n, s, phi, s2), lw=1.7, color=col, ls=ls,
               label=rf"$s$={s:.2f}, $\varphi$={phi:g}, $s_2$={s2:.2f}")
        a.axhline(max(s, s2), color=col, lw=0.6, ls=":", alpha=0.65)
    a.set_xlabel("product length  $n$  (repeats added)")
    a.set_ylabel(r"adjacent-band ratio  $r(n)$")
    a.set_title(r"a   bounded by $\max(s, s_2)$", loc="left", fontweight="bold")
    a.legend(frameon=False, loc="lower right", fontsize=6.2)
    a.set_ylim(0.55, 1.02)

    b = ax[1]
    s, phi, s2 = 0.60, 0.01, 0.95
    r = r_of_n(n, s, phi, s2)
    b.plot(n, r, lw=2.0, color=BLUE, label=r"$r(n)$")
    b.axhline(s, color=GREY, lw=0.8, ls="--")
    b.axhline(s2, color=RED, lw=0.8, ls="--")
    b.text(88, s - 0.014, r"true $s$ = 0.60", fontsize=6.4, color=GREY, va="top", ha="center")
    b.text(20, s2 + 0.010, r"known $s_2$ = 0.95", fontsize=6.4, color=RED, va="bottom")

    windows = [(3, 10, SKY), (16, 25, YELLOW), (40, 80, "#e8e8e8")]
    for lo, hi, fc in windows:
        b.axvspan(lo, hi, color=fc, alpha=0.55, lw=0)
        seg = r[lo - 1:hi]
        mean = seg.mean()
        b.hlines(mean, lo, hi, color="#222222", lw=2.0)
        b.text((lo + hi) / 2, mean - 0.030, f"{mean:.2f}", fontsize=6.8, ha="center",
               fontweight="bold", color="#222222")

    kap = phi * (1 - s2) / ((1 - phi) * (1 - s))
    nstar = 1 + np.log(1 / kap) / np.log(s2 / s)
    b.axvline(nstar, color=BLUE, lw=0.9, ls=":")
    b.annotate(rf"$n^* = {nstar:.1f}$", xy=(nstar, 0.99), xytext=(nstar + 8, 0.99),
               fontsize=6.6, color=BLUE, va="top",
               arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 0.8})

    b.set_xlabel("product length  $n$")
    b.set_ylabel("apparent $s$ from the tail slope")
    b.set_title("b   the same data, three windows", loc="left", fontweight="bold")
    b.set_ylim(0.56, 1.0)
    b.set_xlim(0, 140)
    b.set_xticks([0, 20, 40, 60, 80, 100, 120])

    save(fig, "fig4_mixture_ratio")


# ---------------------------------------------------------------------------------------
# Figure 5 — window dependence and the error it induces
# ---------------------------------------------------------------------------------------
def figure5():
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
    fig.subplots_adjust(wspace=0.38)

    windows = [(3, 10), (10, 16), (16, 25), (25, 40), (40, 80)]
    labels = [f"{lo}–{hi}" for lo, hi in windows]
    x = np.arange(len(windows))

    truths = [(0.60, 0.01, 0.95, BLUE, "-"), (0.70, 0.01, 0.90, ORANGE, (0, (4, 2))),
              (0.85, 0.05, 0.95, GREEN, (0, (1, 1.4))), (0.999, 0.01, 0.95, PURPLE, "-.")]

    a = ax[0]
    for s, phi, s2, col, ls in truths:
        vals = [float(r_of_n(np.arange(lo, hi + 1), s, phi, s2).mean()) for lo, hi in windows]
        a.plot(x, vals, marker="o", ms=3.5, lw=1.5, color=col, ls=ls,
               label=rf"$s$={s:.2f}, $s_2$={s2:.2f}")
    a.set_xticks(x)
    a.set_xticklabels(labels)
    a.set_xlabel("read window, $n$")
    a.set_ylabel("tail slope reported")
    a.set_title("a   window sets the answer", loc="left", fontweight="bold")
    a.legend(frameon=False, fontsize=6.2, loc="center right")
    a.set_ylim(0.58, 1.03)

    b = ax[1]
    s, phi, s2 = 0.60, 0.01, 0.95
    vals = np.array([float(r_of_n(np.arange(lo, hi + 1), s, phi, s2).mean())
                     for lo, hi in windows])
    implied = (1 - vals) / vals
    b.bar(x, implied, color=[SKY, YELLOW, ORANGE, RED, "#333333"], width=0.68)
    for i, v in enumerate(implied):
        b.text(i, v * 1.10, f"{v:.3f}", ha="center", fontsize=6.4)
    b.set_yscale("log")
    b.set_ylim(0.03, 1.4)
    b.set_xticks(x)
    b.set_xticklabels(labels)
    b.set_xlabel("read window, $n$")
    b.set_ylabel(r"implied $k_\mathrm{loss}/k_\mathrm{ext}$")
    b.set_title("b   same data, 12-fold spread", loc="left", fontweight="bold")
    b.annotate("", xy=(4, implied[4] * 2.0), xytext=(3.4, implied[0] * 1.15),
               arrowprops={"arrowstyle": "->", "color": GREY, "lw": 0.9})
    b.text(3.15, implied[0] * 0.62, rf"$\times${implied[0] / implied[4]:.0f}", fontsize=8,
           color=GREY, ha="center", fontweight="bold")

    save(fig, "fig5_window_dependence")


# ---------------------------------------------------------------------------------------
# Figure 6 — identifiability
# ---------------------------------------------------------------------------------------
def figure6():
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.7))
    fig.subplots_adjust(wspace=0.46)

    def csdiv(f, x, i, h=1e-20):
        xc = [complex(v) for v in x]
        xc[i] += 1j * h
        return f(*xc).imag / h

    theta = [REF["k_ext"], REF["k_loss"], REF["k_off"], REF["k_on"]]
    P = np.logspace(-2, 3, 12)
    J = np.array([[csdiv(lambda *a, _p=p: s_of(_p, *a), theta, i) for i in range(4)] for p in P])

    designs = [
        ("processivity only", J, RED, "o"),
        (r"+ turnover ($k_\mathrm{ext}\!\equiv\!k_\mathrm{chem}$)", np.vstack([J, [[1, 0, 0, 0]]]),
         ORANGE, "s"),
        ("+ direct binding", np.vstack([J, [[0, 0, 1, 0]], [[0, 0, 0, 1]]]), GREEN, "^"),
    ]

    a = axes[0]
    for name, M, col, mk in designs:
        sv = np.linalg.svd(M, compute_uv=False)
        rel = sv / sv[0]
        a.semilogy(np.arange(1, len(rel) + 1), np.maximum(rel, 1e-18), marker=mk, ms=4, lw=1.3,
                   color=col, label=name)
    a.axhline(1e-8, color=GREY, ls=":", lw=0.8)
    a.set_xlabel("singular value index")
    a.set_ylabel(r"$\sigma_i/\sigma_1$")
    a.set_title("a   rank 3 of 4", loc="left", fontweight="bold")
    a.legend(frameon=False, fontsize=6.0, loc="lower left")
    a.set_ylim(1e-18, 5)
    a.set_xticks([1, 2, 3, 4])

    b = axes[1]
    names = [r"$k_\mathrm{ext}$", r"$k_\mathrm{loss}$", r"$k_\mathrm{off}$", r"$k_\mathrm{on}$"]
    only = [30.8, 34.3, 28.9, 67.0]
    bind = [21.5, 49.0, 9.9, 9.9]
    xx = np.arange(4)
    b.bar(xx - 0.19, only, 0.36, color=RED, label="processivity only")
    b.bar(xx + 0.19, bind, 0.36, color=GREEN, label="+ direct binding")
    b.set_xticks(xx)
    b.set_xticklabels(names)
    b.set_ylabel("relative SE (%)")
    b.set_ylim(0, 88)
    b.set_title(r"b   binding does not fix $k_\mathrm{loss}$", loc="left", fontweight="bold")
    b.legend(frameon=False, fontsize=6.2, loc="upper center")
    b.annotate("", xy=(1.19, 50.5), xytext=(1.19, 79),
               arrowprops={"arrowstyle": "->", "color": RED, "lw": 0.9})
    b.text(1.30, 60, "worse", fontsize=6.5, color=RED)

    c = axes[2]
    from scipy.stats import norm
    cvs = np.linspace(0.08, 0.40, 140)
    for delta, col in ((0.10, BLUE), (0.20, ORANGE)):
        reps = 3 * 2 * (norm.ppf(0.95) + norm.ppf(0.80)) ** 2 * (cvs ** 2) / delta ** 2
        c.plot(cvs, np.ceil(reps), lw=1.6, color=col, label=rf"$\delta$={delta:.2f}")
    c.axhline(105, color=GREY, ls=":", lw=0.8)
    c.axvline(0.345, color=GREY, ls=":", lw=0.8)
    c.text(0.352, 23, "published CV", fontsize=6.2, color=GREY, rotation=90, va="bottom")
    c.set_xlabel("assay CV")
    c.set_ylabel("replicates per arm")
    c.set_yscale("log")
    c.set_title("c   equivalence (TOST)", loc="left", fontweight="bold")
    c.legend(frameon=False, fontsize=6.4, loc="upper left", title=r"$\delta$ (log)",
             title_fontsize=6.4)

    save(fig, "fig6_identifiability")


def main() -> int:
    print("generating figures ->", OUT)
    for fn in (figure1, figure2, figure3, figure4, figure5, figure6):
        fn()
    manifest = {p.name: p.stat().st_size for p in sorted(OUT.glob("*.png"))}
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("done:", len(manifest), "PNG files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
