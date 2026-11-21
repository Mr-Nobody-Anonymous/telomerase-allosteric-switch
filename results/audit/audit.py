"""Consolidated no-lab audit for the telomerase-allosteric-switch data.

Four analyses that require no laboratory, no GPU and no network:

  A. Structural distance audit of the K757 / Tyr176-Gln794 claim, read from
     the deposited mmCIF files in data/raw/.
  B. Paired re-analysis of the Padmanaban, Tesmer & Nandakumar 2023 (LSA)
     deposited Fig 3B/3C panel, with Holm correction and an unpaired
     sensitivity check.
  C. Power / sample-size calculation from the CVs measured in (B).
  D. Two-species mixture recovery study, i.e. a test of our own estimator.

Every number that appears in results/audit/AUDIT.md is produced by this file.

Usage
-----
    python results/audit/audit.py --selftest
    python results/audit/audit.py            # run from the repository root

Requires numpy, scipy and gemmi (gemmi only for section A).

Scope and limits
----------------
* Section A reports distances between atoms as deposited. It says nothing
  about whether the cryo-EM density supports those side chains; at 3.5 A a
  side-chain pair cannot be placed to better than a few tenths of an
  Angstrom, so a 2.34 A contact is a property of the model, not a
  measurement of the protein.
* Section B re-analyses deposited *summary* numbers with n = 3. It cannot
  overturn the authors, who hold the raw gels. Pairing across rows is
  assumed from the file layout and is tested as a sensitivity check.
* Section D is synthetic. It tests the estimator against a known truth; it
  is not evidence about any real gel.
"""

import json
import math
import sys

import numpy as np
from scipy import stats

ROOT = "."
T2 = stats.t.ppf(0.975, 2)


# --------------------------------------------------------------------------
# A. Structural distance audit
# --------------------------------------------------------------------------
def structure_audit(pdbs=("7TRE", "7TRD", "7TRF"), chain="A"):
    """Distances between the atoms that carry the K757 claim, as deposited."""
    import gemmi

    out = {}
    for pid in pdbs:
        structure = gemmi.read_structure(f"{ROOT}/data/raw/{pid}.cif")
        atoms = {
            (ch.name, res.seqid.num, atom.name): atom.pos
            for ch in structure[0]
            for res in ch
            for atom in res
        }

        def dist(num_a, name_a, num_b, name_b):
            p = atoms.get((chain, num_a, name_a))
            q = atoms.get((chain, num_b, name_b))
            return None if p is None or q is None else p.dist(q)

        lys_nz = atoms.get((chain, 757, "NZ"))
        nearest = []
        if lys_nz is not None:
            nearest = sorted(
                (lys_nz.dist(pos), f"{ch}{num}:{name}")
                for (ch, num, name), pos in atoms.items()
                if not (ch == chain and num == 757)
            )[:3]

        out[pid] = {
            "Y176_OH--Q794_OE1": dist(176, "OH", 794, "OE1"),
            "K757_NZ--Q794_OE1": dist(757, "NZ", 794, "OE1"),
            "K757_NZ--Y176_OH": dist(757, "NZ", 176, "OH"),
            "nearest_to_K757_NZ": [(round(d, 2), label) for d, label in nearest],
        }
    return out


# --------------------------------------------------------------------------
# B. Paired re-analysis of the deposited mutant panel
# --------------------------------------------------------------------------
# Padmanaban S, Tesmer VM, Nandakumar J. Life Sci Alliance 6:e202201727 (2023).
# Deposited source data, Fig 3B (processivity, +PT lanes) and Fig 3C (activity,
# WT-normalised within each experiment). Three rows = three experiments.
PROC = {
    "WT": [0.50879, 0.554671, 0.39286],
    "K749E": [0.369453, 0.325139, 0.142375],
    "R756E": [0.320793, 0.34624, 0.195215],
    "K757E": [0.203158, 0.179349, 0.207296],
    "F759E": [0.199718, 0.094687, 0.056956],
    "K760E": [0.269501, 0.301299, 0.220829],
    "H762E": [0.499718, 0.458955, 0.287844],
    "F776E": [0.229079, 0.146742, -0.33789],
    "L786E": [0.239126, 0.1406, 0.133442],
    "D788R": [0.283557, 0.46741, 0.186042],
}
ACT = {
    "WT": [1, 1, 1],
    "K749E": [1.132561, 0.147964, 0.420697],
    "R756E": [1.004672, 0.123393, 0.270876],
    "K757E": [0.674734, 0.105502, 0.501341],
    "F759E": [0.695901, 0.128053, 0.335026],
    "K760E": [0.539852, 0.274883, 0.182804],
    "H762E": [1.387196, 0.220485, 0.432991],
    "F776E": [0.976458, 0.11072, 0.168121],
    "L786E": [0.573581, 0.042355, 0.128089],
    "D788R": [0.606591, 0.065442, 0.031861],
}


def paired_ratio(x, w):
    """Geometric-mean ratio to WT, paired by experiment; CI on the log scale."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    log_ratio = np.log(x) - np.log(w)
    mean = log_ratio.mean()
    se = log_ratio.std(ddof=1) / math.sqrt(len(log_ratio))
    p = stats.ttest_rel(np.log(x), np.log(w)).pvalue
    return math.exp(mean), math.exp(mean - T2 * se), math.exp(mean + T2 * se), float(p)


def unpaired_ratio(x, w):
    """Sensitivity check: Welch on logs, so pairing is not assumed."""
    lx, lw = np.log(np.asarray(x, float)), np.log(np.asarray(w, float))
    mean = lx.mean() - lw.mean()
    se = math.sqrt(lx.var(ddof=1) / len(lx) + lw.var(ddof=1) / len(lw))
    p = stats.ttest_ind(lx, lw, equal_var=False).pvalue
    return math.exp(mean), math.exp(mean - T2 * se), math.exp(mean + T2 * se), float(p)


def holm(pvals, alpha=0.05):
    """Holm-Bonferroni step-down. Returns [(threshold, rejected), ...]."""
    idx = [i for i, p in enumerate(pvals) if p is not None and not math.isnan(p)]
    order = sorted(idx, key=lambda i: pvals[i])
    m = len(order)
    out = [(None, False)] * len(pvals)
    still_rejecting = True
    for k, i in enumerate(order):
        threshold = alpha / (m - k)
        still_rejecting = still_rejecting and pvals[i] <= threshold
        out[i] = (threshold, still_rejecting)
    return out


def reanalysis():
    res = {"processivity": {}, "activity": {}}

    ps = []
    for mut in PROC:
        if mut == "WT":
            continue
        if min(PROC[mut]) <= 0:
            res["processivity"][mut] = {
                "valid": False,
                "class": "unusable",
                "reason": "non-positive replicate; a fraction of signal cannot be negative",
            }
            ps.append(None)
            continue
        r, lo, hi, p = paired_ratio(PROC[mut], PROC["WT"])
        res["processivity"][mut] = {
            "valid": True,
            "ratio": r,
            "lo": lo,
            "hi": hi,
            "p": p,
            "unpaired": unpaired_ratio(PROC[mut], PROC["WT"]),
        }
        ps.append(p)
    for mut, (thr, rejected) in zip([m for m in PROC if m != "WT"], holm(ps)):
        if res["processivity"][mut]["valid"]:
            res["processivity"][mut].update(holm_threshold=thr, holm=bool(rejected))

    ps2 = []
    for mut in ACT:
        if mut == "WT":
            continue
        logs = np.log(np.array(ACT[mut], float))  # already a ratio to WT
        se = logs.std(ddof=1) / math.sqrt(3)
        res["activity"][mut] = {
            "ratio": math.exp(logs.mean()),
            "lo": math.exp(logs.mean() - T2 * se),
            "hi": math.exp(logs.mean() + T2 * se),
            "p": float(stats.ttest_1samp(logs, 0.0).pvalue),
        }
        ps2.append(res["activity"][mut]["p"])
    for mut, (thr, rejected) in zip([m for m in ACT if m != "WT"], holm(ps2)):
        res["activity"][mut].update(holm_threshold=thr, holm=bool(rejected))

    for endpoint in res:
        for entry in res[endpoint].values():
            if not entry.get("valid", True):
                continue
            if entry["holm"] and entry["hi"] < 1:
                entry["class"] = "reduced (strong)"
            elif entry["hi"] < 1:
                entry["class"] = "reduced (uncorrected only)"
            elif entry["lo"] > 1:
                entry["class"] = "increased"
            else:
                entry["class"] = "not resolvable"
    return res


# --------------------------------------------------------------------------
# C. Power / sample size
# --------------------------------------------------------------------------
CVS = {
    "processivity median": 0.35,
    "processivity K757E": 0.08,
    "activity median": 0.74,
    "activity D788R": 1.37,
}


def power(n, fold, cv, alpha, paired=True, sims=3000, seed=11):
    rng = np.random.default_rng(seed)
    sd = math.sqrt(math.log(1 + cv**2))
    if paired:
        y = rng.normal(math.log(fold), sd / math.sqrt(n), size=sims)
        return float((np.abs(y / (sd / math.sqrt(n))) > stats.t.ppf(1 - alpha / 2, n - 1)).mean())
    a = rng.normal(0, sd, (sims, n))
    b = rng.normal(math.log(fold), sd, (sims, n))
    t = (b.mean(1) - a.mean(1)) / np.sqrt(a.var(1, ddof=1) / n + b.var(1, ddof=1) / n)
    return float((np.abs(t) > stats.t.ppf(1 - alpha / 2, 2 * n - 2)).mean())


def min_n(fold, cv, alpha, paired, target=0.80, nmax=40):
    return next((n for n in range(2, nmax + 1) if power(n, fold, cv, alpha, paired) >= target), None)


def power_table():
    out = {}
    for name, cv in CVS.items():
        out[name] = {
            f"{fold}x": {
                "paired": min_n(fold, cv, 0.05, True),
                "unpaired": min_n(fold, cv, 0.05, False),
                "paired_holm9": min_n(fold, cv, 0.05 / 9, True),
            }
            for fold in (1.25, 1.5, 2.0, 3.0)
        }
        smallest = next((f for f in np.arange(1.05, 6.001, 0.05) if power(3, f, cv, 0.05) >= 0.80), None)
        out[name]["detectable_at_n3"] = None if smallest is None else round(float(smallest), 2)
    return out


def tost_n(cv, margin_log, power_=0.90):
    """Preparations for a one-sample TOST at the given log margin."""
    sd = math.sqrt(math.log(1 + cv**2))
    return math.ceil(((stats.norm.ppf(0.975) + stats.norm.ppf(power_)) * sd / margin_log) ** 2)


# --------------------------------------------------------------------------
# D. Mixture recovery study (synthetic)
# --------------------------------------------------------------------------
def mixture(n, s, s2, phi):
    return (1 - phi) * (1 - s) * s ** (n - 1) + phi * (1 - s2) * s2 ** (n - 1)


def recovery(reps=60, seed=7, n_max=80, true_s=0.60):
    rng = np.random.default_rng(seed)
    n = np.arange(1, n_max + 1)
    rows = []
    for cv in (0.15, 0.35):
        sigma = math.sqrt(math.log(1 + cv**2))
        for phi in (0.005, 0.01, 0.02, 0.05):
            for s2 in (0.90, 0.95, 0.99):
                intensity = mixture(n, true_s, s2, phi)
                keep = intensity > 1e-3 * intensity.max()
                nk, ik = n[keep], intensity[keep]
                obs = ik[None, :] * np.exp(rng.normal(0, sigma, (reps, keep.sum())))
                windows = {}
                for low, high in ((3, 10), (10, 30), (30, 60)):
                    mask = (nk >= low) & (nk <= high)
                    if mask.sum() < 3:
                        windows[f"{low}-{high}"] = None
                        continue
                    x, y = nk[mask], np.log(obs[:, mask])
                    xm = x.mean()
                    slope = ((x - xm) * (y - y.mean(1, keepdims=True))).sum(1) / ((x - xm) ** 2).sum()
                    windows[f"{low}-{high}"] = float(np.exp(slope).mean())
                rows.append(
                    {"cv": cv, "phi": phi, "s2": s2, "windows": windows, "true_majority_s": true_s}
                )
    return rows


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------
def selftest():
    r, lo, hi, _ = paired_ratio([1, 2, 3], [1, 2, 3])
    assert abs(r - 1) < 1e-12 and abs(lo - 1) < 1e-9, "an identity pair must give a ratio of 1"

    res = reanalysis()
    for endpoint in ("processivity", "activity"):
        for mut, entry in res[endpoint].items():
            if entry.get("valid", True):
                assert entry["lo"] > 0, f"a back-transformed log CI must be positive: {endpoint}/{mut}"
            if entry.get("holm"):
                assert entry["p"] <= 0.05

    surviving = [m for m, e in res["processivity"].items() if e.get("holm")]
    assert surviving == ["K760E"], surviving
    assert res["processivity"]["F776E"]["class"] == "unusable"

    for s, s2, phi in ((0.6, 0.95, 0.05), (0.6, 0.9, 0.01), (0.8, 0.99, 0.02)):
        ratios = [mixture(k, s, s2, phi) / mixture(k - 1, s, s2, phi) for k in range(2, 60)]
        assert min(ratios) >= min(s, s2) - 1e-9
        assert max(ratios) <= max(s, s2) + 1e-9

    assert tost_n(0.74, 0.405) > tost_n(0.35, 0.405), "a noisier assay must need more preparations"
    print("selftest: 6/6 invariants passed")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
        sys.exit()

    results = {
        "structure": structure_audit(),
        "reanalysis": reanalysis(),
        "power": power_table(),
        "recovery": recovery(),
    }
    results["tost"] = {
        f"cv{int(c * 100)}_margin{int(100 * (1 - math.exp(-m)))}%": tost_n(c, m)
        for c in (0.35, 0.74)
        for m in (0.405, 0.2231)
    }

    with open("results/audit/audit_results.json", "w") as handle:
        json.dump(results, handle, indent=2, default=str)

    print("=== A. structure ===")
    for pid, entry in results["structure"].items():
        rounded = {k: (round(v, 2) if isinstance(v, float) else v) for k, v in entry.items()}
        print(f"  {pid}: {rounded}")
    print("=== B. reanalysis ===")
    for endpoint in ("processivity", "activity"):
        print(f"  -- {endpoint}")
        for mut, e in results["reanalysis"][endpoint].items():
            if not e.get("valid", True):
                print(f"     {mut:7s} unusable")
                continue
            print(
                f"     {mut:7s} {e['ratio']:5.2f} [{e['lo']:5.2f},{e['hi']:5.2f}] "
                f"p={e['p']:.3f} holm={e['holm']!s:5s} {e['class']}"
            )
    print("=== C. power ===")
    for name, entry in results["power"].items():
        print(
            f"  {name:22s} n(1.5x): paired={entry['1.5x']['paired']} "
            f"unpaired={entry['1.5x']['unpaired']} holm={entry['1.5x']['paired_holm9']} "
            f"| detectable at n=3: {entry['detectable_at_n3']}"
        )
    print(f"  TOST preparations: {results['tost']}")
    print("=== D. recovery ===")
    for row in results["recovery"]:
        w = row["windows"]
        fmt = lambda v: "  n/a" if v is None else f"{v:6.3f}"
        print(
            f"  phi={row['phi']:5.3f} s2={row['s2']:4.2f} cv={row['cv']:4.2f} | "
            f"{fmt(w['3-10'])} {fmt(w['10-30'])} {fmt(w['30-60'])}"
        )
    print("\nwrote results/audit/audit_results.json")
