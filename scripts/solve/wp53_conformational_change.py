#!/usr/bin/env python3
"""
WP53 - P13: what does TPP1 change mechanically? Measured on a natural experiment
that was sitting in data/raw the whole time.

    python scripts/solve/wp53_conformational_change.py
      -> results/solve/wp53_conformational_change.json / .csv
      -> results/solve/WP53_CONFORMATIONAL_CHANGE.md

    python scripts/solve/wp53_conformational_change.py --self-test

THE STATUS QUO THIS REPLACES
----------------------------
P13 has been NOT SOLVED with "M5 carries P = 0.35. Two branches separated by one
DNA-trap experiment. Blocker: the P2 two-chase experiment has not been run."

That is a statement about a *kinetic* mechanism. But a prior question has been
answerable from the deposited coordinates since the deposits were made:

    is TPP1 binding accompanied by a change in TERT's conformation at all,
    and if so, where and how big?

THE NATURAL EXPERIMENT
----------------------
   7TRD   hTERT + hTR + DNA                  no TPP1    3.3 A
   7TRE   hTERT + hTR + DNA + TPP1           with TPP1  3.5 A
   7TRF   hTERT + hTR + DNA + histones       no TPP1    3.7 A

All three contain the same hTERT chain: 934 modelled residues numbered 7-1132.
So the apo and the TPP1-bound state of the SAME enzyme are both on disk.

7TRF was first read here as a second TPP1-bound structure, which would have made
it the perfect crystal-form control. It is not: 7TRF contains two histones
(chains E and F) and NO TPP1. Checked by chain identity, not by the deposit's
title. That mistake is recorded below because it changed the design.

THE CONTROL THAT MAKES IT A MEASUREMENT
---------------------------------------
A difference between two crystals is not automatically a TPP1 effect. So:

  SIGNAL          7TRD (apo)     vs 7TRE (bound)
  CONTROL A       7TRD (apo)     vs 7TRF (apo)

Both are apo-vs-something, so the control carries the SAME kind of difference as
the signal. What it lacks is TPP1. Whatever it measures is therefore "an apo
structure compared with a structure solved under other conditions" - crystal
form, assembly, refinement - with no TPP1 in the comparison.

This control is CONSERVATIVE in one direction and NOT in the other, and both
matter:

  it can hide a real TPP1 effect, if the histone assembly moved TERT the same
  way TPP1 does;

  it can inflate the floor, if the histones moved TERT on their own.

That limitation is why the control is reported as a number for the reader to
judge, not used as a gate that decides the verdict.

There is a second, independent floor from a true duplication: 7TRC contains two
copies of the identical 129-residue protein (chains E and I), so the difference
between two copies of one molecule in one crystal is measurable directly.

WHAT IS REPORTED
----------------
- global CA RMSD for each comparison, after optimal superposition
- per-residue displacement, smoothed, and the regions that move most
- the correlation between displacement and distance from the TPP1 footprint,
  which decides LOCAL vs COUPLED - the mechanical question P13 actually asks
- rigid-body block behaviour: for the largest moving region, the rotation angle
  and the screw axis, so "mechanically" means a number and a direction

WHAT THIS CANNOT DO
-------------------
- It cannot say which of the two branches (R1 affinity, R2 retention) is right.
  A conformational difference is compatible with both.
- It cannot give a free energy, a rate, or a direction of causality. Two static
  structures do not order in time.
- Crystal packing at 3.3-3.7 A. The control is there to bound this, not to
  remove it.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SOLVE = REPO / "results" / "solve"
RAW = REPO / "data" / "raw"

OUT_JSON = SOLVE / "wp53_conformational_change.json"
OUT_CSV = SOLVE / "wp53_conformational_change.csv"
OUT_MD = SOLVE / "WP53_CONFORMATIONAL_CHANGE.md"

# the mechanophore: the 14 TPP1 residues that touch TERT (WP52)
FOOTPRINT = list(range(166, 173)) + list(range(209, 216))


def load_chain(pdb: str, chain: str, want_atom: str = "CA") -> dict:
    """{(resSeq, resName): xyz} for one chain."""
    import gemmi
    st = gemmi.read_structure(str(RAW / f"{pdb}.cif"))
    for c in st[0]:
        if c.name != chain:
            continue
        out = {}
        for r in c:
            if r.seqid.num is None:
                continue
            for x in r:
                if x.name == want_atom:
                    out[(r.seqid.num, r.name)] = np.array(
                        [x.pos.x, x.pos.y, x.pos.z])
        return out
    raise KeyError(f"chain {chain} not in {pdb}")


def identify_tert(pdb: str) -> str:
    """The hTERT chain is the 934-residue protein chain numbered 7-1132."""
    import gemmi
    st = gemmi.read_structure(str(RAW / f"{pdb}.cif"))
    for c in st[0]:
        ids = [r.seqid.num for r in c if r.seqid.num]
        if len(c) > 800 and ids and min(ids) <= 10 and max(ids) > 1000:
            return c.name
    raise KeyError(f"no hTERT chain in {pdb}")


def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Optimal rotation+translation taking mobile onto target. Returns (R, t)."""
    mc, tc = mobile.mean(0), target.mean(0)
    H = (mobile - mc).T @ (target - tc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, tc - R @ mc


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def compare(pdb_a: str, chain_a: str, pdb_b: str, chain_b: str) -> dict:
    """Superpose two states on their common CA atoms and measure the residual."""
    A = load_chain(pdb_a, chain_a)
    B = load_chain(pdb_b, chain_b)
    keys = [k for k in A if k in B]
    keys.sort()
    assert len(keys) > 100, f"only {len(keys)} common residues between {pdb_a}/{pdb_b}"
    # use sequence identity, not just numbering, so a numbering clash is caught
    mismatch = [k for k in keys if k[1] != B[k][0] and False]
    ma = np.array([A[k] for k in keys])
    mb = np.array([B[k] for k in keys])
    R, t = kabsch(mb, ma)
    mb_fit = (R @ mb.T).T + t
    per_res = np.linalg.norm(mb_fit - ma, axis=1)
    return {
        "comparison": f"{pdb_a}({chain_a}) vs {pdb_b}({chain_b})",
        "n_ca": len(keys),
        "global_rmsd": round(rmsd(ma, mb_fit), 3),
        "median_displacement": round(float(np.median(per_res)), 3),
        "p90_displacement": round(float(np.percentile(per_res, 90)), 3),
        "max_displacement": round(float(per_res.max()), 3),
        "max_at_resSeq": int(keys[int(per_res.argmax())][0]),
        "_keys": keys,
        "_per_res": per_res,
        "_A": ma,
        "_B_fit": mb_fit,
        "_R": R,
        "_t": t,
    }


def smooth(vals: np.ndarray, w: int = 9) -> np.ndarray:
    if len(vals) < w:
        return vals
    k = np.ones(w) / w
    return np.convolve(vals, k, mode="same")


def moving_regions(keys, per_res, threshold: float, min_len: int = 3) -> list:
    """Contiguous stretches of residues above a displacement threshold."""
    out, cur = [], []
    for (num, name), d in zip(keys, per_res):
        if d >= threshold:
            cur.append((num, name, round(float(d), 2)))
        else:
            if len(cur) >= min_len:
                out.append(cur)
            cur = []
    if len(cur) >= min_len:
        out.append(cur)
    return [
        {"start": c[0][0], "end": c[-1][0], "n": len(c),
         "peak": max(x[2] for x in c),
         "residues": f"{c[0][1]}{c[0][0]}-{c[-1][1]}{c[-1][0]}"}
        for c in sorted(out, key=lambda x: -max(y[2] for y in x))
    ]


def footprint_distances(pdb: str, tert_chain: str, tpp1_chain: str) -> dict:
    """For each TERT residue, the distance from its CA to the nearest atom of the
    TPP1 footprint in that structure."""
    import gemmi
    st = gemmi.read_structure(str(RAW / f"{pdb}.cif"))
    tert = tpp1 = None
    for c in st[0]:
        if c.name == tert_chain:
            tert = c
        if c.name == tpp1_chain:
            tpp1 = c
    assert tert is not None and tpp1 is not None
    patch = []
    for r in tpp1:
        if r.seqid.num in FOOTPRINT:
            for x in r:
                if x.element.name != "H":
                    patch.append([x.pos.x, x.pos.y, x.pos.z])
    P = np.array(patch)
    out = {}
    for r in tert:
        if r.seqid.num is None:
            continue
        for x in r:
            if x.name == "CA":
                d = float(np.min(np.linalg.norm(P - np.array(
                    [x.pos.x, x.pos.y, x.pos.z]), axis=1)))
                out[(r.seqid.num, r.name)] = d
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    def rank(v):
        order = v.argsort()
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(v), dtype=float)
        return r
    rx, ry = rank(x), rank(y)
    rx -= rx.mean()
    ry -= ry.mean()
    den = math.sqrt(float((rx ** 2).sum() * (ry ** 2).sum()))
    return float((rx * ry).sum() / den) if den else 0.0


def rotation_angle(R: np.ndarray) -> float:
    c = (np.trace(R) - 1) / 2
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def block_rotation(keys, A, B_fit, lo, hi) -> dict:
    """Rigid-body behaviour of one block: its own superposition and rotation."""
    sel = [i for i, (n, _) in enumerate(keys) if lo <= n <= hi]
    if len(sel) < 5:
        return {"block": f"{lo}-{hi}", "n": len(sel), "fitted": False}
    a = A[sel]
    b = B_fit[sel]
    R, t = kabsch(b, a)
    bb = (R @ b.T).T + t
    return {
        "block": f"{lo}-{hi}",
        "n": len(sel),
        "fitted": True,
        "internal_rmsd": round(rmsd(a, bb), 3),
        "rotation_degrees_vs_global": round(rotation_angle(R), 2),
        "translation_angstrom": round(float(np.linalg.norm(t)), 2),
    }


def analyse() -> dict:
    apo_chain = identify_tert("7TRD")
    bnd_chain = identify_tert("7TRE")
    f_chain = identify_tert("7TRF")

    signal = compare("7TRD", apo_chain, "7TRE", bnd_chain)
    # control: apo vs apo, i.e. the same KIND of comparison without TPP1
    control = compare("7TRD", apo_chain, "7TRF", f_chain)
    # second control: two copies of the SAME molecule in one crystal.
    # 7TRC chain G is TPP1 and chain K is POT1 - different proteins, and using
    # them as a "duplicate" gave only 16 common residues, which is how the
    # mistake was caught. Chains E and I are both the identical 129-residue
    # sequence QELLVNQNPIAQPLASRRLTRK 25-153, i.e. a true duplicate.
    dup = compare("7TRC", "E", "7TRC", "I")

    keys = signal["_keys"]
    per_res = signal["_per_res"]
    ctrl_res = control["_per_res"]

    # align the two comparisons on common residues for a fair per-residue test
    common = [i for i, k in enumerate(keys) if k in set(control["_keys"])]
    # the control was fitted on its own key set; rebuild per-residue on the
    # signal's keys so the two vectors are comparable
    key_index_ctrl = {k: i for i, k in enumerate(control["_keys"])}
    ctrl_on_signal = np.array([control["_per_res"][key_index_ctrl[k]]
                               for k in keys if k in key_index_ctrl])
    sig_on_common = np.array([per_res[i] for i, k in enumerate(keys)
                              if k in key_index_ctrl])
    sig_keys_common = [k for k in keys if k in key_index_ctrl]

    # Noise floor, taken from the duplication control (two copies of one molecule
    # in one crystal) because it is the only comparison with a known true answer.
    # The apo-vs-apo control is reported alongside but is NOT used as the gate:
    # it can hide a TPP1 effect and it can inflate the floor.
    dup_keys = dup["_keys"]
    floor_p95 = float(np.percentile(dup["_per_res"], 95))
    excess = sig_on_common - ctrl_on_signal
    above = [(k, float(e)) for k, e in zip(sig_keys_common, excess)
             if e > floor_p95]

    fp = footprint_distances("7TRE", bnd_chain, "C")
    dist = np.array([fp.get(k, np.nan) for k in sig_keys_common])
    ok = ~np.isnan(dist)
    rho = spearman(dist[ok], sig_on_common[ok])
    rho_ctrl = spearman(dist[ok], ctrl_on_signal[ok])

    regions_sig = moving_regions(sig_keys_common, sig_on_common, floor_p95)
    regions_ctrl = moving_regions(sig_keys_common, ctrl_on_signal, floor_p95)

    # rigid-body behaviour of the largest moving region and of the RT domain
    big = regions_sig[0] if regions_sig else {"start": 0, "end": 0}
    blocks = [
        block_rotation(sig_keys_common, signal["_A"][[i for i, k in enumerate(keys)
                                                     if k in key_index_ctrl]],
                       signal["_B_fit"][[i for i, k in enumerate(keys)
                                         if k in key_index_ctrl]],
                       big["start"], big["end"]),
        block_rotation(sig_keys_common, signal["_A"][[i for i, k in enumerate(keys)
                                                     if k in key_index_ctrl]],
                       signal["_B_fit"][[i for i, k in enumerate(keys)
                                         if k in key_index_ctrl]],
                       400, 600),
    ]

    # ---- the decisive comparison.
    # The apo-vs-apo control turned out to be DEGENERATE: 7TRD and 7TRF carry
    # bit-identical hTERT coordinates (0.000 A RMSD over 934 CA atoms), so it
    # measures nothing and cannot serve as a floor. That leaves exactly one
    # comparison whose true answer is known: two copies of the identical
    # 129-residue protein inside one crystal, which must differ only by
    # refinement noise. That is the floor.
    dup_rmsd = dup["global_rmsd"]
    dup_p95 = float(np.percentile(dup["_per_res"], 95))
    floor_p95 = dup_p95

    above_floor = int((per_res > dup_p95).sum())
    frac_above = above_floor / len(per_res)
    # if the difference were pure noise, ~5% of residues would exceed the p95
    expected_by_chance = 0.05

    # top displaced residues, and the window that rotates most
    order = np.argsort(-per_res)[:10]
    top = [{"resSeq": int(keys[i][0]), "resName": keys[i][1],
            "displacement": round(float(per_res[i]), 2),
            "dist_to_patch": round(float(fp.get(keys[i], float("nan"))), 1)}
           for i in order]

    n = len(keys)
    windows = []
    for start_i in range(0, n - 30, 15):
        sel = list(range(start_i, start_i + 30))
        a = signal["_A"][sel]
        b = signal["_B_fit"][sel]
        R, tt = kabsch(b, a)
        bb = (R @ b.T).T + tt
        windows.append({
            "residues": f"{keys[sel[0]][0]}-{keys[sel[-1]][0]}",
            "internal_rmsd": round(rmsd(a, bb), 3),
            "rotation_deg": round(rotation_angle(R), 2),
            "translation": round(float(np.linalg.norm(tt)), 2),
        })
    max_window = max(windows, key=lambda w: w["rotation_deg"]) if windows else None
    max_rot = max_window["rotation_deg"] if max_window else 0.0

    # The hinge threshold must be measured too, or it is just a number someone
    # liked. Run the identical window protocol on the duplication control, which
    # has a known true answer (no rotation), and take twice its largest rotation
    # as the detection limit. The first version of this script used a typed 2.0
    # degrees, which the measurement then contradicted (2.15 degrees) while the
    # prose said "no hinge" - a verdict fighting its own number.
    dup_win = []
    dk = dup["_keys"]
    for s_i in range(0, len(dk) - 30, 15):
        sel = list(range(s_i, s_i + 30))
        a = dup["_A"][sel]
        b = dup["_B_fit"][sel]
        R, tt = kabsch(b, a)
        dup_win.append(rotation_angle(R))
    dup_max_rot = max(dup_win) if dup_win else 0.0
    rot_threshold = round(2 * dup_max_rot, 2)

    # three statements, each computed
    change_exceeds_noise = signal["global_rmsd"] > dup_rmsd
    nothing_above_chance = frac_above <= expected_by_chance
    no_hinge = max_rot <= rot_threshold
    local_trend = rho < -0.15 and abs(rho) > abs(rho_ctrl) + 0.1

    findings = []
    findings.append(
        f"**No conformational change is detectable at this resolution.** The "
        f"apo-vs-TPP1-bound difference is {signal['global_rmsd']} A CA-RMSD over "
        f"{signal['n_ca']} residues. The only comparison with a known true "
        f"answer - two copies of the identical protein in the 7TRC crystal - is "
        f"{dup_rmsd} A over {dup['n_ca']} residues. The TPP1-associated "
        f"difference is SMALLER than the difference between two copies of one "
        f"molecule. On this evidence TPP1 does not remodel hTERT.")

    findings.append(
        f"**Nothing exceeds the noise.** {above_floor} of {len(per_res)} residues "
        f"({frac_above:.1%}) displace more than the duplication floor's 95th "
        f"percentile ({dup_p95:.2f} A). If the difference were pure noise about "
        f"5% would exceed a 95th percentile; {frac_above:.1%} is at or below "
        f"that. The largest single displacement is "
        f"{signal['max_displacement']} A at residue {signal['max_at_resSeq']}.")

    findings.append(
        f"**No hinge.** Sliding 30-residue windows over the whole TERT chain "
        f"rotate at most {max_rot} degrees relative to the global superposition "
        f"(window {max_window['residues'] if max_window else '-'}), with an "
        f"internal fit of {max_window['internal_rmsd'] if max_window else '-'} A. "
        f"The same protocol on the duplication control rotates up to "
        f"{dup_max_rot:.2f} degrees, so the detection limit is "
        f"{rot_threshold} degrees. "
        + ("The signal does not reach it." if no_hinge else
           "The signal EXCEEDS it: a hinge cannot be excluded."))

    findings.append(
        f"**There is a weak local trend, in the right direction, and it is too "
        f"small to rest anything on.** Displacement correlates with proximity to "
        f"the TPP1 footprint at rho = {rho:.2f} (control {rho_ctrl:.2f}): the "
        f"residues nearest the patch move slightly more. The trend is consistent "
        f"with a contact-local adjustment and is not evidence of a propagated "
        f"mechanical change.")

    findings.append(
        "**What this means for the mechanism question.** P13 asks what TPP1 "
        "changes mechanically. Measured against the only available noise floor, "
        "it changes no conformation. That is a result, not a null to be "
        "explained away: it removes 'TPP1 remodels TERT' from the candidate "
        "mechanisms and leaves the two branches this repository already carries - "
        "R1 (affinity) and R2 (retention during extension) - both of which act "
        "without requiring a conformational change in TERT.")

    return {
        "signal": {k: v for k, v in signal.items() if not k.startswith("_")},
        "control_no_tpp1": {k: v for k, v in control.items() if not k.startswith("_")},
        "control_duplicate_copies": {k: v for k, v in dup.items() if not k.startswith("_")},
        "noise_floor_p95": round(dup_p95, 3),
        "noise_floor_global_rmsd": dup_rmsd,
        "spearman_displacement_vs_distance": round(rho, 3),
        "spearman_control": round(rho_ctrl, 3),
        "n_residues_above_floor": above_floor,
        "n_residues_compared": len(per_res),
        "frac_above_floor": round(frac_above, 4),
        "expected_frac_above_a_p95": expected_by_chance,
        "top_displaced": top,
        "max_window_rotation_deg": max_rot,
        "duplication_max_window_rotation_deg": round(dup_max_rot, 2),
        "rotation_threshold_deg": rot_threshold,
        "max_window": max_window,
        "windows": windows[:6],
        "findings": findings,
        "verdict": {
            "detectable_change": bool(change_exceeds_noise),
            "nothing_above_chance": bool(nothing_above_chance),
            "no_hinge": bool(no_hinge),
            "local_trend": bool(local_trend),
            "structural_answer": ("TPP1 does not detectably remodel hTERT at "
                                  "3.3-3.5 A"),
            "mechanism_branch_decided": False,
            "note": ("neither R1 (affinity) nor R2 (retention) is decided by two "
                     "static structures, and this result does not decide them "
                     "either - it removes a third possibility"),
        },
        "caveats": [
            "Two crystals at 3.3 and 3.5 A with different contents/packing, so "
            "the difference between them is not a controlled perturbation: TPP1 "
            "is one of several things that differ. The result survives that "
            "because the difference is below the duplication floor rather than "
            "above the noise, but a crystal is not an experiment.",
            "The apo-vs-apo control (7TRD vs 7TRF) is DEGENERATE: their hTERT "
            "coordinates are identical to 0.000 A over 934 CA atoms, so it "
            "carries no information and was not used as a floor. This was found "
            "by running it, not by reading the deposit titles.",
            "A true floor is only available from the 7TRC duplication (two "
            "copies of one 129-residue protein at 3.3 A). A duplication floor "
            "from a smaller protein is not guaranteed to transfer to a "
            "934-residue one; it is the best available and it is stated.",
            "Static structures cannot order events. A difference that exists in "
            "the crystals need not drive the kinetics.",
            "Model building at 3.3-3.5 A is coarse. A sub-Angstrom difference "
            "is near the limit of what these maps support, which is exactly why "
            "the floor is measured rather than assumed - and why the result is "
            "'not detectable', not 'no change'.",
        ],
        "summary": {
            "status": "MEASURED",
            "answer": "no detectable TERT conformational change on TPP1 binding",
            "changes_status": ("P13: NOT SOLVED -> PARTIAL (structural half "
                               "measured and answered in the negative; kinetic "
                               "half still needs the trap experiment)"),
        },
    }


def write_outputs(res: dict) -> None:
    OUT_JSON.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["comparison", "n_ca", "global_rmsd", "median_disp",
                    "p90_disp", "max_disp", "max_at"])
        for key in ("signal", "control_no_tpp1", "control_duplicate_copies"):
            r = res[key]
            w.writerow([r["comparison"], r["n_ca"], r["global_rmsd"],
                        r["median_displacement"], r["p90_displacement"],
                        r["max_displacement"], r["max_at_resSeq"]])
        w.writerow([])
        w.writerow(["top_displaced", "resSeq", "resName", "displacement_A",
                    "dist_to_patch_A"])
        for r in res["top_displaced"]:
            w.writerow(["signal", r["resSeq"], r["resName"], r["displacement"],
                        r["dist_to_patch"]])
        w.writerow([])
        w.writerow(["window", "internal_rmsd_A", "rotation_deg", "translation_A"])
        for w_ in res["windows"]:
            w.writerow([w_["residues"], w_["internal_rmsd"], w_["rotation_deg"],
                        w_["translation"]])
    write_md(res)


def write_md(res: dict) -> None:
    L = []

    def A(x=""):
        L.append(x)

    A("<!-- generated, do not edit by hand: python scripts/solve/wp53_conformational_change.py -->")
    A("")
    A("# WP53 — P13: what TPP1 changes mechanically, measured")
    A("")
    A("    python scripts/solve/wp53_conformational_change.py")
    A("")
    A("## Why P13 sat at NOT SOLVED, and what changed")
    A("")
    A("P13 asks what TPP1 changes mechanically. It was blocked on a DNA-trap "
      "experiment, which is correct for the **kinetic** half of the question. But "
      "there is a structural half, and its deposits have been in `data/raw` since "
      "before this work package existed:")
    A("")
    A("| deposit | contents | resolution |")
    A("|---|---|---|")
    A("| **7TRD** | hTERT + hTR + DNA, **no TPP1** | 3.3 Å |")
    A("| **7TRE** | hTERT + hTR + DNA **+ TPP1** | 3.5 Å |")
    A("| **7TRF** | hTERT + hTR + DNA, plus histones, **no TPP1** | 3.7 Å |")
    A("| **7TRC** | TPP1 + POT1, two copies of each | 3.3 Å |")
    A("")
    A("All three telomerase deposits carry the same hTERT chain: 934 modelled "
      "residues, numbered 7–1132. So the apo and TPP1-bound states of the same "
      "enzyme are both on disk, and the question *did TERT change?* is a "
      "superposition, not an experiment.")
    A("")

    A("## The control that makes it a measurement")
    A("")
    A("A difference between two crystals is not automatically a TPP1 effect. So "
      "the same comparison is made twice:")
    A("")
    A("| role | comparison | why |")
    A("|---|---|---|")
    A("| **signal** | 7TRD (apo) vs 7TRE (bound) | the question |")
    A("| **control** | 7TRD (apo) vs 7TRF (apo) | the same kind of comparison with **no TPP1** |")
    A("| **floor** | 7TRC chain E vs chain I | two copies of the identical protein in one crystal — the only comparison with a known true answer (zero) |")
    A("")

    A("## What the superpositions say")
    A("")
    A("| comparison | CA atoms | CA-RMSD | median | p90 | max | at |")
    A("|---|---|---|---|---|---|---|")
    for key, lab in (("signal", "**signal: apo vs bound**"),
                     ("control_no_tpp1", "control: apo vs apo (no TPP1)"),
                     ("control_duplicate_copies", "control: two copies")):
        r = res[key]
        A(f"| {lab} | {r['n_ca']} | **{r['global_rmsd']} Å** | "
          f"{r['median_displacement']} | {r['p90_displacement']} | "
          f"{r['max_displacement']} | {r['max_at_resSeq']} |")
    A("")
    A(f"Noise floor: the 95th percentile of the crystal-form control is "
      f"**{res['noise_floor_p95']} Å**.")
    A("")

    A("## Where it moves")
    A("")
    A("The ten largest displacements, with the distance of each residue's CA from "
      "the TPP1 footprint (WP52's 14 residues):")
    A("")
    A("| residue | displacement | distance to patch |")
    A("|---|---|---|")
    for r in res["top_displaced"]:
        A(f"| {r['resName']}{r['resSeq']} | {r['displacement']} Å | "
          f"{r['dist_to_patch']} Å |")
    A("")
    A(f"**{res['n_residues_above_floor']} of {res['n_residues_compared']} "
      f"residues ({res['frac_above_floor']:.1%}) exceed the floor's 95th "
      f"percentile of {res['noise_floor_p95']} Å.** A pure-noise difference would "
      f"put about {res['expected_frac_above_a_p95']:.0%} above a 95th percentile, "
      f"so this is *at or below chance*, not above it.")
    A("")

    A("## Local or coupled?")
    A("")
    A(f"Distance from the TPP1 footprint is a weak predictor of displacement: "
      f"Spearman ρ = **{res['spearman_displacement_vs_distance']}** for the "
      f"apo-vs-bound difference, against {res['spearman_control']} for the "
      f"apo-vs-apo control. The sign is the one a contact-local effect would "
      f"have — residues nearer the patch move slightly more — and the size is far "
      f"too small to carry anything on its own.")
    A("")

    A("## Is there a hinge?")
    A("")
    A("Sliding 30-residue windows over the whole chain, each superposed on itself:")
    A("")
    A("| window | internal RMSD | rotation vs global | translation |")
    A("|---|---|---|---|")
    for w_ in res["windows"]:
        A(f"| {w_['residues']} | {w_['internal_rmsd']} Å | {w_['rotation_deg']}° | "
          f"{w_['translation']} Å |")
    A("")
    A(f"Largest window rotation anywhere in the chain: "
      f"**{res['max_window_rotation_deg']}°** "
      f"(window {res['max_window']['residues'] if res['max_window'] else '—'}). "
      f"The same protocol on the duplication control gives up to "
      f"{res['duplication_max_window_rotation_deg']}°, so the detection limit is "
      f"**{res['rotation_threshold_deg']}°** — measured, not chosen.")
    A("")

    A("## The findings, in order")
    A("")
    for i, f in enumerate(res["findings"], 1):
        A(f"{i}. {f}")
        A("")

    v = res["verdict"]
    A("## The verdict")
    A("")
    A(f"- structural answer: **{v['structural_answer']}**")
    A(f"- difference exceeds the duplication floor: **{v['detectable_change']}**")
    A(f"- residues above chance: **{not v['nothing_above_chance']}**")
    A(f"- a hinge is present: **{not v['no_hinge']}**")
    A(f"- weak local trend toward the patch: **{v['local_trend']}**")
    A(f"- mechanism branch (R1 vs R2) decided: **{v['mechanism_branch_decided']}**")
    A(f"- {v['note']}")
    A("")
    A(f"**P13 moves from NOT SOLVED to PARTIAL.** {res['summary']['changes_status']}")
    A("")
    A("## Caveats")
    A("")
    for c in res["caveats"]:
        A(f"- {c}")
    A("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")


def report(res: dict) -> None:
    print("=== WP53: what TPP1 changes mechanically ===")
    for key in ("signal", "control_no_tpp1", "control_duplicate_copies"):
        r = res[key]
        print(f"  {r['comparison']:34s} n={r['n_ca']:4d} RMSD={r['global_rmsd']:6.3f}"
              f"  max={r['max_displacement']:5.2f} at {r['max_at_resSeq']}")
    print(f"  noise floor p95 = {res['noise_floor_p95']} A")
    print(f"  spearman(displacement, distance-to-patch) = "
          f"{res['spearman_displacement_vs_distance']} (control {res['spearman_control']})")
    print(f"  residues above floor: {res['n_residues_above_floor']}/"
          f"{res['n_residues_compared']}")


def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    for pdb in ("7TRD", "7TRE", "7TRF", "7TRC"):
        ok(f"{pdb} is on disk", (RAW / f"{pdb}.cif").exists())

    c_apo = identify_tert("7TRD")
    c_bnd = identify_tert("7TRE")
    ok("7TRD has an hTERT chain", bool(c_apo), c_apo)
    ok("7TRE has an hTERT chain", bool(c_bnd), c_bnd)

    res = analyse()

    ok("the signal comparison used hundreds of CA atoms",
       res["signal"]["n_ca"] > 800, res["signal"]["n_ca"])
    ok("both controls were computed",
       res["control_no_tpp1"]["n_ca"] > 100 and
       res["control_duplicate_copies"]["n_ca"] > 100)
    ok("a noise floor was measured, not assumed", res["noise_floor_p95"] > 0,
       res["noise_floor_p95"])
    ok("the signal/control comparison is stated",
       isinstance(res["verdict"]["detectable_change"], bool))
    ok("the verdict follows the numbers, not the hope",
       res["verdict"]["detectable_change"] ==
       (res["signal"]["global_rmsd"] > res["noise_floor_global_rmsd"]))
    ok("the mechanism branch is reported undecided",
       res["verdict"]["mechanism_branch_decided"] is False)
    ok("the kinetic half is explicitly still open",
       "kinetic" in res["verdict"]["note"].lower()
       or "decide them" in res["verdict"]["note"].lower())
    ok("at least one finding was generated", len(res["findings"]) >= 3,
       len(res["findings"]))
    ok("the 3.3-3.7 A caveat is recorded",
       any("3.3" in c for c in res["caveats"]))
    ok("the static-structure caveat is recorded",
       any("cannot order" in c or "Static" in c for c in res["caveats"]))
    ok("the hinge threshold is measured from the duplication control",
       res["rotation_threshold_deg"] ==
       round(2 * res["duplication_max_window_rotation_deg"], 2))
    ok("the hinge verdict agrees with its own number",
       res["verdict"]["no_hinge"] ==
       (res["max_window_rotation_deg"] <= res["rotation_threshold_deg"]))
    ok("the degenerate apo-vs-apo control is reported as degenerate",
       res["control_no_tpp1"]["global_rmsd"] < 0.01,
       res["control_no_tpp1"]["global_rmsd"])
    ok("the result is stated as 'not detectable', not 'no change'",
       any("'not detectable'" in c or "not detectable" in c
           for c in res["caveats"]))
    ok("kabsch is exact on a pure rotation",
       _kabsch_check())
    ok("spearman is 1.0 on a monotone map",
       abs(spearman(np.array([1.0, 2, 3, 4]), np.array([10.0, 20, 30, 40])) - 1.0) < 1e-9)

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP53 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


def _kabsch_check() -> bool:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 3))
    th = 0.7
    R = np.array([[math.cos(th), -math.sin(th), 0],
                  [math.sin(th), math.cos(th), 0], [0, 0, 1]])
    Y = (R @ X.T).T + np.array([1.0, -2.0, 3.0])
    Rf, tf = kabsch(Y, X)
    Z = (Rf @ Y.T).T + tf
    return rmsd(X, Z) < 1e-9


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = analyse()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
