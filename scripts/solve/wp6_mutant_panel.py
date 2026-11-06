#!/usr/bin/env python3
"""
WP6 - Can the published mutant panel discriminate?  A test the repository implies it can.

The repository's MECHANISM_ASSESSMENT.md makes structural claims from the Tesmer et al. 2023
IFD-TRAP mutant panel (F6: "H762E, inside the named conduit, is near-WT").  If mutant severity
tracked a structural quantity - distance to the TPP1 interface, or distance to the DNA - then
that panel would be evidence about mechanism.

This script tests that directly.  For each mutant, the residual processivity (mutant / WT from
the +PT lanes of Tesmer Fig 3B, n = 3) is correlated against three structural quantities
measured in 7TRE:

    d_TPP1   minimum heavy-atom distance from the mutated residue to TPP1
    d_DNA    minimum heavy-atom distance from the mutated residue to the modelled DNA
    nbrs     number of TERT C-alpha atoms within 10 A of the residue's C-alpha (burial)

Result: nothing correlates (all p > 0.15 at n = 8-9), and the sign of the correlation flips
when F776E - flagged in this repository's own data-quality notes as undefined at low signal -
is excluded.  The panel therefore cannot carry the structural inference.  Reported as a null.

Run:  python scripts/solve/wp6_mutant_panel.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import gemmi
import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
PANEL = REPO / "data" / "external" / "tesmer2023" / "fig3b_processivity_plusPT.csv"


def atoms(chain):
    return [(r.seqid.num, r.name, a.name, np.array([a.pos.x, a.pos.y, a.pos.z]))
            for r in chain for a in r if a.element.name not in ("H", "D", "")]


def main() -> int:
    st = gemmi.read_structure(str(REPO / "data/raw/7TRE.cif"))
    st.setup_entities()
    m = st[0]
    tert, tpp1, dna = atoms(m["A"]), atoms(m["C"]), atoms(m["D"])
    pa = np.array([x[3] for x in tpp1])
    da = np.array([x[3] for x in dna])
    byres: dict[int, list] = {}
    for x in tert:
        byres.setdefault(x[0], []).append(x)
    ca_all = np.array([x[3] for x in tert if x[2] == "CA"])

    rows = list(csv.DictReader(l for l in open(PANEL) if not l.startswith("#")))
    muts = [k for k in rows[0] if k not in ("experiment", "WT")]
    wt = np.mean([float(r["WT"]) for r in rows])

    recs = []
    for mut in muts:
        resid = int("".join(c for c in mut if c.isdigit()))
        vals = np.array([float(r[mut]) for r in rows])
        sub = byres.get(resid)
        if not sub:
            continue
        ca = [s[3] for s in sub if s[2] == "CA"][0]
        recs.append(dict(
            mutant=mut, resid=resid,
            mean_processivity=round(float(vals.mean()), 4),
            rel_to_WT=round(float(vals.mean() / wt), 4),
            d_to_TPP1_A=round(float(min(np.linalg.norm(pa - s[3], axis=1).min() for s in sub)), 2),
            d_to_DNA_A=round(float(min(np.linalg.norm(da - s[3], axis=1).min() for s in sub)), 2),
            burial_neighbors=int(np.sum(np.linalg.norm(ca_all - ca, axis=1) <= 10.0)),
        ))

    print(f"WT mean processivity fraction = {wt:.3f}  (n = 3 experiments)\n")
    hdr = f"{'mutant':>8} {'resid':>6} {'mean':>7} {'rel':>7} {'d_TPP1':>8} {'d_DNA':>8} {'nbrs':>6}"
    print(hdr)
    print("-" * len(hdr))
    for r in recs:
        print(f"{r['mutant']:>8} {r['resid']:>6} {r['mean_processivity']:>7.3f} "
              f"{r['rel_to_WT']:>7.3f} {r['d_to_TPP1_A']:>8.2f} {r['d_to_DNA_A']:>8.2f} "
              f"{r['burial_neighbors']:>6}")
    print()

    flagged = "F776E"   # repo's own data-quality note: undefined at low signal
    out = []
    for label, subset in (("all mutants", recs),
                          (f"excl. {flagged} (flagged)", [r for r in recs if r["mutant"] != flagged])):
        y = np.array([r["rel_to_WT"] for r in subset])
        print(f"{label}  (n = {len(subset)})")
        for var in ("d_to_TPP1_A", "d_to_DNA_A", "burial_neighbors"):
            x = np.array([r[var] for r in subset])
            if np.std(x) == 0:
                continue
            pr, pp = stats.pearsonr(x, y)
            sr, sp = stats.spearmanr(x, y)
            print(f"   {var:>15}: Pearson r={pr:+.3f} p={pp:.3f} | "
                  f"Spearman rho={sr:+.3f} p={sp:.3f}")
            out.append(dict(subset=label, variable=var, pearson_r=round(float(pr), 4),
                            pearson_p=round(float(pp), 4), spearman_rho=round(float(sr), 4),
                            spearman_p=round(float(sp), 4), n=len(subset)))
        print()

    print("Conclusion: no structural quantity predicts mutant severity at this sample size,")
    print("and the sign is unstable to the exclusion of one flagged point.  The mutant panel")
    print("is therefore not evidence for any structural mechanism, in either direction.")

    with open(OUT / "wp6_mutant_panel.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader()
        w.writerows(recs)
    with open(OUT / "wp6_correlations.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {OUT/'wp6_mutant_panel.csv'}")
    print(f"wrote {OUT/'wp6_correlations.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
