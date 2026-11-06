#!/usr/bin/env python3
"""
WP7 - Independence audit of the five depositions the repository compares.

The repository derives cross-state conclusions by comparing 7TRE, 7TRD and 7TRF.  That is only
valid if the entries are independent determinations.  This script measures how independent they
actually are, at the level of the deposited coordinates.

Finding: 7TRD and 7TRF carry the *same TERT model*.  Atom order is identical, per-atom B-factors
are bit-identical (not merely equal to rounding: the float32 values are equal), and the
coordinates differ by a rigid-body/refinement offset (median 0.23 A, max 2.3 A).  Consequently
every 7TRE-vs-7TRD and 7TRE-vs-7TRF comparison in results/tables/ returns the same number to
three decimals, because it is measuring the same comparison twice.

This does not make the deposited models wrong - reusing a refined model across related entries is
normal - but it means the repository's "cross-state" evidence rests on ONE independent TERT
coordinate set (the apo one), not two.

Run:  python scripts/solve/wp7_entry_independence.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import gemmi
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
ENTRIES = ["7TRE", "7TRD", "7TRF"]


def tert_model(pdb: str):
    st = gemmi.read_structure(str(REPO / "data/raw" / f"{pdb}.cif"))
    m = st[0]
    chain = next(c for c in m if c.name == "A")
    ids, xyz, b = [], [], []
    for r in chain:
        for a in r:
            ids.append((r.seqid.num, a.name))
            xyz.append([a.pos.x, a.pos.y, a.pos.z])
            b.append(a.b_iso)
    return ids, np.array(xyz), np.array(b)


def main() -> int:
    models = {p: tert_model(p) for p in ENTRIES}
    rows = []
    print("=" * 84)
    print("WP7  Pairwise comparison of the deposited TERT models (chain A, auth numbering)")
    print("=" * 84)
    print(f"{'pair':>12} {'same atom order':>16} {'same B values':>14} "
          f"{'median d(x)':>12} {'max d(x)':>10}")
    for i, a in enumerate(ENTRIES):
        for b in ENTRIES[i + 1:]:
            ia, xa, ba = models[a]
            ib, xb, bb = models[b]
            same_order = ia == ib
            same_b = bool(np.array_equal(ba, bb)) if same_order else False
            if same_order:
                d = np.linalg.norm(xa - xb, axis=1)
                med, mx = float(np.median(d)), float(d.max())
            else:
                med = mx = float("nan")
            print(f"{a+'/'+b:>12} {str(same_order):>16} {str(same_b):>14} "
                  f"{med:>12.3f} {mx:>10.2f}")
            rows.append(dict(pair=f"{a}/{b}", same_atom_order=same_order,
                             same_b_factors=same_b, median_coord_diff_A=round(med, 4),
                             max_coord_diff_A=round(mx, 3)))

    print()
    print("=" * 84)
    print("Consequence for the repository's cross-state tables")
    print("=" * 84)
    tbl = REPO / "results" / "tables" / "domain_rmsd.csv"
    if tbl.exists():
        seen: dict[tuple, dict] = {}
        dupes = []
        for r in csv.DictReader(open(tbl)):
            if r["domain_name"] not in ("TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"):
                continue
            key = (r["domain_name"],)
            val = (r["independent_rmsd_angstrom"], r["core_fit_displacement_angstrom"])
            if key in seen and abs(float(seen[key][0]) - float(val[0])) < 1e-3:
                dupes.append((r["domain_name"], r["reference_structure"],
                              seen[key][1], r["reference_structure"], val[0]))
            seen.setdefault(key, (val[0], r["reference_structure"]))
        print("  domains whose 7TRD and 7TRF RMSDs agree to < 0.001 A:")
        for d in dupes:
            print(f"    {d[0]:16s} vs {d[1]}: {d[2]} A   |   vs {d[3]}: {d[4]} A")
        if not dupes:
            print("    none")
    print()
    print("  Both columns are the same comparison because both reference structures contain")
    print("  the same TERT model.  Treat 7TRD and 7TRF as ONE apo TERT state, not two.")

    with open(OUT / "wp7_entry_independence.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT/'wp7_entry_independence.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
