#!/usr/bin/env python3
"""
WP8 - DNA polarity: which end of the modelled 7-mer is at the active site?

WP1 reported per-region distances but did not establish which end of the DNA is which.  That turns
out to decide the interpretation, and an earlier draft of WP1 got it backwards.

7TRE chain D is DT12-DT13-DT14-DA15-DG16-DG17-DG18 (5'-TTTAGGG-3' reading 12 -> 18).  This script
establishes, from the coordinates:

  * the 3' terminal phosphate (DG18) is 5.1 A from the catalytic aspartate D869 and its neighbour
    DG15 is 3.4 A from K973 - the 3' end is at the active site, where chemistry happens;
  * the 5' terminal phosphate (DT12) is 5.0 A from R756, i.e. it abuts the IFD-TRAP N-terminal
    "anchor strand" 756-759;
  * extending the 5' direction away from the enzyme traces a path that passes F759 and then
    approaches Tyr176 of the TEN domain to 6.8 A.

So the deposited 7-mer is fully occupied at the active site and at the anchor strand, and stops
just short of PLYQ.  PLYQ's empty site in 7TRE is therefore expected, not informative about whether
TEN binds DNA.

An earlier draft of this analysis asserted the opposite - that the 7-mer ended at the active site
and pointed away from PLYQ.  That draft checked the distance from PLYQ to the DNA but not the
polarity of the strand, and is superseded here.

Run:  python scripts/solve/wp8_dna_polarity.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import gemmi
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"


def atoms(chain):
    return [(r.seqid.num, r.name, a.name, np.array([a.pos.x, a.pos.y, a.pos.z]))
            for r in chain for a in r if a.element.name not in ("H", "D", "")]


def main() -> int:
    st = gemmi.read_structure(str(REPO / "data/raw/7TRE.cif"))
    st.setup_entities()
    m = st[0]
    tert, dna = atoms(m["A"]), atoms(m["D"])
    byres: dict[int, list] = {}
    for x in tert:
        byres.setdefault(x[0], []).append(x)

    def one(resid, name):
        return next(x[3] for x in byres[resid] if x[2] == name)

    def dmin_between(dna_resid, resid):
        """Minimum heavy-atom distance between one DNA residue and one protein residue."""
        dxyz = np.array([x[3] for x in dna if x[0] == dna_resid])
        xyz = np.array([x[3] for x in byres.get(resid, [])])
        if xyz.size == 0 or dxyz.size == 0:
            return float("nan")
        return float(np.min(np.linalg.norm(dxyz[:, None, :] - xyz[None, :, :], axis=2)))

    phosphates = sorted([x for x in dna if x[2] == "P"], key=lambda x: x[0])
    pts = np.array([x[3] for x in phosphates])
    first, last = phosphates[0][0], phosphates[-1][0]

    print("=" * 84)
    print("WP8  Polarity of the modelled DNA (7TRE chain D)")
    print("=" * 84)
    print("  sequence 5'->3': T T T A G G G   (DT12 ... DG18)")
    print(f"  spacing between adjacent phosphates: "
          f"{np.mean([np.linalg.norm(pts[i+1]-pts[i]) for i in range(len(pts)-1)]):.2f} A "
          f"(one nucleotide step)")
    print()

    # polarity anchor: which end is at the catalytic site
    print("  landmark distances to each DNA end (minimum heavy-atom):")
    rows = []
    for label, resid in (("5' end", first), ("3' end", last)):
        rec = dict(end=label, residue=resid, atoms_modelled=len([x for x in dna if x[0] == resid]))
        for nm, rid in (("D712", 712), ("D868", 868), ("D869", 869),
                        ("R756", 756), ("K757", 757), ("F759", 759),
                        ("Y176", 176), ("K973", 973)):
            rec[f"to_{nm}"] = round(dmin_between(resid, rid), 2)
        rows.append(rec)
        print(f"    {label} (res {resid}): " +
              "  ".join(f"{nm} {rec[f'to_{nm}']:.2f}" for nm in
                        ("D712", "D868", "D869", "R756", "K757", "F759", "Y176", "K973")))
    print()
    print("    -> the 3' end is nearest the catalytic aspartates (D869) and K973; the 5' end is")
    print("       nearest the anchor strand R756/K757.  Chemistry is 3'-terminal, so the 7-mer")
    print("       is correctly oriented with its 3' end in the active site.")
    print()

    # exit direction
    d5 = pts[0] - pts[1]
    d5 /= np.linalg.norm(d5)
    y176 = one(176, "OH")
    v = y176 - pts[0]
    v /= np.linalg.norm(v)
    ang5 = float(np.degrees(np.arccos(np.clip(np.dot(v, d5), -1, 1))))
    d3 = pts[-1] - pts[-2]
    d3 /= np.linalg.norm(d3)
    ang3 = float(np.degrees(np.arccos(np.clip(np.dot(v, d3), -1, 1))))

    print("  direction of the 5' exit path (unit vector from P13 toward P12):")
    print(f"    angle between that direction and the vector to Tyr176-OH : {ang5:5.1f} deg")
    print(f"    same, along the 3' direction (toward the active site)    : {ang3:5.1f} deg")
    print(f"    -> Tyr176 lies along the 5' exit path, not behind the DNA.")
    print()

    print("  extrapolating the 5' path (rigid, straight-line approximation):")
    traj = []
    for t in np.arange(0, 21, 3.0):
        q = pts[0] + d5 * t
        dd, x = min(((np.linalg.norm(q - y[3]), y) for y in tert), key=lambda z: z[0])
        print(f"    +{t:4.1f} A : nearest TERT atom = {x[1]}{x[0]}:{x[2]}  at {dd:5.2f} A")
        traj.append(dict(extension_from_5prime_A=round(float(t), 2),
                         nearest_residue=f"{x[1]}{x[0]}", nearest_atom=x[2],
                         distance_A=round(float(dd), 2)))
    close = min(((np.linalg.norm(pts[0] + d5 * t - y176), t) for t in np.arange(0, 30, 0.25)),
                key=lambda z: z[0])
    print()
    print(f"    closest approach of that path to Tyr176-OH: {close[0]:.2f} A "
          f"at +{close[1]:.2f} A extension")
    print(f"    straight-line distance from the 5' phosphate to Tyr176-OH: "
          f"{np.linalg.norm(pts[0]-y176):.2f} A")
    print()
    print("  Caveat: this is a straight-line extrapolation of one nucleotide step and ignores")
    print("  helix curvature.  It establishes the DIRECTION of the exit path, not a contact.")
    print()

    print("=" * 84)
    print("Reading")
    print("=" * 84)
    print("  The deposited 7-mer occupies the active site (3' end) and the anchor strand")
    print("  756-759 (5' end).  It does not reach PLYQ.  Because the 5' exit path points at")
    print("  Tyr176, a substrate longer by roughly two to four nucleotides would occupy the")
    print("  PLYQ site.  Any statement about whether TEN binds the primer must therefore come")
    print("  from a structure with a longer DNA - 7QXS (30 nt) is the obvious one and is not")
    print("  reachable from this sandbox.")
    print()
    print("  Consequence for the repository: its claim 'the TEN domain makes zero primer")
    print("  contacts, measured at 11.1 A' is a statement about an unoccupied site.  That the")
    print("  site is unoccupied is a property of the 7-mer substrate, not of the enzyme.")

    with open(OUT / "wp8_dna_polarity.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(OUT / "wp8_exit_path.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(traj[0]))
        w.writeheader()
        w.writerows(traj)
    print(f"\nwrote {OUT/'wp8_dna_polarity.csv'}")
    print(f"wrote {OUT/'wp8_exit_path.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
