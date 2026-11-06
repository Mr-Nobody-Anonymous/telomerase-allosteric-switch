#!/usr/bin/env python3
"""
WP9 - Where on TEN does TPP1 act?  The FGF loop vs the PLYQ anchor motif.

WP1-WP8 treated "retention" as a single mechanism with a single structural route (the TEN anchor
site engaging DNA).  Two published results break that assumption, and they disagree with each
other about *where* on TEN the functional element is:

  Liu et al. 2022   PLYQ (174-177) is the DNA anchor; Y176A/Q177A reduces activity and RAP.
  Zaug et al. 2010  G100 is required for POT1-TPP1 stimulation; G100V abolishes it, while the
                    enzyme retains full basal activity and responsiveness to other effectors.

Both cannot be "the TPP1-responsive DNA grip".  This script measures where each sits.

Result: they are on opposite faces of TEN.

    residue   d(TPP1)   d(DNA)     role
    A46         3.35     36.46     closest TPP1 contact on TEN
    F101        5.81     26.77     FGF loop
    G100        7.38     24.34     Zaug 2010: G100V kills stimulation
    F99         8.81     18.66     FGF loop
    Y176       19.21     11.11     Liu 2022: DNA anchor

The FGF loop (99-101) is the one TEN element that is simultaneously (i) in direct contact with
TPP1 and (ii) required for TPP1's effect, while being 19-27 A from the DNA - too far to grip it.
PLYQ is the DNA-proximal element but is 19 A from TPP1, so TPP1 cannot act on it directly.

TEN is also a *tethered* domain: residues 198-321 are unmodelled in all three TERT entries, so TEN
is joined to the ring by a 124-residue linker and has hinge freedom.  Only 36 of 160 TEN C-alpha
atoms have a ring C-alpha within 10 A.

Interpretation this licenses: TPP1's effect on TEN is expressed at the TPP1-TEN interface, at the
FGF loop, and must be transmitted through TEN's own orientation (set by an unmodelled linker)
rather than by TPP1 gripping DNA.  That is a *conformational* mechanism that produces retention.
It is neither "TPP1 holds the DNA" nor "TPP1 activates catalysis".

Run:  python scripts/solve/wp9_ten_face.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import gemmi
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

RESIDUES = [
    (46, "A46", "closest TPP1 contact on TEN"),
    (101, "F101", "FGF loop"),
    (100, "G100", "Zaug 2010: G100V abolishes stimulation"),
    (99, "F99", "FGF loop"),
    (176, "Y176", "Liu 2022: DNA anchor"),
    (177, "Q177", "Liu 2022: DNA anchor"),
    (174, "P174", "PLYQ"),
    (175, "L175", "PLYQ"),
]


def atoms(chain):
    return [(r.seqid.num, r.name, a.name, np.array([a.pos.x, a.pos.y, a.pos.z]), a.b_iso)
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

    def dmin(resid, pts):
        return min(float(np.linalg.norm(pts - x[3], axis=1).min()) for x in byres[resid])

    print("=" * 88)
    print("WP9  The two faces of the TEN domain")
    print("=" * 88)
    print(f"{'residue':>8} {'d(TPP1)':>9} {'d(DNA)':>8} {'med B':>7}  role")
    print("-" * 88)
    recs = []
    for resid, label, role in RESIDUES:
        if resid not in byres:
            print(f"{label:>8}  not modelled")
            continue
        b = float(np.median([x[4] for x in byres[resid]]))
        dt, dd = dmin(resid, pa), dmin(resid, da)
        print(f"{label:>8} {dt:>9.2f} {dd:>8.2f} {b:>7.1f}  {role}")
        recs.append(dict(residue=label, resid=resid, d_TPP1_A=round(dt, 2),
                         d_DNA_A=round(dd, 2), median_B=round(b, 2), role=role))

    print()
    print("  Both cannot be the same functional element:")
    gf = next(r for r in recs if r["residue"] == "G100")
    y176 = next(r for r in recs if r["residue"] == "Y176")
    print(f"    G100  is {gf['d_TPP1_A']:.1f} A from TPP1 and {gf['d_DNA_A']:.1f} A from DNA")
    print(f"    Y176  is {y176['d_TPP1_A']:.1f} A from TPP1 and {y176['d_DNA_A']:.1f} A from DNA")
    print(f"    -> they are {abs(gf['d_DNA_A'] - y176['d_DNA_A']):.1f} A apart in DNA distance and")
    print(f"       {abs(gf['d_TPP1_A'] - y176['d_TPP1_A']):.1f} A apart in TPP1 distance: opposite faces.")
    print()

    # --- TEN is tethered, not docked
    ten = sorted(r for r in byres if r <= 196)
    ring = [x for x in tert if x[0] >= 326]
    ring_ca = np.array([x[3] for x in ring if x[2] == "CA"])
    ten_ca = np.array([x[3] for x in byres[r] if False] if False else
                      [x[3] for r in ten for x in byres[r] if x[2] == "CA"])
    buried = int(np.sum(np.min(np.linalg.norm(ten_ca[:, None, :] - ring_ca[None, :, :], axis=2),
                               axis=1) <= 10.0))
    print("  TEN as a tethered domain:")
    print(f"    TEN residues modelled        : {min(ten)}-{max(ten)}  (n = {len(ten)})")
    print(f"    residues 197-321 modelled    : {[r for r in range(197, 322) if r in byres]} "
          f"-> an unmodelled 124-residue linker")
    print(f"    TEN C-alphas with a ring CA within 10 A : {buried} of {len(ten_ca)}")
    print()
    print("    -> TEN is joined to the TERT ring by a long unmodelled linker.  Its orientation")
    print("       relative to the ring is a hinge degree of freedom, not a fixed pose.")
    print()

    print("=" * 88)
    print("What this changes")
    print("=" * 88)
    print("  'Retention' and 'allostery' are not the only two options, and neither is a")
    print("  mechanism.  There are three distinguishable routes by which TPP1 could raise RAP:")
    print()
    print("    R1  TPP1 changes the enzyme's DNA affinity directly          (TPP1 contacts DNA)")
    print("        -> measured: NO.  TPP1 is 28.7 A from the DNA (WP1).")
    print("    R2  TPP1, via the FGF loop, sets how TEN engages the DNA    (TPP1 -> TEN -> DNA)")
    print("        -> 3.35 A contact, and the element required for stimulation (G100) is")
    print("           24 A from DNA, so it acts on TEN's orientation, not on the DNA itself.")
    print("    R3  TPP1 changes catalysis allosterically                    (TPP1 -> active site)")
    print("        -> 35.8-42.1 A across the RT ring (WP4); no short path; unsupported.")
    print()
    print("  R2 is a conformational mechanism that produces retention.  It is not 'the DNA is")
    print("  held tighter' and not 'the active site is faster'.  It is the TEN/DNA engagement")
    print("  described by Balch et al. 2026 and Zaug et al. 2010, and it is the branch the")
    print("  repository's binary framing cannot represent.")

    with open(OUT / "wp9_ten_face.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader()
        w.writerows(recs)
    print(f"\nwrote {OUT/'wp9_ten_face.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
