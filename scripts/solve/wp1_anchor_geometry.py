#!/usr/bin/env python3
"""
WP1 - Geometry of the TERT DNA anchor site in the 7TRE / 7TRD / 7TRF depositions.

Motivation
----------
`results/mechanism/structural_validation.md` states:

    "The TEN domain contributes 0 primer-contacting residues; Tyr176, the residue at the
     centre of the original bipartite-clamp hypothesis, is 11.1 A from the nearest primer atom."

and uses that measurement to invalidate the TEN/TPP1 clamp model. The present script tests
whether the measurement supports that conclusion, by asking a question the original analysis
did not: *where is the DNA that is being measured against, and is it in a position for the
anchor site to engage it at all?*

Result:

    7TRE models 7 nucleotides of DNA (DT12-DT13-DT14-DA15-DG16-DG17-DG18, i.e. TTTAGGG).
    Polarity matters, and was established separately in WP8:

        the 3' terminal residue DG18 is 2.65 A from catalytic D868 and 4.03 A from D869
        the 5' terminal residue DT12 is 5.02 A from R756 in the IFD-TRAP anchor strand
        the 5' exit direction points at Tyr176 (26.7 deg off-axis, 6.8 A closest approach)

    So the 7-mer fills the active site and the anchor strand and stops just short of PLYQ.
    Per-region distances (7TRE / 7TRD / 7TRF):

        RT channel 940-984           2.6 / 2.5 / 2.5 A   the DNA is docked here
        anchor strand 756-759        5.0 / 4.4 / 4.3 A   contact
        PLYQ 174-177                11.1 / 11.9 / 11.9 A  measured to an unoccupied site
        IFD 766-775                 22.2 / 22.1 / 22.1 A  beyond the 7-mer's reach

    The PLYQ number is therefore a property of the substrate used for the deposition, not
    evidence that TEN does not bind DNA.  Liu et al. 2022 (Nature 604:578) identified PLYQ as
    the anchor and showed Y176A/Q177A reduces activity and RAP.  See WP8 for the polarity
    analysis, which is what licenses this reading.

    NOTE: an earlier draft of this script's header claimed the 7-mer "terminates at the active
    site" and that the distance "measures an unoccupied anchor site" while implying the primed
    end pointed away from PLYQ.  The polarity result in WP8 supersedes that wording: the 5' end
    of the modelled fragment is at the anchor strand and its exit path runs into PLYQ.

Reproduce
---------
    python scripts/solve/wp1_anchor_geometry.py

Writes results/solve/wp1_anchor_geometry.csv and prints the table.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import gemmi

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
OUTDIR = REPO / "results" / "solve"

DNA3 = {"DA", "DT", "DG", "DC", "DU"}
RNA3 = {"A", "U", "G", "C", "I"}

# hTERT (UniProt O14746) regions of interest, author numbering as used in 7TRE/7TRD/7TRF
REGIONS = {
    "PLYQ_anchor_174_177": [174, 175, 176, 177],
    "anchor_strand_756_759": [756, 757, 759],
    "IFD_766_775": list(range(766, 776)),
    "helix_a25_772_780": list(range(772, 781)),
    "TEN_99FGF101": [99, 100, 101],
    "RT_channel_940_984": list(range(940, 985)),
}

ENTRIES = ["7TRE", "7TRD", "7TRF"]


def classify(model: gemmi.Model) -> dict[str, tuple[str, int]]:
    """Return {chain_name: (kind, n_residues)} using residue composition only."""
    out: dict[str, tuple[str, int]] = {}
    for chain in model:
        names = [r.name for r in chain]
        if not names:
            continue
        n = len(names)
        if all(r in DNA3 for r in names):
            kind = "DNA"
        elif all(r in RNA3 for r in names):
            kind = "RNA"
        elif n > 600:
            kind = "TERT"
        elif n >= 100:
            kind = "PROTEIN"
        else:
            kind = "other"
        out[chain.name] = (kind, n)
    return out


def heavy_coords(chain: gemmi.Chain) -> list[tuple[int, str, str, float, float, float]]:
    out = []
    for res in chain:
        for at in res:
            if at.element.name in ("H", "D", ""):
                continue
            out.append(
                (res.seqid.num, res.name, at.name, at.pos.x, at.pos.y, at.pos.z)
            )
    return out


def min_dist(a: list, b: list):
    best = (float("inf"), None, None)
    for x in a:
        x0, x1, x2 = x[3], x[4], x[5]
        for y in b:
            d2 = (x0 - y[3]) ** 2 + (x1 - y[4]) ** 2 + (x2 - y[5]) ** 2
            if d2 < best[0]:
                best = (d2, x, y)
    return math.sqrt(best[0]), best[1], best[2]


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    summary = {}

    for pdb in ENTRIES:
        structure = gemmi.read_structure(str(RAW / f"{pdb}.cif"))
        structure.setup_entities()
        model = structure[0]
        kinds = classify(model)

        tert_chain = max((c for c, (k, _) in kinds.items() if k == "TERT"),
                         key=lambda c: kinds[c][1])
        dna_chain = min((c for c, (k, _) in kinds.items() if k == "DNA"),
                        key=lambda c: kinds[c][1])
        tpp1_chain = None
        prot = [c for c, (k, n) in kinds.items() if k == "PROTEIN" and 120 <= n <= 200]
        if prot:
            tpp1_chain = prot[0]

        tert = heavy_coords(model[tert_chain])
        dna = heavy_coords(model[dna_chain])
        dna_res = sorted({(a[0], a[1]) for a in dna})

        print("=" * 78)
        print(f"{pdb}: TERT={tert_chain} DNA={dna_chain} TPP1={tpp1_chain}")
        print(f"  modelled DNA: {len(dna_res)} nt  "
              f"{'-'.join(r[1] + str(r[0]) for r in dna_res)}")

        summary[pdb] = {
            "tert_chain": tert_chain,
            "dna_chain": dna_chain,
            "tpp1_chain": tpp1_chain,
            "dna_nt_modelled": len(dna_res),
        }

        for label, nums in REGIONS.items():
            sub = [a for a in tert if a[0] in nums]
            if not sub:
                rows.append(dict(pdb=pdb, region=label, min_dist_A=None,
                                 n_modelled=0, n_expected=len(nums),
                                 closest_pair="not modelled"))
                continue
            d, a, b = min_dist(sub, dna)
            pair = f"{a[1]}{a[0]}:{a[2]}<->{b[1]}{b[0]}:{b[2]}"
            n_mod = len({x[0] for x in sub})
            rows.append(dict(pdb=pdb, region=label, min_dist_A=round(d, 2),
                             n_modelled=n_mod, n_expected=len(nums), closest_pair=pair))
            print(f"    {label:22s} {d:6.2f} A  ({n_mod}/{len(nums)} modelled)  {pair}")

        if tpp1_chain:
            tpp1 = heavy_coords(model[tpp1_chain])
            d, a, b = min_dist(tpp1, dna)
            print(f"    {'TPP1_to_DNA':22s} {d:6.2f} A  {a[1]}{a[0]}:{a[2]}<->{b[1]}{b[0]}:{b[2]}")
            rows.append(dict(pdb=pdb, region="TPP1_to_DNA", min_dist_A=round(d, 2),
                             n_modelled=None, n_expected=None,
                             closest_pair=f"{a[1]}{a[0]}:{a[2]}<->{b[1]}{b[0]}:{b[2]}"))
            summary[pdb]["tpp1_dna_min_A"] = round(d, 2)
            d, a, b = min_dist(tpp1, tert)
            print(f"    {'TPP1_to_TERT':22s} {d:6.2f} A  {a[1]}{a[0]}:{a[2]}<->{b[1]}{b[0]}:{b[2]}")
            rows.append(dict(pdb=pdb, region="TPP1_to_TERT", min_dist_A=round(d, 2),
                             n_modelled=None, n_expected=None,
                             closest_pair=f"{a[1]}{a[0]}:{a[2]}<->{b[1]}{b[0]}:{b[2]}"))

    import csv

    with open(OUTDIR / "wp1_anchor_geometry.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb", "region", "min_dist_A",
                                           "n_modelled", "n_expected", "closest_pair"])
        w.writeheader()
        w.writerows(rows)

    (OUTDIR / "wp1_chains.json").write_text(json.dumps(summary, indent=1))
    print(f"\nwrote {OUTDIR/'wp1_anchor_geometry.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
