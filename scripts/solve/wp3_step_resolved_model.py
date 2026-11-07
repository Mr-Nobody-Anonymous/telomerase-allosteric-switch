#!/usr/bin/env python3
"""
WP3 - A step-resolved model that can actually express the measurements, and its predictions.

Why a new model is needed
-------------------------
WP2 shows the repository's eq.1 cannot express the Latrick & Cech (EMBO J 2010) data: it has
one primer k_off and no translocation step, whereas the measurements require a dissociation
rate that depends on which step of the catalytic cycle the enzyme is in, and a separate
translocation step.

The model below is the smallest one that does express them.  Per repeat the enzyme goes through

    extension  ->  translocation  ->  (next repeat)

with translocation either succeeding (carrying the enzyme into the next repeat) or failing
(the primer is lost).  Everything is measured, not fitted:

    phi      translocation efficiency       0.57 -> 0.99        (L&C Fig 2F)
    k_tr     translocation rate constant    0.07 -> 0.10 /min   (L&C Fig 2F)
    k_pol    nucleotide incorporation        unchanged          (L&C Results)
    k_off,rest  resting primer dissociation  0.012 /min, unchanged (L&C Fig 2)

From phi = k_tr/(k_tr + k_dissoc) we get k_dissoc = k_tr (1 - phi)/phi, and therefore a
*dissociation rate during translocation* that POT1-TPP1 lowers ~50-fold, while leaving the
resting-state k_off untouched.

What this model predicts, with no free parameters
-------------------------------------------------
    per-repeat survival contributed by translocation, s_phi = phi
    fold-change in that survival                     = phi_PT / phi_noPT
    if the survival also scales with the rate at which the enzyme commits to the next repeat,
    the combined prediction is (phi_PT/phi_noPT) * (k_tr,PT/k_tr,noPT)

and the latter is compared with the independently measured 2.7 +/- 0.6-fold processivity
increase.

Run:  python scripts/solve/wp3_step_resolved_model.py
"""

from __future__ import annotations

import csv
import math
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# Latrick & Cech, EMBO J 29:924 (2010)
PHI_NO_PT, PHI_PT = 0.57, 0.99
KTR_NO_PT, KTR_PT = 0.07, 0.10          # /min
KOFF_RESTING = 0.012                    # /min
RAP_FOLD_OBS, RAP_FOLD_SD = 2.7, 0.6    # processivity increase, n = 4


def k_dissoc_during_translocation(k_tr: float, phi: float) -> float:
    """From phi = k_tr/(k_tr + k_dissoc)."""
    return k_tr * (1.0 - phi) / phi


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    kd_no = k_dissoc_during_translocation(KTR_NO_PT, PHI_NO_PT)
    kd_pt = k_dissoc_during_translocation(KTR_PT, PHI_PT)

    print("=" * 88)
    print("WP3.1  Step-resolved quantities implied by the published measurements")
    print("=" * 88)
    rows = [
        dict(quantity="translocation efficiency phi", no_PT=PHI_NO_PT, with_PT=PHI_PT,
             fold=round(PHI_PT / PHI_NO_PT, 3), units="fraction"),
        dict(quantity="translocation rate k_tr", no_PT=KTR_NO_PT, with_PT=KTR_PT,
             fold=round(KTR_PT / KTR_NO_PT, 3), units="/min"),
        dict(quantity="k_dissoc during translocation", no_PT=round(kd_no, 4),
             with_PT=round(kd_pt, 4), fold=round(kd_no / kd_pt, 2), units="/min"),
        dict(quantity="resting primer k_off", no_PT=KOFF_RESTING, with_PT=KOFF_RESTING,
             fold=1.0, units="/min"),
    ]
    print(f"  {'quantity':34s} {'no PT':>9} {'+PT':>9} {'fold':>7}  units")
    for r in rows:
        print(f"  {r['quantity']:34s} {r['no_PT']:>9.4f} {r['with_PT']:>9.4f} "
              f"{r['fold']:>7.2f}  {r['units']}")
    print()
    print("  Read the last two rows together.  POT1-TPP1 changes the dissociation rate of the")
    print("  enzyme *during translocation* by ~50-fold and the resting-state k_off by 1.0-fold.")
    print("  Any model that carries a single dissociation rate for the primer cannot hold both.")
    print()

    print("=" * 88)
    print("WP3.2  Prediction of the observed processivity increase - no free parameters")
    print("=" * 88)
    pred_phi = PHI_PT / PHI_NO_PT
    pred_comb = pred_phi * (KTR_PT / KTR_NO_PT)
    print(f"  (a) survival = phi only .................. {pred_phi:.2f}x")
    print(f"  (b) survival = phi * (k_tr ratio) ........ {pred_comb:.2f}x")
    print(f"  observed processivity increase ........... {RAP_FOLD_OBS:.1f} +/- "
          f"{RAP_FOLD_SD:.1f}x  (n=4)")
    for label, pred in (("(a)", pred_phi), ("(b)", pred_comb)):
        z = (pred - RAP_FOLD_OBS) / RAP_FOLD_SD
        verdict = "consistent" if abs(z) <= 2 else "inconsistent at 2 SD"
        print(f"  {label} agreement: {z:+.2f} SD vs observed -> {verdict}")
    print()
    print("  With one free parameter (the weight of the k_tr term, 0..1) the prediction spans")
    print(f"  {pred_phi:.2f}-{pred_comb:.2f}x.  This range does NOT contain the observed")
    print(f"  {RAP_FOLD_OBS}x - the upper end falls {(RAP_FOLD_OBS - pred_comb)/RAP_FOLD_SD:.2f} SD")
    print("  short, though every point in the range is within 2 SD of it.  So the family")
    print("  sits entirely BELOW the central observation.  On a log scale the phi ratio")
    print(f"  accounts for {math.log(pred_phi)/math.log(RAP_FOLD_OBS)*100:.0f}% of the effect, the k_tr ratio")
    print(f"  for a further {math.log(pred_comb/pred_phi)/math.log(RAP_FOLD_OBS)*100:.0f}%, leaving {(1-math.log(pred_comb)/math.log(RAP_FOLD_OBS))*100:.0f}%")
    print(f"  ({(RAP_FOLD_OBS-pred_comb)/RAP_FOLD_SD:.2f} SD) unexplained.")
    print("  (An earlier version of this script claimed the range 'always brackets the")
    print("  observation'; that was arithmetically false and is corrected here - WP18.)")
    print("  Conversely, no parameterisation of the repository's eq.1 reproduces both")
    print("  the unchanged resting k_off and the 2.7-fold increase (WP2.3).")
    print()

    print("=" * 88)
    print("WP3.3  Where this leaves the mechanism question")
    print("=" * 88)
    concl = [
        ("Is there any published quantity that POT1-TPP1 changes which a "
         "DNA-retention model cannot explain?", "no"),
        ("Is there any published quantity that requires a mechanical/conformational "
         "term at the active site?", "no"),
        ("Does the step-resolved retention model predict the observed 2.7x within "
         "experimental error?", "yes - 2.48x, 0.4 SD"),
        ("Is the allosteric alternative excluded, or merely unsupported?", "unsupported"),
        ("What would exclude it?", "a measured PT-dependent change in a rate at a site "
                                   "TPP1 does not contact, on an enzyme with no DNA"),
    ]
    for q, a in concl:
        print(f"  Q: {q}\n     A: {a}\n")

    with open(OUT / "wp3_step_resolved.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["quantity", "no_PT", "with_PT", "fold", "units"])
        w.writeheader()
        w.writerows(rows)

    (OUT / "wp3_predictions.json").write_text(json.dumps({
        "prediction_phi_only": round(pred_phi, 4),
        "prediction_phi_times_ktr": round(pred_comb, 4),
        "observed": RAP_FOLD_OBS,
        "observed_sd": RAP_FOLD_SD,
        "z_phi_only": round((pred_phi - RAP_FOLD_OBS) / RAP_FOLD_SD, 3),
        "z_combined": round((pred_comb - RAP_FOLD_OBS) / RAP_FOLD_SD, 3),
        "k_dissoc_no_PT": round(kd_no, 5),
        "k_dissoc_PT": round(kd_pt, 5),
        "resting_koff_no_PT": KOFF_RESTING,
        "resting_koff_PT": KOFF_RESTING,
    }, indent=1))

    print(f"wrote {OUT/'wp3_step_resolved.csv'}")
    print(f"wrote {OUT/'wp3_predictions.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
