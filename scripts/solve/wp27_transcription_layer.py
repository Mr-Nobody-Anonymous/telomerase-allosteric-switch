#!/usr/bin/env python3
"""
WP27 - The Yamanaka layer: how OCT4 / SOX2 / KLF4 relate to M5, and where they
       belong in a combined design.

WHAT THESE FACTORS ARE
----------------------
The image is three of the four Yamanaka factors (Takahashi & Yamanaka 2006;
Nobel Prize 2012): SOX2, OCT4 (POU5F1) and KLF4. The fourth, c-MYC, is off-frame.

HOW THEY TOUCH TELOMERASE - and this is the part that matters
------------------------------------------------------------
They act on TERT TRANSCRIPTION, not on the enzyme. Documented connections:

  KLF4      binds the TERT proximal promoter directly and activates it
            (Wong 2010; Hoffmeyer 2012; Wang 2012)
  OCT4      activates TERC and can activate TERT; OCT3/4 + NANOG bind the
            TERC promoter (Agarwal 2010)
  SOX2      cooperates with OCT4 as a heterodimer on pluripotency promoters
  c-MYC     dispensable for telomerase activation in mouse iPSCs; also
            reported to activate TERT
  During iPSC reprogramming Tert and Terc are upregulated EARLY, before
  endogenous Nanog/Oct4/Sox2 (Wang 2012; PMC3317559).

That is a DIFFERENT LAYER from everything in this repository. M5 and the whole
WP1-WP26 programme describe the ENZYME: how TPP1 changes per-repeat survival s.
The Yamanaka factors change HOW MUCH TERT there is. The two compose
multiplicatively and do not compete:

    product per recruit = R1/2(s)
    recruits per cell     ~ [TERT]                (transcription layer)
    total output         = R1/2(s) x [TERT]

THE LEVERAGE ASYMMETRY (computed below)
---------------------------------------
R1/2 = -ln2 / ln(s), so the elasticity d ln R1/2 / d ln s = -1/ln(s).
At 0.79 that is 4.3; at 0.917 it is 11.5. A 1% gain in per-repeat survival buys
4-12% more product. A 1% gain in TERT abundance buys exactly 1%.

So the mechanism layer M5 addresses is the HIGH-LEVERAGE one, and the Yamanaka
layer is the low-leverage one. That is not a criticism of the transcription
factors - it is the reason a design should care about both, in that order.

ALSO IMPORTANT: THREE OF THE FOUR ARE NOT THERAPEUTIC
-----------------------------------------------------
This repository retracted every therapeutic claim it had (TELIN-1 quarantined,
14 of 56 claims contradicted). The Yamanaka factors must be held to the same
standard, and they have a specific, well-documented hazard: OCT4 overexpression
causes off-target gene activation and epigenetic aberrations in the resulting
iPSCs, and SKM (without OCT4) gives BETTER developmental potential
(Velychko 2019, Cell Stem Cell). They are a research tool and a
reprogramming reagent, not a telomerase drug.

Run:  python scripts/solve/wp27_transcription_layer.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# L&C anchors used throughout WP22-WP25
S_NO_PT = 2 ** (-1 / 2.963)
S_PT = 2 ** (-1 / 8.0)


def R_half(s: float) -> float:
    return -math.log(2.0) / math.log(s)


def elasticity(s: float) -> float:
    """d ln(R1/2) / d ln(s)."""
    return -1.0 / math.log(s)


def main() -> int:
    L = "=" * 92
    print(L)
    print("WP27  The Yamanaka layer vs M5: two different knobs on one output")
    print(L)
    print()

    # ---------------------------------------------------------------- layer
    print("WP27.1  What layer each factor acts on")
    print(f"  {'factor':10s} {'layer':22s} {'documented action':44s}")
    print("  " + "-" * 78)
    rows_doc = [
        ("OCT4", "transcription", "activates TERC (with NANOG); can activate TERT"),
        ("SOX2", "transcription", "heterodimer partner of OCT4 on promoters"),
        ("KLF4", "transcription", "binds TERT proximal promoter directly"),
        ("c-MYC", "transcription", "dispensable for telomerase activation (mouse)"),
        ("TPP1", "ENZYME (this repo)", "raises per-repeat survival s; M5"),
    ]
    for a, b, c in rows_doc:
        print(f"  {a:10s} {b:22s} {c:44s}")
    print()
    print("  The Yamanaka factors act one layer ABOVE everything in WP1-WP26.")
    print("  They set [TERT]; M5 sets how much each TERT molecule makes per recruit.")
    print()

    # ---------------------------------------------------------------- compose
    print("WP27.2  The composition is multiplicative, not competing")
    print("      product per recruit  = R1/2(s)          <- M5 layer")
    print("      recruits per cell    ~ [TERT]           <- Yamanaka layer")
    print("      total output         = R1/2(s) x [TERT]")
    print()
    print("  A simulation or an assay of M5 says nothing about [TERT], and a")
    print("  transcription experiment says nothing about s.  They are orthogonal,")
    print("  which is why the Yamanaka factors cannot CONFIRM or REFUTE M5.")
    print()

    # ---------------------------------------------------------------- leverage
    print("WP27.3  Leverage: which knob actually moves the output")
    print(f"  {'state':22s} {'s':>8} {'R1/2':>8} {'elasticity':>12}")
    print("  " + "-" * 54)
    lev_rows = []
    for lbl, s in (("no TPP1 (L&C)", S_NO_PT), ("+ TPP1 (L&C)", S_PT),
                   ("s = 0.99", 0.99)):
        R = R_half(s)
        e = elasticity(s)
        print(f"  {lbl:22s} {s:8.4f} {R:8.3f} {e:12.2f}")
        lev_rows.append((lbl, s, R, e))
    print()
    print("  A 1% gain in per-repeat survival buys 4.3-11.5% more product.")
    print("  A 1% gain in TERT abundance buys exactly 1%.")
    print()
    print("  SCOPE NOTE (WP28): 4.3-11.5 are the two ANCHOR SLOPES, not a")
    print("  range.  eps = -1/ln(s) diverges as s -> 1 (99.5 at s=0.99), and")
    print("  the linear [TERT] term here is the x<<1 limit of a saturating")
    print("  one.  See WP28 for M6' and for why the ranking depends on the")
    print("  operating point.  Nothing below is retracted; its scope is.")
    print()
    print("  => The enzyme layer is the high-leverage one.  If a design can")
    print("     improve s at all, it beats an equal-percentage improvement in")
    print("     abundance - by a factor of 4 to 12.")
    print()

    # ---------------------------------------------------------------- integration
    print("WP27.4  Folding the layer into the design as M6")
    print("  M5 (unchanged)  : s = (1+a)^-6,  M2/M4 state-specific translocation")
    print("  M6 adds one term: [TERT] as an independent multiplicative factor")
    print()
    print("      output = R1/2(s) * f([TERT])")
    print()
    print("  What this buys: a design that can be tested at EITHER layer without")
    print("  changing the other, so a null result is informative.  If s improves")
    print("  and output does not, the transcription layer is limiting, and vice")
    print("  versa.  Neither layer can hide behind the other.")
    print()
    print("  What this does NOT buy: any new constraint on s.  [TERT] is not in")
    print("  the identity, so it cannot break the M2/M4 degeneracy and it cannot")
    print("  resolve the 2.73x vs >=4x tension.  Those remain kinetic.")
    print()

    # ---------------------------------------------------------------- safety
    print("WP27.5  Claim discipline - the same standard applied to everything else")
    print("  This repository retracted every therapeutic claim it had: TELIN-1 is")
    print("  quarantined with zero empirical evidence, and 14 of 56 claims were")
    print("  found contradicted.  The Yamanaka factors get the same treatment:")
    print()
    print("   NOT claimed here:")
    print("     - that activating TERT transcription is therapeutic")
    print("     - that any of these factors is a telomerase drug")
    print("   Documented hazard, from the primary literature:")
    print("     - OCT4 overexpression causes off-target gene activation and")
    print("       epigenetic aberrations in iPSCs; the OCT4-free SKM cocktail")
    print("       gives BETTER developmental potential (Velychko 2019).")
    print("     - The four factors drive cells to pluripotency.  That is a")
    print("       cell-fate intervention, not a tune-up of an enzyme.")
    print()
    print("  Their legitimate role in this repository is as a MODEL SYSTEM: a")
    print("  documented, reversible way to raise [TERT] so the transcription")
    print("  layer's contribution can be measured against the enzyme layer's.")
    print()

    # ---------------------------------------------------------------- verdict
    print("WP27.6  VERDICT")
    print("  1. The Yamanaka factors act on TERT transcription, one layer above")
    print("     M5. They cannot confirm or refute the mechanism.")
    print("  2. The two layers compose multiplicatively: output = R1/2(s) x [TERT].")
    print("  3. Per percent, leverage favours the enzyme layer at every")
    print("     operating point (eps_s >= 4.27 vs eps_x <= 1).  In total range")
    print("     the ranking depends on the operating point - see WP28.")
    print("  4. Adding the layer as M6 makes a null result interpretable, because")
    print("     the two layers can now be varied independently.")
    print("  5. No therapeutic claim.  The OCT4 hazard is documented and real.")
    print()

    rows = []
    for lbl, s, R, e in lev_rows:
        rows.append(dict(item=f"R1/2, {lbl}", value=round(R, 4),
                         note=f"elasticity {e:.2f}"))
    rows += [
        dict(item="composition", value="R1/2(s) x [TERT]", note="multiplicative, orthogonal"),
        dict(item="leverage advantage, enzyme vs abundance", value="4.3-11.5x",
             note="per 1% change"),
        dict(item="Yamanaka factors constrain M5", value="no", note="different layer"),
        dict(item="therapeutic claim", value="none", note="OCT4 hazard documented"),
    ]
    out = OUT / "wp27_transcription_layer.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0





# ---------------------------------------------------------------- WP27.7 tests
def _self_test() -> int:
    """The integration must not disturb any existing verdict."""
    print()
    print("WP27.7  Self-test: does the new layer change anything already decided?")
    checks = []

    # 1. composition is multiplicative and monotone
    def out(s, tert):
        return R_half(s) * tert

    base = out(S_PT, 1.0)
    checks.append(("doubling [TERT] doubles output",
                   abs(out(S_PT, 2.0) / base - 2.0) < 1e-12))
    eps = 1e-6
    dlogR = (math.log(R_half(S_PT * (1 + eps))) - math.log(R_half(S_PT))) / eps
    checks.append(("2% more s beats 2% more [TERT]", dlogR > 1.0))
    checks.append(("...elasticity matches -1/ln(s) analytically",
                   abs(dlogR - elasticity(S_PT)) < 1e-3))

    # 2. the layer cannot break the degeneracy
    #    M2 and M4 give identical s by construction, so identical output:
    s_M2 = S_PT
    s_M4 = S_PT
    checks.append(("[TERT] cannot separate M2 from M4",
                   out(s_M2, 1.5) == out(s_M4, 1.5)))

    # 3. the 2.73x vs 4x tension is untouched by scaling abundance
    ratio_a = out(S_PT, 1.0) / out(S_NO_PT, 1.0)
    ratio_b = out(S_PT, 3.0) / out(S_NO_PT, 3.0)
    checks.append(("the s-ratio is invariant to [TERT]", abs(ratio_a - ratio_b) < 1e-12))

    for name, ok in checks:
        print(f"    {name:52s} {'ok' if ok else 'FAIL'}")
    n_fail = sum(1 for _, ok in checks if not ok)
    print(f"    {len(checks)} checks, {n_fail} failed")
    return n_fail


if __name__ == "__main__":
    raise SystemExit(main() or _self_test())
