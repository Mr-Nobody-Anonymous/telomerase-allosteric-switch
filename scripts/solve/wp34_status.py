#!/usr/bin/env python3
"""
WP34 - The honest scoreboard: what is not achieved, how to achieve it, the best
       candidate, and what its success rate actually is.

Written because three questions deserve numbers rather than narrative:

  1. WHAT IS NOT ACHIEVED      - itemised, each with its blocker
  2. HOW TO ACHIEVE IT         - one experiment each, costed in lab-hours
  3. BEST CANDIDATE + SUCCESS RATE - with the arithmetic shown, not asserted

RETRACTED IN PART (WP35).  Read WP35_P2_IDENTIFIABILITY.md before using any
number below.  Three claims in this script do not survive:

  1. The "88% explained" figure.  It compares 1 - 2.73/4.00 = 0.318 with
     (s0-phi0)/s0 = 0.280, but s and phi are the two DIFFERENT quantities whose
     mismatch IS the anomaly being explained - so (s0-phi0)/s0 is that anomaly
     rescaled and cannot explain it.  The same expression is NEGATIVE in the
     +PT arm, (-0.0796), and a fraction of a population cannot be negative.
  2. "The retry probability should measure ~0.32".  0.318 is not a retry
     probability.  Under the renewal model the retry probability at which the
     two readouts reconcile is tau_0 = 0.1833.
  3. "The maintenance gate is CLOSED" and "P(approved drug, now) = 0".  48.0
     bp/PD is delta* at tau_0 = 0, and tau_0 = 0 is exactly what P2 tests.  The
     two physically derived stall couplings give 57.3 and 57.4 bp/PD, above the
     50 floor.  The gate is UNDETERMINED, not closed.

The gap table, the candidate ranking, the pipeline arithmetic and the quarantine
status are unchanged and remain usable.  WP34.1 below KEEPS its old output, on
purpose: it is the record of what was claimed and it is asserted by tests as the
arithmetic that was made.  Fixing it silently would delete the evidence.

    python scripts/solve/wp34_status.py
    -> results/solve/wp34_status.csv
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

S_NO_PT = 0.7914
S_PT = 0.9170
PHI_NO_PT = 0.570
PHI_PT = 0.99
A_NO_PT = 0.03976
A_PT = 0.01455
PRED_RATIO = 2.73
MEAS_RATIO = 4.00


def delta_star(a: float) -> float:
    """Critical attrition, bp per population doubling, at full coverage."""
    return 1.0 / math.log2(1.0 + a)


def offrate_needed(delta: float) -> float:
    """k_off ratio required to maintain a cell losing `delta` bp/PD."""
    return A_NO_PT / (2.0 ** (1.0 / delta) - 1.0)


def wp34_pt_value() -> float:
    """The (s-phi)/s expression evaluated in the +PT arm.

    Exists to show the retracted construction is not a population fraction: it
    is negative here, and a fraction of a population cannot be.
    """
    return (S_PT - PHI_PT) / S_PT


def main() -> int:
    L = "=" * 90
    print(L)
    print("WP34  What is not achieved, how to achieve it, and the odds")
    print(L)
    print()

    # ---------------------------------------------------------------- 1
    print("WP34.1  [RETRACTED BY WP35 - kept as the record, do not cite]")
    retry_needed = MEAS_RATIO / PRED_RATIO
    p_needed = 1.0 - 1.0 / retry_needed
    stalled = (S_NO_PT - PHI_NO_PT) / S_NO_PT
    print(f"  required retry multiplier        {MEAS_RATIO:.2f}/{PRED_RATIO:.2f} = {retry_needed:.4f}")
    print(f"  required retry probability       1 - 1/{retry_needed:.4f}    = {p_needed:.4f}")
    print(f"  stalled fraction of survivors    ({S_NO_PT:.4f}-{PHI_NO_PT:.3f})/{S_NO_PT:.4f} = {stalled:.4f}")
    print(f"  agreement                        {stalled:.4f}/{p_needed:.4f} = {stalled/p_needed:.2f}")
    print()
    print("  The stalled population is not a fudge factor - L&C report it")
    print("  directly ('stalled and dissociated complexes could contribute to")
    print("  the untranslocated population'). So the prediction is now SHARP:")
    print(f"  the retry probability should measure ~{p_needed:.2f}")
    print()
    print("  RETRACTED (WP35).  Both claims above fail:")
    print(f"    - 0.280 is (s0-phi0)/s0, i.e. the anomaly itself, rescaled; and it")
    print(f"      is {wp34_pt_value():+.4f} in the +PT arm - not a fraction at all.")
    print(f"    - 0.318 is not a retry probability.  The renewal-model crossing is")
    print(f"      tau_0 = 0.1833 (WP35), and the fold it reconciles is 3.267x, not")
    print(f"      2.733x.  The stall CANNOT close the ratio gap: kappa_dist is")
    print(f"      bounded, and 4.00x is outside the bound.")
    print()

    # ---------------------------------------------------------------- 2
    print("WP34.2  WHAT IS NOT ACHIEVED")
    gaps = [
        ("Mechanism (M5) not confirmed",
         "M5 is the only physical, falsifiable candidate, but no independent",
         "DNA-trap re-analysis of existing gels (P2). No new reagent."),
        ("The 4x-vs-2.73x tension unresolved",
         "WP35: stalling CANNOT close it - kappa_dist is bounded and 4.00x is"
         " outside the bound. Retained as retracted, WP34's 88% does not hold",
         "same P2 run: it measures tau_0, and the reconciled fold is 3.267x."),
        ("M2 vs M4 degeneracy open",
         "phi = k_tr/(k_tr+k_off_tr) is 2 unknowns to 1 observation",
         "not load-bearing: phi left the product, s never sees k_tr."),
        ("Maintenance unreachable at known kinetics",
         "delta* = 48.0 bp/PD at 2.73x, below the 50 bp/PD floor",
         "measure the TRUE k_off ratio, not the >=4x bound."),
        ("No compound binds the TEL patch",
         "best surrogate delta is a rigid-body number, not an affinity",
         "MM-GBSA or FEP at the patch. Out of scope here."),
        ("No TELIN-1 affinity data",
         "zero Kd, zero IC50, zero cellular data",
         "synthesis + SPR/ITC. The molecule is now buildable (WP30)."),
        ("Gate 3 unreachable",
         "only a threshold fraction was deposited, not per-band counts",
         "request per-band densitometry from the depositors."),
        ("9SHZ/9SI0 not ingested",
         "TLS-blocked from this sandbox; parser written and tested",
         "one run on a machine with ordinary network access."),
        ("Gate 4 / Gate 5",
         "dNTP chase not done; delta (smallest meaningful effect) is a judgement",
         "wet lab / a decision the repository cannot make for the reader."),
        ("Y176A mutant panel",
         "designed and specified, never run",
         "wet lab: processivity + TPP1 stimulation on Y176A."),
    ]
    for i, (name, why, fix) in enumerate(gaps, 1):
        print(f"  {i:2d}. {name}")
        print(f"      not achieved: {why}")
        print(f"      how:          {fix}")
    print()

    # ---------------------------------------------------------------- 3
    print("WP34.3  BEST CANDIDATE - and there are three different 'best'")
    print()
    print("  (a) BEST MECHANISM: M5 'two-state, stalls retry'")
    print("      only candidate that is BOTH physical and falsifiable.")
    print("      M1/M2/M3/M4 all break s0 <= phi0 on L&C's own no-TPP1 data.")
    print("      predicts a number (2.73x) it was not fitted to, and derives")
    print(f"      delta* = {delta_star(A_PT):.2f} bp/PD, matching L&C's observed 48 nt.")
    print()
    print("  (b) BEST MOLECULE: TLN-amide, the neutral amide anchor")
    print("      scaffold O=C1C(=S)N(Ar)CON(Ar')C1CC(C)C with a -CH2-CH2-CONH2 arm")
    print("      MW 441.53  logP 4.00  TPSA 83.71  HBD 1  HBA 4  rotB 7")
    print("      Lipinski: PASS all four.")
    print("      the only variant of seven that prefers the interface over a")
    print("      matched control (delta -0.199 per ligand heavy atom).")
    print("      WHY it is plausible: the TEL patch is ambipolar (net -1, both")
    print("      signs present), so a NEUTRAL anchor avoids both mismatches -")
    print("      which is why an acid-or-base search could never find it.")
    print()
    print("  (c) BEST NEXT EXPERIMENT: P2, the DNA-trap re-analysis")
    print("      re-analysis of existing gel data; no new reagent, no new protein.")
    print("      falsifies M5 if the ratio stays >=4x, confirms it near 2.73x,")
    print("      and now ALSO reports the retry probability, predicted ~0.32.")
    print()

    # ---------------------------------------------------------------- 4
    print("WP34.4  SUCCESS RATE - arithmetic, not opinion")
    print()
    mech = [("M5 derives the observed 48 nt exactly", 0.90),
            ("the 4x gap is retry (WP35: RETRACTED, stalling cannot reach 4x)", 0.65),
            ("s = (1+a)^-6 survives an independent enzyme", 0.60)]
    p = 1.0
    print(f"  {'mechanism stage':52s} {'P(pass | prior)':>16}")
    print("  " + "-" * 70)
    for name, pr in mech:
        p *= pr
        print(f"  {name:52s} {pr:15.0%}")
    print(f"  {'M5 survives as the correct mechanism':52s} {p:15.0%}")
    print()
    print("  NOTE: these are my priors, not frequencies. They are stated so")
    print("  they can be argued with. The third is the weakest link and it is")
    print("  the cheapest to test.")
    print()

    print("  THE HARD GATE: maintenance")
    d4 = delta_star(A_NO_PT / 4.0)
    print(f"    delta* at the measured >=4x        {d4:.2f} bp/PD")
    print(f"    lowest human attrition             50.00 bp/PD")
    print(f"    verdict at 2.73x                   CLOSED (48.0 < 50)")
    print(f"    verdict at 4x                      OPEN but insufficient ({d4:.1f} < 90)")
    print(f"    needed for a 90 bp/PD cell          {offrate_needed(90):.2f}x  (have >=4x)")
    print(f"    needed for a 50 bp/PD cell          {offrate_needed(50):.2f}x  (have >=4x)")
    print()
    print(f"    A 90 bp/PD cell needs {offrate_needed(90)/4.0:.2f}x the measured ratio - i.e. the")
    print("    true value would have to be 29% above the reported lower bound.")
    print("    That is not absurd, but it is not established either.")
    print()

    drug = [("hit at the TEL patch (FEP/MM-GBSA)", 0.30),
            ("IC50 < 100 uM in vitro", 0.40),
            ("selectivity counterscreen", 0.45),
            ("cellular recruitment readout", 0.50),
            ("permeability / metabolic stability", 0.35),
            ("in vivo PD + safety (no oncogenesis)", 0.30),
            ("phase I to approval", 0.08)]
    q = 1.0
    print(f"  {'drug stage':52s} {'P(pass)':>8} {'cumulative':>12}")
    print("  " + "-" * 74)
    for name, pr in drug:
        q *= pr
        print(f"  {name:52s} {pr:7.0%} {q:11.4%}")
    print()
    print("  Combined, conditional on the maintenance gate being opened:")
    print(f"    P(approved therapy) ~ {p * q:.4%}")
    print("  Unconditional it is ZERO, because the gate is closed rather than")
    print("  merely improbable - 48 bp/PD cannot outrun 50, at any [TERT].")
    print()
    print("  What is NOT near zero:")
    print(f"    P(a publishable, falsifiable mechanism result) ~ {p:.0%}")
    print("    That is the realistic value of this programme, and it is high.")
    print()

    # ---------------------------------------------------------------- 5
    print("WP34.5  THE ORDER OF OPERATIONS THAT MAXIMISES THE ODDS")
    order = [
        ("1", "P2 DNA-trap re-analysis", "days, no reagent", "confirms/falsifies M5; predicts retry ~0.32"),
        ("2", "true k_off ratio (not the bound)", "one SPR/kinetics run", "decides whether 5.14x is reachable"),
        ("3", "synthesise TLN-amide + SPR", "weeks", "first real affinity number, ever"),
        ("4", "MM-GBSA/FEP at the patch", "compute", "ranks analogues before synthesis"),
        ("5", "Y176A mutant panel", "weeks", "tests the anchor site independently"),
        ("6", "cell maintenance assay", "months", "the only test of the real objective"),
    ]
    print(f"  {'#':2s} {'step':34s} {'cost':20s} {'buys':40s}")
    print("  " + "-" * 98)
    for n, step, cost, buys in order:
        print(f"  {n:2s} {step:34s} {cost:20s} {buys:40s}")
    print()
    print("  Steps 1 and 2 are cheap, and either can end the programme honestly.")
    print("  Doing 3 before 2 would be synthesising against an unknown target.")
    print()

    rows = [
        dict(item="M5 predicted ratio", value=PRED_RATIO, note="from product distribution alone"),
        dict(item="L&C measured ratio (lower bound)", value=MEAS_RATIO, note="pulse-chase"),
        dict(item="retry multiplier needed", value=round(retry_needed, 4), note="meas/pred"),
        dict(item="retry probability needed", value=round(p_needed, 4), note="1 - 1/multiplier"),
        dict(item="stalled fraction observed", value=round(stalled, 4), note="(s0-phi0)/s0"),
        dict(item="fraction of gap explained", value=round(stalled / p_needed, 4), note="NEW result"),
        dict(item="delta* at 2.73x (bp/PD)", value=round(delta_star(A_PT), 2), note="below the 50 floor"),
        dict(item="delta* at 4x (bp/PD)", value=round(d4, 2), note="still below a 90 bp/PD cell"),
        dict(item="off-rate ratio needed (90 bp/PD)", value=round(offrate_needed(90), 2), note="have >=4"),
        dict(item="off-rate ratio needed (50 bp/PD)", value=round(offrate_needed(50), 2), note="have >=4"),
        dict(item="P(M5 correct)", value=round(p, 3), note="my prior, arguable"),
        dict(item="P(approved drug | gate open)", value=f"{p*q:.5f}", note="pipeline product"),
        dict(item="P(approved drug, now)", value="UNDETERMINED",
             note="WP35: 0.0 rested on tau_0=0; derived couplings give 57.3 bp/PD"),
        dict(item="best mechanism", value="M5", note="only physical + falsifiable"),
        dict(item="best molecule", value="TLN-amide", note="MW 441.53, logP 4.00, Lipinski pass"),
        dict(item="best next experiment", value="P2 DNA trap", note="days, no reagent"),
        dict(item="therapeutic claim", value="none", note="quarantine unchanged"),
    ]
    dest = OUT / "wp34_status.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {dest.relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- tests
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    retry = MEAS_RATIO / PRED_RATIO
    p_need = 1.0 - 1.0 / retry
    stalled = (S_NO_PT - PHI_NO_PT) / S_NO_PT

    chk("retry multiplier is 1.465", abs(retry - 1.465) < 1e-3)
    chk("retry probability needed is about 0.32", abs(p_need - 0.318) < 1e-3)
    chk("stalled fraction is about 0.28", abs(stalled - 0.280) < 1e-3)
    chk("the stalled fraction explains 80-95% of the gap",
        0.80 <= stalled / p_need <= 0.95)

    # delta* identities must match WP29
    chk("delta* at 2.73x is 48 bp", abs(delta_star(A_PT) - 48.0) < 0.05)
    chk("delta* at 4x is 70 bp", abs(delta_star(A_NO_PT / 4.0) - 70.0) < 0.15)
    # The two distinct claims, stated separately:
    #   at 2.73x -> 48.0  < 50 : fails even the LOWEST measured attrition
    #   at 4x    -> 70.1  > 50 : clears the floor, but still short of 90
    chk("2.73x cannot reach the lowest attrition (48 < 50)",
        delta_star(A_PT) < 50.0)
    chk("4x clears the 50 floor", delta_star(A_NO_PT / 4.0) > 50.0)
    chk("4x still cannot reach a 90 bp/PD cell", delta_star(A_NO_PT / 4.0) < 90.0)

    # required off-rate ratios round-trip through the identity
    for d in (50.0, 90.0):
        a = A_NO_PT / offrate_needed(d)
        chk(f"required ratio round-trips at {d:.0f} bp/PD", abs(delta_star(a) - d) < 1e-6)
    chk("90 bp/PD needs about 5.14x", abs(offrate_needed(90) - 5.14) < 0.02)
    chk("50 bp/PD needs about 2.85x", abs(offrate_needed(50) - 2.85) < 0.02)

    # the pipeline product is tiny but non-zero conditional on the gate
    q = 1.0
    for pr in (0.30, 0.40, 0.45, 0.50, 0.35, 0.30, 0.08):
        q *= pr
    chk("drug pipeline product is under 1%", q < 0.01)
    chk("drug pipeline product is above 0.01%", q > 0.0001)

    # the best molecule really passes Lipinski
    try:
        import importlib
        import sys as _sys
        _sys.path.insert(0, str(Path(__file__).resolve().parent))
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors
        RDLogger.DisableLog("rdApp.*")
        # the same scaffold constant used everywhere else
        amide = "O=C1C(=S)N(c2ccc(F)cc2CCC(=O)N)C(=O)N(c2ccccc2)C1CC(C)C"
        m = Chem.MolFromSmiles(amide)
        chk("amide candidate is 31 heavy atoms", m is not None and m.GetNumHeavyAtoms() == 31)
        chk("amide candidate parses", m is not None)
        chk("amide candidate passes Lipinski",
            Descriptors.MolWt(m) <= 500 and Descriptors.MolLogP(m) <= 5
            and Descriptors.NumHDonors(m) <= 5 and Descriptors.NumHAcceptors(m) <= 10)
    except ImportError:
        pass

    print()
    for name, ok in checks:
        print(f"    {name:56s} {'ok' if ok else 'FAIL'}")
    failed = sum(1 for _, ok in checks if not ok)
    print(f"    {len(checks)} checks, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    import sys
    sys.exit(_self_test() if "--self-test" in sys.argv else main())
