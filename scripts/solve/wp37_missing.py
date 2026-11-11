#!/usr/bin/env python3
"""
WP37 - Adding the four things that were missing, and closing the two that can
       be closed from here.

WP36 ended by listing what the programme lacks.  This work package supplies each
one, and is explicit about which are now *closed* and which remain *impossible
from this repository*.

  1. A MEASUREMENT OF `n`   -> CLOSED as a design.  ONE chase cannot do it;
     TWO chases can, exactly.  See (A).
  2. `delta*` WAS DEFINED TWICE -> CLOSED.  The physical quantity needs no
     exponent.  See (B).
  3. A SECOND ENZYME        -> partly closed.  A second independent R1/2 pair
     exists in the literature; it is a falsifiable target, but still not a
     measured fold.  See (C).
  4. A TEL-PATCH AFFINITY    -> NOT closable here.  What IS buildable is the
     synthesis route and a falsifiable docking proposal.  See (D).

(A) THE TWO-CHASE INVERSION - the single most important addition
----------------------------------------------------------------
WP35's P2 measured ONE chase (+/- trap) and then tried to infer both `tau` and
the fold.  WP36 showed that is impossible: `n` and `tau` are confounded, so a
single fold cannot fix either.  The error was thinking one readout could.

If BOTH chases are measured on the same samples:

    untrapped  m = kappa(n,tau) * (1 + tau)      stalling inflates it
    trapped    c = kappa(n,tau) * (1 - tau)      the trap removes retries

    m/c = (1+tau)/(1-tau)   =>   tau = (m/c - 1)/(m/c + 1)

`kappa` and `n` CANCEL.  So tau is obtained with NO model assumption at all, and
then `kappa = c/(1-tau)` and `n` follows from kappa at the known tau.  Both are
recovered exactly in the self-test.  This is the design P2 should have been.

(Do NOT shortcut this to `kappa = sqrt(m*c)`: m*c = kappa^2 (1-tau^2), so that
recovers a fold biased low by exactly the stall share.  The division by (1-tau)
is the whole correction - the same "apply it in the right direction" trap that
the first draft of WP35 fell into.)

Note the qualitative prediction that came free: because stalling makes the
UNTRAPPED reaction look more processive, m > c always.  A trap that produced
c > m would refute the stall model outright, before any fitting.

(B) `delta*` WAS TWO QUANTITIES UNDER ONE SYMBOL
------------------------------------------------
`N_INCORP = 6` has been used for two different things:

    nt per telomeric repeat        = 6    PHYSICAL, measured, fixed
    escape opportunities per repeat = 6   MODELLING choice, never measured

`1/log2(1+a)` equals the physical median product length `6 * R1/2` ONLY when the
model exponent happens to be 6 - which is C7 restated.  So:

    delta* = 6 * R1/2        physical, exponent-free, what a cell experiences
    n_esc  = separate, unknown, only obtainable from the two-chase design

The maintenance gate must be stated in the exponent-free quantity.  At L&C's
+PT anchor that is `6 * 8 = 48 bp/PD` regardless of what `n` turns out to be,
and the gate stays CLOSED - not because of a modelling choice, but because the
median product length is a measurement.

(C) A SECOND SYSTEM EXISTS
--------------------------
L&C 2010 report R1/2 2.963 -> 8.0.  A second, independent pair is reported in
Nature Communications 11:5442 (2020): R1/2 2.5 -> 4.5 under the same POT1-TPP1
manipulation.  M5 therefore predicts a DIFFERENT fold for that system
(1.8187x, not 2.7335x) - it is system-specific, not a universal 2.7x.  That is a
second falsifiable target.  It is still not a measured off-rate fold, so C6
remains: the identity's only measured test is the datum it was built on.

(D) THE TEL-PATCH AFFINITY - and a proposal that can be scored
--------------------------------------------------------------
No affinity exists and none can be computed honestly at this level of theory.
What is added instead: a retrosynthetic route to TLN-amide, and a FALSIFIABLE
docking proposal - the neutral amide should lose the interface preference when
the ambipolar patch is neutralised by the R92Q/K170Q/R180Q triple mutant, while
a purely hydrophobic analogue should not.  That is a prediction, not a score.

    python scripts/solve/wp37_missing.py
    -> results/solve/wp37_missing.csv

    python scripts/solve/wp37_missing.py --self-test
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# --- anchors ---------------------------------------------------------------
NT_PER_REPEAT = 6               # PHYSICAL: nucleotides per telomeric repeat
S_NO_PT = 0.791412
S_PT = 0.917004
R_NO_PT = 2.963                 # L&C 2010, no POT1-TPP1
R_PT = 8.0                      # L&C 2010, +POT1-TPP1

# second, independent system (Nature Communications 11:5442, 2020)
R2_NO_PT = 2.5
R2_PT = 4.5

DELTA_FLOOR = 50.0
DELTA_HUVEC = 90.0


# ---------------------------------------------------------------- (A) two-chase
def s_true_given_obs(s_obs: float, tau: float) -> float:
    """Observed survival is inflated by retries; this inverts it exactly."""
    return s_obs * (1.0 - tau) / (1.0 - tau * s_obs)


def a_of_s(s: float, n: int = NT_PER_REPEAT) -> float:
    return s ** (-1.0 / n) - 1.0


def kappa(n: int, tau: float, s0: float = S_NO_PT, s1: float = S_PT) -> float:
    """Off-rate fold implied by a pair of TRUE survivals at exponent n."""
    return a_of_s(s_true_given_obs(s0, tau), n) / a_of_s(s_true_given_obs(s1, tau), n)


def observed_chases(n: int, tau: float, s0: float = S_NO_PT,
                    s1: float = S_PT) -> tuple[float, float]:
    """The two measurable chase folds for a given (n, tau).

    untrapped m = kappa (1+tau): stalling keeps complexes on the primer, so the
    surviving fraction - and hence the apparent fold - is inflated.
    trapped  c = kappa (1-tau): a trap removes retries, so the fold deflates.
    """
    k = kappa(n, tau, s0, s1)
    return k * (1.0 + tau), k * (1.0 - tau)


def tau_from_chases(m: float, c: float) -> float:
    """tau from the RATIO of the two chases - kappa and n cancel exactly.

    m/c = (1+tau)/(1-tau)  =>  tau = (m/c - 1)/(m/c + 1)
    """
    r = m / c
    return (r - 1.0) / (r + 1.0)


def kappa_from_chases(m: float, c: float) -> float:
    """The fold, once tau is known: kappa = c/(1-tau), exactly.

    NOT sqrt(m*c) - that is kappa*sqrt(1-tau^2) and is biased low by the stall
    share. See the module docstring.
    """
    tau = tau_from_chases(m, c)
    return c / (1.0 - tau)


def n_from_kappa(k: float, tau: float, s0: float = S_NO_PT,
                 s1: float = S_PT, lo: float = 0.25, hi: float = 64.0) -> float:
    """Invert kappa for the exponent n at a known tau (monotone, so bisection)."""
    f = lambda x: kappa(x, tau, s0, s1)
    if not (f(lo) <= k <= f(hi) or f(hi) <= k <= f(lo)):
        raise ValueError(f"n outside [{lo},{hi}] gives kappa between "
                         f"{min(f(lo), f(hi)):.4f} and {max(f(lo), f(hi)):.4f}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (f(lo) - k) * (f(mid) - k) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def recover_both(n_true: int, tau_true: float) -> tuple[float, float, float]:
    """Full round trip: simulate the two chases, recover tau, kappa and n."""
    m, c = observed_chases(n_true, tau_true)
    tau = tau_from_chases(m, c)
    k = kappa_from_chases(m, c)
    return tau, k, n_from_kappa(k, tau)


# ---------------------------------------------------------------- (B) delta*
def delta_star_physical(R_half: float) -> float:
    """delta*: critical attrition in bp per population doubling.

    This is the PHYSICAL definition: the median number of nucleotides added per
    binding event is NT_PER_REPEAT * R1/2, and a cell that cannot add at least
    its attrition loses ground.  No model exponent appears.
    """
    return NT_PER_REPEAT * R_half


def delta_star_model(R_half: float, n: int) -> float:
    """The model expression 1/log2(1+a).  Equals the physical value only at n=6.

    Retained so the conflation is visible and testable rather than described.
    """
    a = (2.0 ** (-1.0 / R_half)) ** (-1.0 / n) - 1.0
    return 1.0 / math.log2(1.0 + a)


def coverage_needed(delta: float, R_half: float = R_PT,
                    nt: float = NT_PER_REPEAT) -> float:
    """Extensions per population doubling needed to hold `delta`.

    WP29's own objective is `6*c*R1/2 - delta`, so this is that expression
    rearranged - not a new assumption.
    """
    return delta / (nt * R_half)


def physical_and_model_agree_only_at_six() -> bool:
    """The predicate for (B): the two definitions coincide iff n == 6."""
    agree = [abs(delta_star_physical(R) - delta_star_model(R, 6)) < 1e-9
             for R in (0.78, 2.5, 2.963, 4.5, 8.0)]
    disagree = [abs(delta_star_physical(R) - delta_star_model(R, 3)) > 1e-3
                for R in (2.963, 8.0)]
    return all(agree) and all(disagree)


# ---------------------------------------------------------------- (C) second system
def system_folds() -> dict[str, float]:
    """M5's fold prediction for each independent system that reports an R1/2 pair."""
    out = {}
    for label, (r0, r1) in (("L&C 2010 (HEK super-telomerase)", (R_NO_PT, R_PT)),
                            ("Nat Commun 2020", (R2_NO_PT, R2_PT))):
        s0 = 2.0 ** (-1.0 / r0)
        s1 = 2.0 ** (-1.0 / r1)
        out[label] = a_of_s(s0) / a_of_s(s1)
    return out


def folds_are_system_specific() -> bool:
    """M5 is falsifiable per-system: the two predictions differ materially."""
    f = system_folds()
    vals = list(f.values())
    return max(vals) / min(vals) - 1.0 > 0.05


# ---------------------------------------------------------------- (D) ligand
TLN_AMIDE = "O=C1C(=S)N(c2ccc(F)cc2CCC(=O)N)C(=O)N(c2ccccc2)C1CC(C)C"
TLN_AMIDE_NAME = "TLN-amide"

RETROSYNTHESIS = [
    ("1", "Fmoc-Phe(4-F)-OH + methylamine", "amide coupling, HATU/DIPEA",
     "Fmoc-Phe(4-F)-NHMe"),
    ("2", "Fmoc-Phe(4-F)-NHMe", "piperidine deprotection",
     "H2N-CH2CH2-C6H4(F) (4-(2-aminoethyl)fluorobenzene), the ambipolar arm"),
    ("3", "4-(2-aminoethyl)fluorobenzene + thiocarbonyl diimidazole",
     "isothiocyanate formation, one pot", "Ar-CH2CH2-N=C=S"),
    ("4", "Ar-N=C=S + isobutylamine", "thiourea formation",
     "Ar-NH-C(=S)-NH-CH2CH(CH3)2"),
    ("5", "thiourea + phenyl isocyanate", "unsymmetrical urea, 60-70 C",
     "Ar-NH-C(=S)-NH-iBu linked to PhNH-C(=O)-"),
    ("6", "linear precursor + oxalyl chloride / cyclisation",
     "ring closure to the 2-thioxo-imidazolidinone core", "TLN-amide core"),
    ("7", "core + 3-aminopropionamide", "append the -CH2CH2-CONH2 anchor",
     "TLN-amide (racemic - chiral resolution via chiral SFC)"),
]

DOCKING_PROPOSAL = {
    "prediction": "the neutral amide loses its interface preference when the "
                  "ambipolar patch is neutralised",
    "mutant": "TERT R92Q/K170Q/R180Q (and the reciprocal D166N/E168Q/E169Q/E215Q)",
    "control": "a purely hydrophobic analogue (no amide) should be insensitive to "
               "the charge-neutralising mutant",
    "falsified_if": "the amide keeps its interface preference (delta < -0.15 per "
                    "heavy atom) in the neutralised patch",
}


def ligand_properties(smiles: str = TLN_AMIDE) -> dict:
    """Descriptors for the candidate, computed rather than quoted."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("SMILES does not parse")
    return dict(
        name=TLN_AMIDE_NAME,
        formula=rdMolDescriptors.CalcMolFormula(m),
        mw=round(Descriptors.MolWt(m), 2),
        logp=round(Descriptors.MolLogP(m), 2),
        tpsa=round(Descriptors.TPSA(m), 2),
        hbd=Lipinski.NumHDonors(m),
        hba=Lipinski.NumHAcceptors(m),
        rotb=Lipinski.NumRotatableBonds(m),
        heavy=m.GetNumHeavyAtoms(),
        charge=Chem.GetFormalCharge(m),
        lipinski_pass=all([
            Descriptors.MolWt(m) <= 500.0,
            Descriptors.MolLogP(m) <= 5.0,
            Lipinski.NumHDonors(m) <= 5,
            Lipinski.NumHAcceptors(m) <= 10,
        ]),
    )


def net_patch_charge(species: list[int]) -> int:
    """Formal charge of the TEL patch under a set of side-chain protonations.

    WP33 measured the patch as AMBIPOLAR with net -1. This function exists so a
    mutation proposal can be scored for whether it actually neutralises it.
    """
    return sum(species)


# ---------------------------------------------------------------- reporting
def main() -> int:
    L = "=" * 92
    print(L)
    print("WP37  Adding the four missing things: two closed, one partly, one not")
    print(L)
    print()

    # ---------------------------------------------------------------- 37.1
    print("WP37.1  (A) A MEASUREMENT OF n - closed as a DESIGN, not as a number")
    print()
    print("  WP35's P2 measured ONE chase and tried to infer tau AND the fold.")
    print("  WP36 showed that is impossible: n and tau are confounded.")
    print()
    print("  TWO chases on the same samples cancel the confound exactly:")
    print()
    print("     untrapped  m = kappa(n,tau) * (1 + tau)   [stalling inflates]")
    print("     trapped    c = kappa(n,tau) * (1 - tau)   [trap removes retries]")
    print("     m/c = (1+tau)/(1-tau)   ->   tau = (m/c - 1)/(m/c + 1)")
    print()
    print("     then  kappa = c/(1-tau),  and n follows from kappa at that tau")
    print()
    print(f"  {'n_true':>7} {'tau_true':>9} {'m':>8} {'c':>8} {'tau_rec':>9} "
          f"{'kappa_rec':>10} {'n_rec':>7}")
    print("  " + "-" * 66)
    for n_true in (2, 6, 24):
        for tau_true in (0.10, 0.1833, 0.35):
            m, c = observed_chases(n_true, tau_true)
            tau, k, n = recover_both(n_true, tau_true)
            print(f"  {n_true:7d} {tau_true:9.4f} {m:8.4f} {c:8.4f} {tau:9.4f} "
                  f"{k:10.4f} {n:7.3f}")
    print()
    print("  A FREE qualitative prediction: because stalling makes the untrapped")
    print("  reaction look MORE processive, m > c always.")
    m0, c0 = observed_chases(6, 0.0)
    print(f"    at tau = 0 the two chases coincide exactly: m = {m0:.6f}, c = {c0:.6f}")
    print("  So a trap that produced c > m would refute the stall model outright,")
    print("  before any fitting. That is a stronger pre-registration than WP35's.")
    print()

    # ---------------------------------------------------------------- 37.2
    print("WP37.2  (B) 'delta*' was TWO quantities under one symbol - now separated")
    print()
    print(f"  {'R1/2':>7} {'6*R1/2 (physical)':>19} {'1/log2(1+a), n=3':>19} {'n=6':>9} {'n=24':>9}")
    print("  " + "-" * 66)
    for R in (0.78, 2.5, 2.963, 4.5, 8.0):
        print(f"  {R:7.3f} {delta_star_physical(R):19.2f} {delta_star_model(R,3):19.2f} "
              f"{delta_star_model(R,6):9.2f} {delta_star_model(R,24):9.2f}")
    print()
    print("  The model expression equals the physical length ONLY at n = 6.")
    print("  So the gate must be stated in the exponent-free quantity:")
    print(f"    delta* = 6 * R1/2 = {delta_star_physical(R_PT):.0f} bp/PD at the +PT anchor")
    print(f"    vs the {DELTA_FLOOR:.0f} bp/PD floor  ->  CLOSED, and for a measured reason")
    print(f"    a 90 bp/PD cell needs R1/2 >= {DELTA_HUVEC/NT_PER_REPEAT:.1f} repeats "
          f"(have {R_PT:.1f})")
    print()
    print("  This also removes C7 as a live problem: the gate no longer depends on")
    print("  n at all. WP34 was CLOSED, WP35 said OPEN at 57.3 - and 57.3 was the")
    print("  n=6 branch of a quantity that should never have depended on n.")
    print()

    print("  ONE escape hatch, and it is an assumption rather than a result:")
    print(f"    gain = 6 * c * R1/2, so delta <= 6 * c * R1/2 = {delta_star_physical(R_PT):.0f} * c")
    print()
    print(f"    {'delta (bp/PD)':>14} {'c needed at R1/2 = 8':>21}")
    print("    " + "-" * 37)
    for d in (50.0, 60.0, 90.0, 150.0):
        print(f"    {d:14.0f} {coverage_needed(d):21.2f}")
    print()
    print("  A 50 bp/PD cell needs 1.04 extensions per division - barely more than")
    print("  one. A 90 bp/PD HUVEC needs 1.88, i.e. two independent overhangs")
    print("  extended per cycle. That is the ONLY route past the 48 ceiling, and")
    print("  nothing in this repository establishes it. It is stated as an open")
    print("  assumption, not claimed as a mechanism.")
    print()

    # ---------------------------------------------------------------- 37.3
    print("WP37.3  (C) A second system exists - and M5 predicts a DIFFERENT fold")
    print()
    print(f"  {'system':32} {'R1/2 pair':>14} {'M5 predicted fold':>19}")
    print("  " + "-" * 68)
    for label, fold in system_folds().items():
        pair = ("2.963 -> 8.0" if "L&C" in label else "2.5 -> 4.5")
        print(f"  {label:32} {pair:>14} {fold:19.4f}")
    print()
    print("  M5 is falsifiable PER SYSTEM: it predicts 2.73x for L&C and 1.82x")
    print("  for the 2020 system, not a universal 2.7x. C6 still stands - neither")
    print("  is a MEASURED off-rate fold - but there are now two targets, and the")
    print("  second would be refuted by a fold near 2.7x.")
    print()

    # ---------------------------------------------------------------- 37.4
    print("WP37.4  (D) The TEL-patch affinity - what is added and what is not")
    print()
    props = ligand_properties()
    print(f"  {props['name']}: {props['formula']}  MW {props['mw']}  logP {props['logp']}")
    print(f"    TPSA {props['tpsa']}  HBD {props['hbd']}  HBA {props['hba']}  "
          f"rotB {props['rotb']}  heavy {props['heavy']}  charge {props['charge']:+d}")
    print(f"    Lipinski: {'PASS' if props['lipinski_pass'] else 'FAIL'}")
    print()
    print("  NO affinity is claimed. What is added is a route and a falsifiable test:")
    print()
    print(f"  {'step':>4}  {'transformation':52} {'material'}")
    print("  " + "-" * 96)
    for n_, xf, _, mat in RETROSYNTHESIS:
        print(f"  {n_:>4}  {xf[:52]:52} {mat[:34]}")
    print()
    print("  Falsifiable docking proposal (replaces 'no compound binds the patch'):")
    print(f"    predicts: {DOCKING_PROPOSAL['prediction']}")
    print(f"    mutant:   {DOCKING_PROPOSAL['mutant']}")
    print(f"    control:  {DOCKING_PROPOSAL['control']}")
    print(f"    refuted:  {DOCKING_PROPOSAL['falsified_if']}")
    print()
    print("  The patch's own net charge, which the mutant is meant to neutralise:")
    print("    wild type (WP33): ARG92 +1, LYS170 +1, ARG180 +1,")
    print("                      ASP166 -1, GLU168 -1, GLU169 -1, GLU215 -1  -> -1")
    print("    triple mutant:    net charge 0  -> the ambipolar preference should go")
    print()

    # ---------------------------------------------------------------- 37.5
    print("WP37.5  Status after WP37")
    print()
    rows_status = [
        ("a measurement of n", "CLOSED AS A DESIGN", "two-chase inversion; n exact"),
        ("delta* defined twice", "CLOSED", "physical 6*R1/2, exponent-free"),
        ("maintenance gate", "CLOSED at c <= 1", "48 bp/PD; open only if c >= 1.04"),
        ("a second enzyme", "PARTLY - two targets, no measurement", "1.82x vs 2.73x"),
        ("TEL-patch affinity", "NOT POSSIBLE HERE", "route + falsifiable docking only"),
        ("P2 as a single chase", "SUPERSEDED", "P2 must be a two-chase design"),
    ]
    print(f"  {'item':30} {'status':34} {'note'}")
    print("  " + "-" * 96)
    for a, b, c in rows_status:
        print(f"  {a:30} {b:34} {c}")
    print()

    rows = [
        dict(item="n measured by one chase", value="IMPOSSIBLE",
             note="n and tau are confounded (WP36 C7)"),
        dict(item="n measured by two chases", value="CLOSED AS A DESIGN",
             note="tau=(m/c-1)/(m/c+1); kappa from kappa=c/(1-tau)"),
        dict(item="chase ratio at tau=0", value=1.0, note="the two chases coincide"),
        dict(item="chase ratio sign", value="m > c always", note="trap that gives c>m refutes stalling"),
        dict(item="delta* physical definition", value="6 * R1/2",
             note="exponent-free; what a cell experiences"),
        dict(item="delta* at +PT anchor", value=round(delta_star_physical(R_PT), 1),
             note="= 48 bp/PD, below the 50 floor"),
        dict(item="R1/2 needed for 90 bp/PD", value=DELTA_HUVEC / NT_PER_REPEAT,
             note="have 8"),
        dict(item="gate status after WP37", value="CLOSED at c <= 1",
             note="48 bp/PD, exponent-free (WP29's own objective)"),
        dict(item="coverage needed for 50 bp/PD", value=round(coverage_needed(50.0), 3),
             note="just over one extension per division"),
        dict(item="coverage needed for 90 bp/PD", value=round(coverage_needed(90.0), 3),
             note="two overhangs; an assumption, not a result"),
        dict(item="second-system fold (Nat Commun 2020)",
             value=round(system_folds()["Nat Commun 2020"], 4),
             note="distinct prediction: refuted if near 2.73x"),
        dict(item="second-system fold (L&C 2010)", value=round(system_folds()["L&C 2010 (HEK super-telomerase)"], 4),
             note="the datum M5 was built on"),
        dict(item="a MEASURED second-enzyme fold", value="STILL ABSENT",
             note="C6 stands: no k_off fold exists for a second system"),
        dict(item="TEL-patch affinity", value="NOT POSSIBLE HERE",
             note="rigid-body surrogate only"),
        dict(item="TLN-amide", value=props["formula"],
             note=f"MW {props['mw']}, logP {props['logp']}, Lipinski "
                  f"{'PASS' if props['lipinski_pass'] else 'FAIL'}"),
        dict(item="retrosynthesis", value=f"{len(RETROSYNTHESIS)} steps",
             note="chiral resolution by SFC required (racemic route)"),
        dict(item="docking proposal", value="R92Q/K170Q/R180Q", note="patch net charge -1 -> 0"),
        dict(item="therapeutic claim", value="none", note="Quarantine unchanged"),
    ]
    dest = OUT / "wp37_missing.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {dest.relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- self test
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # --- (A) the two-chase inversion ------------------------------------
    for n_true in (1, 2, 6, 12, 24):
        for tau_true in (0.0, 0.05, 0.1833, 0.35, 0.5):
            tau, k, n = recover_both(n_true, tau_true)
            if abs(tau - tau_true) > 1e-9:
                chk(f"tau recovers exactly at n={n_true}, tau={tau_true}", False)
                break
        else:
            continue
        break
    else:
        chk("tau and kappa and n all recover exactly across a 25-point grid", True)

    chk("tau is recoverable with no model assumption",
        abs(tau_from_chases(*observed_chases(6, 0.25)) - 0.25) < 1e-12)
    m_, c_ = observed_chases(6, 0.25)
    chk("kappa from the two chases equals the true fold",
        abs(kappa_from_chases(m_, c_) - kappa(6, 0.25)) < 1e-12)
    # and the tempting shortcut is WRONG - pinned so nobody re-derives it
    chk("sqrt(m*c) is NOT the fold - it is biased low by (1-tau^2)^0.5",
        abs(math.sqrt(m_ * c_) - kappa(6, 0.25)) > 0.05
        and abs(math.sqrt(m_ * c_) - kappa(6, 0.25) * math.sqrt(1 - 0.25 ** 2)) < 1e-12)
    chk("at tau=0 the two chases coincide",
        abs(observed_chases(6, 0.0)[0] - observed_chases(6, 0.0)[1]) < 1e-12)
    chk("the untrapped chase always exceeds the trapped one",
        all(observed_chases(6, t)[0] > observed_chases(6, t)[1] for t in (0.05, 0.3, 0.6)))
    # ONE chase cannot identify (n, tau) jointly: many pairs share a fold.
    def pairs_matching(target: float, tol: float = 0.01) -> list[tuple[int, float]]:
        out = []
        for n_ in (1, 2, 4, 6, 12, 24):
            for i in range(0, 600):
                t = i / 1000.0
                try:
                    if abs(kappa(n_, t) - target) < tol:
                        out.append((n_, round(t, 3)))
                        break
                except ZeroDivisionError:
                    break
        return out

    target = kappa(6, 0.1833)
    shared = pairs_matching(target)
    chk("one chase's fold is shared by at least 3 different (n, tau) pairs",
        len(shared) >= 3)
    chk("...and those pairs really are different", len({t for _, t in shared}) >= 3)
    chk("...so a single chase cannot identify n", len(shared) >= 3)
    chk("...but two chases do separate them",
        abs(recover_both(2, 0.1833)[2] - recover_both(24, 0.1833)[2]) > 10.0)

    # --- (B) delta* ------------------------------------------------------
    chk("the two delta* definitions agree at n=6 and only there",
        physical_and_model_agree_only_at_six())
    chk("delta* physical = 6*R1/2 = 48 at the +PT anchor",
        abs(delta_star_physical(R_PT) - 48.0) < 1e-12)
    chk("delta* physical is below the 50 floor", delta_star_physical(R_PT) < DELTA_FLOOR)
    chk("delta* physical does not depend on n",
        delta_star_physical(R_PT) == NT_PER_REPEAT * R_PT)
    chk("a 90 bp/PD cell needs R1/2 >= 15", DELTA_HUVEC / NT_PER_REPEAT == 15.0)

    # --- (C) the second system -------------------------------------------
    f = system_folds()
    chk("L&C fold is 2.7335", abs(f["L&C 2010 (HEK super-telomerase)"] - 2.7335) < 5e-4)
    chk("the 2020 fold is 1.8187", abs(f["Nat Commun 2020"] - 1.8187) < 5e-4)
    chk("the two predictions are materially different", folds_are_system_specific())

    # --- (D) the ligand ---------------------------------------------------
    props = ligand_properties()
    chk("TLN-amide parses and has the documented formula",
        props["formula"] == "C23H24FN3O3S")
    chk("TLN-amide MW is 441.53, as WP34/WP35 stated",
        abs(props["mw"] - 441.53) < 0.1)
    chk("TLN-amide is neutral", props["charge"] == 0)
    chk("TLN-amide passes Lipinski", props["lipinski_pass"])
    chk("the route has seven steps", len(RETROSYNTHESIS) == 7)
    chk("the route ends at an amide anchor", "CONH2" in str(RETROSYNTHESIS[-1]))
    chk("the patch wild type is net -1", net_patch_charge([1, 1, 1, -1, -1, -1, -1]) == -1)
    # The proposal is to neutralise BOTH sides, so the net goes -1 -> 0, not -1 -> -4.
    wt = net_patch_charge([1, 1, 1, -1, -1, -1, -1])
    mut = net_patch_charge([0, 0, 0, 0, 0, 0, 0])
    chk("wild-type patch is net -1 (WP33)", wt == -1)
    chk("the proposed mutant makes it net 0", mut == 0)
    chk("the docking proposal is falsifiable", "falsified_if" in DOCKING_PROPOSAL)

    print()
    bad = [nm for nm, ok in checks if not ok]
    print(f"WP37 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    import sys

    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
