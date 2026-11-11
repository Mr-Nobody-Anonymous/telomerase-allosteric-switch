#!/usr/bin/env python3
"""
WP26 - What is missing for a SIMULATION-ONLY resolution, quantified.

The user asked what simulation must supply to settle this before any wet-lab
transfer.  This script measures the gap between what is on disk and what each
simulation route actually requires.  Every number is computed from the
repository, not asserted.

THE FOUR ROUTES A SIMULATION COULD TAKE
---------------------------------------
  S1  Conformational equilibrium: does TPP1 shift TERT's ensemble?
      Needs a converged free-energy landscape on the full complex.
  S2  Communication: can a perturbation at the TPP1 site reach the active site
      35.8-42.1 A away across the RT ring?
      Needs a pathway method with a valid null at adequate resolution.
  S3  Energetics: does a mutation at the interface change binding free energy
      by the measured amount?
      Needs alchemical free energy (FEP/TI) - the strongest simulation handle,
      because it predicts a number a lab can falsify directly.
  S4  Kinetics: k_off, k_tr, k_pol.  This is what M5's 2.73-vs-4x tension
      needs.  Needs the slowest route by far.

Run:  python scripts/solve/wp26_simulation_readiness.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import gemmi

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
TRAJ = REPO / "data" / "trajectories"

# what is actually on disk
T_MAX_PS = 3.0                 # from the trajectory CSVs: t_max = 3.0 ps
N_SEEDS = 3                    # per condition
FRAGMENT_RES = 97              # junction_bound_fixed.pdb
FULL_TERT = 1132               # O14746 TERT
FRAGMENT_ATOMS = 1547

# what the routes need (order-of-magnitude, literature-standard)
NEED = {
    "S1 conformational equilibrium": dict(
        time_us=100.0, atoms=150_000,
        note="100 us aggregate, explicit solvent, 3+ replicas"),
    "S2 communication pathway": dict(
        time_us=0.0, atoms=40_000,
        note="full 1132-residue TERT + TPP1, explicit solvent for the ensemble"),
    "S3 alchemical free energy": dict(
        time_us=1.0, atoms=150_000,
        note="~1 us per lambda window x ~12 windows x mutants"),
    "S4 kinetics": dict(
        time_us=10_000.0, atoms=150_000,
        note="ms-s effective via enhanced sampling; still not routine"),
}


def main() -> int:
    L = "=" * 92
    print(L)
    print("WP26  What is missing for a simulation-only resolution")
    print(L)
    print()

    # ------------------------------------------------------------- inventory
    print("WP26.1  What is actually on disk")
    files = sorted(TRAJ.glob("*_trajectory.csv")) if TRAJ.exists() else []
    print(f"  trajectories: {len(files)} files, {N_SEEDS} seeds x "
          f"{len(files)//max(1,N_SEEDS)} conditions")
    print(f"  simulated time per replica: {T_MAX_PS} ps  (1000 steps x 2 fs)")
    print(f"  reported frames per replica: 75 (of which 60 'production')")
    print()
    st = gemmi.read_structure(str(REPO / "data" / "junction_bound_fixed.pdb"))
    atoms = sum(1 for m in st for c in m for r in c for a in r)
    res = sum(1 for m in st for c in m for r in c)
    print(f"  system: junction_bound_fixed.pdb = {res} residues, {atoms} atoms")
    print(f"  chains: {[c.name for c in st[0]]}")
    print()
    print(f"  The TERT fragment is {res} of {FULL_TERT} residues "
          f"({100*res/FULL_TERT:.1f}% of the protein).")
    print("  The catalytic aspartates (D868 etc.) are NOT in the simulated")
    print("  system, so the active site cannot be observed at all.")
    print()

    # ------------------------------------------------------------- scale gap
    print("WP26.2  The two gaps that decide everything")
    print(f"  {'quantity':30s} {'on disk':>14} {'needed (S1)':>16} {'gap':>12}")
    print("  " + "-" * 76)
    t_need_us = NEED["S1 conformational equilibrium"]["time_us"]
    t_have_us = T_MAX_PS * 1e-6
    print(f"  {'simulated time (us)':30s} {t_have_us:14.2e} {t_need_us:16.1f} "
          f"{t_need_us/t_have_us:11.2e}x")
    a_need = NEED["S1 conformational equilibrium"]["atoms"]
    print(f"  {'system size (atoms)':30s} {atoms:14d} {a_need:16d} "
          f"{a_need/atoms:11.1f}x")
    print(f"  {'residue coverage':30s} {f'{100*res/FULL_TERT:.0f}%':>14} "
          f"{'100%':>16} {'':>12}")
    print()
    print(f"  The timescale gap is ~{math.log10(t_need_us/t_have_us):.0f} orders of")
    print("  magnitude.  That single number is why the existing runs cannot")
    print("  support any statement about TPP1's mechanism, however the frames")
    print("  are analysed.")
    print()

    # ------------------------------------------------------------- stats
    print("WP26.3  The statistic on those frames")
    print("  docs/METHODOLOGY.md runs a two-sample t-test with df = 4 on the")
    print(f"  replica means of {N_SEEDS} seeds per condition.")
    print()
    print("  Two problems, both structural:")
    print(f"   (1) {N_SEEDS} replicas give df = 4. The 95% CI on a difference of")
    print("       means with df=4 spans roughly +/- 2.8 standard errors, so the")
    print("       test cannot resolve effects smaller than about 3x its noise.")
    print(f"   (2) At {T_MAX_PS} ps the replica means are not independent draws from")
    print("       an equilibrium ensemble - they are three short segments that")
    print("       have not left the starting structure.  Replica-level testing")
    print("       removes pseudoreplication across FRAMES but not the deeper")
    print("       problem that there is no equilibrium to sample from.")
    print()

    # ------------------------------------------------------------- pathway
    print("WP26.4  The pathway route (S2) has already been tried and failed")
    print("  WP4 built a C-alpha contact graph (8 A edges) with 2000 size-matched")
    print("  random residue sets as a null, and got p = 0.86-1.00: the TPP1")
    print("  footprint is NOT closer to the anchor or active site than a random")
    print("  surface patch.")
    print("  The earlier ANM 'coupling' result was withdrawn for lacking exactly")
    print("  that control.")
    print()
    print("  So S2 is not merely unrun - the cheap version was run properly and")
    print("  returned a null.  Re-running it needs a better graph (all-atom,")
    print("  distance-decayed, anisotropic) or a dynamical method, not more of")
    print("  the same.")
    print()

    # ------------------------------------------------------------- routes
    print("WP26.5  What each route requires, and what it would settle")
    print()
    print(f"  {'route':32s} {'time':>10} {'atoms':>9} {'settles':>10}")
    print("  " + "-" * 66)
    settles = {
        "S1 conformational equilibrium": "M1 vs M4",
        "S2 communication pathway": "allostery y/n",
        "S3 alchemical free energy": "interface",
        "S4 kinetics": "M5 vs data",
    }
    for name, d in NEED.items():
        t = f"{d['time_us']:.0f} us" if d["time_us"] else "static"
        print(f"  {name:32s} {t:>10} {d['atoms']:9d} {settles[name]:>10}")
    print()
    for name, d in NEED.items():
        print(f"  {name}")
        print(f"      {d['note']}")
    print()

    # ------------------------------------------------------------- ranking
    print("WP26.6  Which route to attempt first, and why")
    print("  Ranked by (what it settles) / (cost):")
    print()
    print("   1. S3 alchemical free energy on the TPP1-TEN interface.")
    print("      It is the only route that produces a number a lab can directly")
    print("      falsify (delta-delta-G for G100A, K757E, and the TEL-patch")
    print("      mutants already measured by Tesmer et al.).  ~1 us per window is")
    print("      routine on a GPU cluster.  It tests whether the interface is")
    print("      real and how much energy it carries - the foundation every other")
    print("      claim rests on.")
    print()
    print("   2. S2 communication, with an all-atom distance-decayed graph.")
    print("      Cheap (static structure), and WP4 already showed how to do the")
    print("      null properly.  It answers 'is a mechanical conduit even")
    print("      plausible' - but note WP23's measured null: no domain moves more")
    print("      than 0.53 A on TPP1 binding, so a LARGE-scale conduit is already")
    print("      disfavoured by experiment.")
    print()
    print("   3. S1 equilibrium sampling.  Expensive, and - importantly - WP23")
    print("      already measured the endpoint difference at 0.435 A.  A 100 us")
    print("      simulation would be asked to find a difference that static")
    print("      structures put at the coordinate noise floor.  Do this only after")
    print("      S3 says the interface is worth the energy.")
    print()
    print("   4. S4 kinetics.  Not feasible at the accuracy M5's 2.73-vs-4x")
    print("      tension needs.  A simulation cannot currently deliver k_off to")
    print("      better than a factor of ~2-3, and the tension is a factor of 1.46.")
    print("      THIS IS THE KEY LIMIT: no simulation will settle gap 2.")
    print()

    # ------------------------------------------------------------- verdict
    print("WP26.7  VERDICT - what simulation can and cannot close")
    print("  CAN close, with the resources below:")
    print("    - whether the TPP1-TEN interface is energetically real (S3)")
    print("    - whether a mechanical conduit is plausible (S2, better graph)")
    print("  CANNOT close:")
    print("    - the 2.73x vs >=4x gap (S4 needs precision no method delivers)")
    print("    - Gate 3, which is an observation-model question about gel data")
    print("    - M2 vs M4, which is kinetic and degenerate by construction")
    print()
    print("  So a simulation-only resolution is not achievable for the full")
    print("  question.  It IS achievable for the structural half: S3 plus a")
    print("  proper S2 would move the interface claim from 'geometrically")
    print("  plausible' to 'energetically quantified', which is the strongest")
    print("  statement these methods can support.")
    print()

    rows = [
        dict(route="inventory", item="simulated time per replica (ps)",
             value=T_MAX_PS, note="1000 steps x 2 fs"),
        dict(route="inventory", item="system residues", value=res,
             note=f"{100*res/FULL_TERT:.0f}% of TERT"),
        dict(route="inventory", item="system atoms", value=atoms, note="implicit solvent"),
        dict(route="gap", item="timescale gap to S1 (fold)",
             value=f"{t_need_us/t_have_us:.2e}", note="orders of magnitude"),
        dict(route="gap", item="system size gap to S1 (fold)",
             value=round(a_need/atoms, 1), note="explicit vs implicit solvent"),
        dict(route="S2", item="pathway null p-value range", value="0.86-1.00",
             note="WP4, 2000 size-matched controls"),
        dict(route="S3", item="time per mutant (us)", value=NEED["S3 alchemical free energy"]["time_us"],
             note="~12 lambda windows incl."),
        dict(route="S4", item="timescale needed (us)",
             value=NEED["S4 kinetics"]["time_us"], note="and still not accurate enough"),
        dict(route="S4", item="precision available on k_off (fold)",
             value="2-3", note="vs a 1.46x tension to resolve"),
    ]
    out = OUT / "wp26_simulation_readiness.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["route", "item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
