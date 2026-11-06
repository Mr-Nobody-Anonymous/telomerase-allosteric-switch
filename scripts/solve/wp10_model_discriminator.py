#!/usr/bin/env python3
"""
WP10 - Making the four models produce different measurements.

Requested: stop asking which mechanism "looks" right and design a test where retention, allostery
and translocation predict DIFFERENT observables.  This script derives that test exactly, and it
turns out to depend on two measurements that are already standard.

Model class (minimal, six nucleotide additions then one translocation per repeat):

    s  = [k_pol/(k_pol + k_off)]^6 * k_tr/(k_tr + k_off)        per-repeat survival
    phi = k_tr/(k_tr + k_off)                                   translocation efficiency

Both are measurable in the same reaction: s from the product-length distribution, phi from the
single-round translocation assay.

EXACT IDENTITY (sympy-verified, parameter-free):

    log(s) - log(phi) = 6*log(k_pol) - 6*log(k_pol + k_off) = -6*log(1 + k_off/k_pol)

so s/phi depends ONLY on the ratio a = k_off/k_pol.  In particular

    d/dlog(k_tr)  [ log s - log phi ] = 0        exactly
    d/dlog(k_off) [ log s - log phi ] = -6k_off/(k_off+k_pol)
    d/dlog(k_pol) [ log s - log phi ] = +6k_off/(k_off+k_pol)

That yields a discrimination table in which each model makes a DIFFERENT prediction:

    model            TPP1 changes    phi          s            s/phi
    M1 retention     k_off down       rises        rises        RISES  (phi saturates at 1)
    M2 translocation k_tr up          rises        rises        UNCHANGED  <- unique
    M3 catalytic     k_pol up         UNCHANGED    rises        RISES      <- unique
    M4 combined      any of them      rises        rises        rises

Two numbers from one sample identify the class:
  * phi unchanged                  -> M3 (catalytic only)
  * phi changed, s/phi unchanged   -> M2 (translocation only)
  * both changed                   -> M1 or M4

and the closed-form estimator a = (s/phi)^(-1/6) - 1 recovers k_off/k_pol with no fitting.

WHAT THIS DOES AND DOES NOT SETTLE
----------------------------------
It settles the DESIGN question: the discriminating pair is (s, phi) measured in the same
samples with and without TPP1.  It does not settle the mechanism here, because the published
record does not contain a consistent (s, phi) pair in matched units - Latrick & Cech report phi,
and report processivity as a fraction of product above a length cut-off, which is not per-repeat
survival s.  Running the table on mismatched units would be exactly the error the repository's
preregistration forbids (holding a cutoff-dependent estimator to be a rate).

Run:  python scripts/solve/wp10_model_discriminator.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import sympy as sp

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"


def survival(k_pol, k_tr, k_off):
    return (k_pol / (k_pol + k_off)) ** 6 * (k_tr / (k_tr + k_off))


def phi_of(k_tr, k_off):
    return k_tr / (k_tr + k_off)


def main() -> int:
    # ---------------------------------------------------------------- exact identity
    kpol, ktr, koff = sp.symbols("k_pol k_tr k_off", positive=True)
    s = (kpol / (kpol + koff)) ** 6 * (ktr / (ktr + koff))
    phi = ktr / (ktr + koff)
    D = sp.simplify(sp.log(s) - sp.log(phi))

    print("=" * 90)
    print("WP10.1  The exact identity")
    print("=" * 90)
    print(f"  log(s) - log(phi) = {D}")
    print(f"  d/dlog(k_tr)  = {sp.simplify(sp.diff(D, ktr) * ktr)}   <- zero for EVERY parameter set")
    print(f"  d/dlog(k_off) = {sp.simplify(sp.diff(D, koff) * koff)}")
    print(f"  d/dlog(k_pol) = {sp.simplify(sp.diff(D, kpol) * kpol)}")
    assert sp.simplify(sp.diff(D, ktr)) == 0
    print()
    print("  A change in k_tr therefore cannot move s/phi.  Every other rate can.")
    print("  Corollary: if TPP1 changed only k_tr, s would rise by exactly the phi ratio.")
    print()

    # ---------------------------------------------------------------- discrimination table
    base = dict(k_pol=2.0, k_tr=0.07, k_off=0.012)
    s0, p0 = survival(**base), phi_of(base["k_tr"], base["k_off"])

    pert = {
        "M1 retention      (k_off 4x lower)": dict(base, k_off=base["k_off"] / 4),
        "M2 translocation  (k_tr 1.74x up)": dict(base, k_tr=base["k_tr"] * 1.74),
        "M3 catalytic      (k_pol 2.4x up)": dict(base, k_pol=base["k_pol"] * 2.4),
        "M4 combined       (both)": dict(base, k_off=base["k_off"] / 4, k_tr=base["k_tr"] * 1.74),
    }

    print("=" * 90)
    print("WP10.2  What each model predicts, from one (s, phi) pair")
    print("=" * 90)
    print(f"  baseline: k_pol={base['k_pol']}, k_tr={base['k_tr']}, k_off={base['k_off']}")
    print(f"            s = {s0:.5f}   phi = {p0:.5f}   s/phi = {s0/p0:.5f}")
    print()
    hdr = (f"  {'model':36s} {'phi ratio':>10} {'s ratio':>9} {'(s/phi) ratio':>14} "
           f"{'verdict':>22}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    rows = []
    for name, kw in pert.items():
        s1, p1 = survival(**kw), phi_of(kw["k_tr"], kw["k_off"])
        pr, sr, dr = p1 / p0, s1 / s0, (s1 / p1) / (s0 / p0)
        verdict = ("phi flat -> CATALYTIC" if abs(pr - 1) < 1e-3 else
                   "s/phi flat -> TRANSLOCATION" if abs(dr - 1) < 1e-3 else
                   "both move -> RETENTION/COMBINED")
        print(f"  {name:36s} {pr:>10.3f} {sr:>9.3f} {dr:>14.3f} {verdict:>22}")
        rows.append(dict(model=name.strip(), phi_ratio=round(pr, 4), s_ratio=round(sr, 4),
                         s_over_phi_ratio=round(dr, 4), verdict=verdict))
    print()
    print("  Each model owns a signature no other model reproduces:")
    print("    phi unchanged            -> catalytic only")
    print("    s/phi unchanged          -> translocation only")
    print("    both change              -> retention and/or combined")
    print()

    # ---------------------------------------------------------------- recovery estimator
    print("=" * 90)
    print("WP10.3  Recovering k_off/k_pol from one measured (s, phi) pair, with no fitting")
    print("=" * 90)
    print("  a = k_off/k_pol = (s/phi)^(-1/6) - 1")
    print()
    print("  *** CAVEAT ADDED AFTER AUDIT (wp11): this is an AGGREGATE, not a rate ratio. ***")
    print("  With per-step k_pol rates the identity becomes log(s/phi) = -sum_i log(1+k_off/k_pol_i)")
    print("  and the inversion above returns the geometric mean of the per-step factors.  Bias")
    print("  reaches ~2x when one of the six steps is 10x slower, and a further ~2x if the")
    print("  effective number of dissociation-competent steps is not 6.  Use the CHANNEL for")
    print("  classification; do not report the number as a rate.")
    print()
    print(f"  {'s/phi':>10} {'implied a':>12} {'comment':>44}")
    for r in (0.5, 0.8, 0.9, 0.95, 0.99, 1.0):
        a = r ** (-1 / 6) - 1
        print(f"  {r:>10.3f} {a:>12.5f}")
    print()
    print("  Note the direction: s/phi close to 1 forces a small k_off/k_pol, i.e. the enzyme")
    print("  loses the primer far more often during translocation than during extension, or")
    print("  polymerises much faster than it dissociates.  That ordering is itself measurable.")
    print()

    # ---------------------------------------------------------------- what is missing
    print("=" * 90)
    print("WP10.4  Why the table cannot be run on the published numbers yet")
    print("=" * 90)
    print("  The discriminating pair needs s and phi in matched units from the SAME samples.")
    print()
    print("  Published (Latrick & Cech 2010):")
    print("    phi                     0.57 -> 0.99      (single-round translocation assay)")
    print("    k_tr                    0.07 -> 0.10 /min (same assay)")
    print("    processivity            ~2.7x            (product-length distribution)")
    print()
    print("  Processivity reported as a fraction of product above a length cut-off is s^n for")
    print("  some window n, not s.  Deriving s from it requires knowing n, and the repository's")
    print("  own manuscript (section 3 of MANUSCRIPT.md) shows that estimator is window-dependent")
    print("  in a mixture - a 12-fold swing in implied k_off/k_pol from identical data.")
    print("  So the table is a DESIGN, not a result.  What it specifies is small:")
    print()
    print("    1. one reaction per condition, split two ways;")
    print("    2. full product-length distribution -> s by likelihood, not by a cut-off;")
    print("    3. single-round translocation on the same material -> phi;")
    print("    4. report the pair.  The model class then reads off from the table above.")
    print()
    print("  No new assay technology and no supercomputer.  Both measurements already exist in")
    print("  the literature separately; they have simply never been reported as a matched pair.")

    with open(OUT / "wp10_discriminator.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT/'wp10_discriminator.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
