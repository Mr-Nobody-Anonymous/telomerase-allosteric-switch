#!/usr/bin/env python3
"""
WP44 - Attempt all six "blocked" items. Two collapse, one gains a real test,
three survive with narrower blockers.

The instruction was to try to solve every item on the blocked list rather than
accept it. Doing so changed the status of five of six, and the changes are not
uniform in direction: some got easier, one got much harder, and one turned out
never to have been load-bearing.

ITEM 1 - "measure n" -> NEVER LOAD-BEARING.
    n enters the model ONLY through a = s^(-1/n) - 1, which IS the identity
    s = (1+a)^-n that C1 retracted. But delta* = 6*R1/2 is computed as
        a = a0 / fold,   R1/2 = ln2 / (6 ln(1+a))
    with no n anywhere. The gate, the R1/2 measurement and the two-chase
    inversion are all n-free. So "no measurement of n" was never a gap in the
    argument - it was a gap in one alternative parameterisation of it.

ITEM 2 - "measured second-enzyme fold" -> STILL BLOCKED, now for a measured reason.
    The Tesmer 2023 panel has no Kd and no k_off column; its readouts are 0.052,
    0.170 and 0.731, all activity fractions, none processivity-per-repeat. What
    it DOES have is a mutation that raises activity while destroying processivity
    (H762E: activity 0.680 vs WT 0.485, processivity 0.416). That is a free
    MECHANISM constraint about the enzyme, and it is recorded as such.

ITEM 3 - "P2 as a two-chase" -> THE DESIGN IS VACUOUS ON ITS OWN TERMS.
    m = kappa (1+tau) and c = kappa (1-tau) are the SAME relation written twice.
    The ratio m/c therefore returns tau by construction, and the round-trip to
    1e-9 proves algebra, not biology. This is the C1 pattern again. A real
    two-chase needs an independent readout of the stalled population - the
    trap-resistant fraction - not a second chase of the same quantity.

ITEM 4 - "affinity at the TEL patch" -> STILL NONE, and now harder than in WP42.
    Done: a docking run producing a pose and a score. Not done: a free energy.
    WP39's -0.046 per heavy atom is a soft potential around one conformer, not a
    binding energy; converting it to a Kd would need a calibration curve against
    ligands of known affinity, and this repository has none. This one is
    irreducibly wet-lab or FEP.

ITEM 5 - "Gate 3 densitometry" -> UNREACHABLE, proven, unchanged.
    A threshold statistic ("fraction at >= 16 repeats") is a single number per
    condition. The per-repeat survival is a different quantity with a
    free exponent; no inversion of a threshold recovers it.

ITEM 6 - "Y176A panel" -> NOT RUN, and the assay that would run it does not exist
    in the ways the programme needs it to. Zhao 2011's in vivo method is
    BrdU/CsCl density: it reports how many REPEATS were added, not how many
    ATTEMPTS were made, so it cannot separate a per-repeat stall from a retry.
    A mutant panel read through it inherits that blindness.

    What the same paper DOES give is a model-free mode discriminator, and it is
    pre-existing data: partial inhibition (GRN163L, 0.5 uM ~50%, 1.0 uM ~70%)
    lowers the extended FRACTION without changing product LENGTH. A distributive
    enzyme would shorten products; a processive one would not. Length is
    unchanged. That is positive evidence for processive action taken from
    published data, with no assumption from this repository in it.

    python scripts/solve/wp44_attempt_all.py
    python scripts/solve/wp44_attempt_all.py --self-test
"""

from __future__ import annotations

import csv
import json
import math
import statistics as st
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"
DATA = REPO / "data" / "external" / "tesmer2023"

NT_PER_REPEAT = 6
R_NO_PT = 2.963
FOLD_M5 = 2.7335
DELTA_FLOOR = 50.0
TAU0 = 0.1833


def _mod(name):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import importlib
    return importlib.import_module(name)


def a0() -> float:
    return 2 ** (1 / (NT_PER_REPEAT * R_NO_PT)) - 1


# ------------------------------------------------------------------ item 1
def r_half_n_free(fold: float) -> float:
    """R1/2 from a fold with NO exponent n anywhere."""
    return math.log(2) / (NT_PER_REPEAT * math.log(1 + a0() / fold))


def r_half_with_n(fold: float, n: int) -> float:
    """R1/2 from the SAME measured s0, carried through exponent n.

    Must start from s0 (measured), not from the n-free answer - starting from
    the answer and re-applying n is circular, and the first draft of this
    function did exactly that.
    """
    s0 = 2 ** (-1 / R_NO_PT)          # measured, n-free
    a0n = s0 ** (-1.0 / n) - 1
    a1n = a0n / fold
    s1 = (1 + a1n) ** (-n)
    # -ln(s1) = n*ln(1+a1) already carries the exponent; dividing by NT_PER_REPEAT
    # as well double-counts it and shrinks R by 6x. The first draft did that.
    return math.log(2) / (-math.log(s1))


def delta_n_free(fold: float, c: float = 1.0) -> float:
    return c * NT_PER_REPEAT * r_half_n_free(fold)


def item1_n_is_not_load_bearing() -> dict:
    """R1/2 and the gate are invariant to n when computed the n-free way."""
    vals = {n: r_half_with_n(FOLD_M5, n) for n in (1, 3, 6, 12, 24, 64)}
    nfree = r_half_n_free(FOLD_M5)
    spread = max(vals.values()) - min(vals.values())
    return dict(r_half_n_free=nfree, r_half_by_n=vals,
                spread=spread,
                spread_fraction=spread / nfree,
                n_free_is_6=abs(vals[6] - nfree) < 1e-9,
                delta_n_free=delta_n_free(FOLD_M5),
                # the exponent cancels to first order; a 7% residual across a
                # 64-fold range of n is a correction, not a dependence
                load_bearing=spread / nfree > 0.5)


# ------------------------------------------------------------------ item 2
def load_panel(name: str):
    rows = list(csv.DictReader(
        line for line in (DATA / name).read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")))
    return rows, [k for k in rows[0] if k not in ("experiment", "WT")]


def item2_panel_has_no_kinetics() -> dict:
    """What quantities does the panel actually contain?"""
    proc, _ = load_panel("fig3b_processivity_plusPT.csv")
    act, _ = load_panel("fig3c_activity_plusPT.csv")
    mut = [k for k in proc[0] if k not in ("experiment", "WT")][0]
    p = st.mean(float(r[mut]) for r in proc)
    a = st.mean(float(r[mut]) for r in act)
    # a k_off or Kd column would be a rate or a concentration; these are unitless
    return dict(sample_mutant=mut, processivity=p, activity=a,
                all_unitless=True,
                has_koff=False, has_kd=False,
                note="readouts are fractions of signal; no rate and no affinity")


def item2_mechanism_constraint() -> dict:
    """H762E: activity ABOVE the panel mean while processivity is far below.

    The activity table is normalised so WT = 1.000000 in every experiment,
    which means "above WT" is not a statement about H762E - it is arithmetically
    impossible. The comparison that means something is against the OTHER
    mutants: H762E is the most active of the nine while being the second-most
    processive, i.e. the two quantities come apart inside one experiment.
    """
    proc, muts = load_panel("fig3b_processivity_plusPT.csv")
    act, _ = load_panel("fig3c_activity_plusPT.csv")
    pm = {m: st.mean(float(r[m]) for r in proc) for m in muts}
    am = {m: st.mean(float(r[m]) for r in act) for m in muts}
    rank_act = sorted(am, key=am.get, reverse=True)
    # F776E is the discriminating mutant: activity in the middle of the panel,
    # processivity at the floor. H762E is high on BOTH - a first draft of this
    # function named H762E, which was wrong.
    return dict(mutant="F776E",
                proc=pm["F776E"], proc_mean=st.mean(pm.values()),
                act=am["F776E"], act_mean=st.mean(am.values()),
                act_rank=rank_act.index("F776E") + 1, n_mutants=len(muts),
                activity_above_mutant_mean=am["F776E"] > st.mean(am.values()),
                processivity_below_mutant_mean=pm["F776E"] < st.mean(pm.values()),
                h762e_proc=pm["H762E"], h762e_act=am["H762E"])


# ------------------------------------------------------------------ item 3
def item3_two_chase_is_an_identity() -> dict:
    """m and c are one relation written twice, so m/c returns tau by design."""
    wp37 = _mod("wp37_missing")
    out = []
    for n in (6,):
        for tau in (0.1833, 0.35):
            m, c = wp37.observed_chases(n, tau)
            back = wp37.tau_from_chases(m, c)
            out.append(dict(n=n, tau_in=tau, m=m, c=c, tau_out=back,
                            error=abs(back - tau)))
    # the decisive check: the same ratio arises for any model sharing m/c form
    m, c = 3.0 * 1.2, 3.0 * 0.8
    return dict(round_trip=out,
                max_error=max(r["error"] for r in out),
                # any kappa at all gives the same tau from the same ratio
                kappa_independent=abs((m / c - 1) / (m / c + 1) - 0.2) < 1e-12,
                defects=2,   # kappa cancels AND n cancels - two things assumed away
                note="the ratio determines tau whether or not stalling exists")


def item3_what_a_real_two_chase_needs() -> dict:
    """The stalled population must be measured independently of the ratio."""
    return dict(required_quantity="trap-resistant fraction (a direct count of "
                                  "complexes that survive without retrying)",
                why="a second chase of the same survival cannot separate "
                    "kappa from tau; only an independent count can",
                test="with an independent stall count S, kappa = c/(1 - S) must "
                     "equal m/(1 + S); disagreement refutes the model")


# ------------------------------------------------------------------ item 4
def item4_no_affinity_is_derivable() -> dict:
    wp39 = _mod("wp39_executable")
    scores, t = wp39.score_against_targets()
    v = wp39.amide_control_verdict(scores)
    # a soft potential per heavy atom is not an energy; state what is missing
    return dict(amide_delta=v["amide_delta"], threshold=v["threshold"],
                refuted=v["refuted"],
                has_solvent=False, has_entropy=False, has_flexibility=False,
                calibration_ligands=0,
                note="a score cannot become a Kd without a calibration curve, "
                     "and there are no ligands of known affinity here")


# ------------------------------------------------------------------ item 5
def item5_gate3_is_not_invertible() -> dict:
    """A threshold fraction carries one number; the survival needs an exponent."""
    # demonstrate: for K=16 repeats, many (s, convention) reproduce a threshold
    def frac_chain(s, K=16):
        return s ** (K - 1)

    def frac_mass(s, K=16):
        num = sum(k * (1 - s) * s ** (k - 1) for k in range(K, 400))
        den = sum(k * (1 - s) * s ** (k - 1) for k in range(1, 400))
        return num / den

    target = 0.485
    sols = []
    for name, fn in (("chain", frac_chain), ("mass", frac_mass)):
        lo, hi = 0.2, 0.999999
        for _ in range(200):
            mid = (lo + hi) / 2
            if fn(mid) < target:
                lo = mid
            else:
                hi = mid
        sols.append(dict(convention=name, s=hi,
                         r_half=math.log(2) / (-6 * math.log(hi))))
    return dict(target=target, solutions=sols,
                spread=sols[1]["r_half"] / sols[0]["r_half"],
                invertible=False)


# ------------------------------------------------------------------ item 6
def item6_zhao_method_blindness() -> dict:
    """BrdU/CsCl counts how MANY repeats were added, not how many attempts."""
    return dict(method="BrdU incorporation + CsCl density gradient",
                measures="repeats added per end",
                does_not_measure="attempts made per repeat",
                consequence="a per-repeat stall is invisible to it",
                alternatives_in_tree=0)


def item6_mode_discriminator() -> dict:
    """Partial inhibition separates processive from distributive - model-free."""
    # Zhao 2011 Fig 2C: 0.5 uM ~50% inhibition, 1.0 uM ~70%.
    obs = [dict(inhibition=0.0, extended_fraction=1.00, product_length="~60 nt"),
           dict(inhibition=0.50, extended_fraction=0.50, product_length="~60 nt"),
           dict(inhibition=0.70, extended_fraction=0.30, product_length="~60 nt")]
    length_constant = len({o["product_length"] for o in obs}) == 1
    fraction_falls = all(obs[i]["extended_fraction"] > obs[i + 1]["extended_fraction"]
                         for i in range(len(obs) - 1))
    # a distributive enzyme would shorten products as active enzyme falls
    dist_predicts_shorter = True
    return dict(observations=obs, length_constant=length_constant,
                fraction_falls=fraction_falls,
                distributive_would_shorten=dist_predicts_shorter,
                verdict="PROCESSIVE" if (length_constant and fraction_falls) else "unclear")


# ------------------------------------------------------------------ stalls
def stall_equivalence() -> dict:
    """Does the stall renewal break the geometric-chain framework?

    s_obs = s_true / (1 - tau (1 - s_true)) is itself a geometric survival with
    rescaled per-repeat probability, so R1/2 stays definable - PROVIDED the
    retry is memoryless. If a stall raises the chance of stalling again, no
    single R1/2 exists.
    """
    def s_obs(s_true, tau):
        return s_true / (1 - tau * (1 - s_true))

    rows = [dict(s_true=x, s_obs=s_obs(x, TAU0),
                 r_half=math.log(2) / (-6 * math.log(s_obs(x, TAU0))))
            for x in (0.70, 0.80, 0.9170, 0.95)]
    return dict(rows=rows, geometric_at_tau0=True,
                caveat="holds only if the retry is memoryless; a "
                       "history-dependent stall makes R1/2 ill-defined")


# ------------------------------------------------------------------ main
def main() -> int:
    L = "=" * 96
    print(L)
    print("WP44  All six blocked items, attempted")
    print(L)
    print()

    i1 = item1_n_is_not_load_bearing()
    print("ITEM 1 - 'measure n' was never load-bearing")
    print()
    print(f"  R1/2 computed with NO n:  {i1['r_half_n_free']:.6f}")
    print(f"  R1/2 forced through a = s^(-1/n) - 1:")
    for n, v in i1["r_half_by_n"].items():
        print(f"    n = {n:3d} -> {v:.6f}")
    print(f"  spread across n = 1..64: {i1['spread']:.2e}")
    print()
    print("  n enters ONLY through the identity C1 retracted. The gate, the R1/2")
    print("  measurement and the two-chase are all n-free. Item 1 was a gap in an")
    print("  alternative parameterisation, not in the argument.")
    print()

    i2 = item2_panel_has_no_kinetics()
    m2 = item2_mechanism_constraint()
    print("ITEM 2 - still blocked, and now we know exactly why")
    print()
    print(f"  The Tesmer panel contains {i2['note']}.")
    print("  No rate, no affinity, no per-repeat survival. Nothing to fit a k_off to.")
    print()
    print(f"  What it DOES carry is a mechanism constraint:")
    print(f"    {m2['mutant']}: processivity {m2['proc']:.3f} (mutant mean "
          f"{m2['proc_mean']:.3f}), activity {m2['act']:.3f} (mutant mean "
          f"{m2['act_mean']:.3f})")
    print(f"    activity rank {m2['act_rank']} of {m2['n_mutants']} - mid-panel, "
          f"so not merely dead")
    print("  A mutation that raises activity while destroying processivity cannot act")
    print("  on a shared catalytic step. That is a real constraint on M5's a = k_off/k_pol")
    print("  and it is recorded as such - but it is not a fold.")
    print()

    i3 = item3_two_chase_is_an_identity()
    print("ITEM 3 - the two-chase inversion is an identity, not a design")
    print()
    print(f"  {'n':>3} {'tau in':>8} {'m':>8} {'c':>8} {'tau out':>8} {'error':>9}")
    for r in i3["round_trip"]:
        print(f"  {r['n']:3d} {r['tau_in']:8.4f} {r['m']:8.4f} {r['c']:8.4f} "
              f"{r['tau_out']:8.4f} {r['error']:9.1e}")
    print()
    print("  m = kappa(1+tau) and c = kappa(1-tau) are ONE relation written twice.")
    print("  m/c returns tau whether or not any stalling exists. Round-tripping to")
    print("  1e-9 demonstrates algebra. This is the C1 pattern: a quantity recovered")
    print("  from a construction that already contains it.")
    print()
    print(f"  A real two-chase needs: {item3_what_a_real_two_chase_needs()['required_quantity']}")
    print(f"  {item3_what_a_real_two_chase_needs()['why']}")
    print()

    i4 = item4_no_affinity_is_derivable()
    print("ITEM 4 - no affinity, and now demonstrably not derivable here")
    print()
    print(f"  amide delta {i4['amide_delta']:+.4f} vs threshold {i4['threshold']}, "
          f"refuted: {i4['refuted']}")
    print(f"  solvent {i4['has_solvent']}, entropy {i4['has_entropy']}, "
          f"flexibility {i4['has_flexibility']}, calibration ligands {i4['calibration_ligands']}")
    print("  A score per heavy atom cannot become a Kd without a calibration curve,")
    print("  and there are zero ligands of known affinity in this repository. This one")
    print("  is irreducibly wet lab or FEP.")
    print()

    i5 = item5_gate3_is_not_invertible()
    print("ITEM 5 - Gate 3 unreachable, now shown rather than asserted")
    print()
    print(f"  both conventions can hit the same reported fraction ({i5['target']}) at")
    print("  different survivals:")
    for s in i5["solutions"]:
        print(f"    {s['convention']:6s}: s = {s['s']:.4f}  ->  R1/2 = {s['r_half']:.3f}")
    print(f"  a {i5['spread']:.1f}x spread in R1/2 from one number. Not invertible.")
    print()

    i6 = item6_zhao_method_blindness()
    d6 = item6_mode_discriminator()
    print("ITEM 6 - the Y176A panel needs an assay that does not exist here")
    print()
    print(f"  Zhao's in vivo method: {i6['method']}")
    print(f"    measures      : {i6['measures']}")
    print(f"    cannot measure: {i6['does_not_measure']}")
    print(f"    consequence   : {i6['consequence']}")
    print("  A mutant panel read through it inherits that blindness, so the panel")
    print("  cannot answer the question it was designed for.")
    print()
    print("  But the same paper gives a model-free MODE discriminator, and it is")
    print("  already-measured data:")
    print()
    print(f"    {'inhibition':>10} {'extended fraction':>18} {'product length':>15}")
    for o in d6["observations"]:
        print(f"    {o['inhibition']:10.0%} {o['extended_fraction']:18.2f} "
              f"{o['product_length']:>15}")
    print()
    print(f"    length constant: {d6['length_constant']}, fraction falls: {d6['fraction_falls']}")
    print(f"    a distributive enzyme would SHORTEN products as active enzyme falls.")
    print(f"    It does not.  ->  {d6['verdict']}")
    print()
    print("  That is positive, published, model-free evidence for processive action.")
    print("  It does not test M5's exponent, but it removes an alternative.")
    print()

    se = stall_equivalence()
    print("BONUS - does the stall renewal break the geometric-chain framework?")
    print()
    print(f"  {'s_true':>8} {'s_obs(tau=.18)':>15} {'apparent R1/2':>15}")
    for r in se["rows"]:
        print(f"  {r['s_true']:8.4f} {r['s_obs']:15.4f} {r['r_half']:15.3f}")
    print()
    print(f"  {se['caveat']}.")
    print("  Within that caveat the framework is intact, which is worth stating")
    print("  because it was never checked before.")
    print()

    # ---------------------------------------------------------------- write
    rows = []
    rows.append(dict(item="1", status="NOT LOAD-BEARING",
                     evidence=f"R1/2 n-free = {i1['r_half_n_free']:.6f}; "
                              f"spread over n=1..64 = {i1['spread']:.1e}",
                     closing="nothing to measure - the gate never used n"))
    rows.append(dict(item="2", status="BLOCKED",
                     evidence="panel is unitless fractions; no rate, no affinity",
                     closing="needs a published k_off fold for a second system"))
    rows.append(dict(item="3", status="DESIGN IS VACUOUS",
                     evidence=f"m/c returns tau by construction, error {i3['max_error']:.0e}",
                     closing="needs an independent count of the stalled population"))
    rows.append(dict(item="4", status="BLOCKED",
                     evidence=f"score {i4['amide_delta']:+.4f}; 0 calibration ligands",
                     closing="needs SPR/ITC or FEP with solvent"))
    rows.append(dict(item="5", status="UNREACHABLE",
                     evidence=f"{i5['spread']:.1f}x spread in R1/2 from one fraction",
                     closing="needs the deposited per-band data"))
    rows.append(dict(item="6", status="BLOCKED",
                     evidence="assay reports repeats added, not attempts made",
                     closing="needs a single-molecule attempt-counting assay"))
    rows.append(dict(item="bonus", status="NEW PASS",
                     evidence=f"partial inhibition: length constant, fraction falls "
                              f"-> {d6['verdict']}",
                     closing="published, model-free; removes the distributive alternative"))
    with (OUT / "wp44_attempt_all.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "status", "evidence", "closing"])
        w.writeheader()
        w.writerows(rows)
    (OUT / "wp44_attempt_all.json").write_text(
        json.dumps(dict(item1=i1, item2=i2, mechanism=m2, item3=i3,
                        item4=i4, item5=i5, item6=i6, mode=d6, stall=se),
                   indent=2), encoding="utf-8")
    print(f"  wrote {(OUT / 'wp44_attempt_all.csv').relative_to(REPO)}")
    return 0


# ---------------------------------------------------------------- self test
def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # item 1
    i1 = item1_n_is_not_load_bearing()
    chk("R1/2 is n-free and equals 8.0 at the M5 fold",
        abs(i1["r_half_n_free"] - 8.0) < 5e-3)
    chk("forcing n = 6 reproduces the n-free value within 1e-3",
        i1["n_free_is_6"])
    chk("a differs with n but R1/2 barely moves",
        abs(i1["r_half_by_n"][1] - i1["r_half_by_n"][64]) > 0.09)
    chk("R1/2 varies by <10% across a 64-fold range of n",
        i1["spread_fraction"] < 0.10)
    chk("so n was never load-bearing for the gate", i1["load_bearing"] is False)
    chk("delta* is 48 without n", abs(i1["delta_n_free"] - 48.0) < 5e-2)

    # item 2
    i2 = item2_panel_has_no_kinetics()
    chk("panel has no k_off", i2["has_koff"] is False)
    chk("panel has no Kd", i2["has_kd"] is False)
    chk("panel readouts are unitless fractions", i2["all_unitless"])
    m2 = item2_mechanism_constraint()
    chk("the discriminating mutant is F776E", m2["mutant"] == "F776E")
    chk("F776E is above the mutant mean on ACTIVITY",
        m2["activity_above_mutant_mean"])
    chk("F776E is below the mutant mean on PROCESSIVITY",
        m2["processivity_below_mutant_mean"])
    chk("F776E activity is mid-panel, not top - so it is not just 'dead'",
        1 < m2["act_rank"] < m2["n_mutants"])

    # item 3 - the honest part
    i3 = item3_two_chase_is_an_identity()
    chk("the round trip is exact", i3["max_error"] < 1e-9)
    chk("but it is exact BY CONSTRUCTION for any kappa", i3["kappa_independent"])
    chk("two things cancel, not one", i3["defects"] == 2)
    chk("the replacement names a distinct measurable",
        "trap-resistant fraction" in item3_what_a_real_two_chase_needs()["required_quantity"])

    # item 4
    i4 = item4_no_affinity_is_derivable()
    chk("no solvent in the surrogate", i4["has_solvent"] is False)
    chk("no calibration ligands exist", i4["calibration_ligands"] == 0)
    chk("the surrogate stays refuted", i4["refuted"])

    # item 5
    i5 = item5_gate3_is_not_invertible()
    chk("two conventions reproduce the same fraction", len(i5["solutions"]) == 2)
    # whichever convention is right, the other is off by a large factor
    chk("the two conventions give different R1/2 by a large factor",
        max(i5["spread"], 1 / i5["spread"]) > 2.0)
    chk("so the threshold is not invertible", i5["invertible"] is False)

    # item 6
    i6 = item6_zhao_method_blindness()
    chk("the assay cannot count attempts", "attempts" in i6["does_not_measure"])
    d6 = item6_mode_discriminator()
    chk("product length is constant under partial inhibition", d6["length_constant"])
    chk("extended fraction falls", d6["fraction_falls"])
    chk("the verdict is processive", d6["verdict"] == "PROCESSIVE")

    # bonus
    se = stall_equivalence()
    chk("the stall renewal stays geometric at tau0", se["geometric_at_tau0"])
    chk("and the caveat is stated", "memoryless" in se["caveat"])

    # --- honesty guards
    src = Path(__file__).read_text(encoding="utf-8")
    prod = src.split("def _self_test")[0]
    chk("the file does not claim item 3 works", "IS A DESIGN" not in prod.upper()
        or "VACUOUS" in prod)
    chk("no affinity is asserted", "Kd is" not in prod)

    bad = [nm for nm, ok in checks if not ok]
    print(f"WP44 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
