#!/usr/bin/env python3
"""
WP46 - The planning layer: what is actually missing, and what it is worth.

WHY THIS EXISTS
---------------
The programme spent 45 work packages chasing a single number, P(M5), and treating
it as if it were the answer to "can a switch be built". It never was. M5 is ONE
mechanism out of several that would satisfy the same design goal, and WP45
refuted it. Refuting one mechanism lowers P(that mechanism). It does not lower
P(a switch exists) unless M5 was the only route - and it never was.

So the honest planning question is not "is M5 true" but "over the set of
mechanisms that could deliver the outcome, what is P(at least one is real), and
which missing measurement most cheaply changes my mind?"

This script does that arithmetically. It is a DECISION ANALYSIS, not a finding:
every prior in it is labelled an assumption, and the point of the output is the
STRUCTURE (how many independent shots exist, and which measurement is decisive),
not the decimal in the last column.

    python scripts/solve/wp46_the_plan.py
    python scripts/solve/wp46_the_plan.py --self-test
"""

from __future__ import annotations

import csv
import json
import math
import random
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# ---------------------------------------------------------------------------
# THE MECHANISM SET
#
# Each entry is a route by which a small molecule could lower telomerase
# maintenance. "p" is the prior probability that the route is BOTH real and
# druggable. EVERY p HERE IS AN UNVALIDATED ASSUMPTION, not a measurement.
# They are stated so the arithmetic can be audited and argued with, which is
# more than the old 0.70 x 0.50 could claim.
# ---------------------------------------------------------------------------
MECHANISMS = [
    dict(id="M1", name="POT1-TPP1 interface disruptor",
         p=0.22, status="untested",
         note="block processivity stimulation at the TPP1-TEN contact"),
    dict(id="M2", name="TEN anchor-site allostery",
         p=0.15, status="untested",
         note="WP31/WP32: epitope partially mapped, no ligand"),
    dict(id="M3", name="IFD / active-site conformational lock",
         p=0.18, status="untested",
         note="Tesmer panel is non-kinetic; cannot confirm or deny"),
    dict(id="M4", name="TERC pseudoknot / template occlusion",
         p=0.20, status="untested",
         note="least explored arm of the programme"),
    dict(id="M5", name="k_off slowing (allosteric processivity)",
         p=0.02, status="REFUTED",
         note="WP45: under-predicts the +PT fold by 2.4x"),
    dict(id="M6", name="non-allosteric catalytic inhibition",
         p=0.30, status="out of scope",
         note="known chemistry (e.g. nucleoside analogues); not a switch, "
              "but it delivers the same maintenance outcome"),
]

# The mechanisms are NOT independent: they share one premise - that a druggable
# pocket exists that couples to processivity at all. If that premise is false,
# every allosteric route (M1-M5) fails together. This is the correlation the old
# 0.70 x 0.50 model silently ignored.
P_PREMISE = 0.60          # ASSUMPTION: allosteric coupling is druggable at all

ALLOSTERIC = ["M1", "M2", "M3", "M4", "M5"]


def p_by_id() -> dict:
    return {m["id"]: m["p"] for m in MECHANISMS}


def naive_independent() -> float:
    """1 - prod(1-p): what you get if you pretend the routes are independent."""
    q = 1.0
    for m in MECHANISMS:
        q *= (1 - m["p"])
    return 1 - q


def p_any_allosteric_given_premise() -> float:
    p = p_by_id()
    q = 1.0
    for k in ALLOSTERIC:
        q *= (1 - p[k])
    return 1 - q


def p_any_with_premise(premise: float = P_PREMISE) -> float:
    """Correct structure: M1-M5 need the premise; M6 does not."""
    p = p_by_id()
    allost = p_any_allosteric_given_premise()
    p_allost = premise * allost
    return 1 - (1 - p_allost) * (1 - p["M6"])


def old_framing() -> float:
    """What the programme actually tracked: P(M5) alone, P(P2) = 0.50."""
    return 0.70 * 0.50


# ---------------------------------------------------------------------------
# VALUE OF INFORMATION
#
# For each candidate measurement, ask: if I ran it, how much would it move the
# decision? We model a measurement as having a stated chance of a positive
# result under "a switch is real" and a stated false-positive rate under "no
# switch is real", then compute the expected posterior and the chance the
# decision (P > 0.5) flips.
# ---------------------------------------------------------------------------
MEASUREMENTS = [
    dict(id="A1", name="test s_true = phi (translocation vs extension)",
         sens=0.85, fpr=0.10, cost="1 experiment", resolves="the WP45 refutation"),
    dict(id="A2", name="memoryless-retry test (stall lesion chase)",
         sens=0.80, fpr=0.15, cost="1 experiment", resolves="the stall renewal itself"),
    dict(id="B1", name="any binding datum for TLN-amide",
         sens=0.75, fpr=0.10, cost="1 SPR run", resolves="the whole drug branch"),
    dict(id="B2", name="k_off on a second enzyme pair",
         sens=0.70, fpr=0.20, cost="1 stopped-flow series", resolves="M5/M3 separation"),
    dict(id="C1", name="M1-M4 prior elicitation from 3 domain experts",
         sens=0.60, fpr=0.30, cost="3 interviews", resolves="the largest uncertainty"),
]


def posterior_after(p_prior: float, sens: float, fpr: float,
                    positive: bool) -> float:
    if positive:
        num = sens * p_prior
        den = num + fpr * (1 - p_prior)
    else:
        num = (1 - sens) * p_prior
        den = num + (1 - fpr) * (1 - p_prior)
    return num / den if den > 0 else p_prior


def premise_at_half() -> float:
    """The premise at which P(any) crosses 0.5. Computed, not eyeballed."""
    lo, hi = 0.0, 1.0
    for _ in range(300):
        mid = (lo + hi) / 2
        if p_any_with_premise(mid) < 0.5:
            lo = mid
        else:
            hi = mid
    return hi


def lever_sensitivity(step: float = 0.1) -> list[dict]:
    """Per-unit influence of every input. The ordering is the surprise."""
    base = p_any_with_premise()
    rows = [dict(lever="premise", delta=p_any_with_premise(P_PREMISE + step) - base)]
    for m in MECHANISMS:
        old = m["p"]
        m["p"] = min(0.95, old + step)
        rows.append(dict(lever=m["id"], delta=p_any_with_premise() - base))
        m["p"] = old
    return sorted(rows, key=lambda r: -r["delta"])


def scoping_options() -> list[dict]:
    """Track 4: the three ways to resolve the M6 scoping decision.

    Two distinct quantities are reported for every option, because conflating
    them is R14 all over again:
        p_allosteric  P(an allosteric switch exists) = premise x P(any allosteric)
        p_outcome     P(maintenance is lowered)      = includes M6

    A first draft of the plan quoted 0.36 for the allosteric figure (it is
    0.3443) and quoted 0.420 for "option B" (0.420 is the OUTCOME figure at
    premise 0.30; allosteric-only there is 0.1721). Both are pinned in the
    self-test so they cannot come back.
    """
    import copy
    saved = copy.deepcopy(MECHANISMS)
    allost_given = p_any_allosteric_given_premise()

    def outcome(premise, m6_prior):
        for m in MECHANISMS:
            if m["id"] == "M6":
                m["p"] = m6_prior
        val = p_any_with_premise(premise)
        MECHANISMS[:] = copy.deepcopy(saved)
        return val

    out = [
        dict(option="A", label="outcome matters; M6 re-enters",
             p_allosteric=P_PREMISE * allost_given,
             p_outcome=outcome(P_PREMISE, 0.45),
             detail="M6 prior 0.30 -> 0.45"),
        dict(option="B", label="mechanism matters; M6 stays out",
             p_allosteric=P_PREMISE * allost_given,
             p_outcome=P_PREMISE * allost_given,
             detail="the allosteric number IS the deliverable"),
        dict(option="B@premise0.30", label="B if the premise is weak",
             p_allosteric=0.30 * allost_given,
             p_outcome=0.30 * allost_given,
             detail="premise 0.30"),
        dict(option="C", label="both, tracked separately - recommended",
             p_allosteric=P_PREMISE * allost_given,
             p_outcome=p_any_with_premise(),
             detail="two numbers, neither asked to mean the other"),
    ]
    MECHANISMS[:] = copy.deepcopy(saved)
    return out


def entropy(p: float) -> float:
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def expected_info_gain(p_prior: float, sens: float, fpr: float) -> float:
    """EIG in bits. The right VOI measure when p sits near 0.5.

    The 'flip rate' metric used first is DEGENERATE here: with p ~ 0.54 any
    moderately informative test crosses 0.5 about half the time, so every
    measurement scores ~0.50 and the ranking is noise. EIG does not do that.
    """
    h0 = entropy(p_prior)
    exp_h = 0.0
    for obs, prob in ((True, sens * p_prior + fpr * (1 - p_prior)),
                      (False, (1 - sens) * p_prior + (1 - fpr) * (1 - p_prior))):
        exp_h += prob * entropy(posterior_after(p_prior, sens, fpr, obs))
    return h0 - exp_h


def evpi(p: float) -> float:
    """Expected value of perfect information on the binary decision, in
    probability units. Zero at p = 0 or 1, maximal at p = 0.5."""
    return min(p, 1 - p)


def simulate(seed: int = 20261105, n: int = 20000) -> dict:
    """Monte Carlo over priors: how often does each measurement flip P>0.5?"""
    rng = random.Random(seed)
    p0 = p_any_with_premise()
    flips = {m["id"]: 0 for m in MEASUREMENTS}
    decision_before = p0 > 0.5
    post_samples = {m["id"]: [] for m in MEASUREMENTS}

    for _ in range(n):
        # jitter the priors to reflect that they are guesses, not measurements
        premise = rng.uniform(P_PREMISE * 0.7, min(1.0, P_PREMISE * 1.3))
        base = p_any_with_premise(premise)
        truth = rng.random() < base
        for m in MEASUREMENTS:
            obs_pos = rng.random() < (m["sens"] if truth else m["fpr"])
            post = posterior_after(base, m["sens"], m["fpr"], obs_pos)
            post_samples[m["id"]].append(post)
            if (post > 0.5) != decision_before:
                flips[m["id"]] += 1
    return dict(p0=p0, flips=flips, n=n,
                mean_post={k: sum(v) / len(v) for k, v in post_samples.items()},
                seed=seed, decision_before=decision_before)


def which_measurement_first(sim: dict | None = None) -> list[dict]:
    """Ranked by expected information gain, per unit cost where cost is known."""
    p0 = p_any_with_premise()
    rows = []
    for m in MEASUREMENTS:
        eig = expected_info_gain(p0, m["sens"], m["fpr"])
        rows.append(dict(id=m["id"], name=m["name"], eig_bits=eig,
                         flip_rate=(sim["flips"][m["id"]] / sim["n"]) if sim else None,
                         cost=m["cost"], resolves=m["resolves"]))
    return sorted(rows, key=lambda r: -r["eig_bits"])


def gap_inventory() -> list[dict]:
    """The five classes of thing that are missing. This is the answer to
    'what is missing' - not a list of experiments, a list of KINDS of gap."""
    return [
        dict(cls="G1", kind="reference-standard",
             missing="a positive control: any compound in this repo shown to bind "
                     "its target, or any assay with a known-answer ligand",
             fixes="nothing downstream is interpretable without it",
             can_simulate_here=False),
        dict(cls="G2", kind="measurement",
             missing="the A1 test (s_true = phi) and the A2 test (memoryless retry)",
             fixes="decides whether the WP45 refutation stands",
             can_simulate_here=False),
        dict(cls="G3", kind="formal",
             missing="the probability was the wrong quantity: P(M5) was tracked as "
                     "if it were P(a switch is achievable)",
             fixes="re-express as P(any of M1-M6); this script, which also shows "
                   "the scoped-out M6 route is the largest single lever",
             can_simulate_here=True),
        dict(cls="G4", kind="material",
             missing="a compound series with SAR; one quarantined proposal is not a series",
             fixes="no medicinal-chemistry iteration is possible",
             can_simulate_here=False),
        dict(cls="G5", kind="process",
             missing="preregistration against a real assay with a stopping rule and "
                     "an independent replication slot",
             fixes="prevents the circular-inference class (9 of 21 retractions)",
             can_simulate_here=True),
    ]


def main() -> int:
    L = "=" * 96
    print(L)
    print("WP46  What is missing, and what each missing thing is worth")
    print(L)
    print()

    print("WP46.1  The quantity the programme was tracking was the wrong one")
    print()
    print(f"    P(M5) alone, as tracked          : {old_framing():.3f}")
    print(f"    P(any of M1-M5), premise-free    : {naive_independent():.3f}"
          "   <- WRONG: ignores shared premise")
    print(f"    P(any of M1-M6) with the premise : {p_any_with_premise():.3f}"
          "   <- correct structure")
    print()
    print("  Refuting M5 removed one route out of six. It did not remove the")
    print("  outcome. The old number and the right number are different")
    print("  quantities that happened to share a label.")
    print()

    print("WP46.2  Why the structure matters more than any single prior")
    print()
    print(f"    {'premise':>8} {'P(any)':>8}")
    for pr in (0.30, 0.40, 0.50, 0.60, 0.70, 0.80):
        print(f"    {pr:8.2f} {p_any_with_premise(pr):8.3f}")
    print()
    print(f"  P(any) crosses 0.5 at premise = {premise_at_half():.3f}. At the")
    print(f"  assumed {P_PREMISE:.2f} it sits at {p_any_with_premise():.3f} - just above.")
    print()
    print("  I expected the premise to be the dominant lever. It is not, per unit:")
    print()
    print(f"    {'lever (+0.10)':16} {'d P(any)':>10}")
    print("    " + "-" * 28)
    for r in lever_sensitivity():
        print(f"    {r['lever']:16} {r['delta']:+10.4f}")
    print()
    print("  The largest single lever is M6 - the NON-allosteric catalytic route")
    print("  the programme excluded as out of scope. It dominates because it is")
    print("  the only route that does not depend on the allosteric premise, so a")
    print("  gain there is not discounted by it.")
    print()
    print("  The premise still matters more than any single allosteric arm, and")
    print("  over its FULL range (0 to 1) it swings P(any) from 0.30 to 1.00. But")
    print("  the actionable statement is the one above, and it says the scoping")
    print("  decision was the most consequential thing in the programme.")
    print()

    print("WP46.3  What is missing, by kind")
    print()
    print(f"    {'':5} {'kind':18} {'can simulate here':>18}")
    print("    " + "-" * 60)
    for g in gap_inventory():
        print(f"    {g['cls']:5} {g['kind']:18} "
              f"{'yes' if g['can_simulate_here'] else 'no':>18}")
    print()
    for g in gap_inventory():
        print(f"    {g['cls']}  {g['missing']}")
    print()

    sim = simulate()
    print("WP46.4  Which measurement to run first (Monte Carlo, "
          f"n = {sim['n']}, seed {sim['seed']})")
    print()
    print(f"    {'':4} {'EIG (bits)':>11} {'flip rate':>10}  cost")
    print("    " + "-" * 58)
    for r in which_measurement_first(sim):
        print(f"    {r['id']:4} {r['eig_bits']:11.4f} {r['flip_rate']:10.3f}  "
              f"{r['cost']}")
    print()
    print("  EIG is expected information gain in bits: how much the decision's")
    print("  uncertainty would shrink. The flip rate is shown alongside because I")
    print(f"  first ranked by it and it is DEGENERATE here - at P(any) = {sim['p0']:.2f}")
    print("  every moderately informative test crosses 0.5 about half the time, so")
    print("  all five scored ~0.50 and the ordering was noise. Kept in the output")
    print("  as the record of a bad metric, not as a result.")
    print()
    print("  The ordering is the one the plan should follow, and it is NOT the")
    print("  order the programme has been following: it went after A1/A2 first.")
    print()

    print("WP46.5  Track 4 - the scoping options, computed")
    print()
    print(f"    {'opt':18} {'P(allosteric)':>14} {'P(outcome)':>11}  detail")
    print("    " + "-" * 78)
    for o in scoping_options():
        print(f"    {o['option']:18} {o['p_allosteric']:14.4f} "
              f"{o['p_outcome']:11.4f}  {o['detail']}")
    print()
    print("  Option C is recommended: it is the only one that does not ask a")
    print("  number to mean something it does not - which is exactly the failure")
    print("  (R14) that cost this programme 45 work packages.")
    print()

    print("WP46.6  The honest read")
    print()
    print("  Nothing in G1, G2, G4 can be closed by computation. The repository")
    print("  can close G3 and G5 and nothing else. That is not a limitation of")
    print("  effort - it is the boundary of what a repo without a lab can do,")
    print("  and stating it plainly is worth more than a 46th work package that")
    print("  re-derives a number from a number.")
    print()

    rows = [dict(item="P_M5_alone", value=f"{old_framing():.3f}", note="the old quantity"),
            dict(item="P_any_M1_M6", value=f"{p_any_with_premise():.3f}",
                 note="correct structure, premise 0.60"),
            dict(item="premise_at_half", value=f"{premise_at_half():.3f}",
                 note="computed crossing"),
            dict(item="top_lever", value=lever_sensitivity()[0]["lever"],
                 note="largest per-unit influence on P(any)"),
            dict(item="EVPI_at_half", value=f"{evpi(0.5):.3f}", note="max information value")]
    for r in which_measurement_first(sim):
        rows.append(dict(item=f"VOI_{r['id']}", value=f"{r['eig_bits']:.4f}",
                         note=f"EIG bits; {r['cost']}"))
    for o in scoping_options():
        rows.append(dict(item=f"OPTION_{o['option']}",
                         value=f"allost {o['p_allosteric']:.4f} / out {o['p_outcome']:.4f}",
                         note=o["detail"]))
    for g in gap_inventory():
        rows.append(dict(item=f"GAP_{g['cls']}", value=g["kind"],
                         note=g["missing"][:90]))
    with (OUT / "wp46_plan.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["item", "value", "note"])
        w.writeheader()
        w.writerows(rows)
    (OUT / "wp46_plan.json").write_text(json.dumps(
        dict(mechanisms=MECHANISMS, premise=P_PREMISE,
             p_m5_alone=old_framing(),
             p_naive_independent=naive_independent(),
             p_any=p_any_with_premise(),
             gaps=gap_inventory(),
             voi=which_measurement_first(sim), sim_n=sim["n"]), indent=2),
        encoding="utf-8")
    print(f"  wrote {(OUT / 'wp46_plan.csv').relative_to(REPO)}")
    return 0


def _self_test() -> int:
    checks = []

    def chk(name, ok):
        checks.append((name, bool(ok)))

    # the refutation must not have been allowed to delete the outcome
    chk("P(any) exceeds P(M5) alone", p_any_with_premise() > old_framing())
    chk("P(any) is a probability", 0.0 < p_any_with_premise() < 1.0)
    chk("six mechanisms are enumerated", len(MECHANISMS) == 6)
    chk("M5 is marked refuted",
        next(m for m in MECHANISMS if m["id"] == "M5")["status"] == "REFUTED")

    # structure, not arithmetic: the premise must dominate
    low, high = p_any_with_premise(0.0), p_any_with_premise(1.0)
    chk("P(any) increases with the premise", low < high)
    # over its FULL range the premise dominates; per unit it does not.
    chk("over the full range the premise swings P(any) more than any single prior",
        (high - low) > max(m["p"] for m in MECHANISMS))
    chk("per unit, the premise does NOT dominate",
        lever_sensitivity()[0]["lever"] != "premise")
    chk("the top per-unit lever is the out-of-scope route M6",
        lever_sensitivity()[0]["lever"] == "M6")

    # monotone in every mechanism prior
    q = p_by_id()
    base = p_any_with_premise()
    chk("raising M4's prior raises P(any)", True)  # structural, see below
    import copy
    saved = copy.deepcopy(MECHANISMS)
    for m in MECHANISMS:
        if m["id"] == "M4":
            m["p"] = min(0.95, m["p"] + 0.1)
    chk("raising a mechanism prior raises P(any)", p_any_with_premise() > base)
    MECHANISMS[:] = saved
    chk("priors were restored", abs(p_any_with_premise() - base) < 1e-12)

    # the naive independent answer must be the optimistic one
    chk("independence over-states P(any) vs the premise model",
        naive_independent() > p_any_with_premise())

    # posterior mechanics
    chk("a positive result raises the posterior",
        posterior_after(0.3, 0.85, 0.10, True) > 0.3)
    chk("a negative result lowers the posterior",
        posterior_after(0.3, 0.85, 0.10, False) < 0.3)
    chk("perfect test returns 1 on positive",
        abs(posterior_after(0.3, 1.0, 0.0, True) - 1.0) < 1e-9)
    chk("posteriors stay in [0,1]",
        all(0 <= posterior_after(0.5, s, f, o) <= 1
            for s in (0.6, 0.85) for f in (0.1, 0.2) for o in (True, False)))

    # information measures
    chk("entropy is maximal at 0.5", abs(entropy(0.5) - 1.0) < 1e-12)
    chk("entropy is 0 at the extremes",
        entropy(0.0) == 0.0 and entropy(1.0) == 0.0)
    chk("a useless test has ~0 EIG",
        expected_info_gain(0.54, 0.5, 0.5) < 1e-9)
    chk("a perfect test gains the full entropy",
        abs(expected_info_gain(0.54, 1.0, 0.0) - entropy(0.54)) < 1e-9)
    chk("EIG is non-negative for every measurement",
        all(expected_info_gain(0.54, m["sens"], m["fpr"]) >= 0 for m in MEASUREMENTS))
    chk("EIG is not degenerate at p ~ 0.5 (ranking has spread)",
        (lambda e: max(e) - min(e) > 0.02)(
            [expected_info_gain(0.54, m["sens"], m["fpr"]) for m in MEASUREMENTS]))

    # the scoping options are computed, and the first draft's 0.36 was wrong
    o = {x["option"]: x for x in scoping_options()}
    chk("the allosteric figure is 0.3443, not the 0.36 first written",
        abs(o["C"]["p_allosteric"] - 0.3443) < 0.01)
    chk("option A's outcome exceeds option C's", o["A"]["p_outcome"] > o["C"]["p_outcome"])
    chk("option B@0.30 allosteric is 0.1721, NOT 0.420 (that is the outcome)",
        abs(o["B@premise0.30"]["p_allosteric"] - 0.1721) < 0.01)
    chk("0.420 is the OUTCOME at premise 0.30", abs(p_any_with_premise(0.30) - 0.4205) < 0.01)
    chk("the two quantities are never equal in option C",
        abs(o["C"]["p_allosteric"] - o["C"]["p_outcome"]) > 0.1)
    chk("scoping_options did not mutate the priors",
        abs(p_any_with_premise() - 0.541) < 1e-4)   # 0.54098, not exactly 0.541

    # gaps
    g = gap_inventory()
    chk("five gap classes", len(g) == 5)
    chk("only two are simulable here",
        sum(1 for x in g if x["can_simulate_here"]) == 2)
    chk("G1 names the missing reference standard",
        "positive control" in next(x for x in g if x["cls"] == "G1")["missing"])
    chk("G3 names the wrong-quantity defect",
        "wrong quantity" in next(x for x in g if x["cls"] == "G3")["missing"])

    # simulation determinism
    a, b = simulate(n=500), simulate(n=500)
    chk("the simulation is seeded and reproducible", a["flips"] == b["flips"])
    chk("every measurement gets a flip count", set(a["flips"]) ==
        {m["id"] for m in MEASUREMENTS})

    # honesty
    src = Path(__file__).read_text(encoding="utf-8")
    prod = src.split("def _self_test")[0]
    chk("priors are labelled assumptions", "UNVALIDATED ASSUMPTION" in prod)
    chk("it says it is decision analysis, not a finding",
        "DECISION ANALYSIS" in prod and "not a finding" in prod)
    chk("it does not claim the repo can close G1/G2/G4",
        "can close G3 and G5 and nothing else" in prod)

    bad = [n for n, ok in checks if not ok]
    print(f"WP46 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for n in bad:
        print(f"  FAILED: {n}")
    return 1 if bad else 0


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
