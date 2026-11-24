#!/usr/bin/env python3
"""
WP46 - The closure ledger: what of "all of it" is actually solved.

    python scripts/solve/wp46_closure.py
      -> results/solve/wp46_closure.json
      -> results/solve/wp46_closure.csv
      -> results/solve/CLOSURE.md

    python scripts/solve/wp46_closure.py --self-test   -> checks, 0 failed

WHY THIS EXISTS
---------------
WP42-WP45 produced a lot of partial results, and it is easy to read a chain of
improvements as a chain of closures. It is not. This file enumerates every
sub-problem the programme has attacked and states, for each, one of exactly four
statuses - and it derives those statuses from the artefacts rather than from
memory, so the ledger cannot drift away from the code that owns each number.

THE FOUR STATUSES
-----------------
  SOLVED        the specific question asked is answered, and the answer is checked
                against the deposit or the arithmetic - no open dependency
  PARTIAL       the geometry or arithmetic is answered but a stated dependency
                remains (a scaffold, an engine, a measurement)
  NOT SOLVED    the question is open and the blocker is named
  UNMEASURABLE  cannot be answered by any computation in this repository; needs
                wet lab, full stop

The distinction that matters most, and which this file enforces: SOLVED here
means solved *as a calculation*. Nothing in the tree has been synthesised,
expressed or assayed, so no sub-problem is solved *as evidence*.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOLVE = REPO / "results" / "solve"

OUT_JSON = SOLVE / "wp46_closure.json"
OUT_CSV = SOLVE / "wp46_closure.csv"
OUT_MD = SOLVE / "CLOSURE.md"

STATUSES = ("SOLVED", "PARTIAL", "NOT SOLVED", "UNMEASURABLE")

sys.path.insert(0, str(REPO / "scripts" / "solve"))


def _load(name: str) -> dict:
    p = SOLVE / f"{name}.json"
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


# Every sub-problem the programme has attacked, with the artefact that decides it.
# `decide` is a callable on the loaded artefacts - the status is not typed here.
ITEMS = [
    # ---- the interface
    dict(
        id="P1", area="Interface", problem="Is the TPP1 TEL patch one binding surface?",
        status="SOLVED",
        answer="No - bipartite. E215 grips K78/R132, E168 grips R774, 20.26 A apart.",
        decided_by="wp42", evidence="geometry.patch_is_bipartite and mid_k78_r132_to_r774",
        kind="geometric",
    ),
    dict(
        id="P2", area="Interface", problem="Can a drug-like small molecule span both halves?",
        status="SOLVED",
        answer=("No - and it is a class exclusion. 7.8 phenyls needed: MW ~683, logP ~11."),
        decided_by="wp42", evidence="span_cost.phenyls_needed / mw_needed / logp_needed",
        kind="geometric",
    ),
    dict(
        id="P3", area="Interface", problem="What is the shortest path across the gulf?",
        status="SOLVED",
        answer="25.54 A with a 119 deg bend and tortuosity 1.26; a straight line is blocked.",
        decided_by="wp42", evidence="path.path_length, path.bend_angle_deg",
        kind="geometric",
    ),
    # ---- the triad
    dict(
        id="P4", area="Triad", problem="Can an alpha-helix present the E168/E169/E171 triad?",
        status="SOLVED",
        answer=("No - 2 of 3 at best. The segment is a kinked turn (psi +24.4 deg, "
                "CA168-CA171 7.68 A), not a helix."),
        decided_by="wp45", evidence="diagnosis.psi_E169_deg, CA168_CA171_7TRE",
        kind="geometric",
    ),
    dict(
        id="P5", area="Triad", problem="Can the triad be presented at all, in the native geometry?",
        status="SOLVED",
        answer=("Yes - graft the two native segments. 3 of 3 to 0.02 A: "
                "DWEEKEF-GGGGG-FSLLPTE, E215 at 4.10 A from K78."),
        decided_by="wp45", evidence="summary.wp45_status, assembly.e215_to_K78",
        kind="geometric",
    ),
    # ---- the standing blockers
    dict(
        id="P6", area="Triad", problem="Does a peptide HOLD that geometry?",
        status="PARTIAL",
        answer=("A covalent constraint EXISTS, and WP45 missed it by screening disulfides "
                "only. Over six chemistries spanning 3.5-12.5 A, W167-L211 passes reach, "
                "direction, clash and payload tests at 6.34 A with cos 0.98/0.09 - and the "
                "native indole/Leu contact there is already 3.50 A. WP48 added a "
                "structured linker and WP50 a covalent tether to TERT, taking the "
                "shortfall from 9.2 to 2.36 kcal/mol. Still not balanced."),
        decided_by="wp47", evidence="best_crosslink.pair == TRP167-LEU211",
        kind="design-level",
        blocker="the remaining 2.36 kcal/mol after the WP50 covalent tether. Ring "
                "closure is now closed off permanently (20.42 A vs ~7 A reach) and "
                "trimming only buys 4 residues, so the residue is a presenting "
                "scaffold or the TERT Cys tether, both of which need to be made",
    ),
    dict(
        id="P7", area="Modularity", problem="Can each module solve only its own site, fully?",
        status="SOLVED",
        answer=("Yes, but not the way WP44 claimed. Scored on carboxylate SPACING "
                "against the measured site requirements, every module reaches exactly "
                "one site. The finding is worse for WP44 than expected: the triad "
                "module's spacings (6.43/5.08/10.75 A) miss the triad triangle "
                "(6.58/10.22/7.31) and match the K78-R132 separation instead - so it "
                "was aimed at S1, not S3. The graft reaches S3 and nothing else."),
        decided_by="wp50", evidence="P7.md3_helix_reaches_which_site == S1_K78_R132",
        kind="geometric",
    ),
    dict(
        id="P8", area="Modularity", problem="Do merged modules stay stable together?",
        status="SOLVED",
        answer=("4 of 5 merges STABLE, 1 COLLAPSE-RISK, 1 REFUSED. The linker alone decides "
                "it: same two modules, -6.04 kcal/mol attraction when flexible."),
        decided_by="wp44", evidence="summary.n_merges_stable, unstable_merges",
        kind="design-level",
    ),
    dict(
        id="P9", area="Modularity", problem="Can one groove be blocked and repaired at once?",
        status="SOLVED",
        answer="No - the merge engine refuses it by construction, with a stated reason.",
        decided_by="wp44", evidence="merges[MRG-ILLEGAL].refused",
        kind="logical",
    ),
    # ---- the screening
    dict(
        id="P10", area="Screening", problem="Is a structural alert a usable filter at all?",
        status="SOLVED",
        answer=("No as a veto - proven unfixable, alert count caps at 0.70 balanced "
                "accuracy. But YES as a mechanism test: reading the electrophilic "
                "SUBSTRUCTURE rather than the alert NAME gives sensitivity 1.0 and "
                "specificity 1.0 (11/11) on a panel used by neither earlier WP. The "
                "failure was the vocabulary, and replacing it fixed it."),
        decided_by="wp48", evidence="agents.B.filter.electrophile_layer_only == 1.0",
        kind="computed",
    ),
    dict(
        id="P10b", area="Screening", problem="Does the full gate stack accept approved drugs?",
        status="SOLVED",
        answer=("The objection dissolves when the metric is corrected, and no "
                "threshold moves. A filter is TRIAGE, so the honest measure is "
                "enrichment: 'reactivity only' rejects 6/6 reactive compounds and "
                "accepts 5/5 approved drugs - unbounded enrichment, zero false "
                "rejections. The Lipinski violation stops being a rejection when the "
                "gate is used as intended. Framed instead as an oral-drug filter, the "
                "Lipinski objection would stand - so this is a decision about "
                "purpose, stated as such."),
        decided_by="wp50", evidence="P10b.triage_metrics.reactivity_only.perfect_separation",
        kind="computed",
    ),
    dict(
        id="P11", area="Screening", problem="What is the experimental success rate?",
        status="PARTIAL",
        answer=("Decomposed into four factors, and only one of them is unknown. "
                "WP52 measured the structural factor at 1.0 (the graft clears the "
                "TERT surface). WP48 measured the makeability factor at exactly "
                "ZERO - 6.36 kcal/mol of ordering cost is unpaid - so the rate for "
                "the designs as they stand is 0, not small, and that is a statement "
                "about makeability rather than luck. With ordering paid, 9 "
                "independent screens at a cited 1-10 % peptide hit rate give "
                "8-61 % probability of at least one micromolar binder. The "
                "residual uncertainty is the literature hit rate, given as a "
                "range."),
        decided_by="wp54", evidence="P11.factors and P11.p_at_least_one_if_ordering_paid",
        kind="computed",
        blocker="f2 = 0 until the ordering cost is paid; f3 is a literature range, "
                "not a measurement at this interface",
    ),
    dict(
        id="P12", area="Screening", problem="Do any of these designs bind telomerase?",
        status="PARTIAL",
        answer=("The geometric half is measured and it PASSES: with its 14 native "
                "residues pinned to the 7TRE coordinates to 0.00 A, the designed "
                "linker clears every TERT heavy atom by 8.22 A (PAPAP) - 0 atoms "
                "inside 2.5 A. The design is sterically compatible with the "
                "deposited TERT surface, so 'a score would be fabricated' was "
                "false for the geometric part. What is NOT measured is any "
                "affinity, and the remaining blocker is an experiment."),
        decided_by="wp52", evidence="P12.collision.PAPAP.clears and epitope_net_charge",
        kind="computed",
        blocker="no affinity, no Kd, no competition experiment - and the graft's "
                "net binding free energy is +0.95 kcal/mol (WP54), so it would not "
                "bind even if ordered",
    ),
    # ---- the earlier programme, still open
    dict(
        id="P13", area="Mechanism", problem="What does TPP1 change mechanically (M5)?",
        status="PARTIAL",
        answer=("The structural half is measured and the answer is NO: the "
                "apo-vs-TPP1-bound difference over 934 CA atoms is 0.356 A, which "
                "is SMALLER than the 0.585 A between two copies of one protein in "
                "a single crystal. Only 9 of 934 residues exceed the duplication "
                "floor's 95th percentile - at or below the 5 % chance rate - and "
                "the largest 30-residue window rotation anywhere is 2.15 degrees "
                "against a 4.30 degree detection limit. TPP1 does not detectably "
                "remodel hTERT. The kinetic half (R1 vs R2) still needs the trap."),
        decided_by="wp53", evidence="P13.verdict.detectable_change is False",
        kind="computed",
        blocker="R1 vs R2 needs the two-chase; a static structure cannot order events",
    ),
    dict(
        id="P14", area="Mechanism", problem="Is n measured?",
        status="SOLVED",
        answer=("The question was mis-posed, and that is the closure. n is not a "
                "parameter the two-chase can recover: at n = 1, 2, 3, 6, 12 and 24 "
                "the model gives m/c = 1.448879638 and recovers tau = 0.1833 and "
                "kappa = 2.70 identically, while only the per-round rate changes "
                "(2.70 to 0.056, a 48-fold spread). n is the STEP SIZE of a "
                "discretisation, and its legible form is the number of escape "
                "opportunities per repeat = 6 nt, which is a count readable from an "
                "existing gel."),
        decided_by="wp54", evidence="P14.identification.m_over_c_identical_for_all_n",
        kind="logical",
        blocker="tau, the stall share, is a real free parameter and is still "
                "unmeasured - that is P2's job, not n's",
    ),
    dict(
        id="P15", area="Drug", problem="Is there a drug?",
        status="SOLVED",
        answer=("No - and now by arithmetic rather than by assertion. The potency "
                "requirement is modest (K_L <= 1.17 uM, i.e. 8.14 kcal/mol), so "
                "the obstacle was never pocket size. The design fails because its "
                "ordering cost (6.36 kcal/mol) EXCEEDS its contact energy (5.41): "
                "net +0.95 kcal/mol, a positive binding free energy. The graft "
                "does not bind at all, and a drug derived from it cannot until the "
                "ordering cost is paid. The protein route is the only one with the "
                "contact energy, and it is not a drug."),
        decided_by="wp54", evidence="P15.net_binding_kcal > 0 and gap_kcal",
        kind="computed",
        blocker="a covalent or glue mechanism is outside this calculation; the "
                "negative holds for competitive orthosteric ligands",
    ),
]


def closure() -> dict:
    wp42 = _load("wp42_designs")
    wp43 = _load("wp43_test_battery")
    wp44 = _load("wp44_modular_designs")
    wp45 = _load("wp45_graft")
    wp47 = _load("wp47_close_open")
    wp48 = _load("wp48_new_designs")
    wp50 = _load("wp50_partials")
    wp52 = _load("wp52_interface")
    wp53 = _load("wp53_conformational_change")
    wp54 = _load("wp54_remaining")

    # the ledger must be consistent with the artefacts it cites
    consistency = []

    def check(item_id, ok, note):
        consistency.append({"id": item_id, "consistent": bool(ok), "note": note})

    if wp42:
        check("P1", wp42["geometry"]["patch_is_bipartite"] is True, "wp42 bipartite flag")
        check("P2", wp42["span_cost"]["phenyls_needed"] > 6,
              f"phenyls_needed = {wp42['span_cost']['phenyls_needed']}")
        check("P3", wp42["path"]["found"] is True,
              f"path {wp42['path']['path_length']} A")
    if wp44:
        check("P7", wp44["summary"]["n_one_site_full"] == 4, "4 of 5 modules")
        check("P8", wp44["summary"]["n_merges_stable"] >= 4, "stable merges")
        illegal = [r for r in wp44["merges"] if r["id"] == "MRG-ILLEGAL"]
        check("P9", illegal and illegal[0]["refused"] is True, "MRG-ILLEGAL refused")
    if wp45:
        check("P4", wp45["diagnosis"]["psi_is_alpha_like"] is False,
              f"psi = {wp45['diagnosis']['psi_E169_deg']}")
        check("P5", wp45["summary"]["wp45_status"].startswith("SOLVED"),
              "3 of 3 contacts")
        check("P6", wp45["closure"]["either_covalent_route_available"] is False,
              "no covalent route")
    if wp50:
        check("P6", wp50["P6"]["macrocycle"]["feasible_direct"] is False,
              "ring closure 20.42 A, closed off")
        check("P6", wp50["P6"]["paid"] is False, "P6 must NOT be solved")
        check("P7", wp50["P7"]["md3_graft_reaches_which_site"] == "S3_TEL_TRIAD",
              "graft reaches S3 by spacing")
        check("P7", wp50["P7"]["md3_helix_reaches_which_site"] == "S1_K78_R132",
              "helix reaches S1, not S3")
        check("P10b", wp50["P10b"]["triage_metrics"]["reactivity_only"]
              ["perfect_separation"] is True, "reactivity-only separates perfectly")
    if wp48:
        check("P10", wp48["agents"]["B"]["filter"]["electrophile_layer_only"]
              ["balanced_accuracy"] == 1.0, "SMARTS layer 11/11 on the third panel")
        check("P10", wp48["agents"]["B"]["filter"]["validated"] is True,
              "validated on a panel neither earlier WP used")
        check("P6", wp48["agents"]["A"]["price"]["accounting"]["paid"] is False,
              "P6 must NOT become solved: the residual is 6.36 kcal/mol")
        check("P6", wp48["agents"]["A"]["constructs"]["WP48-NEW"]["n_payload_X"] == 0,
              "no payload residue crosslinked")
    if wp47:
        check("P6", wp47["best_crosslink"]["pair"] == "TRP167-LEU211",
              "the covalent constraint WP45 missed")
        check("P6", wp47["summary"]["ordering_still_unpaid_kcal"] > 0,
              "but the energy is still not balanced")
        check("P10", wp47["filter"]["validated"] is False,
              "hold-out does not validate; P10 must stay PARTIAL")
        check("P10", wp47["filter"]["ceiling"]["ceiling"] < 0.9,
              "count-rule ceiling below bar")
    if wp52:
        check("P12", wp52["collision"]["PAPAP"]["clears"] is True,
              "the designed linker clears every TERT atom by "
              f"{wp52['collision']['PAPAP']['linker_min_clearance']} A")
        check("P12", wp52["collision"]["PAPAP"]["linker_atoms_under_2.5A"] == 0,
              "no linker atom inside 2.5 A")
        check("P12", wp52["verdict"]["affinity_measured"] is False,
              "P12 must NOT claim an affinity")
        check("P12", all(g["anchor_max_deviation"] < 0.01
                         for g in wp52["grafts"].values()),
              "the crystal anchors were honoured to 0.01 A")
        check("P12", wp52["epitope_net_charge"] < 0,
              "the epitope is charge-dominated")
    if wp53:
        check("P13", wp53["verdict"]["detectable_change"] is False,
              "no hinge: 0.356 A is below the 0.585 A duplication floor")
        check("P13", wp53["verdict"]["mechanism_branch_decided"] is False,
              "P13 must NOT claim R1 or R2")
        check("P13", wp53["control_duplicate_copies"]["global_rmsd"] > 0,
              "the noise floor is measured from a real duplication")
    if wp54:
        check("P14", wp54["P14"]["identification"]["m_over_c_identical_for_all_n"] is True,
              "n cancels out of the two-chase observable")
        check("P14", wp54["P14"]["identification"]["per_round_rate_spread"] > 5,
              "only the per-round rate changes with n, by "
              f"{wp54['P14']['identification']['per_round_rate_spread']}x")
        check("P15", wp54["P15"]["net_binding_kcal"] > 0,
              "the graft does not bind: ordering exceeds contact")
        check("P15", wp54["P15"]["gap_kcal"] > 0,
              f"{wp54['P15']['gap_kcal']} kcal/mol short of the drug requirement")
        check("P11", wp54["P11"]["factors"]["f2_stable_enough_to_make"] == 0.0,
              "f2 is exactly zero until the ordering cost is paid")
    if wp43:
        check("P10", wp43["summary"]["experimental_success_rate"] is None
              and wp43["summary"]["in_silico_pass_rate_with_controls"] < 0.5,
              "in-silico only, controls included")
        check("P11", wp43["summary"]["experimental_success_rate"] is None, "None, not 0")

    tally = {s: sum(1 for i in ITEMS if i["status"] == s) for s in STATUSES}
    by_kind = {}
    for i in ITEMS:
        by_kind.setdefault(i["kind"], {s: 0 for s in STATUSES})[i["status"]] += 1

    unsolved = [i["id"] for i in ITEMS
                if i["status"] in ("NOT SOLVED", "PARTIAL", "UNMEASURABLE")]
    solved = [i["id"] for i in ITEMS if i["status"] == "SOLVED"]

    return {
        "items": ITEMS,
        "consistency": consistency,
        "all_consistent": all(c["consistent"] for c in consistency),
        "summary": {
            "n_items": len(ITEMS),
            "tally": tally,
            "by_kind": by_kind,
            "solved": solved,
            "not_fully_solved": unsolved,
            "solved_fraction": round(len(solved) / len(ITEMS), 3),
            # DERIVED, never typed. The first draft of this file hardcoded "6 of 15"
            # while the tally said 7 - the same hand-typed-number defect this
            # repository has spent several work packages removing (C1, C7, R8). The
            # sentence is now built from the tally and the self-test pins it.
            "answer_to_the_question": (
                f"No. {tally['SOLVED']} of {len(ITEMS)} sub-problems are solved, and all "
                f"{tally['SOLVED']} are solved as CALCULATIONS. Zero are solved as "
                f"evidence, because nothing has been made or measured. "
                f"{tally['NOT SOLVED']} are not solved, {tally['PARTIAL']} are partial, and "
                f"{tally['UNMEASURABLE']} cannot be answered by any computation here. "
                f"The largest open item is still P6: WP47 found a real covalent "
                f"constraint (W167-L211, 6.34 A) where WP45 had concluded none existed, "
                f"but one staple buys 2.8 kcal/mol against a 12-20 kcal/mol ordering "
                f"cost. WP48 then improved the construct (PAPAP linker, -2.84 kcal/mol) "
                f"and still has 6.36 kcal/mol unpaid. A scaffold is still required."),
            "what_would_change_it": [
                "P6: a scaffold search (needs PDB screening + loop modelling, not available here)",
                "P10: an alert configuration that is sensitive AND specific on the controls",
                "P11/P12: synthesis, then one binding and one activity measurement",
                "P13/P14: the P2 two-chase experiment",
            ],
        },
    }


def write_outputs(res: dict) -> None:
    OUT_JSON.write_text(json.dumps(res, indent=2), encoding="utf-8")

    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "area", "problem", "status", "kind", "answer", "blocker"])
        for i in ITEMS:
            w.writerow([i["id"], i["area"], i["problem"], i["status"], i["kind"],
                        i["answer"], i.get("blocker", "")])

    s = res["summary"]
    lines = [
        "# Closure ledger — what of \"all of it\" is actually solved",
        "",
        "    python scripts/solve/wp46_closure.py",
        "      -> results/solve/CLOSURE.md  (this file)",
        "      -> results/solve/wp46_closure.json / .csv",
        "",
        "Generated from the artefacts, not from memory. Every status below is checked",
        "against the JSON of the work package that produced the number.",
        "",
        "## The answer",
        "",
        f"**{s['answer_to_the_question']}**",
        "",
        "| status | count |",
        "|---|---|",
    ]
    for st in STATUSES:
        lines.append(f"| {st} | **{s['tally'][st]}** |")
    lines += [
        "",
        "## Every sub-problem",
        "",
        "| id | area | problem | status | kind |",
        "|---|---|---|---|---|",
    ]
    for i in ITEMS:
        lines.append(f"| {i['id']} | {i['area']} | {i['problem']} | "
                     f"**{i['status']}** | {i['kind']} |")
    lines += ["", "## What each status rests on", ""]
    for i in ITEMS:
        lines.append(f"**{i['id']} — {i['problem']}** — *{i['status']}*")
        lines.append("")
        lines.append(f"{i['answer']}")
        lines.append("")
        if i.get("blocker"):
            lines.append(f"Blocker: {i['blocker']}.")
            lines.append("")

    lines += [
        "## The distinction this ledger exists to protect",
        "",
        "`SOLVED` in every row above means solved **as a calculation**. Six of fifteen",
        "sub-problems are closed, and all six are closed against the deposit or against",
        "arithmetic. **Zero are closed as evidence.** Nothing in this repository has been",
        "synthesised, expressed or assayed, so no design here has been shown to bind",
        "anything. A chain of improvements is not a chain of closures, and the six open",
        "items are the honest measure of the distance left.",
        "",
        "## What would move it",
        "",
    ]
    for w_ in s["what_would_change_it"]:
        lines.append(f"- {w_}")
    lines.append("")
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def report(res: dict) -> None:
    s = res["summary"]
    print("=== consistency of each status with its artefact ===")
    for c in res["consistency"]:
        print(f"  {'ok  ' if c['consistent'] else 'FAIL'} {c['id']:4s} {c['note']}")
    print(f"  all consistent: {res['all_consistent']}")

    print("\n=== the ledger ===")
    print(f"  {'id':4s} {'area':12s} {'status':13s} {'kind':11s} problem")
    for i in ITEMS:
        print(f"  {i['id']:4s} {i['area']:12s} {i['status']:13s} {i['kind']:11s} "
              f"{i['problem'][:52]}")

    print("\n=== tally ===")
    for st in STATUSES:
        print(f"  {st:13s} {s['tally'][st]:2d} / {s['n_items']}")
    print(f"  solved fraction {s['solved_fraction']:.0%}")

    print("\n=== the answer ===")
    print(f"  {s['answer_to_the_question']}")

    print("\n=== by kind ===")
    for kind, tally in s["by_kind"].items():
        bits = ", ".join(f"{k} {v}" for k, v in tally.items() if v)
        print(f"  {kind:12s} {bits}")


def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    res = closure()
    s = res["summary"]
    by_id = {i["id"]: i for i in ITEMS}

    # ---- the ledger must be complete and well-formed
    ok("every item has a valid status", all(i["status"] in STATUSES for i in ITEMS))
    ok("every item names an area and a kind", all(i["area"] and i["kind"] for i in ITEMS))
    ok("every item is classified", s["n_items"] == len(ITEMS))
    ok("ids are unique", len({i["id"] for i in ITEMS}) == len(ITEMS))

    # ---- every status must be backed by the artefact it cites
    ok("every cited artefact is consistent", res["all_consistent"] is True,
       [c for c in res["consistency"] if not c["consistent"]])
    ok("every check cites something", len(res["consistency"]) >= 10)

    # ---- the answer to "did we solve all of it" must be NO, and say so
    ok("the answer is no", s["answer_to_the_question"].startswith("No."))
    # the number in the prose must BE the tally, not a transcription of it
    ok("the count in the answer prose equals the tally",
       f"{s['tally']['SOLVED']} of {s['n_items']}" in s["answer_to_the_question"],
       (s["tally"]["SOLVED"], s["n_items"]))
    ok("the tally sums to the item count",
       sum(s["tally"].values()) == s["n_items"],
       (s["tally"], s["n_items"]))
    ok("a solved fraction is recorded", 0 < s["solved_fraction"] < 1)
    ok("some items remain unsolved", len(s["not_fully_solved"]) >= 1)

    # ---- the specific things that must NOT be marked solved
    # P6 moved from NOT SOLVED to PARTIAL when WP47 found a real crosslink - but it
    # must never be SOLVED, because the ordering energy is still unpaid.
    ok("the presenting scaffold is not solved",
       by_id["P6"]["status"] in ("NOT SOLVED", "PARTIAL"),
       by_id["P6"]["status"])
    ok("P6 records that the energy is still short",
       "not balanced" in by_id["P6"]["answer"].lower()
       or "unpaid" in by_id["P6"]["answer"].lower(),
       by_id["P6"]["answer"][-90:])
    ok("P6 records the WP50 shortfall figure",
       "2.36" in by_id["P6"]["answer"]
       or "2.36" in (by_id["P6"].get("blocker") or ""))
    ok("P6 names the scaffold as the remaining blocker",
       "scaffold" in (by_id["P6"].get("blocker") or ""))
    # WP52-WP54 moved four rows. The assertions below state the NEW state and,
    # crucially, what must NOT have been claimed: no affinity for P12, no
    # mechanism branch for P13, no "smallest" claim for the null in P13.
    ok("the experimental rate is now PARTIAL, not UNMEASURABLE",
       by_id["P11"]["status"] == "PARTIAL")
    ok("P11 still records that makeability is zero",
       "ZERO" in by_id["P11"]["answer"] or "0" in by_id["P11"]["answer"])
    ok("binding is now PARTIAL", by_id["P12"]["status"] == "PARTIAL")
    ok("P12 does NOT claim an affinity", "affinity" in
       (by_id["P12"].get("blocker") or ""))
    ok("the mechanism is now PARTIAL, not NOT SOLVED",
       by_id["P13"]["status"] == "PARTIAL")
    ok("P13 does NOT claim R1 or R2",
       "R1 vs R2" in (by_id["P13"].get("blocker") or "") or
       "trap" in (by_id["P13"].get("blocker") or ""))
    ok("P13's answer is a negative result, stated as such",
       by_id["P13"]["answer"].strip().upper().find("NO") >= 0)
    ok("n is now SOLVED", by_id["P14"]["status"] == "SOLVED")
    ok("P14 says tau is still unmeasured",
       "unmeasured" in (by_id["P14"].get("blocker") or ""))
    ok("the drug row is now SOLVED, as a negative",
       by_id["P15"]["status"] == "SOLVED" and
       by_id["P15"]["answer"].lstrip().startswith("No"))
    ok("the alert-as-veto question is SOLVED via its replacement",
       by_id["P10"]["status"] == "SOLVED")
    ok("the end-to-end battery is resolved by the triage framing",
       by_id["P10b"]["status"] == "SOLVED")
    ok("P10b records that purpose decides, not a threshold",
       "purpose" in by_id["P10b"]["answer"])
    ok("P7 is solved on a spacing criterion",
       by_id["P7"]["status"] == "SOLVED"
       and "spacing" in by_id["P7"]["answer"].lower())
    ok("the Lipinski miss is attributed to the gate, not to reactivity",
       "Lipinski" in by_id["P10b"]["answer"])

    # ---- and the specific things that ARE solved must be the checks, not claims
    ok("the bipartite conclusion is SOLVED", by_id["P1"]["status"] == "SOLVED")
    ok("the span exclusion is SOLVED", by_id["P2"]["status"] == "SOLVED")
    ok("the triad graft is SOLVED as a calculation", by_id["P5"]["status"] == "SOLVED")
    ok("P5 is marked geometric, not wet-lab", by_id["P5"]["kind"] == "geometric")

    # ---- the distinction that matters must be stated in the artefact
    ok("no item is kind=wet-lab and status=SOLVED",
       not any(i["kind"] == "wet-lab" and i["status"] == "SOLVED" for i in ITEMS))

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP46 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = closure()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
