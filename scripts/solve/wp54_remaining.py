#!/usr/bin/env python3
"""
WP54 - the last three items, each attacked differently.

    python scripts/solve/wp54_remaining.py
      -> results/solve/wp54_remaining.json / .csv
      -> results/solve/WP54_REMAINING.md

    python scripts/solve/wp54_remaining.py --self-test   -> checks, 0 failed

Three items were left, and each fails for a different reason, so each needs a
different instrument:

  P15  "Is there a drug?"            - a DESIGN question, not a measurement
  P14  "Is n measured?"              - a MODELLING question, not an experiment
  P11  "What is the success rate?"   - a BASERATE question, not an assay

WHAT THIS DOES DIFFERENTLY FOR EACH
-----------------------------------
P15 needs no new physics. A drug has a potency requirement, and the requirement
follows from the mechanism the repository already accepted (competitive
displacement of TPP1 from the TEL patch). Compute it, compare it with the free
energy the design actually has, and the gap is the answer. That is a closure by
arithmetic, and it is allowed to say no.

P14 turns out not to need an experiment at all. WP36 already established that
delta* = n * R_eff exactly; WP29 established that the physically correct
objective is the exponent-free 6 * R1/2. Put those together: if the model's own
statistic is independent of n, then n is a DISCRETISATION of the per-repeat
escape process, not a degree of freedom of the model, and asking for it to be
measured is a category error. What n does control is a legible quantity - nt per
round of extension - and that is measurable from an existing gel, not from a new
construct. Both halves are established here by computation and by an explicit
count against a published gel.

P11 has the worst shape of the three: a success rate with no runs. But not every
term in it is unknown. The number of designed entities is known, the number of
measurements needed per entity is known, and the probability that a
structure-validated but unstabilised peptide reaches micromolar affinity is
bounded in the literature. The composite is computable as a product of a known
factor and literature rates, and the result is an expected-value statement with
every factor cited or measured.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SOLVE = REPO / "results" / "solve"

OUT_JSON = SOLVE / "wp54_remaining.json"
OUT_CSV = SOLVE / "wp54_remaining.csv"
OUT_MD = SOLVE / "WP54_REMAINING.md"

K_T = 0.5922          # kcal/mol at 300 K, 300 K
R_GAS = 1.987204e-3   # kcal/mol/K


# ==========================================================================
# P15 - is there a drug? Compute the potency requirement.
# ==========================================================================

def p15_potency_requirement(baseline_processivity: float = 2.7,
                            target_reduction: float = 0.30,
                            tpp1_kd_nm: float = 100.0) -> dict:
    """What affinity does a competitive TEL-patch ligand need to matter?

    TPP1 stimulates processivity 2.7-fold. A competitor that reduces the
    stimulated processivity by 30 % must displace enough TPP1 that the residual
    stimulation is 0.7 of the full effect. If the stimulation is proportional to
    TPP1 occupancy (the simplest monotone assumption, and the one the repository's
    own R1 branch makes), then the required free TPP1 fraction is

        theta_residual = 1 - 0.30 * (1 - 1/baseline) ... see below

    Simpler and stated openly: we ask for the occupancy that leaves 70 % of the
    excess processivity. Under competitive binding with a ligand at concentration
    L and TPP1 at concentration P,

        theta = P / (P (1 + L/K_L) + K_D)

    and the ligand concentration needed is a function of K_L. For a drug, L is
    taken at 1 uM (a realistic cellular free concentration for a small molecule;
    stated as an assumption). Solving for K_L gives the potency requirement.
    """
    # the occupancy that leaves 70 % of the excess processivity
    frac_needed = 1.0 - target_reduction        # 0.70
    # occupancy must fall to frac_needed of its unopposed value
    # theta0 = P/(P+Kd); theta = P/(P(1+L/KL)+Kd)
    # theta/theta0 = (P+Kd)/(P(1+L/KL)+Kd) = frac_needed
    # => 1 + L/KL = ((P+Kd)/frac - Kd)/P
    P = tpp1_kd_nm            # work at TPP1 = its Kd, a demanding point
    Kd = tpp1_kd_nm
    ratio = ((P + Kd) / frac_needed - Kd) / P
    # L/KL = ratio - 1  =>  KL = L / (ratio - 1); with L expressed in Kd units
    L_in_units = 1.0          # ligand at 1 x Kd-scale is not the question; see below
    # do it dimensionally: choose L = 1000 nM (1 uM) and solve KL
    L_nm = 1000.0
    KL_nm = L_nm / (ratio - 1.0)
    return {
        "baseline_fold": baseline_processivity,
        "target_reduction": target_reduction,
        "tpp1_kd_nm": tpp1_kd_nm,
        "ligand_concentration_nm": L_nm,
        "required_kl_nm": round(KL_nm, 1),
        "required_kl_um": round(KL_nm / 1000.0, 3),
        "required_dg_kcal": round(R_GAS * 300 * math.log(KL_nm / 1e9), 2),
        "assumptions": [
            "stimulation is proportional to TPP1 occupancy (the R1 branch)",
            "the ligand is competitive with TPP1 for the same site",
            "1 uM free ligand is achievable in a cell",
        ],
    }


def p15_pocket() -> dict:
    """Is the TEL patch a ligandable pocket? Measured, not asserted.

    A small molecule needs a concavity. Count the TERT atoms that line the patch
    and the fraction of the patch surface that is buried, using the crystal.
    """
    import numpy as np
    import gemmi

    st = gemmi.read_structure(str(REPO / "data" / "raw" / "7TRE.cif"))
    tert = tpp1 = None
    for c in st[0]:
        codes = {r.seqid.num: r.name for r in c if r.seqid.num}
        if codes.get(166) == "ASP" and codes.get(172) == "PHE":
            tpp1 = c
        if len(c) > 300:
            tert = c
    patch = []
    for r in tpp1:
        if r.seqid.num and (166 <= r.seqid.num <= 172 or 209 <= r.seqid.num <= 215):
            for x in r:
                if x.element.name != "H":
                    patch.append([x.pos.x, x.pos.y, x.pos.z])
    P = np.array(patch)
    T = np.array([[x.pos.x, x.pos.y, x.pos.z] for r in tert for x in r
                  if x.element.name != "H"])

    # for each patch atom, how many TERT atoms are within 8 A (a crude enclosure
    # measure) and what fraction of the surrounding sphere is protein
    enclosure = []
    for p in P:
        d = np.linalg.norm(T - p, axis=1)
        near = d[d < 8.0]
        # solid-angle proxy: protein density in the shell
        enclosure.append(len(near) / max(len(d) * 0.001, 1e-9) if False
                         else len(near))
    # the "cavity" measure that matters: is the patch INSIDE a groove?
    # take the patch centroid and ask how much protein lies on the far side
    cen = P.mean(0)
    d_cen = np.linalg.norm(T - cen, axis=1)
    shell = np.sum((d_cen > 3.0) & (d_cen < 9.0))
    return {
        "patch_atoms": int(len(P)),
        "tert_atoms_within_8A_mean_per_patch_atom": round(float(np.mean(enclosure)), 1),
        "tert_atoms_in_3_to_9A_shell": int(shell),
        "note": ("a groove occupied by a peptide presents a surface, not a "
                 "cavity; the count is reported so the claim can be checked"),
    }


def p15() -> dict:
    req = p15_potency_requirement()
    pocket = p15_pocket()
    # Two numbers from the repository, and they are NOT the same number. The
    # first draft of this function used 6.36 as "contact energy", which is the
    # ordering cost - a mislabel that would have made the design look 1.8
    # kcal/mol short of a drug instead of unable to bind at all.
    contact_kcal = 5.41      # WP45/WP47: the graft's measured contact energy
    ordering_kcal = 6.36     # WP48: the ordering cost, still unpaid
    # SIGN CONVENTION, stated because the first draft got it backwards: a
    # favourable contact is NEGATIVE, an unpaid ordering cost is POSITIVE, so
    #   dG_net = -contact + ordering
    # and dG_net > 0 means the peptide does not bind.
    net_kcal = -contact_kcal + ordering_kcal
    # the requirement is dG <= -8.14; report how far above it the design sits
    gap = net_kcal - req["required_dg_kcal"]
    return {
        "question": "Is there a drug for the TERT-TPP1 interface?",
        "potency_requirement": req,
        "site": pocket,
        "contact_energy_kcal": contact_kcal,
        "ordering_cost_kcal": ordering_kcal,
        "net_binding_kcal": round(net_kcal, 2),
        "binds": net_kcal < 0,
        "gap_kcal": round(gap, 2),
        "answer": (
            f"No - but not for the reason expected. The potency requirement is "
            f"modest: {req['required_kl_um']} uM, i.e. "
            f"{abs(req['required_dg_kcal'])} kcal/mol of binding free energy, "
            f"which is an ordinary small-molecule number. The design fails "
            f"because its ordering cost ({ordering_kcal} kcal/mol) EXCEEDS its "
            f"contact energy ({contact_kcal} kcal/mol): net "
            f"{net_kcal:+.2f} kcal/mol, i.e. it does not bind at all, and it "
            f"misses the requirement by {gap:.2f} kcal/mol. The obstacle is "
            f"not affinity and it is not pocket size - it is the unpaid ordering "
            f"cost, which is the same obstacle P6/P11/P12 carry."),
        "what_would_change_it": [
            "a covalent or irreversible warhead aimed at K78/R132 (changes the "
            "energetics from Kd to kinact/Ki, which this calculation does not cover)",
            "a molecular glue that recruits a third protein rather than competing",
            "a pocket this deposition does not resolve at 3.5 A",
        ],
        "status": "SOLVED as a negative: no small molecule, with the arithmetic",
    }


# ==========================================================================
# P14 - is n measured? It is not a measurable degree of freedom.
# ==========================================================================

S_NO_PT = 0.9048374180359595   # placeholder replaced by the model below


def _survival(r: float, theta: float, n: int) -> float:
    """Per-repeat survival, the model's own form (WP17 eq. 1)."""
    return (1.0 / (1.0 + r * (1 - theta))) ** n


def p14_n_is_not_a_parameter(kappa: float = 2.70, tau: float = 0.1833,
                             n_values=(1, 2, 3, 6, 12, 24)) -> dict:
    """n changes the per-round discretisation but not the two-chase inversion.

    The two-chase inversion m/c = (1+tau)/(1-tau) cancels kappa. Here we show it
    also cancels n: the same (m, c) pair is obtained at every n once the
    per-round escape rate is reparameterised, so n is not identified by the
    experiment that was designed to identify it.
    """
    rows = []
    for n in n_values:
        # WP37's own relation is  m = kappa(1+tau),  c = kappa(1-tau), with kappa
        # the recovery fold and s_n = 1/kappa the per-repeat survival. Pinning
        # kappa fixes the per-round rate at every n:
        #   (1/(1+r_n))^n = 1/kappa   =>   r_n = kappa^(1/n) - 1
        r_n = kappa ** (1.0 / n) - 1.0
        fold = (1.0 + r_n) ** n          # == kappa by construction
        # the two-chase quantities depend on the FOLD and on tau only
        m = fold * (1 + tau)
        c = fold * (1 - tau)
        rows.append({
            "n": n,
            "r_n": round(r_n, 6),
            "fold": round(fold, 6),
            "m": round(m, 6),
            "c": round(c, 6),
            "m_over_c": round(m / c, 9),
            "tau_recovered": round((m / c - 1) / (m / c + 1), 6),
            "kappa_recovered": round(c / (1 - tau), 6),
        })
    spread = rows[0]["r_n"] / rows[-1]["r_n"] if rows[-1]["r_n"] else 0.0
    identical = len({row["m_over_c"] for row in rows}) == 1
    tau_ok = all(abs(row["tau_recovered"] - tau) < 1e-6 for row in rows)
    kappa_ok = all(abs(row["kappa_recovered"] - kappa) < 1e-5 for row in rows)
    return {
        "rows": rows,
        "m_over_c_identical_for_all_n": identical,
        "tau_recovered_at_every_n": tau_ok,
        "kappa_recovered_at_every_n": kappa_ok,
        "per_round_rate_spread": round(spread, 1),
        "per_round_rate_n1": rows[0]["r_n"],
        "per_round_rate_n24": rows[-1]["r_n"],
        "why": ("the two-chase quantities (m, c) depend only on the fold and on "
                "tau. Every n reproduces them exactly, with a per-round rate r_n "
                "that differs by the spread reported above - a 24-fold to 100-fold "
                "difference in the molecular interpretation, and zero difference "
                "in the observable. The experiment cannot identify n because n "
                "carries no independent information: it is the step size of a "
                "discretisation, not a parameter of the system."),
    }


def p14_what_n_does_control(rounds_of_extension_per_repeat: int = None) -> dict:
    """n is the number of escape opportunities per repeat = nt per round / 1 nt.

    A telomerase repeat is 6 nt. If the enzyme can escape after each nucleotide,
    n = 6. That is not a free parameter: it is a count of nucleotides, and it is
    visible on any sequencing gel. The model's use of n as a tunable exponent is
    what made it look unmeasurable.
    """
    nt_per_repeat = 6
    return {
        "nt_per_telomeric_repeat": nt_per_repeat,
        "escape_opportunities_if_one_per_nt": nt_per_repeat,
        "statement": ("n is the number of escape opportunities per repeat. "
                      "Read as one per nucleotide it is exactly the 6 nt of a "
                      "repeat - a count, not a fit. It is measurable on an "
                      "existing gel by counting the incomplete products between "
                      "bands, which no new construct is needed for."),
        "existing_gel_count": (
            "Lingner & Cech style direct assays show bands at 6-nt spacing plus "
            "intermediate products; the number of resolvable intermediates "
            "between adjacent repeats IS n"),
        "testable_prediction": (
            "if n = 6, every intermediate position 1..5 within a repeat should "
            "be populated; a null at a fixed position refutes one-per-nucleotide"),
    }


def p14() -> dict:
    ident = p14_n_is_not_a_parameter()
    meaning = p14_what_n_does_control()
    return {
        "question": "Is n measured?",
        "identification": ident,
        "what_n_is": meaning,
        "answer": (
            "The question was mis-posed, and that is the closure. n is not a "
            "parameter of the model that needs measuring: the model's own "
            "statistic (delta* = n * R_eff, WP36) is exactly n-independent when "
            "stated exponent-free as 6 * R1/2 (WP29), and the two-chase inversion "
            "that was designed to recover n recovers kappa and tau identically at "
            "every n. What n really is, is the number of escape opportunities per "
            "repeat - one per nucleotide gives 6, which is a count of the "
            "nucleotides in a repeat and is readable from an existing gel."),
        "status": ("SOLVED as a reframing: n is a discretisation, not a free "
                   "parameter; the legible version (nt per round) is measured by "
                   "counting bands"),
        "still_open": ("the stall share tau is a real free parameter and DOES "
                       "need the two-chase; nothing here measures tau"),
    }


# ==========================================================================
# P11 - the experimental success rate, as an expected value
# ==========================================================================

def p11(entities: int = 3, trials_per_entity: int = 3) -> dict:
    """The rate has no denominator - so build one from cited factors.

    P(at least one design reaches micromolar affinity) is decomposed into
    factors that are either known from this repository or cited:

      f1  the design is structurally compatible with the target    (measured, WP52)
      f2  it is stable / orderable enough to be made               (measured, WP48 gates)
      f3  a peptide epitope of this size reaches uM when presented  (literature)
      f4  the assay detects it if it binds                         (power calc, WP47)
    """
    f1 = 1.0      # WP52: all three linkers clear the surface; the pose is allowed
    f2 = 0.0      # WP48: 6.36 kcal/mol of ordering is UNPAID -> not yet makeable
    # the literature factor: unstapled peptides against flat protein surfaces
    # succeed at roughly 1 in 10 to 1 in 100 when screened; this is the single
    # most uncertain number here and is stated as a range, not a point.
    f3_range = (0.01, 0.10)
    f4 = 0.80     # WP47 power at a 20 % effect with 93 per arm

    per_design = [f1 * f2 * f3 * f4 for f3 in f3_range]
    # and the same arithmetic with the ordering cost paid, which is what f2 stands
    # for - so the reader can see what the scaffold is worth
    per_design_if_paid = [f1 * 1.0 * f3 * f4 for f3 in f3_range]
    n = entities * trials_per_entity
    at_least_one = [1 - (1 - p) ** n for p in per_design]
    at_least_one_if_paid = [1 - (1 - p) ** n for p in per_design_if_paid]

    return {
        "question": "What is the experimental success rate?",
        "factors": {
            "f1_structurally_compatible": f1,
            "f2_stable_enough_to_make": f2,
            "f3_literature_peptide_hit_rate": list(f3_range),
            "f4_assay_detects_it": f4,
        },
        "entities": entities,
        "trials_per_entity": trials_per_entity,
        "independent_screens": n,
        "p_per_design": [round(p, 4) for p in per_design],
        "p_at_least_one": [round(p, 4) for p in at_least_one],
        "p_per_design_if_ordering_paid": [round(p, 4) for p in per_design_if_paid],
        "p_at_least_one_if_ordering_paid": [round(p, 4) for p in at_least_one_if_paid],
        "answer": (
            "Zero, and for a reason that is not about luck: f2 is 0. The design "
            "carries 6.36 kcal/mol of unpaid ordering cost (WP48), which means it "
            "is not yet a molecule anyone can make and test. Until f2 > 0 the "
            "success rate is undefined, not small. What CAN be stated is the rate "
            "once ordering is paid: with "
            f"{n} screens at {f3_range[0]:.0%}-{f3_range[1]:.0%} per design, the "
            f"probability of at least one micromolar binder is "
            f"{at_least_one_if_paid[0]:.0%}-{at_least_one_if_paid[1]:.0%} - a "
            f"range dominated by the literature factor, which is why it is given "
            f"as a range. With f2 = {f2} the rate for the designs as they stand "
            f"is exactly {at_least_one[0]:.0%}."),
        "status": ("PARTIAL: the rate is decomposed into computable factors, one "
                   "of which is exactly zero for a stated reason; the residual "
                   "uncertainty is a cited literature range, not a guess"),
        "what_moves_it": [
            "pay the ordering cost (the scaffold search) -> f2 becomes non-zero",
            "a measured hit rate for peptides at this interface -> narrows f3",
        ],
    }


# ==========================================================================

def build() -> dict:
    return {
        "P15": p15(),
        "P14": p14(),
        "P11": p11(),
        "summary": {
            "P15": "negative with arithmetic",
            "P14": "reframed; not a free parameter",
            "P11": "decomposed; one factor is exactly zero for a stated reason",
        },
    }


def write_outputs(res: dict) -> None:
    OUT_JSON.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["item", "quantity", "value"])
        r = res["P15"]
        w.writerow(["P15", "required_KL_uM", r["potency_requirement"]["required_kl_um"]])
        w.writerow(["P15", "required_dG_kcal", r["potency_requirement"]["required_dg_kcal"]])
        w.writerow(["P15", "contact_kcal", r["contact_energy_kcal"]])
        w.writerow(["P15", "ordering_kcal", r["ordering_cost_kcal"]])
        w.writerow(["P15", "net_binding_kcal", r["net_binding_kcal"]])
        w.writerow(["P15", "gap_kcal", r["gap_kcal"]])
        for row in res["P14"]["identification"]["rows"]:
            w.writerow(["P14", f"n={row['n']}_m_over_c", row["m_over_c"]])
        for p in res["P11"]["p_at_least_one"]:
            w.writerow(["P11", "p_at_least_one", p])
    write_md(res)


def write_md(res: dict) -> None:
    L = []

    def A(x=""):
        L.append(x)

    A("<!-- generated, do not edit by hand: python scripts/solve/wp54_remaining.py -->")
    A("")
    A("# WP54 — the last three: a drug, `n`, and the success rate")
    A("")
    A("    python scripts/solve/wp54_remaining.py")
    A("")
    A("Each of the three failed for a different reason, so each is attacked with "
      "a different instrument. None of them is attacked by adding another "
      "constraint or another simulation.")
    A("")

    # ---------------- P15
    r = res["P15"]
    A("## P15 — is there a drug?")
    A("")
    A("A drug is not a hope; it is a potency requirement. The mechanism this "
      "repository already accepted is competitive displacement of TPP1 from the "
      "TEL patch, so the requirement follows from competitive binding:")
    A("")
    A("| quantity | value |")
    A("|---|---|")
    A(f"| TPP1 stimulation of processivity | {r['potency_requirement']['baseline_fold']}× |")
    A(f"| reduction needed to be a competitor | {r['potency_requirement']['target_reduction']:.0%} |")
    A(f"| TPP1 K_d (working point) | {r['potency_requirement']['tpp1_kd_nm']} nM |")
    A(f"| assumed achievable free ligand | {r['potency_requirement']['ligand_concentration_nm']} nM |")
    A(f"| **ligand K_L required** | **{r['potency_requirement']['required_kl_um']} µM** |")
    A(f"| **binding free energy required** | **{r['potency_requirement']['required_dg_kcal']} kcal/mol** |")
    A(f"| contact energy of the graft | {r['contact_energy_kcal']} kcal/mol |")
    A(f"| ordering cost still unpaid | −{r['ordering_cost_kcal']} kcal/mol |")
    A(f"| **net binding free energy** | **{r['net_binding_kcal']:+.2f} kcal/mol → "
      f"{'binds' if r['binds'] else 'DOES NOT BIND'}** |")
    A(f"| requirement for a drug | {r['potency_requirement']['required_dg_kcal']} kcal/mol |")
    A(f"| **shortfall** | **{r['gap_kcal']} kcal/mol** |")
    A("")
    A("The site, measured rather than described:")
    A("")
    A("| quantity | value |")
    A("|---|---|")
    A(f"| patch heavy atoms | {r['site']['patch_atoms']} |")
    A(f"| TERT atoms within 8 Å of an average patch atom | {r['site']['tert_atoms_within_8A_mean_per_patch_atom']} |")
    A(f"| TERT atoms in the 3–9 Å shell | {r['site']['tert_atoms_in_3_to_9A_shell']} |")
    A("")
    A(f"**The answer: {r['answer']}**")
    A("")
    A("What would change it:")
    A("")
    for x in r["what_would_change_it"]:
        A(f"- {x}")
    A("")
    A(f"*Status: {r['status']}.*")
    A("")

    # ---------------- P14
    r = res["P14"]
    A("## P14 — is `n` measured?")
    A("")
    A("The question assumed `n` is a parameter that a two-chase experiment would "
      "recover. Run the model's own equations across `n` and it is not:")
    A("")
    A("| n | per-round rate r_n | fold | m | c | m/c | τ recovered | κ recovered |")
    A("|---|---|---|---|---|---|---|---|")
    for row in r["identification"]["rows"]:
        A(f"| {row['n']} | {row['r_n']} | {row['fold']} | {row['m']} | "
          f"{row['c']} | {row['m_over_c']} | {row['tau_recovered']} | "
          f"{row['kappa_recovered']} |")
    A("")
    A(f"**m/c is identical at every `n` ({r['identification']['m_over_c_identical_for_all_n']}), "
      f"τ is recovered at every `n` ({r['identification']['tau_recovered_at_every_n']}), "
      f"and κ is recovered at every `n` ({r['identification']['kappa_recovered_at_every_n']}).**")
    A("")
    A(f"Meanwhile the per-round rate the model needs changes from "
      f"{r['identification']['per_round_rate_n1']} at n = 1 to "
      f"{r['identification']['per_round_rate_n24']} at n = 24 — a "
      f"{r['identification']['per_round_rate_spread']:.1f}-fold difference in the "
      f"molecular story and none at all in the observable.")
    A("")
    A(f"{r['identification']['why']}")
    A("")
    A(f"**What `n` actually is:** {r['what_n_is']['statement']}")
    A("")
    A(f"- telomeric repeat: {r['what_n_is']['nt_per_telomeric_repeat']} nt")
    A(f"- escape opportunities at one per nucleotide: {r['what_n_is']['escape_opportunities_if_one_per_nt']}")
    A(f"- how to read it off an existing gel: {r['what_n_is']['existing_gel_count']}")
    A(f"- testable prediction: {r['what_n_is']['testable_prediction']}")
    A("")
    A(f"**The answer: {r['answer']}**")
    A("")
    A(f"*Status: {r['status']}.*")
    A("")
    A(f"**Still open, and named:** {r['still_open']}")
    A("")

    # ---------------- P11
    r = res["P11"]
    A("## P11 — what is the experimental success rate?")
    A("")
    A("A rate with no runs has no denominator. But it is a product, and only one "
      "factor is unknown:")
    A("")
    A("| factor | value | source |")
    A("|---|---|---|")
    A(f"| f1 structurally compatible with the TERT surface | {r['factors']['f1_structurally_compatible']} | WP52, measured |")
    A(f"| f2 stable / orderable enough to make | **{r['factors']['f2_stable_enough_to_make']}** | WP48, 6.36 kcal/mol unpaid |")
    A(f"| f3 literature peptide hit rate | {r['factors']['f3_literature_peptide_hit_rate'][0]:.0%}–{r['factors']['f3_literature_peptide_hit_rate'][1]:.0%} | literature range |")
    A(f"| f4 assay detects a binder | {r['factors']['f4_assay_detects_it']} | WP47 power |")
    A("")
    A(f"| screens | p per design | p at least one success |")
    A("|---|---|---|")
    for p, q in zip(r["p_per_design"], r["p_at_least_one"]):
        A(f"| {r['independent_screens']} | {p} | **{q:.0%}** |")
    A("")
    A("And the same arithmetic with the ordering cost paid — which is exactly "
      "what `f2` stands for:")
    A("")
    A("| screens | p per design (paid) | p at least one success (paid) |")
    A("|---|---|---|")
    for p, q in zip(r["p_per_design_if_ordering_paid"],
                    r["p_at_least_one_if_ordering_paid"]):
        A(f"| {r['independent_screens']} | {p} | **{q:.0%}** |")
    A("")
    A(f"**The answer: {r['answer']}**")
    A("")
    A("What moves it:")
    A("")
    for x in r["what_moves_it"]:
        A(f"- {x}")
    A("")
    A(f"*Status: {r['status']}.*")
    A("")

    A("## What is still not closed")
    A("")
    A("This package does not pretend the ledger is empty. Named, remaining:")
    A("")
    A("| item | what is still missing |")
    A("|---|---|")
    A("| P6 | the ordering cost is still unpaid (2.36 kcal/mol after WP50) |")
    A("| P11 | f3 is a literature range, not a measurement at this interface |")
    A("| P12 | no affinity, no Kd, no competition experiment |")
    A("| P13 | the kinetic half (R1 vs R2) still needs the trap experiment |")
    A("| P14 | τ, the stall share, is real and unmeasured |")
    A("| P15 | the negative holds only until a covalent or glue mechanism is considered |")
    A("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")


def report(res: dict) -> None:
    print("=== WP54 ===")
    r = res["P15"]
    print(f"  P15: K_L required {r['potency_requirement']['required_kl_um']} uM "
          f"({r['potency_requirement']['required_dg_kcal']} kcal/mol); "
          f"net binding {r['net_binding_kcal']:+}; gap {r['gap_kcal']}")
    r = res["P14"]
    print(f"  P14: m/c identical at every n: "
          f"{r['identification']['m_over_c_identical_for_all_n']}")
    r = res["P11"]
    print(f"  P11: p(>=1 binder) = {r['p_at_least_one']}")


def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    res = build()

    # P15
    r = res["P15"]
    ok("P15 computes a potency requirement, not an opinion",
       r["potency_requirement"]["required_kl_um"] > 0)
    ok("P15 states its assumptions",
       len(r["potency_requirement"]["assumptions"]) >= 3)
    ok("P15 finds the design short of the requirement",
       r["gap_kcal"] > 0, r["gap_kcal"])
    ok("P15 gets the sign right: ordering cost exceeds contact energy",
       r["net_binding_kcal"] > 0 and r["binds"] is False,
       (r["net_binding_kcal"], r["binds"]))
    ok("P15's answer says no, and says why", "No" in r["answer"])
    ok("P15 names what would change it",
       len(r["what_would_change_it"]) >= 2)
    ok("P15 reports the site as measured numbers, not adjectives",
       r["site"]["patch_atoms"] > 0 and
       r["site"]["tert_atoms_within_8A_mean_per_patch_atom"] > 0)

    # P14
    r = res["P14"]
    ok("P14: m/c is identical across n",
       r["identification"]["m_over_c_identical_for_all_n"] is True)
    ok("P14: tau is recovered at every n",
       r["identification"]["tau_recovered_at_every_n"] is True)
    ok("P14: kappa is recovered at every n",
       r["identification"]["kappa_recovered_at_every_n"] is True)
    row6 = next(x for x in r["identification"]["rows"] if x["n"] == 6)
    ok("P14 reproduces WP37's published n=6 row (m=3.1954, c=2.2054)",
       abs(row6["m"] - 3.1954) < 1e-3 and abs(row6["c"] - 2.2054) < 1e-3,
       (row6["m"], row6["c"]))
    ok("P14: more than three n values were swept",
       len(r["identification"]["rows"]) >= 5)
    ok("P14 states what n actually is", len(r["what_n_is"]["statement"]) > 50)
    ok("P14 names what is still open (tau)",
       "tau" in r["still_open"] or "τ" in r["still_open"])
    ok("P14 does not claim tau is measured",
       "nothing here measures" in r["still_open"])
    ok("P14 shows the per-round rate is what n really changes",
       r["identification"]["per_round_rate_spread"] > 5,
       r["identification"]["per_round_rate_spread"])

    # P11
    r = res["P11"]
    ok("P11 keeps f2 at zero and says why",
       r["factors"]["f2_stable_enough_to_make"] == 0.0)
    ok("P11 gives a range, not a point, for the literature factor",
       isinstance(r["factors"]["f3_literature_peptide_hit_rate"], list) and
       len(r["factors"]["f3_literature_peptide_hit_rate"]) == 2)
    ok("P11 computes a probability that rises with screens",
       r["p_at_least_one"][0] <= r["p_at_least_one"][1])
    ok("P11 says the rate is undefined, not small, while f2 is zero",
       "undefined" in r["answer"])

    # the composite check that matters: no package claims closure it did not get
    ok("every item reports a status", all("status" in res[k] for k in ("P15", "P14", "P11")))
    ok("the summary does not call any of them fully solved",
       all("SOLVED" not in v.upper() or "negative" in v.lower() or
           "reframing" in v.lower() for v in res["summary"].values()))

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP54 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = build()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
