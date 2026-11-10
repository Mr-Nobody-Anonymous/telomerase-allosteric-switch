"""
WP35 - Implement P2: the DNA trap, made executable, and three of WP34's
       numbers corrected in the process.

WP34 left the programme at step 1 of its own order of operations:

    "1 | P2 DNA-trap re-analysis | days, no reagent | confirms/falsifies M5"

This work package implements that step.  Building it exposed three defects in
the step as written, and fixing them produced one result that changes the
programme's status.

WHAT WAS WRONG
--------------
(A) P2 IS NOT A RE-ANALYSIS.  It needs a new measurement.

    WP34, WP25 and WP22 all call P2 a re-analysis of existing data.  It is
    not, and this repository had already proved it could not be:

      * WP24 showed the deposited Tesmer 2023 numbers are ONE SCALAR PER
        CONDITION (fraction of signal at >=16 repeats).  WP24 built two
        distributions sharing that scalar to ten decimal places that differ
        4.98x below repeat 16.  A per-repeat survival cannot be recovered from
        a threshold fraction at any level of care.
      * `data/benchmarks/rap_independent_benchmark.csv` DOES have per-repeat
        bands, but it is generated, not measured: `generate_independent_
        benchmarks.py` calls `np.random.seed(42)` and synthesises intensities
        from chosen k_on / k_off / R_max.  "Confirming" a model on it would be
        recovering the numbers used to write it.

    So P2 is a NEW experiment: the existing pulse-chase run +/- a DNA trap.
    WP35.6 powers it and WP35.5 pre-registers both outcomes.

(B) THE STALL CORRECTION HAD THE WRONG SIGN.  This is the substantive error,
    and it was in the first draft of THIS script.

    WP34/WP25 treat stalling as inflating the un-trapped chase, which is right,
    but then apply it to the chase alone.  The product distribution is a
    stall-inflated observation too.  Since a retry means the complex keeps
    elongating, the observed per-repeat survival OVERSTATES the true one:

        s_obs  = s_true / (1 - tau (1 - s_true))      [retries inflate survival]
        s_true = s_obs (1 - tau) / (1 - tau s_obs)    [exact inverse]

    Correcting in that direction makes the distribution fold RISE with tau
    (2.7335 -> 3.27 -> 4.97 at tau = 0, 0.18, 0.50) while the trap-diluted
    chase FALLS (4.00 -> 2.00).  The two readouts CONVERGE on one tau, from
    opposite sides.  The first draft had this running forward, which made the
    fold barely move (2.7335 -> 2.9121) and produced an invented and false
    identity "delta* = 1 + tau".  The self-test caught both; the sign is the
    entire result.

(C) WP34's "88% EXPLAINED" DOES NOT SURVIVE, IN BOTH OF ITS TERMS.

    WP34 compared 1 - 2.73/4.00 = 0.318 with (0.7914 - 0.570)/0.7914 = 0.280.

      i.  0.280 is a function of (s, phi).  s and phi are the two DIFFERENT
          quantities whose mismatch is the anomaly M5 exists to explain, so
          (s0-phi0)/s0 IS that anomaly rescaled.  It cannot explain itself.
          In the +PT arm the same expression is -0.0796, and a fraction of a
          population cannot be negative.
      ii. 0.318 is not the physical retry probability.  tau is defined by the
          renewal model above, and the crossing is at tau_0 = 0.1833.

    What replaces it is stronger: stalling RECONCILES the two observations at
    one tau with no freedom left over, and that tau is measurable.

THE ASSUMPTION THAT DECIDES THE SIGN - derived, not chosen
----------------------------------------------------------
A retry can only follow a FAILURE, and failures per repeat are (1-phi):
0.430 without POT1-TPP1 and 0.010 with it.  Saturating POT1-TPP1 therefore
leaves almost nothing to retry, and the retry probability should scale with the
failure rate:

    tau_1 = tau_0 * (1-phi_PT)/(1-phi_noPT) = tau_0 * 0.0233    [default]

Three couplings are solved and tabulated so no single one is reported alone:

    coupling          tau_0    fold     delta*      gate vs 50 bp/PD
    failure-scaled    0.1833   3.2670   57.30       OPEN
    none              0.1815   3.2741   57.43       OPEN
    same tau          0.3345   2.6621   46.76       CLOSED

"same tau in both arms" is the only variant that keeps the gate shut, and the
phi data contradict it.  But note what that also means: the spread from 46.8 to
57.4 bp/PD is driven by the coupling choice, not by data.  Reporting one number
without the other two would be the same error WP33 caught in WP32.

THE RESULT THAT CHANGES THE PROGRAMME'S STATUS
----------------------------------------------
WP34 concluded: maintenance gate CLOSED at 48.0 < 50 bp/PD, therefore
P(approved drug) = 0%, "because the gate is closed, not merely improbable".

That is a conditional presented as a closure.  48.0 bp/PD is delta* at
tau_0 = 0, and tau_0 = 0 is exactly the value P2 exists to test.  Since the two
physically derived couplings both put delta* at ~57 bp/PD - above the 50 floor -
the honest status is:

    the gate is UNDETERMINED, and one measurement decides it.

    tau_0 = 0        -> 48.0 bp/PD -> CLOSED
    tau_0 = 0.18     -> 57.3 bp/PD -> matches the 50 floor, still short of 90
    fold >= 5.14x    -> 90 bp/PD   -> needs tau_0 >= 0.51, beyond the crossing

P(approved drug) is therefore not 0%.  It is WP34's own pipeline product,
conditional on the gate opening - small, and not a finding of zero.

WHAT P2 COSTS - and why the word "days" was wrong
-------------------------------------------------
The primary contrast is 3.27x against 4.00x, a ratio of 1.22, measured against
the only experiment-level CV this repository holds from a real gel (17.2%,
Tesmer Fig 3B WT):

    independent arms, rho = 0    n = 12 per arm
    paired within reaction, rho = 0.8    n = 3 per arm

L&C report 3 per arm.  So P2 is one extra paired experiment or a campaign of
twelve, and which one it is depends on a DESIGN choice, not on anything in this
script.  Quoting "days, no reagent" without rho = 0.8 was not supportable.

    python scripts/solve/wp35_p2_identifiability.py
    -> results/solve/wp35_p2.csv
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "solve"

# --- anchors, identical to WP22 / WP29 / WP34 so the four are comparable ----
S_NO_PT = 0.791412          # L&C 2010, R1/2 = 2.963 without POT1-TPP1
S_PT = 0.917004             # L&C 2010, R1/2 = 8.0 with POT1-TPP1
PHI_NO_PT = 0.570           # L&C window-limited translocation efficiency
PHI_PT = 0.990
KAPPA_CHASE = 4.00          # L&C pulse-chase, a LOWER BOUND ("at least 4-fold")
N_INCORP = 6                # nucleotides per repeat

# --- Wang 2007, the independent less-processive enzyme ----------------------
S_WANG_0 = 2.0 ** (-1.0 / 0.78)
S_WANG_1 = 2.0 ** (-1.0 / 3.30)

# --- real assay noise: the one MEASURED experiment-level error in the repo ---
# data/external/tesmer2023/fig3b_processivity_plusPT.csv, WT row, 3 independent
# experiments: 0.508790 / 0.554671 / 0.392860 -> CV 17.2%.
# (data/benchmarks/*.csv is generated by np.random.seed(42) and must not be
# used for a power calculation about a real gel.)
FIG3B_WT = (0.508790, 0.554671, 0.392860)

DELTA_FLOOR = 50.0          # lowest measured human fibroblast attrition, bp/PD
DELTA_HUVEC = 90.0          # measured HUVEC attrition, bp/PD


# ---------------------------------------------------------------- core model
def a_of_s(s: float, n: int = N_INCORP) -> float:
    """Off-rate-to-polymerase ratio implied by per-repeat survival: s=(1+a)^-n."""
    return s ** (-1.0 / n) - 1.0


def s_of_a(a: float, n: int = N_INCORP) -> float:
    return (1.0 + a) ** (-n)


def s_obs_given_true(s_true: float, tau: float) -> float:
    """Observed per-repeat survival when a failed step retries with prob tau.

    Renewal: success outright with s_true, or failure then a retry with prob
    tau.  Summing the geometric series over retries:
        s_obs = s_true / (1 - tau (1 - s_true))
    """
    d = 1.0 - tau * (1.0 - s_true)
    if d <= 0.0:
        raise ZeroDivisionError("tau=1: no release at all; survival is undefined")
    return s_true / d


def s_true_given_obs(s_obs: float, tau: float) -> float:
    """Exact inverse of s_obs_given_true.  This is the direction P2 needs.

    Solve s_obs (1 - tau + tau s_true) = s_true:
        s_true = s_obs (1 - tau) / (1 - tau s_obs)
    At tau=0 it returns s_obs unchanged, as it must.
    """
    return s_obs * (1.0 - tau) / (1.0 - tau * s_obs)


#: How the retry probability in the +POT1-TPP1 arm relates to the no-PT arm.
#: This choice decides the sign of the answer, so it is made explicit here and
#: every variant is reported in main().  It is NOT a free parameter to pick the
#: nicest answer from - WP35.2 derives the default from phi.
#:
#:   "failure-scaled"  tau_1 = tau_0 * (1-phi_1)/(1-phi_0)     <- default
#:   "none"            tau_1 = 0                                <- M5 as written
#:   "same"            tau_1 = tau_0                            <- strongest correction
COUPLINGS = ("failure-scaled", "none", "same")

#: Ratio of per-repeat FAILURE probabilities, (1-phi_PT)/(1-phi_noPT).
#: L&C: phi 0.57 -> 0.99, so this is 0.010/0.430 = 0.0233.  A retry can only
#: happen after a failure, so a retry probability should scale with the failure
#: rate; this factor is what makes "stalling matters" and "stalling is gone"
#: two different physical statements rather than two fitted options.
FAILURE_RATIO = (1.0 - PHI_PT) / (1.0 - PHI_NO_PT)


def tau_pt(tau: float, coupling: str = "failure-scaled") -> float:
    """Retry probability in the +PT arm, given the no-PT value, under `coupling`."""
    if coupling == "failure-scaled":
        return tau * FAILURE_RATIO
    if coupling == "none":
        return 0.0
    if coupling == "same":
        return tau
    raise ValueError(f"unknown coupling {coupling!r}; expected one of {COUPLINGS}")


def kappa_dist(tau: float, coupling: str = "failure-scaled",
               s_obs0: float = S_NO_PT, s1: float = S_PT) -> float:
    """Off-rate fold from the product distribution, with stall inflation removed.

    Both arms are stall-corrected, but the +PT correction is governed by
    `coupling` because the two arms need not share a retry probability.

    This is the function whose SIGN was wrong in the first draft.  It must
    increase with tau for the reconciliation to work: a retrying complex keeps
    elongating, so the observed survival overstates the true one, and removing
    the inflation with tau=0 leaves the largest off-rate in the numerator.
    """
    return (a_of_s(s_true_given_obs(s_obs0, tau))
            / a_of_s(s_true_given_obs(s1, tau_pt(tau, coupling))))


def kappa_chase(tau: float, measured: float = KAPPA_CHASE) -> float:
    """Measured chase fold after a DNA trap removes retrying complexes.

    A trapped complex that would have retried is lost instead, so the observed
    bound population decays faster by exactly the share that stalling had been
    sustaining: kappa_chase = measured (1 - tau).
    """
    return measured * (1.0 - tau)


def solve_selfconsistent(measured: float = KAPPA_CHASE,
                         coupling: str = "failure-scaled",
                         tau_max: float = 0.60) -> tuple[float, float]:
    """The tau at which the corrected distribution and the trap chase agree.

    kappa_chase always falls with tau.  kappa_dist rises with tau under every
    coupling, but by a very different amount - 6.3x (failure-scaled and none)
    against 1.03x (same).  So on (0, tau_max) the difference has exactly one
    root when the endpoints straddle zero, and under "same" the root sits near
    the far end.  Bisection is safe and needs no scipy.
    """
    f = lambda t: kappa_dist(t, coupling) - kappa_chase(t, measured)
    lo, hi = 0.0, tau_max
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError(f"no crossing on (0,{tau_max}) for {coupling!r}: "
                         f"f(0)={flo:+.4f} f({tau_max})={fhi:+.4f}")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    tau = 0.5 * (lo + hi)
    return tau, kappa_dist(tau, coupling)


def delta_star(kappa: float, a0: float = None) -> float:
    """Critical attrition, bp per population doubling, at full coverage (WP29)."""
    a0 = a_of_s(S_NO_PT) if a0 is None else a0
    return 1.0 / math.log2(1.0 + a0 / kappa)


def kappa_needed(delta: float, a0: float = None) -> float:
    """Off-rate fold required to maintain a cell losing `delta` bp/PD (WP29)."""
    a0 = a_of_s(S_NO_PT) if a0 is None else a0
    return a0 / (2.0 ** (1.0 / delta) - 1.0)


def tau_for_delta(delta: float, tau_max: float = 0.60) -> float | None:
    """Stall share at which delta* reaches `delta`, or None if out of reach."""
    if delta <= delta_star(kappa_dist(0.0)):
        return 0.0
    if delta >= delta_star(kappa_dist(tau_max)):
        return None
    lo, hi = 0.0, tau_max
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if delta_star(kappa_dist(mid)) < delta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------- the experiment
def assay_cv(values: tuple[float, ...] = FIG3B_WT) -> float:
    """Between-experiment CV of a real gel readout (Tesmer Fig 3B, WT)."""
    n = len(values)
    m = sum(values) / n
    var = sum((v - m) ** 2 for v in values) / (n - 1)
    return math.sqrt(var) / m


def n_per_arm(cv: float, ratio: float, rho: float = 0.0, power: float = 0.80) -> int:
    """Independent experiments per arm to resolve a fold of `ratio` on a log scale.

    For a ratio of two means each measured with CV c and n per arm,
        SE(ln ratio) ~ sqrt(2(1-rho)) c / sqrt(n)
    where rho is the correlation between the paired measurements (rho=0 is two
    independent arms; rho=1 is a within-experiment paired comparison with no
    extra variance).  Requiring SE <= ln(ratio)/(z_{a/2}+z_power) gives
        n >= 2(1-rho) c^2 (z_{a/2}+z_power)^2 / ln(ratio)^2
    z values are hard-coded (two-sided alpha 0.05, 80 % power) to keep this
    script free of scipy.
    """
    z = 1.959964 + 0.841621
    return max(2, math.ceil(2.0 * (1.0 - rho) * cv ** 2 * z ** 2 / math.log(ratio) ** 2))


def cv_tolerated(ratio: float, n: int, rho: float = 0.0) -> float:
    """Largest between-experiment CV that `n` per arm can still resolve `ratio`."""
    z = 1.959964 + 0.841621
    return math.log(ratio) / z * math.sqrt(n / (2.0 * (1.0 - rho)))


def pred_interval(cv: float, n: int, centre: float) -> tuple[float, float]:
    """95% prediction interval on a fold measured at n per arm."""
    half = 1.959964 * math.sqrt(2.0) * cv / math.sqrt(n)
    return centre * math.exp(-half), centre * math.exp(+half)


def ladder(s: float, n_bands: int = 8) -> list[float]:
    """SIMULATED per-band fractions for a geometric ladder with survival s.

    Used only to show that a per-band ladder pins s and a threshold fraction
    does not.  Never reported as data.
    """
    return [s ** (i - 1) * (1.0 - s) for i in range(1, n_bands + 1)]


def threshold_fraction(s: float, at: int = 16) -> float:
    """P(L >= at) for a geometric ladder - the ONLY thing Tesmer 2023 deposited."""
    return s ** (at - 1)


def s_from_ladder_demo(cv: float, n: int, true_s: float = S_NO_PT,
                       seed: int = 20261105) -> float:
    """Recover s from a noisy per-band ladder - the quantity Gate 3 needs.

    Deterministic pseudo-noise (a fixed low-discrepancy sequence) rather than a
    random draw, so the doc's number is reproducible without a seed.
    """
    bands = ladder(true_s, 12)
    noisy = []
    for i, b in enumerate(bands):
        z = math.sin((i + 1) * 12.9898 + seed) * 43758.5453
        z = (z - math.floor(z)) * 2.0 - 1.0
        noisy.append(max(b * (1.0 + cv * z / math.sqrt(n)), 1e-12))
    xs = list(range(1, len(noisy) + 1))
    mx = sum(xs) / len(xs)
    my = sum(math.log(v) for v in noisy) / len(noisy)
    num = sum((x - mx) * (math.log(v) - my) for x, v in zip(xs, noisy))
    den = sum((x - mx) ** 2 for x in xs)
    return math.exp(num / den)


# ---------------------------------------------------------------- reporting
def main() -> int:
    L = "=" * 92
    print(L)
    print("WP35  P2 implemented: the sign that was wrong, and the gate it moves")
    print(L)
    print()

    a0, a1 = a_of_s(S_NO_PT), a_of_s(S_PT)
    tau, kappa = solve_selfconsistent()
    tau_max = 0.60
    # every coupling solved and tabulated - the sign of the conclusion depends
    # on this choice, so no single one is reported without the other two.
    grid = {c: solve_selfconsistent(coupling=c) for c in COUPLINGS}

    # ---------------------------------------------------------------- 35.1
    print("WP35.1  CORRECTION 1 - P2 is not a re-analysis, and this repo proved it")
    print()
    print("  WP24 already showed the deposited numbers are one scalar per condition:")
    print("    two product-length distributions, identical at >=16 repeats to 10 dp,")
    print("    differ 4.98x in the fraction ending at <=3 repeats.")
    print("  The one per-band dataset here is generated, not measured:")
    print("    scripts/generate_independent_benchmarks.py -> np.random.seed(42),")
    print("    intensities synthesised from chosen k_on / k_off / R_max.")
    print()
    print("  => P2 requires a NEW run: the existing pulse-chase +/- a DNA trap.")
    print("     One experiment per arm beyond L&C's three.  Specified in WP35.6.")
    print()

    # ---------------------------------------------------------------- 35.2
    print("WP35.2  CORRECTION 2 - the stall correction had the wrong sign")
    print()
    print("  WP34/WP25 applied stalling to the chase only.  Both observations are")
    print("  stall-inflated, so both must be put on the same tau:")
    print()
    print("    s_obs = s_true / (1 - tau(1 - s_true))       [retries inflate survival]")
    print("    s_true = s_obs (1 - tau) / (1 - tau s_obs)    [exact inverse]")
    print("    kappa_dist(tau) = a(s_true0) / a(s_true1)     [distribution fold]")
    print("    kappa_chase(tau) = 4.00 (1 - tau)             [trap dilutes the chase]")
    print()
    print("  The first draft of this script ran the correction FORWARD instead of")
    print("  inverted, which gave a fold that barely moved (2.7335 -> 2.9121) and")
    print("  an invented 'delta* = 1 + tau' identity that was simply false.  The")
    print("  self-test caught both.  The sign is the entire result.")
    print()
    print("  THE ASSUMPTION THAT DECIDES THE SIGN - and why it is derived, not chosen")
    print()
    print(f"  A retry can only follow a FAILURE, and failures per repeat are")
    print(f"  (1-phi): {1-PHI_NO_PT:.3f} without POT1-TPP1, {1-PHI_PT:.3f} with it.")
    print(f"  That is a ratio of {FAILURE_RATIO:.4f} - saturating POT1-TPP1 leaves")
    print("  almost nothing to retry.  So the default coupling is")
    print()
    print(f"      tau_1 = tau_0 x {FAILURE_RATIO:.4f}")
    print()
    print("  and 'same tau in both arms' is a strictly stronger - and physically")
    print("  unsupported - assumption.  Both are reported below, so the reader")
    print("  does not have to trust the default.")
    print()
    print(f"  {'tau_0':>6} {'s_true0':>8} {'s_true1':>8} {'kappa_dist':>11} {'kappa_chase':>12} {'delta*':>8}")
    print("  " + "-" * 58)
    for t in (0.0, 0.05, 0.10, 0.15, tau, 0.25, 0.30, 0.40, 0.50):
        d = delta_star(kappa_dist(t))
        print(f"  {t:6.3f} {s_true_given_obs(S_NO_PT,t):8.4f} "
              f"{s_true_given_obs(S_PT, tau_pt(t)):8.4f} {kappa_dist(t):11.4f} "
              f"{kappa_chase(t):12.4f} {d:8.2f}")
    print()
    print("  kappa_dist RISES with tau and kappa_chase FALLS: the two readouts")
    print("  converge from opposite sides on one tau, which is why P2 is")
    print("  identifiable at all.")
    print()

    # ---------------------------------------------------------------- 35.3
    print("WP35.3  CORRECTION 3 - WP34's '88%' does not survive either")
    print()
    p34 = (S_NO_PT - PHI_NO_PT) / S_NO_PT
    p34_pt = (S_PT - PHI_PT) / S_PT
    print("  WP34 compared 1 - 2.73/4.00 = 0.318 with (s0-phi0)/s0 = 0.280.")
    print()
    print(f"    0.280 is a function of (s, phi).  s and phi are the two DIFFERENT")
    print(f"    quantities whose mismatch is the anomaly M5 exists to explain.")
    print(f"    (s0-phi0)/s0 IS that anomaly, rescaled - so it cannot explain it.")
    print(f"    In the +PT arm the same expression gives {p34_pt:+.4f}, and a")
    print(f"    fraction of a population cannot be negative.")
    print()
    print(f"    0.318 is not the physical tau either.  The physical tau is defined")
    print(f"    by the renewal model, and the crossing is at tau = {tau:.4f}.")
    print()
    print("  What replaces it is stronger, not weaker: stalling RECONCILES the")
    print("  two observations at one tau with no freedom left, and that tau is")
    print("  measurable.  'Reconciled at a measurable tau' beats '88% explained")
    print("  by a quantity that is the gap itself'.")
    print()

    # ---------------------------------------------------------------- 35.4
    print("WP35.4  THE RESULT THAT CHANGES THE PROGRAMME - the gate moves")
    print()
    print(f"  {'coupling':16} {'tau_0':>7} {'kappa':>8} {'delta*':>9} {'vs 50':>8} {'vs 90':>8}")
    print("  " + "-" * 62)
    for c in COUPLINGS:
        tc, kc_ = grid[c]
        d = delta_star(kc_)
        print(f"  {c:16} {tc:7.4f} {kc_:8.4f} {d:9.2f} "
              f"{'OPEN' if d > DELTA_FLOOR else 'closed':>8} "
              f"{'OPEN' if d > DELTA_HUVEC else 'closed':>8}")
    print()
    print(f"  kappa needed for the 50 bp/PD floor   {kappa_needed(DELTA_FLOOR):.4f}x")
    print(f"  kappa needed for a 90 bp/PD cell      {kappa_needed(DELTA_HUVEC):.4f}x")
    print()
    print("  The two derived couplings agree: correcting only where failures")
    print("  are common gives 57.4 bp/PD, correcting nowhere gives 48.0, and the")
    print("  failure-scaled correction lands on 57.3 because there is almost")
    print(f"  nothing to retry in the +PT arm ({FAILURE_RATIO:.3f} of the failures).")
    print("  Only the 'same tau in both arms' assumption - which the phi data")
    print("  contradict - keeps the gate shut.")
    print()
    print("  So the honest statement is NOT 'the gate is open'.  It is:")
    print()
    print("    the gate's status is set by tau_0, and the two physically derived")
    print("    couplings both put delta* at ~57 bp/PD, above the 50 floor.")
    print("    WP34's 48.0 and its 'P(approved drug) = 0%, because the gate is")
    print("    closed' both rest on M5 as written, i.e. on tau_0 = 0, which is")
    print("    the one value P2 exists to test.  That is a conditional, not a")
    print("    closure, and WP35 replaces it with the conditional.")
    print()

    # ---------------------------------------------------------------- 35.5
    print("WP35.5  FALSIFIABLE IN BOTH DIRECTIONS - what makes it a real test")
    print()
    print(f"  {'trap chase returns':>22} {'tau':>8} {'true fold':>10} {'delta*':>8}  verdict")
    print("  " + "-" * 78)
    for kc in (4.00, 3.60, 3.40, 3.27, 3.00, 2.80):
        tt = 1.0 - kc / KAPPA_CHASE
        try:
            kd = kappa_dist(tt)
            d = delta_star(kd)
            if kc >= 3.99:
                v = "tau=0: M5 unchanged, gate CLOSED"
            elif d > DELTA_FLOOR:
                v = "gate OPEN, maintenance plausible"
            else:
                v = "gate closed, stalling insufficient"
            print(f"  {kc:22.2f} {tt:8.3f} {kd:10.4f} {d:8.2f}  {v}")
        except ZeroDivisionError:
            print(f"  {kc:22.2f} {tt:8.3f} {'degenerate':>10}")
    print()
    print("  Pre-registered: the trap chase should return 3.27x; >= 3.99x means")
    print("  no stall contribution and M5 keeps 2.73x with the gate closed.")
    print("  Both outcomes are publishable and neither is a rescue.")
    print()

    # ---------------------------------------------------------------- 35.6
    print("WP35.6  POWERING P2 from the only real experiment-level noise here")
    print()
    cv = assay_cv()
    print(f"  Tesmer 2023 Fig 3B, WT, 3 independent experiments: {FIG3B_WT}")
    print(f"  mean {sum(FIG3B_WT)/3:.4f}   CV {cv:.3f}")
    print()
    print(f"  {'contrast':38} {'ratio':>8} {'n/arm':>7}")
    print("  " + "-" * 56)
    for label, r in (("3.27x prediction vs 4.00x bound", KAPPA_CHASE / kappa),
                     ("2.73x (tau=0) vs 4.00x bound", KAPPA_CHASE / (a0 / a1)),
                     ("3.27x vs 2.73x  (direct tau estimate)", kappa / (a0 / a1))):
        print(f"  {label:38} {r:8.4f} {n_per_arm(cv, r):7d}")
    n = n_per_arm(cv, KAPPA_CHASE / kappa)
    lo, hi = pred_interval(cv, n, kappa)
    print()
    print(f"  Independent arms at this CV need n = {n} per arm.  That is not")
    print("  'days, no reagent' - it is a campaign, and it is the honest cost of")
    print(f"  a {KAPPA_CHASE/kappa:.2f}x contrast against a {cv:.3f} CV.")
    print()
    print("  A PAIRED design is much cheaper, because the trap is applied to the")
    print("  same reaction.  With correlation rho between paired readouts:")
    print()
    print(f"  {'rho':>6} {'n per arm':>11} {'CV tolerated at that n':>24}")
    print("  " + "-" * 46)
    for rho in (0.0, 0.5, 0.8, 0.9):
        nn = n_per_arm(cv, KAPPA_CHASE / kappa, rho)
        print(f"  {rho:6.1f} {nn:11d} {cv_tolerated(KAPPA_CHASE/kappa, nn, rho):24.4f}")
    print()
    n_paired = n_per_arm(cv, KAPPA_CHASE / kappa, 0.8)
    lo, hi = pred_interval(cv, n_paired, kappa)
    print(f"  At rho = 0.8 (a realistic within-reaction pairing) n = {n_paired} per arm,")
    print(f"  and L&C already have 3.  The requirement is then ONE more.")
    print(f"  95% prediction interval at n={n_paired}: [{lo:.2f}, {hi:.2f}]x")
    print()
    print("  This is the single most decision-relevant number in WP35: whether P2")
    print("  is cheap or expensive depends entirely on pairing, which is a design")
    print("  choice, not a discovery.  Do not quote 'days' without rho = 0.8.")
    print()

    ts, ks = grid["same"]

    # ---------------------------------------------------------------- 35.7
    print("WP35.7  WHAT THIS CHANGES IN WP34")
    print()
    print(f"  {'WP34 said':34}  WP35 says")
    print("  " + "-" * 86)
    print(f"  {'P2 = re-analysis, days, no reagent':34}  P2 = chase +/- trap; n=12, or n=3 paired")
    print(f"  {'stalling inflates the chase only':34}  stalling inflates BOTH readouts")
    print(f"  {'4x gap 88% explained by 0.280':34}  that quantity IS the gap, rescaled")
    print(f"  {'retry probability needed ~0.32':34}  tau_0 = {tau:.4f}; 0.318 was not a tau")
    print(f"  {'M5 prediction stays 2.73x':34}  fold {kappa:.4f}x, or {ks:.4f}x under 'same'")
    print(f"  {'gate CLOSED at 48.0 < 50':34}  conditional: {delta_star(kappa):.1f} / {delta_star(ks):.1f} bp/PD")
    print(f"  {'P(approved drug) = 0%':34}  0% is not a finding, it is the assumption tau_0=0")
    print(f"  {'P(M5) ~35%':34}  ~35%; weakest term is now (A), not (C)")
    print()

    rows = [
        dict(item="P2 is a re-analysis", value=0, note="FALSE - needs chase +/- trap"),
        dict(item="why re-analysis cannot work", value="threshold fraction", note="WP24 proof; >=16 is one scalar"),
        dict(item="per-band dataset in repo", value=0, note="generated (seed 42), not measured"),
        dict(item="kappa at tau=0", value=round(a0 / a1, 4), note="M5 as written"),
        dict(item="physical tau at the crossing", value=round(tau, 4), note="WP34 said 0.318 (non-physical)"),
        dict(item="tau_0 with failure-scaled coupling", value=round(grid["failure-scaled"][0], 4), note="default; tau_1 = tau_0 x 0.0233"),
        dict(item="tau_0 with no +PT correction", value=round(grid["none"][0], 4), note="M5 as written downstream"),
        dict(item="tau_0 with same tau both arms", value=round(ts, 4), note="strongest, and contradicted by phi"),
        dict(item="failure ratio (1-phi_PT)/(1-phi_noPT)", value=round(FAILURE_RATIO, 4), note="why tau_1 ~ 0"),
        dict(item="reconciled true fold", value=round(kappa, 4), note="failure-scaled coupling"),
        dict(item="reconciled fold, same-tau coupling", value=round(ks, 4), note="the pessimistic branch"),
        dict(item="delta* at tau=0", value=round(delta_star(a0 / a1), 2), note="WP29 value, CLOSED vs 50"),
        dict(item="delta* at reconciled tau", value=round(delta_star(kappa), 2), note="failure-scaled: OPEN vs 50"),
        dict(item="delta* under same-tau coupling", value=round(delta_star(ks), 2), note="pessimistic: CLOSED vs 50"),
        dict(item="kappa for 50 bp/PD", value=round(kappa_needed(DELTA_FLOOR), 4), note="below reconciled fold"),
        dict(item="kappa for 90 bp/PD", value=round(kappa_needed(DELTA_HUVEC), 4), note="needs tau ~0.42"),
        dict(item="maintenance gate", value="UNDETERMINED", note="decided by P2, not by a bound"),
        dict(item="WP34 '88% explained'", value=0, note="quantity is the anomaly itself, rescaled"),
        dict(item="(s0-phi0)/s0 in +PT arm", value=round(p34_pt, 4), note="negative: not a fraction"),
        dict(item="assay CV (Tesmer Fig 3B WT)", value=round(cv, 4), note="only measured noise in repo"),
        dict(item="n per arm, independent", value=n, note="at CV 0.172, rho=0"),
        dict(item="n per arm, paired rho=0.8", value=n_paired, note="L&C have 3 -> one more"),
        dict(item="pre-registered trap fold", value=round(kappa, 4), note="falsified if >= 3.99x"),
        dict(item="P(M5 correct)", value=0.35, note="0.90 x 0.65 x 0.60, priors, arguable"),
        dict(item="therapeutic claim", value="none", note="Quarantine unchanged"),
    ]
    dest = OUT / "wp35_p2.csv"
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

    a0, a1 = a_of_s(S_NO_PT), a_of_s(S_PT)
    tau, kappa = solve_selfconsistent()

    # --- anchors round-trip exactly -------------------------------------
    chk("s -> a -> s round-trips", abs(s_of_a(a0) - S_NO_PT) < 1e-12)
    chk("WP29 delta* is reproduced at tau=0", abs(delta_star(kappa_dist(0.0)) - 48.0) < 0.02)
    chk("WP22 fold is reproduced at tau=0", abs(kappa_dist(0.0) - 2.7335) < 5e-4)

    # --- the inverse correction, which is where the sign lives ----------
    chk("tau=0 leaves the observed survival unchanged",
        abs(s_true_given_obs(S_NO_PT, 0.0) - S_NO_PT) < 1e-15)
    chk("true survival is always below the observed one",
        all(s_true_given_obs(S_NO_PT, t) < S_NO_PT for t in (0.05, 0.2, 0.4)))
    # round-trip: squeeze the true value through the forward map, get it back
    for t in (0.05, 0.1815, 0.35):
        st = s_true_given_obs(S_NO_PT, t)
        chk(f"forward/inverse round-trips at tau={t}", abs(s_obs_given_true(st, t) - S_NO_PT) < 1e-12)

    # --- the direction of the correction (the bug the self-test caught) --
    grid = [i / 400 for i in range(240)]           # tau in [0, 0.6)
    diffs = [kappa_dist(grid[i + 1]) - kappa_dist(grid[i]) for i in range(len(grid) - 1)]
    chk("kappa_dist INCREASES with tau (correct sign)", all(d > 0 for d in diffs))
    chk("the fold moves a lot, not a little",
        kappa_dist(0.30) / kappa_dist(0.0) > 1.25)
    chk("the first draft's flat curve is NOT what this computes",
        abs(kappa_dist(0.30) - 2.9121) > 0.5)
    chk("kappa_chase falls by exactly 4.0 per unit tau",
        abs((kappa_chase(0.0) - kappa_chase(1.0)) - KAPPA_CHASE) < 1e-12)

    # --- the crossing ---------------------------------------------------
    chk("crossing tau_0 is about 0.1833", abs(tau - 0.1833) < 2e-3)
    chk("reconciled fold is about 3.267", abs(kappa - 3.2670) < 5e-3)
    chk("the two readouts agree at the crossing",
        abs(kappa_dist(tau) - kappa_chase(tau)) < 1e-9)
    chk("the crossing is unique on (0, 0.6)", all(
        (kappa_dist(t) - kappa_chase(t)) * (t - tau) >= -1e-9
        for t in (i / 200 for i in range(120))))

    # --- WP34's two numbers, both shown to be wrong on their own terms ---
    chk("WP34's 0.318 is not the physical tau", abs(tau - 0.318) > 0.05)
    p34 = (S_NO_PT - PHI_NO_PT) / S_NO_PT
    p34_pt = (S_PT - PHI_PT) / S_PT
    chk("WP34's 0.280 is reproducible", abs(p34 - 0.2798) < 5e-4)
    chk("that expression is NEGATIVE in the +PT arm", p34_pt < 0.0)
    chk("a fraction of a population cannot be negative -> not a fraction", p34_pt < 0.0)

    # --- the gate -------------------------------------------------------
    chk("delta* at tau=0 is below the 50 floor", delta_star(kappa_dist(0.0)) < DELTA_FLOOR)
    chk("delta* at the crossing is ABOVE the 50 floor", delta_star(kappa) > DELTA_FLOOR)
    chk("delta* at the crossing is still below 90", delta_star(kappa) < DELTA_HUVEC)
    chk("the 50 floor needs a smaller fold than the crossing gives",
        kappa_needed(DELTA_FLOOR) < kappa)
    chk("the 90 bp/PD cell needs 5.1428x", abs(kappa_needed(DELTA_HUVEC) - 5.1428) < 1e-3)
    chk("the 90 cell needs a tau above the crossing", tau_for_delta(DELTA_HUVEC) is None
        or tau_for_delta(DELTA_HUVEC) > tau)

    # --- powering -------------------------------------------------------
    cv = assay_cv()
    chk("the real assay CV is about 0.172", abs(cv - 0.172) < 2e-3)
    r = KAPPA_CHASE / kappa
    n_ind = n_per_arm(cv, r)
    n_paired = n_per_arm(cv, r, 0.8)
    chk("independent arms need 12 per arm at rho=0", n_ind == 12)
    chk("a paired design at rho=0.8 needs 3", n_paired == 3)
    chk("pairing is never worse than independence", n_paired <= n_ind)
    lo, hi = pred_interval(cv, n_ind, kappa)
    chk("the prediction interval sits below 4.00x", hi < KAPPA_CHASE)
    chk("the prediction interval excludes no-stall", lo > kappa_dist(0.0))
    chk("a tighter CV cannot need MORE replicates", n_per_arm(cv / 2, r) <= n_ind)
    # n_per_arm ceils, so the exact inverse statement is that n_ind is the
    # SMALLEST arm size whose tolerated CV still covers the observed one.
    chk("cv_tolerated(n) covers the observed CV", cv_tolerated(r, n_ind) >= cv)
    chk("cv_tolerated(n-1) does not", cv_tolerated(r, n_ind - 1) < cv)
    chk("tau itself is only weakly determined by the primary contrast",
        (kappa / kappa_dist(0.0) - 1) < 2.0 * 1.959964 * math.sqrt(2) * cv / math.sqrt(n_ind))

    # --- the ladder demo ------------------------------------------------
    chk("ladder fractions are a proper distribution",
        abs(sum(ladder(0.8, 12)) - (1 - 0.8 ** 12)) < 1e-9)
    chk("the threshold fraction is one number per condition",
        abs(threshold_fraction(0.5, 16) - 0.5 ** 15) < 1e-15)
    chk("a noisy ladder recovers s to within 5% at n=3",
        abs(s_from_ladder_demo(cv, 3) / S_NO_PT - 1.0) < 0.05)
    chk("a ladder at n=1 is still within 10%",
        abs(s_from_ladder_demo(cv, 1) / S_NO_PT - 1.0) < 0.10)

    # --- Wang cross-check ------------------------------------------------
    chk("the less-processive enzyme needs a bigger stall to reconcile",
        (a_of_s(S_WANG_0) / a_of_s(S_WANG_1)) > (a0 / a1))

    # --- the coupling assumption, which decides the sign ---------------
    chk("failure ratio is (1-phi_p)/(1-phi_0) = 0.0233", abs(FAILURE_RATIO - 0.02326) < 1e-4)
    chk("a retry is 43x rarer with POT1-TPP1 present", 1.0 / FAILURE_RATIO > 40.0)
    chk("tau_pt is below tau_0 under failure-scaling",
        tau_pt(0.3, "failure-scaled") < 0.3)
    chk("tau_pt is zero under 'none'", tau_pt(0.3, "none") == 0.0)
    chk("tau_pt is identity under 'same'", abs(tau_pt(0.3, "same") - 0.3) < 1e-15)
    try:
        tau_pt(0.3, "nonsense")
        chk("an unknown coupling raises", False)
    except ValueError:
        chk("an unknown coupling raises", True)

    # every coupling crosses, and each gives a different answer
    solved = {c: solve_selfconsistent(coupling=c) for c in COUPLINGS}
    chk("all three couplings reconcile", len(solved) == 3)
    chk("failure-scaled and 'none' agree closely",
        abs(solved["failure-scaled"][1] - solved["none"][1]) < 0.05)
    chk("'same' differs materially from the derived couplings",
        abs(solved["same"][1] - solved["failure-scaled"][1]) > 0.5)
    chk("'same' is the pessimistic branch (lowers the fold)",
        solved["same"][1] < solved["failure-scaled"][1])

    # --- the guard against silently reverting to a one-arm correction ---
    # a fold that rises by only ~6% over the whole tau range is the bug in
    # the first draft; a correct correction raises it by >50%.
    chk("kappa_dist moves by more than 30% across tau (catches the old bug)",
        kappa_dist(0.30) / kappa_dist(0.0) > 1.30)
    chk("kappa_dist('same') barely moves - which is why it is the weak branch",
        kappa_dist(0.30, "same") / kappa_dist(0.0, "same") < 1.10)
    chk("correcting nothing equals correcting only where failures are common",
        abs(solved["failure-scaled"][1] - solved["none"][1]) < 0.05)

    print()
    bad = [nm for nm, ok in checks if not ok]
    print(f"WP35 self-test: {len(checks) - len(bad)} passed, {len(bad)} failed")
    for nm in bad:
        print(f"  FAILED: {nm}")
    return 1 if bad else 0


if __name__ == "__main__":
    import sys

    if "--self-test" in sys.argv:
        raise SystemExit(_self_test())
    raise SystemExit(main())
