# Open gaps — audited state

Every item below was checked against the tree, not recalled. Each is either
**closed**, **closed to a decision**, or **blocked with the blocker named**.

| # | gap | state now | where |
|---|---|---|---|
| 1 | Gate 3 — per-band densitometry | **blocked, and provably so** | `WP24_GATE3_STATUS.md` |
| 2 | 2.73× vs ≥4× tension (P2/P3) | **re-scoped by WP35 (§13 withdrawn, §14); related C1/C2 in §15** | `WP35_P2_IDENTIFIABILITY.md`, `WP36_CIRCULARITY.md` |
| 3 | M2 vs M4 (`k_tr` vs `k_off_tr`) | **open by construction, not load-bearing** | `WP25_RESOLVING_TESTS.md` |
| 4 | 9SHZ / 9SI0 ingestion | **blocked: channel limitation** | `fetch_rap_fixture.py` |
| 5 | `domain_motion()` dead code | **closed — removed (it was also broken)** | `colab/tpp1_colab.py` |
| 6 | Gate 4 (dNTP chase), Gate 5 (δ) | **experimental / a judgement, unchanged** | `paper/MANUSCRIPT.md` §11 |
| 7 | No TPP1-free matched series | **closed — measured** | `WP23_TPP1_FREE_SERIES.md` |
| 8 | TELIN-1 is two fragments, not a molecule | **repaired; 3 sub-items remain open** | `WP30_LIGAND_COMPLETION.md` |
| 9 | No defined TPP1 epitope for docking | **CLOSED — TPP1 was on disk; assembled at 1.97 Å RMSD** | `WP31_EPITOPE.md` |
| 10 | Gate 3 criterion not operational | **stated, and deliberately not fitted** | this file, §1 |
| 11 | Y176A mutant panel (WP16 follow-up) | **designed, not run** — needs a wet lab | `PLAN.md` |
| 12 | 9SHZ/9SI0 ingestion | **still blocked** — re-tested this session | `fetch_rap_fixture.py` |

---

## 1. Gate 3 — now a proof, not a complaint

`data/external/tesmer2023/` records one scalar per condition: the fraction of
signal at ≥16 repeats. WP24 builds two product-length distributions that share
that number **exactly** (0.3500000000 both) while differing **4.98×** in the
fraction of early products. They describe different enzymes and produce the same
recorded row.

A threshold fraction is therefore invariant to the shape below repeat 16 — which
is exactly where per-repeat survival is decided. **No re-analysis can satisfy
Gate 3.** It needs per-band counts, which this repository does not hold.

Separately, and honestly: the criterion itself ("a majority `s` below the tail
slope") is not yet operational — the tail window and the `s` extraction are
unspecified. Fixing those after seeing which synthetic cases pass would be
fitting the criterion to the result. So the gap is recorded, not papered over.

## 2. The 2.73× vs ≥4× tension

Three explanations exist and they are exhaustive. **E1 (`k_tr` falls) is excluded
numerically** — it needs `k_tr` to drop ~300×, while L&C fit it rising. Two
survive:

| branch | prediction under a DNA trap |
|---|---|
| E2 — stalling inflates the measurement | **falls to 2.73×** |
| E3 — ≥4× is a lower bound | **stays ≥4×** |

One experiment, two distinguishable outcomes, and the assay already exists.

## 3. M2 vs M4

`φ = k_tr/(k_tr + k_off_tr)` is two unknowns against one observation; every
`k_tr` value has a matching `k_off_tr` (0.070/0.00071, 0.100/0.00101,
0.200/0.00202 all give φ = 0.99). Open by construction.

**But M5 does not need it**: φ left the product, so `s = (1+a)⁻⁶` never sees
`k_tr`. The gap is real and not load-bearing.

## 4. 9SHZ / 9SI0

Three channels were tried this session:

| channel | result |
|---|---|
| Python `urllib` / TLS | **blocked** — socket closed |
| `models.rcsb.org` per-atom | works, but ~36 kB of metadata per residue |
| `files.rcsb.org` mmCIF | works — 197 chunks of mixed content |

The parser is written and verified offline against a live single-atom pull
(Y176 CA = 133.669, 153.363, 143.095). `fetch_rap_fixture.py` needs one run on a
machine with ordinary network access. **The code is not the blocker; the sandbox
is.**

## 5. `domain_motion()` — removed

It was dead *and* broken: `if i in set(rt_range) or True` made the filter a
no-op, and it returned a 3-tuple under a `-> float | None` annotation. The
correct implementation already exists inline in `analyse_entry`. Removed (32
lines); colab selftest unchanged at 53/0.

## Simulation-only closure (WP26)

Requested: what simulation must supply to settle this before any lab transfer.
Measured state: the repo's MD is **3.0 ps** per replica on a **97-residue,
1,547-atom** fragment — **~8 orders of magnitude** short of the ~100 µs needed
for an equilibrium answer, and the **catalytic aspartates are not in the system**.

| route | settles | closable? |
|---|---|---|
| **S3** alchemical free energy | is the TPP1–TEN interface energetically real? | **yes** — 1st priority |
| **S2** pathway, all-atom graph | is a mechanical conduit plausible? | yes, but WP23 already measured ≤0.53 Å |
| **S1** equilibrium | M1 vs M4 | expensive and pre-empted by WP23's endpoint null |
| **S4** kinetics | the 2.73× vs ≥4× tension | **no** — needs 1.46× precision on k_off; methods give 2–3× |

So the **structural half is closable in simulation; the kinetic half is not** —
and that is a property of the degeneracy, not of available compute.

## 6. Gates 4 and 5

Gate 4 needs a dNTP-chase experiment. Gate 5 needs a judgement on the smallest
meaningful effect δ — a decision, not a computation, and not one this repository
can make on the reader's behalf.

## 7. TPP1-free matched series — closed

7TRD is TPP1-free and matches 7TRE at **934/934 residues (100 %)** in every
domain. Result: the largest domain displacement is **0.53 Å**, TEN is
**0.435 Å**. That is at or below the coordinate uncertainty at 3.2–3.3 Å.

So: *no displacement larger than the uncertainty is detectable on TPP1 binding.*
This is a **measured null**, and it constrains the mechanism — a large-scale
mechanical switch is not supported. Caveat: 7TRD and 7TRF share one TERT model
(0.0005 Å apart), so there is one apo model, not two.


## 8. TELIN-1's missing part — repaired, with what stays open

`WP30_LIGAND_COMPLETION.md` closed the *chemical* half of gap 8 and left three
things open, named rather than papered over:

| | state |
|---|---|
| one connected molecule, correct formula | **fixed** — `data/ligands/telin1_repaired.pdb` |
| impossible geometry (1.208 Å ring bond, 1.715 Å C–C, 158.9°) | **fixed** |
| which aryl position carries the arm | **open** — three positions tie on formula |
| docked pose (160 clashes, C19 0.17 Å from a backbone N) | **open** — must be redone elsewhere |
| TPP1 epitope to dock against | **fixed in WP31** — TPP1 is in 7TRE chain C and 2I46 chain A |

The last one is the gap that matters. The complex holds TERT, TERC and DNA and
`no TPP1 chain`, so the pose never tested the claimed mechanism. **The part
TELIN-1 is really missing is a target**, and no ligand edit substitutes for it.

## 9. Re-tested this session, still blocked

RCSB is reached over TLS only, and every endpoint fails the handshake from this
sandbox:

| endpoint | result |
|---|---|
| `files.rcsb.org/download/9SHY.cif` | `SSL_ERROR_SYSCALL` |
| `models.rcsb.org/9SHY.pdb` | `SSL_ERROR_SYSCALL` |
| `www.ebi.ac.uk/pdbe/entry-files/download/9shz.cif` | connection refused |
| `github.com` (control) | **200** |

So the block is specific to those hosts, not a general egress failure — gap 12
stays blocked for a reason, and `fetch_rap_fixture.py` still needs one ordinary
run elsewhere.


## 10. Gap 9 retracted, and how it was wrong

Gap 9 said the TPP1 epitope was "blocked: TPP1 absent from every file on
disk". The second clause was false. `data/raw/7TRE.cif` chain C and
`data/raw/2I46.cif` chain A are both TPP1 TEN, 152 residues each with identical
sequence; only `data/ligands/7tre_telin1_complex.pdb` had the chain stripped.
The check that produced the wrong answer looked at the *complex* file and
generalised to the repository.

`WP31_EPITOPE.md` assembles TERT+TPP1 (superposition RMSD **1.966 Å** over
152/152 Cα), defines the recruitment epitope as **26 TERT + 21 TPP1 residues
at 5 Å**, and runs every ligand against it. No network access was needed.

What remains open is smaller and more specific: a real free-energy calculation
(MM-GBSA / FEP) at the TEL patch. That is gap 13.

| # | gap | state now | where |
|---|---|---|---|
| 13 | Affinity at the TEL patch | **open — needs FEP/MM-GBSA; now with a control showing why** | `WP32_PHYSICS.md` |
| 14 | No ligand prefers the epitope | **superseded — the neutral amide does prefer it (WP33); the anchor's charge sign is not the whole story** | `WP33_CORRECTED_SCORING.md` |


## 11. Gap 14 opened by closing a capability gap

WP32 replaced WP31's contact counter with a function that has partial charges,
ionisation, directional hydrogen bonds and burial — and gave it a control patch
it must not prefer. It preferred the control patch for **every** ligand.

Two consequences, both recorded rather than smoothed over:

- **WP31's "convergence" claim is retracted.** TLN-3 led on contact counts;
  with charges it does not lead, and no ligand is interface-preferring.
- **The cause is identified and it is chemistry, not code.** The TERT–TPP1
  interface is net **−17.9**; every ligand in the series is anionic. Swapping
  the anchor from carboxylate to amine cuts the anti-selectivity **6.2×**
  (+20.5 → +3.3), which confirms the mechanism and still does not produce a
  selective molecule.

Gap 13 is therefore better specified than it was: what is needed is not "a
better scorer" but a calculation with flexibility, solvent and entropy. Gap 14
is the design problem that survives it.


## 12. Gap 14 re-scoped, because WP32's evidence was confounded

WP33 found four defects in WP32's comparison and fixed each with a test:

| defect | fix | effect |
|---|---|---|
| control differed by 1.68× in atom count | matched control (0.91 ratio) + per-atom energies | the size confound is gone |
| prose said the patch was both acidic and basic | measured it: **ambipolar, net −1** | "the anchor must be basic" **withdrawn** |
| five of six variants changed the scaffold | matched protonation pairs | charge is now the only variable |
| H-bond test varied distance and angle together | constant 3.00 Å, angle verified | the 120° boundary is actually tested |

The corrected result is more interesting than the one it replaces: the only
variant preferring the interface is the **neutral amide**, and charge does not
order the series monotonically — which is the expected signature of an
**ambipolar** patch. So gap 14 is no longer "these ligands are the wrong charge";
it is "the anchor should probably not be charged at all, and that needs a real
free-energy calculation to test."


## 13. [WITHDRAWN IN PART - see §14] The 1.46x tension is 88% explained

> **This section is withdrawn on two counts, both shown in §14.** The "88 %"
> compares a probability to `(s-phi)/s`, which *is* the anomaly M5 exists to
> explain - so it explains itself. It is also negative in the +PT arm. And the
> consequence drawn from it ("the retry probability should measure about 0.32")
> does not follow: 0.318 is not a retry probability at all. Retained for the
> record, not for use.

WP34 found that the gap between M5's 2.73x and L&C's measured >=4x is almost
entirely accounted for by the stalled population L&C report directly:

    required retry probability   1 - 2.73/4.00        = 0.318
    stalled fraction observed    (0.7914-0.570)/0.7914 = 0.280   -> 88% of the gap

So gap 2 is no longer "an open discrepancy with a plausible explanation". It is
a **falsifiable prediction with a number attached**: the retry probability
should measure about **0.32**, and the P2 DNA-trap re-analysis measures that
quantity directly. This upgrades P2 from "the cheap decisive test" to "the test
that tests a specific number".

Full status, candidate and odds: `WP34_STATUS.md`.



## 14. P2 implemented - and the gate is undetermined, not closed

WP35 built step 1 of WP34's order of operations. It does not re-analyse
anything, because the repository has nothing to re-analyse (WP24's proof), and
doing it correctly changed two conclusions.

**The stall correction had the wrong sign.** WP34 applied stalling to the
chase only. The product distribution is stall-inflated too, and a retry makes
the observed per-repeat survival *overstate* the true one:

    s_obs  = s_true / (1 - tau(1-s_true))
    s_true = s_obs (1 - tau) / (1 - tau s_obs)

Corrected that way, the distribution fold **rises** with tau while the
trap-diluted chase **falls**. They converge on one tau:

    tau_0 = 0.1833     fold = 3.2670     delta* = 57.30 bp/PD

**So the gate is UNDETERMINED.** WP34's 48.0 bp/PD is delta* at tau_0 = 0, and
tau_0 = 0 is exactly what P2 tests. The two physically derived couplings give
57.3 and 57.4 bp/PD - above the 50 floor. The pessimistic "same tau in both
arms" coupling gives 46.8, and the phi data contradict it. The spread from 46.8
to 57.4 is assumption-driven, so all three are reported and none alone.

**P(approved drug) is therefore not 0 %.** It is WP34's own pipeline product,
conditional on the gate opening: ~0.008 %. "0 % because the gate is closed" was
a conditional presented as a closure.

What does not change: 57.3 bp/PD is still short of a 90 bp/PD cell, which needs
5.14x - so the maintenance objective stays out of reach without a real change
in the true off-rate, which is step 2 and remains unmeasured.

**Cost.** P2 is 3.27x against 4.00x, a ratio of 1.224, against a real
experiment-level CV of 17.2 % (Tesmer Fig 3B WT). That is n = 12 per arm
independent, or n = 3 paired (rho = 0.8). "Days, no reagent" assumed the
pairing silently.



## 15. The error CLASS, not just the instances (WP36)

Four failures so far were the same mistake, so fixing them one at a time was
never going to stop the next. WP36 names the class and makes it testable.

**Circular inference:** a quantity is computed from a construction that already
contains the conclusion, then reported as evidence for the conclusion. Every
instance passes its arithmetic. That is why it survived review four times.

| id | claim | what it really is | predicate |
|---|---|---|---|
| **C1** | M5 independently reproduces the observed 48 nt | `δ* = n·R½` **identically**; 47.98 *is* 6×8 | holds for every `n` (swept 1–24) |
| **C2** | the stalled fraction explains 88 % of the gap | `(s₀−φ₀)/s₀` is the anomaly itself, rescaled | negative where φ→1 |
| **C3** | synchronised confirmation vs chance | the expected value came from the null being tested | equals the null exactly |
| **C4** | removing all 12 loops changes output 4.4 % | detector and remover share a definition | invariant for any input |
| **C5** | the "physical" SPR/RAP benchmarks | `np.random.seed(42)` synthetic | recovers its own parameters |
| **C6** | M5 passes an out-of-sample test on Wang 2007 | no measured off-rate exists; 4.48× is `a₀/a₁` from Wang's own R½ | recovers its own input (2e-15) |
| **C7** | δ\* decides the maintenance gate | `δ* = n·R_eff`; at a fixed fold it runs 8.9 → 231 bp/PD over `n = 1…24` | depends on `n` alone |

**C1 is the expensive one.** `δ*` tracks `n` over a 24× span while the fold `κ`
moves only 7.5 %, so `δ*` carries no information — `κ` is the only informative
statistic. WP34 scored "M5 derives the observed 48 nt exactly" at **90 %** and it
was the largest single term in `P(M5)`.

**C6 and C7 are the more serious pair.** C6 means `s = (1+a)⁻⁶` has *never* been
tested on a second enzyme — the only measured fold is the L&C datum it was built
on. C7 means the maintenance gate, the thing WP34 declared closed and WP35
reopened, is a statement about an **assumed exponent** rather than about the
enzyme: at a fixed 3.267× fold, δ\* is 8.9 bp/PD at `n = 1` and 231 bp/PD at
`n = 24`, straddling the 50 floor between `n = 4` and `n = 6`.

**P(M5) is ≈ 35 %** — the same figure WP34 quoted, with none of the same terms.
WP34's was `0.90 × 0.65 × 0.60` (two circular); this is `0.70 × 0.50` (neither
is). The one genuinely independent result in the programme remains `s₀ ≤ φ₀`
holding on L&C's own data.

Guards now in the suite: `δ* == n·R½` for all `n`, the +PT value is negative, the
synthetic fixture reconstructs its inputs, and every retracted document carries a
machine-checked marker. Full text: `WP36_CIRCULARITY.md`.



## 16. The four missing things, and what WP37 did with each

| # | missing | after WP37 |
|---|---|---|
| 1 | a measurement of `n` | **CLOSED as a design** — one chase cannot, two chases can, exactly |
| 2 | `δ*` defined twice | **CLOSED** — physical `6·R½`; the model form equals it only at `n = 6` |
| 3 | a second enzyme | **PARTLY** — a second R½ pair exists (1.82× prediction); still no measured fold |
| 4 | a TEL-patch affinity | **NOT POSSIBLE HERE** — 7-step route + a docking proposal with a negative control |

> **CORRECTED BY WP40 (§18 below).** The conclusion of this section is wrong and
> is retained only for the record. "CLOSED at 48" is a statement about **M5's own
> prediction** (8.0/2.963 = 2.70×), not about L&C's measurement, which is ≥4×
> and gives **70.1 bp/PD**. The gate is also a function of the unmeasured stall
> share and opens at **τ ≥ 0.042**. Read §18 before using anything below.

**The gate is CLOSED, and now for a measured reason** *(superseded — see §18)*. WP29's objective was
already exponent-free: `6·c·R½ − δ`. At the +PT anchor that is **48 bp/PD**,
below the 50 bp/PD floor, and it does not depend on `n`. WP35's "OPEN at 57.3"
substituted a model expression into an objective that had already been fixed
correctly. The only escape hatch is `c > 1` — a 50 bp/PD cell needs 1.04
extensions per division, a HUVEC-like 90 bp/PD cell needs 1.88 — and nothing in
this repository establishes that. It is recorded as an assumption.

**The two-chase inversion** is the substantive addition. `m = κ(1+τ)` and
`c = κ(1−τ)`, so `m/c` gives τ with **no model assumption**, then `κ = c/(1−τ)`
and `n` follows. Exact across a 25-point grid. It also yields a free
pre-registration: `m > c` always, so a trap that *raised* processivity would
refute stalling outright.

**What is still impossible here:** a measured second-enzyme off-rate fold, any
affinity at the TEL patch, `c > 1`, and Gate 3 per-band densitometry. Four
things, all wet lab or deposition, none of them closable by more analysis.



## 17. Two items were misclassified as blocked (WP39)

WP37's closing table said six things could not be closed here. Two of them could,
and both were classified as blocked **before the literature and the code on disk
were checked**.

**Gap 5 — `c > 1` extensions/division. CLOSED, and it closes the gate harder.**
`c` is the fraction of telomere ends extended per cycle, and it is measured per
cell line: H1299 = **1.00** (Zhao 2011, *Mol Cell* 43:251), HeLa = 0.70, hESC =
0.70, HeLa at pH 7.4 = 0.77, at pH 6.8 = 0.48, *S. cerevisiae* = 0.07. At the
best measured value the ceiling is exactly 48 bp/PD — which **misses the 50 floor
by 2 bp**. `c_required(50) = 1.042`, above every measured value, so the escape
hatch is not open anywhere in the data. The gate now depends only on R½, which
would need 8.33 (floor) or 15.0 (HUVEC) against 8.0 measured.

**Gap 4 — the TEL-patch docking criterion. RUN, and the proposal is refuted.**
WP37 wrote a refutation criterion with a threshold, a mutant and a negative
control, and never evaluated it. WP39 evaluated it against the seven-residue TEL
patch (not the 47-residue interface, which is the confound WP33 retired). The
neutral amide moves **−0.046** per heavy atom across a full neutralisation — the
wrong direction, below the 0.15 threshold, and *smaller* than both controls. No
affinity is claimed for anything; the surrogate that was supposed to substitute
for one is now refuted.

**Still genuinely blocked:** a measured `n`, a measured second-enzyme fold, P2 as
a two-chase, any affinity at the patch, and Gate 3 densitometry.



## 18. The 50 bp/PD gate was read off the wrong input (WP40)

`R½ = 8.0` has been used as "the +PT measurement" since WP29. It is not one:

    8.0 / 2.963 = 2.7000x        <- M5's PREDICTED fold

WP29 labels that row *"M5 prediction (2.73×)"*. L&C's measured statement is
**≥4×**, i.e. `R½ ≥ 11.85`, i.e. **70.1 bp/PD**. So:

| anchor | R½ | bp/PD | provenance |
|---|---|---|---|
| no TPP1 | 2.963 | 17.8 | measured |
| "8.0, +PT" | 8.00 | 47.4 | **M5's prediction, fed back** |
| L&C ≥4× | 11.85 | **70.1** | measured lower bound |

**"The gate is CLOSED at 48" is a statement about M5's prediction.** On the
measurement the gate is open, and by a wide margin. The closed gate is therefore
evidence against **M5**, not against the target: the mechanism predicts less
processivity than has already been seen.

**Second correction, same section.** The gate is not a fixed number at all — it
is a function of the stall share τ, and it opens at **τ ≥ 0.042**, not at
WP35's default of 0.1833. "CLOSED at 48" (τ = 0) and "OPEN at 57.3" (τ = 0.1833)
are one function at two modelling choices. The honest statement is that the gate
is open for any stall share above 4.2 %, and P2 exists to measure whether it is
that big.

**And E2 is closer to dead than it looked.** To make an observed 4.00× into a
true 2.7335×, `(1+τ₋)/(1+τ₊) = 1.4633`, which needs `τ₊ < 0` — impossible —
for any `τ₋` below **0.4633**. E2 as written contradicts the stall model it is
meant to rescue.

**Unchanged:** the 90 bp/PD HUVEC ceiling still needs τ ≥ 0.88 or R½ ≥ 15.0, and
nothing here supplies that. No affinity exists. TELIN-1 remains quarantined.



## 19. All six blocked items were attempted (WP44)

They were not all closable, but five changed status and **one was never a gap**.

**Item 1 — "measure `n`" — never load-bearing.** `n` enters only through
`a = s^(−1/n) − 1`, which *is* the identity C1 retracted. Recomputed from the
measured `s₀`, `R½` runs 7.53 → 8.09 across `n = 1…64`: **7.4 % over a 64-fold
range**. `R½`, `δ*` and the gate are `n`-free. Removed from the ledger.

**Item 3 — "P2 as a two-chase" — the design was vacuous.** `m = κ(1+τ)` and
`c = κ(1−τ)` are one relation written twice, so `m/c` returns τ by construction;
the 1e-17 round trip proves algebra (C1 pattern). A real two-chase needs an
**independent count** of the stalled population.

**Item 5 — Gate 3 — unreachable, now shown.** Both conventions reproduce the
reported fraction 0.485 at survivals implying R½ = 2.395 and 1.036 — a **2.3×
spread from one number**.

**Item 6 — the Y176A assay is blind to its own question.** BrdU/CsCl counts
repeats added, not attempts made. But the same paper's partial-inhibition series
is a model-free **mode** test: at 0 / 50 / 70 % inhibition the extended fraction
falls 1.00 → 0.50 → 0.30 while **product length stays ~60 nt**. A distributive
enzyme would shorten products. **Positive, published evidence for processive
action**, with no assumption from this repository in it.

**Item 2 — still blocked**, and now for a measured reason: the Tesmer readouts
are unitless fractions with no rate and no affinity. It adds a constraint —
**F776E** is mid-panel on activity (0.418, rank 5/9) while processivity is at the
floor (0.013), so the two are separable inside one experiment.

**Item 4 — still blocked, demonstrably**: the surrogate has no solvent, no
entropy, no flexibility and **zero ligands of known affinity** to calibrate
against. A score cannot become a `Kd`.

**Bonus.** `s_obs = s_true/(1−τ(1−s_true))` is itself geometric with a rescaled
probability, so the renewal does not break the `R½` framework — **provided the
retry is memoryless**. Nobody had checked; if a stall raises the chance of
stalling again, every `R½` here is an effective value.


---

## §20 — WP45: the creative route, and what it cost

**Route taken.** Instead of chasing a seventh experiment, τ was recovered from
two quantities L&C already published from the same preparation
(`s_obs = 0.7914`, `φ = 0.570`). Both are functions of τ under the stall
renewal, so the pair fixes it: **τ₀ = 0.6506**.

**Result — the programme's leading mechanism is refuted.**
With τ fixed, M5 predicts a +PT fold of **2.734** against **6.654** implied by
L&C's own data — a **2.4× under-prediction**. Registered as **R12**.

**Result — the gate number does not survive.** The stored `R½(+PT) = 8.0` is
recovered only by assuming the +PT stall *vanishes*, while τ(−PT) is
simultaneously large enough to move the fold. Those cannot both hold. Applying τ
consistently gives **δ\* = 7.9 bp/PD**, six times lower. Registered as **R13**.

**Load-bearing assumption: `s_true = φ`.** Everything above rests on reading φ as
the un-retried per-repeat success. It is a modelling step, not a theorem, and it
is the single item that decides whether R12 and R13 stand. Named as **A1**.

**Does this close any gate?** No. It moves them further away. The 90 bp/PD
ceiling is unchanged and is now further out on the true reading.

**Status of the six blocked items:** unchanged from §19 — none converted.
