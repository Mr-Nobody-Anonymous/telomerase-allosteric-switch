# Open gaps — audited state

Every item below was checked against the tree, not recalled. Each is either
**closed**, **closed to a decision**, or **blocked with the blocker named**.

| # | gap | state now | where |
|---|---|---|---|
| 1 | Gate 3 — per-band densitometry | **blocked, and provably so** | `WP24_GATE3_STATUS.md` |
| 2 | 2.73× vs ≥4× tension (P2/P3) | **re-scoped by WP35: §13 is withdrawn, see §14** | `WP35_P2_IDENTIFIABILITY.md` |
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
