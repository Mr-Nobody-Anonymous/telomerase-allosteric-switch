# Running the structural test on Google Colab (or anywhere with internet)

## TL;DR for an agent

```bash
git clone https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch.git
cd telomerase-allosteric-switch
bash colab/run_colab.sh
```

That script picks a working interpreter, installs `numpy scipy gemmi` if needed, runs the
**offline self-test** (and aborts if it fails), then runs the full download + analysis and prints
`colab_out/report.md`.

Two artefacts matter:

* `colab_out/report.md` — human-readable verdict + distance table
* `colab_out/results.json` — machine-readable, everything

**Self-test alone (no internet needed):**

```bash
python colab/tpp1_colab.py --selftest     # or: bash colab/run_colab.sh --selftest
```

**No internet at all, use the three bundled CIFs:**

```bash
python colab/tpp1_colab.py --offline
```

## Files

| file | role |
|---|---|
| `colab/tpp1_colab.py` | **the whole test** — one self-contained script, no notebook dependency |
| `colab/run_colab.sh` | one-command runner with dependency detection and gating |
| `colab/TPP1_structural_test.ipynb` | thin wrapper; calls the script. Contains no logic of its own, so it cannot drift |

## The question it settles

Does the TEN-domain anchor motif **PLYQ (174–177)** contact the DNA once the substrate is long
enough to reach it?

All local structures (7TRE / 7TRD / 7TRF) carry only a **7-nucleotide** DNA. Its 3′ end sits on the
catalytic aspartates (D868 2.65 Å) and its 5′ end on the anchor strand R756 (5.02 Å) — stopping two
nucleotides short of PLYQ. Extrapolating the 5′ exit direction put Tyr176 at **6.8 Å**, but that was
a straight line through one nucleotide step, not an observation. **7QXS / 7QXA (30 nt)** settle it.

### The verdict rule is fixed in advance, in the script

| closest PLYQ–DNA distance in a ≥15 nt structure | verdict |
|---|---|
| ≤ 5.0 Å | **ANCHOR SITE ENGAGES DNA** — the 11.1 Å figure measured an unoccupied site |
| ≥ 10.0 Å | **ANCHOR SITE DOES NOT ENGAGE DNA** — the repository was right |
| 5–10 Å | **INTERMEDIATE — approach, not contact**; report as such |
| no ≥15 nt structure available | **TEST DID NOT RUN** — explicitly *not* a negative result |

All four branches are exercised by the self-test with fabricated records, so the rule is verified to
be reachable rather than assumed.

## What the self-test actually checks

1. **Known values reproduce** against the bundled CIFs: PLYQ = 11.11 Å (7TRE), 11.94 Å (7TRD/7TRF);
   anchor strand = 5.02 / 4.35 / 4.34 Å. Tolerance 0.05 Å.
2. **TPP1 is found by contact** — 7TRE chain C, and **ABSENT** in 7TRD and 7TRF.
   This catches the bug that a chain-length rule causes: 7TRF carries histone H2A/H2B at 90/82
   residues, which a length heuristic mislabels as TPP1. That bug was live in an earlier draft and
   produced a fake "TPP1→DNA = 34.5 Å" for a structure with no TPP1.
3. **Every verdict branch is reachable** (5 fabricated cases).
4. **The RAP-state comparison recovers injected motion**: a known 1.5 Å / 3.0 Å shift in TEN is
   recovered as 1.50 / 3.00 Å while the RT ring reads 0.00 Å.

Item 4 exists because the first version of that comparison was **silently wrong**. It superposed
each domain independently, so a rigid-body motion was absorbed by the fit and reported as ~0.0 Å —
exactly the motion it was meant to detect. The test injects a known answer and fails if it vanishes.

## Interpreting the RAP-state section

The three 2026 states (9SHY / 9SHZ / 9SI0) have **TPP1 present in all three**. So any motion they
show is **cycle motion, not TPP1-driven motion** — the script prints this explicitly next to the
numbers. Do not read it as evidence about TPP1.

## What it cannot do

It cannot settle the mechanism. Structures show what is reachable, not which route TPP1 uses. That
needs the bench measurements in `../results/solve/PREREGISTRATION_WP14.md`:

1. **A on untreated enzyme** — one reaction; sets the operating point.
2. **w with a 95% CI** — active fraction; decides whether the feasibility gate is interpretable.
3. **tau across ≥3 purifications** — decides whether the study can be sized at all.

## On the GPU

Colab Pro VRAM does not help this project.

| task | GPU value | should you? |
|---|---|---|
| this test | none | CPU runtime, ~2 min |
| the kinetic algebra (wp10–wp15) | none | already done |
| long-timescale MD of TEN/TPP1 | real, days–weeks | **no** — cannot identify a rate the observables don't constrain |
| ColabFold on a TERT mutant | modest, minutes | maybe, if a mutant has no structure |
| cryo-EM processing | large | only if you collect data |

The repository already tried 3.0 ps MD and retired an ANM result. The open question is
**experimental identifiability**, not sampling.

## If a download fails

The script reports the failure and the verdict becomes **TEST DID NOT RUN**. It does not fall back
to a weaker claim. Three errors in this project came from reading a missing measurement as a
negative one, so the failure mode is deliberately loud.
