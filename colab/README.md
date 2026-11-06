# Running the structural test on Google Colab (or anywhere with internet)

## TL;DR for an agent

**Read `colab/AGENT_START_HERE.md` first** — that is the document written for you, and
`colab/test_colab.py` verifies it is still accurate.

```bash
git clone https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch.git
cd telomerase-allosteric-switch
bash colab/run_colab.sh
```

That script picks a working interpreter, installs `numpy scipy gemmi` if needed, runs **three
offline self-test suites** (and aborts if any fails), then runs the full download + analysis and
prints `colab_out/report.md`.

**Self-test only (no internet):**

```bash
bash colab/run_colab.sh --selftest
python colab/tpp1_colab.py --selftest   # geometry against the bundled CIFs
python colab/verdict.py --selftest      # the verdict rule against real coordinate subsets
python colab/test_colab.py              # does AGENT_START_HERE.md match the code?
```

## The answer, already measured

The question — does the TEN anchor motif PLYQ (174–177) contact the DNA once the substrate is
long enough? — has now been answered from deposited coordinates. Three independent entries put the same residue on the DNA backbone:


| entry | modelled DNA | Y176 → DNA | catalytic control | verdict |
|---|---|---|---|---|
| **7QXS** | 22 of 30 nt | **2.75 Å** (OH···O3′) | 2.65 Å ✓ | **CONTACT** |
| **7QXA** | 12 of 30 nt | **4.35 Å** (OH···C5′) | 2.93 Å ✓ | **CONTACT** |
| **9SHY** | 10 of 34 nt | **2.78 Å** (OH···OP1 of nt 25) | 2.78 Å ✓ | **CONTACT** |
| 7TRE / 7TRD / 7TRF | 7 nt | 11.11 / 11.94 / 11.94 Å | — | never reaches the site |

**Yes — the anchor site contacts the DNA**, via **Tyr176**, against the **sugar–phosphate
backbone** rather than a base. Full write-up: `results/solve/WP16_ANCHOR_SITE_RESULT.md`.

What that does **not** establish: any mechanism. A backbone contact is equally consistent with a
positioning element and with a steric barrier, and TPP1-stimulated RAP (mechanical vs retention)
remains experimentally open. Do not quote 2.75 Å as precise — these are 3.2–3.8 Å maps; the robust
statement is the factor-of-four ordering, ~3 Å versus ~11 Å.

## Why the first version of this test could not have answered it

The original runner gated on **"at least one entry with ≥15 modelled nucleotides"**. 7QXA models
12 and 9SHY models 10, so both were discarded — and 7QXA is one of the three entries that settle
the question. The runner would then have printed `TEST DID NOT RUN`.

The gate was testing the wrong thing. Usability is now decided by a **positive control**: catalytic
Asp868 must lie within 6 Å of the DNA 3′ end. If the active site is seated, the coordinates mean
something; size is only used to decide whether a *negative* is admissible, never to discard a
structure. Three outcomes are now distinguished, and never conflated:

- `CONTACT` / `NO CONTACT` — measured
- `NOT RESOLVED` — the relevant nucleotides are unmodelled (missing data)
- `TEST DID NOT RUN` — no structure passed the control (failure, **not** a negative result)

## Files

| file | role |
|---|---|
| `colab/AGENT_START_HERE.md` | **the document an agent reads**; verified by `test_colab.py` |
| `colab/tpp1_colab.py` | downloads and analyses full structures |
| `colab/verdict.py` | the verdict rule and its calibration |
| `colab/test_colab.py` | asserts the document's claims against the code |
| `colab/run_colab.sh` | one-command runner, three self-tests, exit codes propagate |
| `colab/TPP1_structural_test.ipynb` | thin wrapper; no logic of its own, so it cannot drift |
| `colab/fixtures/` | real coordinate subsets making every number above reproducible offline |

## Two bugs this caught in itself

1. **Histones mislabelled as TPP1.** A chain-length heuristic picked up 7TRF's H2A/H2B (90/82
   residues) and reported a fake `TPP1→DNA = 34.5 Å` for a structure with no TPP1. TPP1 is now
   identified by contact to TERT.
2. **Self-defeating superposition.** The RAP-state comparison superposed each domain
   independently, absorbing rigid-body motion and reporting **0.0 Å** for exactly the displacement
   it was built to detect. Found only by injecting a known 1.5 Å shift. It now superposes on the
   RT ring and measures the residual.

Both are the same failure mode as the ≥15-nt gate: a measurement that quietly returns a
plausible-looking wrong answer. The self-tests now inject known values rather than checking that
the code merely runs.

## Independent checks used

- **Catalytic control**: D868 → DNA 3′ end is 2.65 / 2.93 / 2.78 Å across three entries. A
  correct active site, verified per structure, not assumed.
- **Two entries, same nucleotide**: 7QXS and 7QXA are independently refined models; both put Y176
  on auth 19.
- **Sanity of scale**: 11.11 Å ÷ 2.75 Å = 4.0×.
- **Pro174** is 7.49 Å from the DNA, so the contact belongs to the tyrosine, not to the whole motif.
