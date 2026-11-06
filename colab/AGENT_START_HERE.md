# START HERE — telomerase TEN anchor-site test

You are an automated agent. Read this file, run one command, report the verdict.
Nothing below requires judgement on your part.

---

## Run this

```bash
bash colab/run_colab.sh
```

That is the whole task. It exits **0** on success and non-zero on failure.

If you prefer Python:

```bash
python colab/tpp1_colab.py            # full run, needs internet
python colab/tpp1_colab.py --selftest # offline check, no internet needed
python colab/verdict.py --selftest    # the geometry + verdict rule, offline
python colab/test_colab.py            # verifies every claim in THIS file
```

Runtime: about two minutes. CPU only. No GPU needed and none is used.

## What it does

1. Downloads telomerase structures from RCSB (7QXS, 7QXA, 9SHY, 9SHZ, 9SI0, 7TRE, 7TRD, 7TRF).
2. Checks each structure's **active site is seated on the DNA** — catalytic Asp868 must be
   within 6 Å of the DNA 3′ end. This decides whether an entry is usable.
3. Measures the TEN anchor motif PLYQ (174–177) against the DNA.
4. Applies the fixed verdict rule and writes `colab_out/report.md`.

## What the answer already is

The measurement has been made and cross-checked. Two independent entries agree:

| entry | Tyr176 to DNA | catalytic control | verdict |
|---|---|---|---|
| **7QXS** | **2.75 Å** (H-bond to a phosphate O3′) | 2.65 Å ✓ | **CONTACT** |
| **7QXA** | **4.35 Å** (vdW to a sugar C5′) | 2.93 Å ✓ | **CONTACT** |
| 9SHY | 30.23 Å, but the relevant nucleotide is unmodelled | 2.78 Å ✓ | NOT RESOLVED |
| 7TRE / 7TRD / 7TRF | 11.11 / 11.94 / 11.94 Å (7-nt DNA, never reaches the site) | — | out of reach |

**Verdict: the TEN PLYQ anchor site does physically contact the DNA** once the substrate is
long enough. The contact belongs to **Tyr176**, and it is to the **sugar–phosphate backbone**,
not to a base — so it cannot be reading sequence.

Full write-up: `results/solve/WP16_ANCHOR_SITE_RESULT.md`.

## What this does NOT show

Read this part before reporting anything.

- **A contact is not a mechanism.** It does not distinguish DNA retention from catalytic
  or allosteric acceleration. The standing question (TPP1-stimulated RAP: mechanical or
  retention?) remains **open**.
- **Backbone contact, not base contact.** Consistent with a positioning element and with a
  steric barrier equally.
- **Do not quote 2.75 Å as a precise number.** These are 3.2–3.8 Å cryo-EM maps. What is
  robust is the ordering: ~3 Å vs ~11 Å, a factor of four.
- **Do not report "no contact" from 9SHY.** Its anchor-contacting nucleotide is not modelled.
  That is missing data, and the rule returns NOT RESOLVED precisely so it cannot be misread.

## Output files

| file | contents |
|---|---|
| `colab_out/report.md` | human-readable verdict, distance table, per-nucleotide profile |
| `colab_out/results.json` | everything above, machine-readable |
| `colab_out/distance_table.csv` | one row per structure |

## If it fails

- **Exit 2 — no Python with numpy/gemmi/scipy.** The runner tries `.venv`, `python3`,
  `python`, and a `--break-system-packages` install fallback. If all fail, install
  manually: `pip install numpy scipy gemmi`.
- **Verdict reads "TEST DID NOT RUN"** — no structure passed the active-site control. That is
  a **failure to measure, not a negative result**. Report it as such. Do not write
  "no contact".
- **Downloads fail but the self-test passed** — the local fixtures and bundled CIFs are fine;
  the network is the problem. Re-run with internet access.

## Verify the claims in this file

```bash
python colab/test_colab.py
```

This does not trust the text above. It re-derives every number in the table from
`colab/fixtures/`, checks the controls pass, checks every verdict branch is reachable,
and fails loudly if this document and the code disagree.

## Files

| file | role |
|---|---|
| `colab/run_colab.sh` | one-command entry point |
| `colab/tpp1_colab.py` | downloads and analyses full structures |
| `colab/verdict.py` | the verdict rule and its calibration |
| `colab/test_colab.py` | verifies this document |
| `colab/fixtures/` | real coordinate subsets that make the numbers reproducible offline |
