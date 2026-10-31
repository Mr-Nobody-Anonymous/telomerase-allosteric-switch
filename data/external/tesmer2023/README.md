# External data: Tesmer et al. 2023, Fig 3 source data

## What is here

Transcribed numeric source data for the IFD-TRAP mutant panel, copied from the publisher's
source-data file for Fig 3. These CSVs are the **primary quantitative dataset** for the
affinity–processivity causality investigation.

| file | upstream | quantity | semantics |
|---|---|---|---|
| `fig3b_processivity_plusPT.csv` | Fig 3B, +PT lanes | processivity (fraction of signal in ≥16-repeat products) | 3 rows = 3 independent experiments |
| `fig3c_activity_plusPT.csv` | Fig 3C, +PT lanes | relative total activity, normalised to WT **within each experiment** | 3 rows = 3 independent experiments |
| `fig3e_motif3n_activity.csv` | Fig 3E | normalised activity, motif 3N (Δ641–651) ± PT | 3 rows = 3 independent experiments |
| `SOURCE_MANIFEST.csv` | — | SHA-256 of each local file, replicate semantics | — |

## Provenance status — read before using these numbers

**Verified:**

* Upstream location: `https://pmc.ncbi.nlm.nih.gov/articles/instance/10055720/bin/LSA-2022-01727_SdataF3.1.pdf`
  (PMC10055720, Tesmer et al., *Life Sci Alliance* 6:e202201727).
* The transcription reproduces the published summary statistics. Independent recomputation of
  column means gives H762E = 85.6% processivity (text: 86%), D788R = 64.3% (text: 64%),
  H762E = 68.0% activity (text: 68%), D788R = 23.5% (text: 24%). Four independently stated
  numbers are recovered from the transcription, which is the check that the cells were copied
  correctly. This check is automated in `tests/test_causality_preregistration.py`.
* The files are three independent experiments per the paper, not technical replicates. No
  technical replicate values are present in the source file.

**Not verified:**

* **No upstream SHA-256.** The sandbox can reach `pmc.ncbi.nlm.nih.gov` only through a text
  extractor; direct download fails (`OpenSSL SSL_connect: SSL_ERROR_SYSCALL`), so the upstream
  file's own hash could not be taken. The manifest records this as `NO` rather than leaving the
  field blank. The transcription is therefore verified *by its ability to reproduce published
  statistics*, not by byte-level comparison with the source. If anyone can download the PDF,
  take its SHA-256 and replace `sha256_upstream`; nothing else needs to change.

## What is **not** in this dataset — the critical limitation

**There is no product-length distribution here.** Each replicate contributes a *single scalar*
per mutant: the fraction of signal in products of ≥16 repeats. The five distribution models that
this investigation was designed to compare (geometric, mixtures of geometric, repeat-dependent
extension, heterogeneous population, conformational) require **per-band intensities as a function
of product length**, which the numeric source file does not contain. Those data exist only inside
the gel images (`LSA-2022-01727_SdataF3.2.tif`, `SdataF3.3.tif`, 8.1 and 7.5 MB) and in the
figures' scanned gels, neither of which can be retrieved in this environment.

Consequence: distribution-shape analysis is **blocked on data, not on method**. See
`results/causality/PREREGISTRATION.md` §5 for how the plan handles this.

## Known data-quality flags

* **`F776E`, experiment 3 = −0.337890`.** A negative processivity fraction is physically
  impossible; it arises when the ≥16-repeat signal is at or below the background estimate. The
  value is retained verbatim rather than clipped, and the mutant is flagged
  `processivity_undefined_at_low_signal` and excluded from quantitative comparison. It is not
  excluded from the qualitative statement "F776E is essentially dead", which is supported by
  activity.
* **`L786E` is labelled `L786` in the upstream table.** Transcribed as `L786E` with the typo noted
  here.
* **The published percentages are approximately, not exactly, the mean of these rows.** The mean
  of the D788R activity rows is 23.5% while the text says 24%; H762E reproduces exactly (68.0% vs
  "68%"). The half-point gap is left visible rather than tuned away, and the preregistration
  therefore treats published percentages as ordinal.
* **Large inter-experiment spread in activity.** Experiment 2 reports mutant activities ~10× lower
  relative to WT than experiment 1 (e.g. K749E: 1.133 vs 0.148). The paper attributes imprecision
  to the low signal-to-noise of *trans*-complementation and states it could not establish
  significance for activity differences. This spread is the dominant uncertainty in everything
  downstream and is quantified in `PREREGISTRATION.md` §4.
