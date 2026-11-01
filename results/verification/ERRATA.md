# ERRATA — independent verification pass of 2026-10-30

This file lists every change made to the repository as a result of re-deriving its headline claims
from the deposited coordinates. The mechanism, tolerances and full results are in
[`CLAIM_VERIFICATION_REPORT.md`](CLAIM_VERIFICATION_REPORT.md); machine-readable results in
[`claim_recheck.csv`](claim_recheck.csv); the corrected scientific claim in
[`corrected_breakthrough_claim.md`](corrected_breakthrough_claim.md).

Reproduce everything with:

```bash
python scripts/verify_all_claims.py     # ~2 min, writes claim_recheck.csv + verification_log.txt
```

---

## A. Changed in place

| File | Was | Is now | Why |
|---|---|---|---|
| `README.md` badge | `pytest 342 passed` (not reproducible on a clean Linux clone) | `pytest 342 passed` (**reproducible** after three fixes) | The documented count was unreachable because `data/manifest/data_manifest.csv` was missing (7 failures / 9 errors), one test opened `data/raw/7tre.cif` (lower-case) on a case-sensitive filesystem (3 errors), and the Biopython fallback in `compute_bsa` ignored `use_auth=False` (1 failure). Fixed in this pass; `PYTHONPATH=src pytest -q` now reports 342 passed. |
| `README.md` header | — | Verification-status banner + prior-art caution on the K757A section | The "separation-of-function" proposal is contradicted by Sekne 2022 and Tesmer 2023. |
| `report.md` §3 citations | Liu 604:578–585; Sekne 375:1173–1179; Wang *NSMB* 14:1011–1017 | Liu 604:578–583; Sekne 375:1173–1176; Wang *Nature* 445:506–510 | Citation metadata corrected against the publisher records. |
| `report.md` Table 5.1 (2I46) | 1.85 / 2.70 Å; chain A 90–241 | 2.70 Å; chain A 90–301, chain B 90–287 | `_reflns.d_resolution_high = 2.7`; model ranges read from the CIF. |
| `report.md` Table 3.1 (Lys78–Glu215) | 2.89 Å / 3.58 Å | 3.584 Å (NZ–OE2) | No atom pair in that residue pair is 2.89 Å apart. |
| `report.md` Table 6.1 (Lys78, Leu769, Pro771) | 2.89, 3.52 (CD2–CD1), 3.86 (CG–CD1) | 3.584; 3.665 (CD2–CD2); 3.524 (CG–CD2) | Atom-pair labels were wrong; the distances belong to different atom pairs. |
| `report.md` §8 (Glu168 shift) | 2.748 Å | 1.253 Å | Re-derived after sequence-verified superposition of 152 common Cα (1.966 Å global RMSD). Matches the project's own `evidence_ledger.csv` CLM-018. |
| `report.md` §9 | "Jansson et al., 2015" | Wang 2007; Latrick & Cech 2010 | Unverifiable citation replaced with the primary sources for TPP1-stimulated RAP. |
| `results/breakthrough_investigation/FINAL_RESEARCH_REPORT.md` §3.1 | four SHA-256 values; 2I46 = 1.50 Å | real `sha256sum` values; 2I46 = 2.70 Å | The printed hashes did not match the shipped files. |
| `results/.../FINAL_RESEARCH_REPORT.md` §4.1 | "Latrick & Cech (*NSMB*, 2010)" | *EMBO J.* 29:924–933 (2010) | Journal corrected. |
| `tests/test_adversarial_b2_challenger2.py` | `data/raw/7tre.cif` | `data/raw/7TRE.cif` | Case-sensitivity portability bug. |
| `data/manifest/data_manifest.csv` | absent | added (identical to root `data_manifest.csv`) | `tests/test_validate.py` and `tests/test_adversarial_m1.py` require both copies; a clean clone failed 7 tests and 9 fixtures. |
| `src/telotrade/interface_geometry.py` (`compute_bsa`) | Biopython fallback indexed chains by author id even when `use_auth=False` | resolves the selected chains to author ids first | `test_compute_bsa_label_chains` returned 331.65 Å² instead of 856.82 Å² whenever freeSASA was unavailable. |
| `tests/test_claim_verification.py` | — | added (10 tests) | Pins the SHA-256s, the 7TRD ≡ 7TRF identity, the K757A-file identity, the TELIN-1 fragment count, the two anchor distances and the deposited resolutions, so the correction documents cannot silently rot. |
| `README.md` test badge | 342 passed (unreachable) | 352 passed (measured) | 342 + the 10 new guard tests, all passing. |

## B. Flagged, not changed

| Item | Finding |
|---|---|
| `results/tables/tel_patch_contacts.csv` | Row order is non-deterministic (unstable sort over tied distances): regenerating it produces a spurious git diff. |
| `src/telotrade/interface_geometry.py::compute_bsa` | With `use_auth=False` the freeSASA branch is label-correct but the Biopython fallback indexes by author chain, returning 331.65 Å² instead of 856.82 Å². `test_compute_bsa_label_chains` only passes when freeSASA is installed. |
| `requirements.txt` → `freesasa` | Listed as a requirement but fails to build in a bare environment; several BSA tests then silently exercise a different code path. |
| `data/raw/*.cif` provenance | Cannot be confirmed offline (RCSB unreachable from this sandbox). Per-residue ADP grouping is consistent with group-B refinement and is not evidence of tampering, but authenticity must be re-established by a fresh download + hash check before publication. |
| `results/experimental_validation/spr_benchmark_results.csv`, `rap_benchmark_results.csv` | Generated from `scripts/generate_independent_benchmarks.py` (`np.random.seed(42)`, hand-set k_on/k_off). They validate the fitting pipeline, not any biology. The "WT reference Kd = 150 ± 25 nM" has no experimental source. |
| `data/trajectories/*.csv` | 2.36 ps production windows on a 97-residue, four-chain fragment with no hTR and no DNA. `data/md_replicate_summary.txt` (3.496/5.296 Å) and `RECALCULATED_REPLICA_STATISTICS.csv` (frame means 3.820/5.620 Å) disagree with each other. |
| `data/junction_bound_k757a_fixed.pdb` | Byte-identical to `data/junction_bound_fixed.pdb` (MD5 `1b284b2f…`); residue 757 is LYS. No K757A system was ever built. |
| `data/ligands/telin1_inhibitor.pdb` | Two disconnected fragments; not a single chemical entity. |
| `results/tables/` "shape complementarity" | No Sc value exists anywhere in `results/`; the claim in the aims was never computed. |
| `scripts/audit_restriction_sites.py` | Fetches NM_198253.3 from NCBI at run time (not reproducible offline), assumes codon 757 = AAG instead of verifying it, and tests a different enzyme panel (BssHII, NgoMIV, NaeI, …) than the README's claimed diagnostics (BstUI, HinP1I, MspI/HpaII). The diagnostic-site claims are therefore unverified. |

## C. Retracted claims (do not cite)

1. "TERT Lys757A is a separation-of-function mutant that preserves basal catalysis" — contradicted by
   Sekne et al., *Science* 375:1173–1176 (2022), who assayed K757A/F759A, and by Tesmer et al.,
   *Life Sci. Alliance* 6:e202201727 (2023), who assayed K757E/F759E.
2. "Tyr176–Gln794 is a state-dependent allosteric switch" — these are the paper's DNA-anchor
   residues (Sekne 2022); the B-factor evidence is invalid across entries (7TRE mean B 16.5 Å² vs
   7TRD 85.9 Å²; Lys757 ratio 0.24 ≈ its domain's 0.26).
3. "TPP1-induced 67–82 % B-factor quench at the DNA exit pore" — global ADP-scale artefact (§above).
4. "MD validates the bipartite clamp / the Y176–Q794 latch" — 2.36 ps on a 97-residue fragment
   without nucleic acid; no K757A system exists.
5. "7TRD and 7TRF independently corroborate the unbound state" — they are the same coordinates
   (0.0005 Å).
6. "TELIN-1 is a computational lead inhibitor" — not a single molecule; 160 clashes down to 0.17 Å;
   no affinity, no synthesis, no literature.
7. "SPR/RAP benchmarks validate the mutant phenotypes" — the data are synthetic.
8. "0.338 Å RMSD proves the absence of allostery" — already retracted by `report.md` itself; restated
   here so it is not resurrected: the value is below the 0.35–0.5 Å coordinate uncertainty of these
   maps and can neither demonstrate nor exclude allostery.

---

## Second pass — 2026-11-01: the identifiability phase

This second errata covers the affinity–processivity investigation
([`results/causality/`](../causality/README.md)) and the manuscript
([`paper/MANUSCRIPT.md`](../../paper/MANUSCRIPT.md)). As with the first pass, the list is of
**errors in our own work**, found by writing the check for each claim rather than by re-reading the
prose.

| # | claim as first written | why it was wrong | corrected to |
|---|---|---|---|
| 9 | "C1a (*k*_off ↓) and C1b (*k*_on ↑) are **identical** at matched *K*_d" | the check written to verify it returned `false`. The loss term is *R* = *k*_off²/(*k*_off + [P]*k*_on), so at matched *K*_d the two differ by an **exact factor λ** | they are not equivalent; the factor is algebraically confounded with *k*_loss |
| 10 | the displayed closed form for *P*\* | the numerator had the **opposite sign**, turning every real crossing into a non-crossing | sign corrected; now checked against `sympy.solve` and a 20 001-point scan (residual 1.1×10⁻¹⁶) |
| 11 | "the mixture bias at saturation is −0.040" | not reproducible at any band window; the true value at *s* = 0.999 is −0.0014 | saturation bias is ≈ 0; the real objection is the **window dependence** (12-fold in implied *k*_loss/*k*_ext) |
| 12 | Δ*s* ≈ φ(1−*s*) as the mixture bias | wrong functional form: omits (1−*s*₂), linear in φ where the truth is not, and cannot change sign | retired; use *r*(*n*) = *s*(1+κρⁿ)/(1+κρⁿ⁻¹) |
| 13 | "direct binding improves practical identifiability to ~10% for all parameters" (cond. no. 25) | neither the condition number nor the *k*_loss figure was reproducible | cond. no. **7.5**; *k*_off/*k*_on 9.9%, *k*_loss gets **worse** (34.3% → 49.0%) |
| 14 | "the mixture mean manufactures a C2 signal" | the mixture *mean* stays at 1.0 at saturation; only the *tail-slope estimator* fabricates C2 | estimator-specific, not a property of the distribution |
| 15 | *k*_off/*k*_ext relative SE quoted as 13% | 12.4932% rounds to **12%** | 12% (also corrected in the manuscript) |

Items 9–12 and 14 are also recorded in the manuscript's §13 and in
[`results/causality/MECHANISM_CLASSIFICATION.md`](../causality/MECHANISM_CLASSIFICATION.md) §0.1, so
that a reader who arrives by any route meets the correction rather than the original claim.

**Not an error, but withdrawn as evidence:** the "three independent scales" argument for
single-turnover incorporation. Adding that measurement to the design adds rank **only** under the
assumption *k*_chem ≡ *k*_ext; with a free coupling parameter the design gains one observable and
one parameter and the rank gain is zero.
