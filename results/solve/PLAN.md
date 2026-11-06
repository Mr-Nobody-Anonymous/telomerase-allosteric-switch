# Plan — resolving the TPP1 processivity mechanism with what is computable here

**Goal (fixed, not to be changed):** determine whether TPP1's stimulation of telomerase
repeat-addition processivity is *mechanical/allosteric* or *DNA retention*.

**Constraint discovered:** this sandbox has PyPI + GitHub egress only. `files.rcsb.org`,
`rest.uniprot.org`, `eutils.ncbi.nlm.nih.gov` are **unreachable from bash** (OpenSSL
`SSL_ERROR_SYSCALL`). The `fetch_page` tool *can* reach RCSB, so individual records are
retrievable but multi-megabyte coordinate files are not practical to pull. Everything below is
therefore designed to be answerable from **local data + published numbers + algebra**.

---

## Work packages, in dependency order

| # | WP | Question | Method | Needs new data? |
|---|---|---|---|---|
| **1** | Geometry of the anchor site | Is the repo's "TEN makes zero primer contacts (Y176 = 11.1 Å)" a real finding or an artefact of the 7-nt DNA in 7TRE? | Measure the DNA exit path, per-residue distances to PLYQ / K757 / IFD β17, in all three DNA-containing depositions | No — local mmCIF |
| **2** | Model-class audit | Is the repo's kinetic model the one that describes the assay? | Re-derive `s`; test the model's own prediction against the *saturating-primer* conditions of Latrick & Cech 2010 | No — algebra + published values |
| **3** | Correct identifiability | Given a correctly specified model, what **is** identifiable, and what is the minimal decisive experiment? | Symbolic rank, Fisher information, profile likelihood, optimal design | No |
| **4** | Structural coupling | Is there a *short* through-protein path from the TPP1 footprint to the DNA anchor site, and is it shorter than chance? | Contact-graph shortest paths + ANM with permutation nulls + region-matched controls | No — local 7TRE/7TRD |
| **5** | Literature discrimination | Do already-published measurements already falsify one of the models? | Formal test of model B against the measured k_off / k_pol / translocation numbers | No |
| **6** | Synthesis | What is resolved, what is still open, what truly needs a lab or a supercomputer? | — | No |

## Falsification discipline (applies to every WP)

* Every number gets a null or a control. The repo's `CLAIM_VERIFICATION_REPORT.md` found 14
  contradicted claims, most from cross-entry comparisons without controls. Do not repeat that.
* A result that only holds under one arbitrary parameter choice is reported as sensitivity, not as
  a finding.
* If a WP is refuted, the WP is rewritten, not reinterpreted.

## Status log

| WP | state | outcome |
|---|---|---|
| 1 | **done** | the "zero TEN primer contacts / 11.1 A" claim measured a 7-nt DNA that ends at the active site; PLYQ is the closest TEN element to it. Corrected. |
| 2 | **done** | eq.1 falsified twice: affinity channel closed at [P]/K_d=100 (0.99% of loss), and no state variable can express the unchanged-resting / 4x-slower-cycling k_off. |
| 3 | **done** | step-resolved retention model predicts 2.48x vs observed 2.7 +/- 0.6x (-0.37 SD), no free parameters. |
| 4 | **done** | TPP1 footprint 35-42 A from catalytic site; contact-graph paths not discriminating (reported as null). |
| 5 | **done** | decisive experiment specified (single-turnover k_chem, no primer) with power and a pre-declared rule. |
| 6 | **done** | mutant panel does not discriminate - null, sign unstable. |
| 7 | **done (added)** | 7TRD and 7TRF share one TERT model; cross-state evidence has one independent apo set. |
| 8 | **done (added, corrects WP1)** | DNA polarity: the 3' end is on the catalytic aspartates and the 5' end on the anchor strand R756; the 5' exit path runs into Tyr176 at 6.8 A. The 11.1 A PLYQ distance measures a site the 7-mer stops short of. |
| 15 | **done** | design validation: heterogeneity bounds PROVEN not sampled (R in [1.2903,1.3041] at A0=0.70); w not identifiable from A and the hazard is a passing gate over an unmeasurable subpopulation; trap needs 2 methods x 5 readouts; hierarchical power (tau=0.10 -> P 7->20 at A0=0.90). |
| 14 | **done** | robustness: z-formula sample sizes corrected by exact-t (6->8 at A=0.90); A<=0.70 gate SURVIVES position-dependent k_pol but FAILS for mixed preparations (needs active fraction w>=0.8); trap moves both baseline and effect size and must be validated; preregistration written. |
| 13 | **done (sample sizes corrected by WP14)** | identifiability + operating range: sensitivity map (A depends only on k_off_ext/k_pol), sample sizes (A=0.99 needs ~700-2700 replicates), Clopper-Pearson intervals, decision rule with INCONCLUSIVE branch, and the first experiment = measure A untreated. A near-1 value is uninformative, not exclusionary. |
| 11 | **done** | assumption audit: k_tr-invariance survives A1-A5; the identity SPLITS into extension (A=s/phi) and translocation (B=phi/(1-phi)) channels; the inversion is an aggregate (+/-2x), not a rate ratio. |
| 12 | **done** | power + uncertainty: design has power only if k_pol is slow (unknown); phi->1 costs 50-101% error at phi=0.99; one extra k_off_ext measurement resolves the extension phase exactly but leaves translocation unresolved. |
| 10 | **done (added, the decisive design)** | exact identity log(s)-log(phi) = -6log(1+k_off/k_pol), invariant to k_tr. Each of M1-M4 gets a unique (s, phi) signature; k_off/k_pol recoverable with no fitting. Needs a matched (s, phi) pair, which no published dataset reports. |
| 9 | **done (added, reframes the question)** | G100 (Zaug 2010: required for stimulation) is 7.4 A from TPP1 but 24.3 A from DNA; Y176 (DNA anchor) is 19.2 A from TPP1. Opposite faces of TEN. TEN is tethered on an unmodelled 124-residue linker. "Retention" splits into R1 (excluded: TPP1 makes no DNA contact) and R2 (supported: TPP1 sets TEN's orientation). |

**Outcome (revised).** R1 — TPP1 changing the enzyme's DNA affinity directly — is excluded: TPP1 is
28.7 A from the DNA and makes no contact. R3 — TPP1 changing catalysis — is unsupported: the active
site is 35.8-42.1 A away across the RT ring. **R2 — TPP1 acting on TEN's orientation, which changes
retention — is supported and is the live hypothesis.** Allostery is not excluded by any published
number, but no number requires it. See RESULTS.md and CORRECTIONS.md.

**Milestone:** WP11-WP15 complete.

**Closed thread — anchor site (WP16).** The PLYQ–DNA contact has been measured in the Sekne
long-DNA depositions: **Tyr176 is 2.75 Å (7QXS) and 4.35 Å (7QXA) from auth nt 19**, both entries
passing a catalytic control. Against 11.11–11.94 Å in the 7-nt structures this is a fourfold
change. It is a backbone contact, not base-specific, and it settles geometry only — the mechanism
remains open. See `WP16_ANCHOR_SITE_RESULT.md`.

**Next structural experiment it raises, not yet done:** PLYQ cavity mutant **Y176A**, assayed for
processivity and for TPP1 stimulation. If the anchor is load-bearing, processivity should fall
while recruitment stays intact. This is falsifiable with the readouts already specified in
`WP15_DESIGN_VALIDATION.md` and needs no new structures. Mechanism unresolved, and now with a measured reason: the
discriminator works only inside an unmeasured operating regime, and the first experiment should
measure that regime (A untreated, one reaction, no TPP1) rather than the mechanism.

**Corrections made after external input** (recorded in CORRECTIONS.md, not silently edited):
(i) the proposed decisive assay (no-primer single-turnover) is physically impossible and is
withdrawn; (ii) the WP3 2.48x is calibration on the same data, not confirmation; (iii) "retention"
was used for two mechanisms and is now split; (iv) the remaining work is not one experiment.

**Plan changes during execution.** (i) A WP7 was added when the entry-independence check became
necessary; (ii) the first affinity-ceiling formulation in WP2 was wrong (it invented a ceiling
rather than deriving one) and was replaced with the exact closed form R = k_off^2/(k_off+k_on[P]);
(iii) a preliminary B-factor comparison between 7TRE and 7TRD was abandoned as unsafe once the
whole-chain medians were found to differ 20-fold, and the reason (shared model / differing
refinement) was traced instead; (iv) **WP1's original reading was wrong and was replaced.** The
first pass checked the distance from PLYQ to the DNA but never checked which end of the DNA was
which; it concluded the fragment pointed away from PLYQ. WP8 established the polarity — 3' end at
the active site, 5' end at the anchor strand — and showed the exit path runs into Tyr176. The
goal was unchanged; the intermediate conclusion was corrected.
