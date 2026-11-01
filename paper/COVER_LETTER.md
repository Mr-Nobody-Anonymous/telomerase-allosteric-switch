# Covering letter — template, not sent

> **This letter has not been sent, and no address in this file is a verified contact address.**
> It is a template. Before sending, the author must (i) confirm each recipient's current
> institutional address from the paper or the group's own website, and (ii) decide whether the
> recipients below are the right ones. Nothing here should be pasted into a mail client unread.

---

## Why this note exists

The manuscript describes a negative result about inference: telomerase processivity distributions
cannot separate the mechanisms currently proposed for TPP1, and a widely used summary statistic — the
tail slope — is window-dependent in a way that can fabricate the very effect a mechanism claim would
rest on. Negative methodological results are only useful if the people who generate the data see them
before the data are interpreted, which is why this is being circulated before publication rather than
after.

It also retracts three errors from its own earlier drafts, listed in §13. That is deliberate: a
methods paper that hides its own corrections is not worth the reader's attention.

## Draft text

> **Subject:** Preprint — what processivity distributions cannot identify about TPP1 (and a caveat on tail-slope estimates)
>
> Dear Professor [Name],
>
> I am writing to share a short analysis that may be relevant to your work on telomerase
> processivity, in case it is useful before it goes any further.
>
> The question it addresses is narrow: given a product-length distribution, what can be inferred
> about whether TPP1 changes primer retention, extension rate, or both? The answer is less than is
> usually assumed. Within a standard per-repeat kinetic model, an extension-rate mechanism and an
> affinity-retention mechanism are the same function of the observable, so under matched effective
> survival they produce identical distributions; and no absolute rate is recoverable from
> processivity data at any precision, only ratios.
>
> The part that may be of more practical use concerns heterogeneous populations. For a two-species
> mixture the adjacent-band ratio has a closed form, and three things follow: the ratio is bounded by
> the most extreme subpopulation; its tail converges to the longest-lived species rather than the
> majority one, so the bias is positive for a faster contaminant and negative for a slower one; and a
> minority component takes over the signal beyond a computable length. In a worked case a true
> processivity of 0.60 reads as 0.61 or 0.95 depending only on which band window is used — a 12-fold
> difference in the implied concentration-independent loss, from the same data. If that holds for
> real gels, tail-slope estimates are not comparable across papers that read different windows.
>
> I want to be clear about what this is not. No mechanism is demonstrated, and I do not claim that
> TPP1 acts through retention, through catalysis, or through anything else. The results are exact
> within a stated model and say nothing about biology until tested against real densitometry, which
> is the next step and needs no new experiment. The paper also corrects three errors from its own
> earlier drafts, including a claim that the two affinity variants are equivalent at matched K_d —
> they are not.
>
> Everything is reproducible from the repository, including the figure generation and a test suite
> that fails if any quoted number drifts from the code:
>
> https://github.com/Mr-Nobody-Anonymous/telomerase-allosteric-switch
> (paper: `paper/MANUSCRIPT.md`; analysis: `results/causality/MECHANISM_CLASSIFICATION.md`)
>
> If you think the mixture result is already established somewhere in the assay literature, I would
> genuinely like to know — a prior derivation would improve the paper and would not weaken it. And if
> you have a gel with recoverable per-band intensities, the falsification test in §10 is a
> densitometry exercise, not a new experiment.
>
> With thanks for your time,
>
> Beamlak Sintayehu

## Who this is relevant to, and why

Identified from the cited literature. **These are not verified contact addresses** — they are the
people whose published work the argument bears on, so that the author can decide whom to approach.

| researcher | group / affiliation as given in the cited work | why the note is relevant |
|---|---|---|
| Jayakrishnan Nandakumar | University of Michigan | corresponding author of the IFD–TRAP / TEN–TPP1 interaction work (ref. 2) whose processivity and activity values are deposited in this repository |
| Ci Ji Lim | University of Wisconsin–Madison | corresponding author of the RPA processivity-factor work (ref. 5); the identifiability argument applies to that system as well |
| Thomas R. Cech | University of Colorado Boulder / HHMI | senior author of the kinetic study arguing for slowed primer dissociation plus aided translocation (ref. 6) — the paper's Result 1 is consistent with both halves of that account |
| Thi Hoang Duong (T.H.D.) Nguyen | MRC Laboratory of Molecular Biology | senior author of the TPP1–POT1 telomerase recruitment structure (ref. 3) |
| Juli Feigon | UCLA | senior author of the active telomerase–TPP1 structure (ref. 4) |
| Kathleen Collins | UC Berkeley | telomerase biochemistry; relevant to the single-turnover versus repeat-cycle distinction (§7) |
| Titia de Lange | Rockefeller University | shelterin and telomerase recruitment; relevant to the assay-design consequences |

The **two most likely to find it useful** are the Nandakumar and Lim groups, because both measure
processivity distributions in the regime the mixture result describes, and both have deposited or
could readily produce per-band data. The Cech group is the most likely to know whether Eq. 3 is
already standard in the assay literature, which is the one question that would most improve the paper.

## Note on tone, if the author edits this

Two things should stay if this is sent, because they are what makes it credible rather than
promotional: the explicit statement that no mechanism is demonstrated, and the invitation to be shown
that the mixture result is not new. The letter should also keep the repository link — with the
command to regenerate every number, the reader can check the central claim in about two seconds,
which is worth more than any amount of framing.
