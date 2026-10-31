# The Bottleneck: Does TPP1 Act at a Distance, or Where It Binds?

**Date:** 2026-10-30 · **Status:** negative result + reframed hypothesis + decisive experiment
**Reproduce:** `python scripts/allosteric_coupling.py` (≈40 s)
**Outputs:** `results/allostery/` — `anm_summary.json`, `system_comparison.csv`,
`cutoff_robustness.csv`, `tpp1_distance_gradient.csv`, `residue_coupling_ranked.csv`,
`analysis_log.txt`

---

## 1. What the bottleneck actually was

The project's whole mechanistic story rested on one sentence: *TPP1 docks 25 Å from the DNA exit
channel, so its effect on processivity must be allosteric.* Every downstream claim — the
"Tyr176–Gln794 switch", "Lys757A separation-of-function", the GNM "entropic clamp" — is a
consequence of that premise.

So the bottleneck is not "which residue is the hinge". It is: **is the premise true?**

Three things had to be checked, in order:

| # | Question | Answer |
|---|---|---|
| 1 | Is TPP1 really far from the functional site? | **No — the premise is misleading.** TPP1 binds the same helix that carries the DNA anchor, with a continuous 2.3 → 26 Å gradient of contacts along it (§2). |
| 2 | Does the mechanical coupling survive a proper null test? | **No.** A distance-matched null shows no coupling beyond what proximity predicts (z = +0.34), with the negative control behaving correctly (§3). |
| 3 | Is there a "two-hand" clamp spanning two DNA-binding modules? | **No, not in these coordinates.** In 7TRE the TEN domain makes *no* contact with the modelled DNA (Tyr176 is 11.1 Å from the primer) — there is only one hand on the DNA (§4). |

---

## 2. Structural fact: TPP1 sits on the anchor helix, not across a void

Heavy-atom distances from TPP1 to each TERT residue (`tpp1_distance_gradient.csv`) give a
monotonic gradient along the IFD-TRAP helix:

| Residue | TERT 756 Arg | 757 Lys | 758 Ala | 759 Phe | 760 Lys | 761 Ser | 762 His |
|---|---|---|---|---|---|---|---|
| d to TPP1 (Å) | 26.0 | **25.4** | 22.0 | 20.0 | 17.6 | 15.2 | 13.0 |

| Residue | 763 Val | 764 Ser | 765 Thr | **766 Leu** | 767 Thr | 769 Leu | 771 Pro | 772 Tyr | **774 Arg** |
|---|---|---|---|---|---|---|---|---|---|
| d to TPP1 (Å) | 13.6 | 9.9 | 7.2 | **3.6** | 4.9 | 3.7 | 3.5 | 3.4 | **2.3** |

* The 21-residue footprint (≤4.5 Å) is 13 residues on the TEN surface (43–46, 78, 121–136) and
  8 on the IFD-TRAP helix (766, 769, 771, 772, 774, 775, 798, 799).
* TPP1 → primer DNA minimum = **28.7 Å** (SER90:OG–DT12:OP1); TPP1 → hTR minimum = **14.2 Å**.
* Backbone path 757 → 766 = **34.2 Å**; through-space distance from TPP1 to 766 = **3.6 Å**.

**Interpretation.** The 25 Å figure is real but it measures the wrong thing. It is the distance
from TPP1 to the *distal end* of a continuous helix that TPP1 is sitting on. There is no void for
a signal to cross: the anchor residues (756/757/759) and the TPP1 footprint (766–799) are two
segments of one helix bundle. Whatever TPP1 does to the anchor site, it does **locally**, by
binding one face of the helix — the class of mechanism that helix-propagation and side-chain
repacking describe, not the class that "long-range allostery" describes.

## 3. Mechanical test: an elastic network with a proper null says no

Design (`scripts/allosteric_coupling.py`): anisotropic network model (Tirion springs, 15 Å
cutoff, 40 lowest non-trivial modes) on **identical** TERT/hTR/DNA coordinates, with three
systems that differ only in TPP1:

* **APO** — TPP1 deleted (1158 beads)
* **BOUND** — TPP1 at its real site (1310 beads)
* **DECOY** — TPP1 translated 250 Å away: same bead count, same internal elasticity, no contact

Observables: coupling index *K* between module centroid motions (1 = rigid staple, 0 =
independent); distance-normalised correlation residual; site-to-site coupling against a
**distance-matched permutation null**; per-residue mobility.

| Observable | APO | BOUND | DECOY | Verdict |
|---|---|---|---|---|
| *K*(TEN 7–196, IFD-TRAP 740–800) | +0.521 | +0.395 | +0.436 | ligation effect **−0.04** (non-specific component −0.09) |
| *K*(docking surface, primer channel) | +0.591 | −0.054 | +0.190 | no clamp |
| Coupling dock → primer channel (residual) | +0.0335 | +0.0489 | +0.0105 | z = **+0.34** vs null (+0.044 ± 0.015) |
| **Positive control:** pairs < 8 Å apart | +0.048 | +0.041 | +0.056 | metric works locally |
| **Positive control:** catalytic site ↔ DNA 3′ end, \|corr\| | 0.333 | 0.103 | 0.193 | barely above inter-chain background (0.203 / 0.087 / 0.102) |

Cutoff robustness of the staple metric: **18 Å → −0.051**, **15 Å → −0.041** (consistent);
12 Å is unusable (11–17 null modes because the bead model is under-constrained there).

**Reading it honestly, in both directions:**

1. The **H_switch** prediction — "docking at 25 Å changes the mechanics of the primer channel
   more than proximity alone allows" — **fails**: the observed coupling is indistinguishable from
   random residue pairs at the same distance, and the decoy is *below* its own null.
2. The **H_clamp** prediction — "TPP1 rigidly staples two DNA-binding modules" — **fails**: *K*
   changes by −0.04 once the non-specific effect of adding a body is subtracted.
3. **But** the same model can barely distinguish an essential link (catalytic site ↔ primer,
   \|corr\| = 0.10 vs 0.09 background). So this is *not* evidence that no coupling exists — it is
   evidence that **this class of model cannot adjudicate the question**, which is exactly what
   the repository's earlier GNM-based claims ("entropic clamp", "vibrational anchor") failed to
   acknowledge.
4. The mobility result is also a caution: TPP1 binding reduces ANM fluctuations most in
   RT 541–739 (ratio 2.40) and least in TEN (1.42) and CTE (1.38) — a **global** stiffening, not
   a localised one, which is the same failure mode the B-factor analysis had.

Validation of the residue ranking is weak and should not be over-sold: the published
loss-of-function positions score AUC 0.748 against all other TERT residues, but several of them
were chosen *because* they contact DNA, so a positive AUC is expected.

## 4. The clamp model also fails on geometry alone

The "two-hand retention" idea requires TERT to hold the primer with two modules whose relative
motion TPP1 locks. In the deposited coordinates:

* TERT residues ≤ 6 Å from the modelled primer: **37 residues, all within the RT/IFD-TRAP core**
  (329, 500–502, 570–571, 631, 681–683, 712, 717, 756–757, 839, 866–869, 930–932, 944–949,
  957–960, 972–981, 1007, 1011).
* TEN-domain residues ≤ 6 Å from the primer: **none.** Tyr176 is 11.1 Å away.

So the TEN domain does not hold this DNA in 7TRE; only the RT/IFD-TRAP side does. The two-hand
clamp is not visible in the structure this repository is built on. (It may still operate on a
longer 5′ overhang, which 7TRE does not model — that is a hypothesis for a new construct, not a
finding.)

## 5. The reframed hypothesis: the helix conduit

> **H_conduit.** TPP1 stimulates processivity by binding one face of the IFD-TRAP helix that
> carries the DNA anchor (residues 756–759), propagating a local conformational/electrostatic
> change along residues 760–766 into the primer exit channel. It does not require long-range
> propagation across the molecule.

This is a different claim from the one the repository published, and it makes a different
prediction:

* **If H_conduit is right:** breaking the continuity between the TPP1 footprint and the anchor
  (residues 760–766) should reduce χ_RAP (TPP1 fold-stimulation) with **recruitment (K_d) intact**,
  because the TPP1 docking surface (766–799 + TEN face) is untouched.
* **If H_conduit is wrong** (recruitment-only model): every mutant in 760–766 that preserves
  structure will keep both K_d and χ_RAP, and χ_RAP will track K_d across the whole mutant panel.

**The discriminating mutant panel** (all single point mutants, TERT background, same
reconstitution):

| Mutant | Position relative to TPP1 footprint | Prediction under H_conduit |
|---|---|---|
| K760A / S761A | 17.6 / 15.2 Å, at the boundary | partial loss of χ_RAP |
| H762A / V763A | 13.0 / 13.6 Å, mid-conduit | partial, additive with above |
| S764A / T765A | 9.9 / 7.2 Å, proximal | weak effect |
| **L766P** (helix break) | 3.6 Å, at the footprint edge | largest χ_RAP loss; helix-breaking control |
| R774A | 2.3 Å, docking hotspot | **should** lose recruitment (control for the assay) |
| K757A/F759A | 25.4 Å, anchor itself | loses basal activity too — *published* (Sekne 2022) |

The powerful internal control is that K757A/F759A (already published) sits at the *far* end of
the helix and destroys basal catalysis, whereas the 760–766 mutants sit in the *middle* and
should spare it. If both classes behave identically, H_conduit is dead and the effect is simply
"anything that perturbs the exit channel".

## 6. What this does and does not achieve

**Achieved**
1. Killed the "25 Å allostery" premise with a measured distance gradient — the single most load-
   bearing assumption in the repository.
2. Showed the repository's GNM claims do not survive a distance-matched null (z = +0.34; the
   decoy is *below* its own null), and stated why a negative from this model class cannot settle
   the biology.
3. Showed the two-hand clamp is absent from the deposited coordinates (TEN has no DNA contact).
4. Turned "find the allosteric pathway" into a bounded, testable helix-conduit hypothesis with a
   discriminating mutant panel whose internal control is an already-published result.

**Not achieved**
* No energetics, no dynamics, no experimental data. The ANM is harmonic and single-structure;
  it cannot see anharmonic/entropic coupling, explicit-solvent effects, or nucleic-acid-specific
  flexibility.
* The ranking of individual residues is not reliable enough to nominate a single "hinge".
* Nothing here establishes a therapeutic window; that question is untouched.

## 7. What would actually settle it (in priority order)

1. **Re-analyse the published gels** (Sekne 2022 Fig. 4E–G; Tesmer 2023 Fig. 3) as χ_RAP at
   matched basal activity. Free, one day, and it either supports or kills the conduit premise
   before any construct is made.
2. **The 760–766 panel above** in the reconstitution system already specified in
   `results/experimental_validation/` — 7 mutants, one SPR run, one RAP gel series.
3. **All-atom MD that is actually an MD** — 100 ns+, explicit solvent, full holoenzyme with hTR
   and a 12-nt primer; the 2.4 ps fragment trajectories in this repository cannot be salvaged.
   Only then can dynamic coupling (or its absence) be measured rather than extrapolated.
4. **A longer-primer construct** (≥18 nt with 5′ overhang) to test the two-hand clamp that 7TRE
   cannot show.
