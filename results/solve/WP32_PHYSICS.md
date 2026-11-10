# WP32 — Editing the part that lacked the capability. It failed, and that is the result

WP31 said, in its own words, that its scoring function was blind to the things
that decide binding. This work package edits that part. It adds the missing
physics, and then it adds something the old function never had: **a control it
can fail.**

**It failed.** The main result here is negative, and it retires WP31's ranking.

    python scripts/solve/wp32_physics_scoring.py
    -> results/solve/wp32_physics_scores.csv

## 1. What was edited

| capability | WP31 | WP32 |
|---|---|---|
| partial charges | **none** | Gasteiger, plus formal charges |
| ionisation state | **none** | pH 7.4 microspecies rules |
| hydrogen bonds | **none** | distance **and** D–H···A angle ≥ 120° |
| desolvation | **none** | burial term (≥ 8 neighbours within 5 Å) |
| clash handling | −50, i.e. "less negative is worse" | **+40, genuinely repulsive** |
| control target | **none** | 3 patches, including one it must not prefer |

```
E = E_elec + E_hbond + E_burial + E_clash
```

## 2. Three real bugs were found by writing the tests first

These are worth recording because each one silently produced a wrong answer:

1. **Clashes were attractive.** The old penalty was *negative*, so a 1.5 Å
   overlap still totalled **−48**, better than a good pose at −7.3, because a
   bare Coulomb sum rewards getting closer without limit. Fixed by making the
   clash term positive **and** killing electrostatics inside the clash radius
   (at overlap the interaction is Pauli repulsion, not charge–charge).
2. **Hydrogen bonds had no direction.** A contact at 2.83 Å but 116.5° was
   being counted. Now rejected; only ≥ 120° counts.
3. **The ligands were never charged.** Gasteiger charges on a *neutral* molecule
   sum to ~0, so the carboxylate that is the entire proposed mechanism was **not
   −1 in the model at all**. Added ionisation: acid → carboxylate (−1), amine →
   ammonium (+1), guanidine → guanidinium (+1).

A fourth thing the tests caught is chemistry, not a bug: **deprotonating the
acid removes its only H-bond donor.** An anionic anchor can accept but can never
donate.

## 3. The control test — and the answer

Each ligand was placed against three patches: the true interface, a **control
patch on TERT ≥ 25 Å away**, and the DNA.

| ligand | E(interface) | E(control) | Δ | prefers |
|---|---|---|---|---|
| TELIN-1 core | −2.4 | −11.3 | **+8.9** | control |
| TLN-1 repaired | −1.7 | −16.0 | **+14.3** | control |
| TLN-3 one anchor | −1.5 | −20.1 | **+18.6** | control |
| TLN-4 no tail | −2.9 | −18.4 | **+15.5** | control |

**All four prefer the control patch.** The function is not resolving the
epitope at all.

Consequence: **WP31's ranking was an artefact.** Contact counting put TLN-3
first by the widest margin (+28); with charges in the model it is not first, and
every ligand prefers a control patch. "Two independent methods converge" was
coincidence. It is retracted here.

## 4. Why it fails — and the control is what makes this answerable

| patch | atoms | net charge | positive | negative |
|---|---|---|---|---|
| interface | 377 | **−17.89** | +14.40 | −32.29 |
| control | 224 | **−1.87** | +12.70 | −14.57 |

The interface carries **16 more negative charge** than the control. Every ligand
in the series is **anionic** (−1.8 to −3.0). An anionic ligand meeting an
anionic patch is repelled, and the model is right to say so.

This is not a modelling accident — it follows from the biology. The TERT–TPP1
recruitment face presents the TEL patch, and the pathogenic TEL-patch mutations
that cause dyskeratosis congenita map onto a **basic** surface. The interface is
acidic because the recruitment mechanism is.

## 5. The design consequence, tested

If the epitope is acidic, the anchor should be **basic**. Three anchors, same
scaffold and linker:

| anchor | formal q | E(interface) | E(control) | Δ | H-bonds |
|---|---|---|---|---|---|
| **TLN-3 carboxylate** | **−1** | 2.15 | −18.33 | **+20.48** | 0 |
| **TLN-5 primary amine** | **+1** | −3.37 | −6.69 | **+3.31** | 0 |
| **TLN-6 guanidinium** | **+1** | −5.37 | −13.58 | **+8.21** | 0 |

Switching the anchor from acid to amine cuts the anti-selectivity from **+20.5
to +3.3 — a 6.2× reduction**, in exactly the direction the electrostatics
predict. The guanidinium lands in between, consistent with its delocalised
charge.

**But none reaches a negative Δ.** No variant *prefers* the epitope. The model
can reject the wrong sign of anchor; it cannot design the right molecule.

That is the honest ceiling of this approach, and it is a real ceiling rather
than a failure of effort: reaching past it needs flexibility, explicit solvent
and configurational entropy — an FEP/MM-GBSA treatment — not a better
hand-rolled function.

## 6. What the licence is on WP31's "convergence"

WP31 reported that TLN-3 was the only ligand to pass the physchem screen *and*
the only one to prefer the interface, and called it the most interesting result
in the programme. That claim is now **withdrawn**. The physchem screen is
unchanged (MW 442.51, logP 4.60, one carboxylate — still passes all four
Lipinski rules), but the structural half of the argument does not survive the
better function.

What replaced it is weaker in one sense and stronger in another: weaker because
nothing is selective, stronger because the **sign of the anchor's charge is now
a testable design variable** with a measured effect.

## 7. Scoreboard

| | state |
|---|---|
| scorer has charges, H-bonds, burial | **fixed** |
| scorer has a control | **fixed** |
| scorer passes its control | **NO — all four prefer the control patch** |
| WP31 "convergence" claim | **retracted** |
| cause of the failure | **identified: the epitope is acidic, the ligands are anionic** |
| design consequence | **tested: acid → amine cuts anti-selectivity 6.2×** |
| any variant prefers the epitope | **no** |
| affinity | **none, and none obtainable here** |

## 8. Tests

26 checks, 0 failed. The ones that matter: a 3.0 Å donor–acceptor pair is a
contact at *exactly* 3.0 Å and not at 8.0 Å; a well-angled H-bond is counted and
a bent one at a similar distance is **not**; a 1.5 Å overlap is a hard clash with
**positive** energy and zero electrostatics; opposite charges attract and like
charges repel; a buried atom is detected; the carboxylic acid ionises to −1, the
amine to +1, the guanidine to +1, without changing the heavy-atom count; and the
control patch is ≥ 20 Å from the interface centre.

**Quarantine unchanged.** No affinity has been measured, none is measurable
here, and a lower Δ is not a binding measurement.
