# WP31 — The missing part was on disk. Assembled it, ran every ligand in all of them

WP30 reported the TPP1 epitope as **blocked**. That was wrong, and this corrects
it: TPP1 was in the repository the whole time.

| file | TPP1 | residues |
|---|---|---|
| `data/raw/7TRE.cif` | **chain C** | 152 (auth 90–241) |
| `data/raw/2I46.cif` | **chain A** | 152, TEL-patch fully resolved |
| `data/ligands/7tre_telin1_complex.pdb` | **absent** | — |

Only the *docked complex* had TPP1 stripped out. So the gap was an
**ingestion** gap, not a data gap, and no network fetch was ever needed.

    python scripts/solve/wp31_epitope.py
    -> data/structures/telomerase_tpp1_assembled.pdb
    -> results/solve/wp31_epitope_scores.csv

## 1. The assembly, and why it is not invented

7TRE's chain C and 2I46's chain A are the **same 152 residues with identical
sequence**. Superposing 2I46's TPP1 onto 7TRE's TPP1 gives

**RMSD 1.966 Å over 152/152 Cα.**

That number is the whole justification for the merge. Two independent crystals
place TPP1 in the same position to within ~2 Å, so dropping 2I46's better-resolved
copy into 7TRE's TERT frame is a rigid transfer, not a guess. If it had been
high, the merge would have been refused.

> **A bug worth recording.** The first implementation of the Kabsch fit mapped
> 7TRE → 2I46, the wrong direction, which silently put TPP1 in the wrong frame
> with **no error at all**. It was caught only because the interface came out
> empty. The self-test now asserts that the transformed TPP1 actually sits near
> TERT (Cα–Cα < 8 Å), so that failure mode cannot pass again.

## 2. The epitope, now defined

At a 5 Å heavy-atom cutoff: **26 TERT residues** and **21 TPP1 residues**.

TPP1 interface residues (auth): `90-95, 112, 114, 166-170, 180, 183, 185, 211-215`

That set contains the **TEL patch** — the region the ligand claims to disrupt —
and it is the difference between docking against a defined interface and docking
against a protein surface, which is what WP30 objected to.

## 3. Every ligand run against every target

Targets: whole TERT (the WP30 target, 7,488 atoms), the TERT–TPP1 interface
(205 atoms), and the TPP1 side alone (172 atoms).

| ligand | target | contacts | clashes | d_min | score |
|---|---|---|---|---|---|
| TELIN-1 core | whole TERT | 74 | 0 | 2.64 | 74.0 |
| TELIN-1 core | **TERT–TPP1 interface** | 64 | 0 | 2.93 | **64.0** |
| TELIN-1 core | TPP1 side | 54 | 0 | 2.69 | 54.0 |
| TLN-1 repaired | whole TERT | 47 | 0 | 2.78 | 47.0 |
| TLN-1 repaired | **TERT–TPP1 interface** | 64 | 0 | 2.67 | **64.0** |
| TLN-1 repaired | TPP1 side | 58 | 0 | 2.84 | 58.0 |
| TLN-3 one anchor | whole TERT | 60 | 0 | — | 60.0 |
| **TLN-3 one anchor** | **TERT–TPP1 interface** | **88** | 0 | — | **88.0** |
| TLN-3 one anchor | TPP1 side | 65 | 0 | — | 65.0 |
| TLN-4 no tail | whole TERT | 71 | 0 | — | 71.0 |
| TLN-4 no tail | TERT–TPP1 interface | 61 | 0 | — | 61.0 |
| TLN-4 no tail | TPP1 side | 56 | 0 | — | 56.0 |

Every pose is **clash-free** — poses with any heavy-atom pair under 2.6 Å are
rejected outright rather than merely penalised, which is the fix that made the
docking behave.

## 4. What the deltas actually say

| ligand | whole-TERT | interface | Δ |
|---|---|---|---|
| TELIN-1 core | 74.0 | 64.0 | **−10.0** |
| TLN-1 repaired | 47.0 | 64.0 | **+17.0** |
| TLN-3 one anchor | 60.0 | 88.0 | **+28.0** |
| TLN-4 no tail | 71.0 | 61.0 | **−10.0** |

A **positive** Δ means the molecule prefers the interface over the rest of the
protein surface — which is the shape of a selective molecule. A **negative** Δ
means it would rather stick somewhere else, and selectivity is not there.

> **RETRACTED (WP32).** This claim does not survive a scorer with charges in
> it. See `WP32_PHYSICS.md`: all four ligands prefer a control patch, TLN-3 is
> no longer first, and the "convergence" was an artefact of counting contacts.
> The paragraph is kept so the retraction has something to point at.

**TLN-3 is the only ligand in the series with a large positive Δ.** Two
independent changes pointed the same way: the physchem screen in WP30 picked it
out (MW 442.51, logP 4.60, one carboxylate, passes all four Lipinski rules), and
the structural run independently prefers it at the interface by the widest
margin. That convergence is the most interesting result in the programme so far,
and it is worth exactly what a crude rigid-body surrogate is worth:

- **No hydrogen bonds.** A salt bridge between the anchor carboxylate and a
  TPP1 lysine is the entire proposed mechanism, and this scorer cannot see it.
- **No electrostatics, no desolvation, no flexibility, no water.** TPP1's TEL
  patch is a shallow, flat surface; real affinity there is dominated by exactly
  the terms that are absent.
- **Low resolution.** A contact count resolves "interface vs not interface",
  not "TLN-1 vs TLN-3".

So the honest statement is: **the design that passed the physchem screen also
fails to prefer anywhere else, and that is consistent with — not evidence for —
selectivity.** It is a reason to run a real calculation, not a finding.

## 5. Scoreboard

| | state |
|---|---|
| TPP1 present in the model | **fixed** — 152 residues, assembled, RMSD 1.97 Å |
| epitope defined | **fixed** — 26 TERT + 21 TPP1 residues at 5 Å |
| every ligand run in all targets | **done** — 12 poses, 0 clashes |
| selectivity signal | **TLN-3, Δ +28** — surrogate only |
| any affinity | **none, and none obtainable here** |
| real free-energy calculation | **not done** — needs an all-atom FEP/MM-GBSA run |

The remaining gap is now **one specific calculation**, not a missing structure:
an MM-GBSA or FEP treatment of TLN-3 at the TEL patch. That is the natural next
step and it is out of scope here.

**Quarantine unchanged.** TELIN-1 remains an unvalidated in-silico proposal
with zero empirical evidence. A clash-free pose is not a binding measurement,
and a contact-count Δ is not an affinity.

## 6. Tests

22 checks, 0 failed — including that Kabsch recovers a known rigid transform
exactly, that both files really carry ≥100 TPP1 residues, that the two copies
superpose below 3 Å, that the interface is a *minority* of TERT (i.e. it is
actually an interface and not the whole protein), that **the assembled TPP1 sits
in the TERT frame**, that the scorer counts a 3.9 Å contact and penalises a
0.5 Å clash with a negative score, and that the placed poses are clash-free.
