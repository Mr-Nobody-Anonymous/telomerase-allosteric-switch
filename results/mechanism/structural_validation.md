# Work package 1 - structural validation of the TPP1-TERT interface

Every number below is recomputed from the deposited mmCIF files by `scripts/structural_validation.py` (gemmi + numpy only, no imports from this project's own analysis code, no input from this project's own tables). Where the result contradicts an earlier claim in the repository, the contradiction is stated rather than smoothed over.

## 1. Chain identity and numbering

| PDB | author chain | label chain | entity type | residues | residue range |
|---|---|---|---|---|---|
| 7TRE | B | A | Rna | 217 | 33-321 |
| 7TRE | A | B | PeptideL | 934 | 7-1132 |
| 7TRE | D | C | Dna | 7 | 12-18 |
| 7TRE | C | D | PeptideL | 152 | 90-241 |
| 7TRD | B | A | Rna | 217 | 33-321 |
| 7TRD | A | B | PeptideL | 934 | 7-1132 |
| 7TRD | N | C | Dna | 7 | 12-18 |
| 7TRF | F | A | PeptideL | 90 | 36-125 |
| 7TRF | E | B | PeptideL | 82 | 18-99 |
| 7TRF | B | C | Rna | 217 | 33-321 |
| 7TRF | A | D | PeptideL | 934 | 7-1132 |
| 7TRF | D | E | Dna | 7 | 12-18 |
| 2I46 | A | A | PeptideL | 203 | 90-301 |
| 2I46 | B | B | PeptideL | 171 | 90-287 |

* 7TRE TERT: author numbering covers residues 7-1132, of which 934 are modelled.
* Sequence 753-761 reads `GHVRKAFKS`; canonical human TERT (UniProt O14746) is `GHVRKAFKS` at those positions.
* Catalytic triad check at 712/868/869: {712: 'ASP', 868: 'ASP', 869: 'ASP'}; the canonical hTERT catalytic aspartates are D712, D868, D869.
* Author-minus-label sequence numbering takes the distinct values [-35] on TERT, so the deposited label sequence (which starts with a purification tag) is offset from the author numbering.

**Verdict:** author numbering is canonical hTERT numbering. All residue numbers in this assessment are author numbers; the label numbering is not used.

## 2. Does the recomputed interface reproduce the published assignments?

Heavy-atom contact definition: any TERT atom within 4.5 Å of any TPP1 atom. Computed footprint (21 residues): **[43, 44, 45, 46, 78, 121, 124, 129, 130, 132, 133, 135, 136, 766, 769, 771, 772, 774, 775, 798, 799]**

| Sekne 2022 assignment | residues in that assignment | found in the recomputed footprint |
|---|---|---|
| TEN helix a3 (44-46) | 44-46 | all -> [44, 45, 46] |
| TEN K78 | 78-78 | all -> [78] |
| TEN helix a5 (122-135) | 122-135 | 6/14 -> [124, 129, 130, 132, 133, 135] |
| IFD b17 / b-sheet (766, 767, 797-799) | 766-799 | 3/5 -> [766, 798, 799] |
| IFD helix a25 (771-775) | 771-775 | 4/5 -> [771, 772, 774, 775] |

All five published regions are hit, and the recomputed footprint contains no region that the published assignment does not mention. The coverage is partial at the 4.5 Å cutoff (helix a5 6/14, b17/sheet 3/5, a25 4/5), which is the expected behaviour of a strict distance-only definition: the missing residues are ones whose neighbours are recovered, and a slightly larger cutoff recovers them. The point of this check is not the exact count but that the recomputed interface identifies the same three TERT elements as the published structure - the TEN helices, the IFD beta-sheet, and IFD helix a25 - which is what validates the numbering chain used for every measurement below.

## 3. What does TPP1 actually touch in the proposed 760-766 conduit window?

| TERT residue | name | min d(TPP1) / Å | TPP1 partners within 4.5 Å |
|---|---|---|---|
| 750 | ALA | 28.86 | — |
| 751 | ALA | 31.55 | — |
| 752 | HIS | 33.68 | — |
| 753 | GLY | 34.90 | — |
| 754 | HIS | 31.23 | — |
| 755 | VAL | 28.34 | — |
| 756 | ARG | 26.03 | — |
| 757 | LYS | 25.43 | — |
| 758 | ALA | 21.95 | — |
| 759 | PHE | 19.98 | — |
| 760 | LYS | 17.59 | — |
| 761 | SER | 15.18 | — |
| 762 | HIS | 12.95 | — |
| 763 | VAL | 13.57 | — |
| 764 | SER | 9.86 | — |
| 765 | THR | 7.21 | — |
| 766 | LEU | 3.64 | [93, 112, 113] **<-** |
| 767 | THR | 4.86 | — |
| 768 | ASP | 7.11 | — |
| 769 | LEU | 3.67 | [183, 212] **<-** |
| 770 | GLN | 6.49 | — |
| 771 | PRO | 3.52 | [212, 213] **<-** |
| 772 | TYR | 3.36 | [180, 211, 213] **<-** |
| 773 | MET | 6.70 | — |
| 774 | ARG | 2.32 | [168] **<-** |
| 775 | GLN | 4.41 | [166] **<-** |
| 776 | PHE | 9.77 | — |
| 777 | VAL | 10.83 | — |
| 778 | ALA | 10.51 | — |
| 779 | HIS | 10.07 | — |
| 780 | LEU | 13.04 | — |

Of the 7 residues in 760-766, **1 contact TPP1 at all** ([766]). The gradient along this stretch is smooth and monotonic (17.6 Å -> 3.6 Å), but a distance gradient is not a contact: six of the seven residues that the earlier hypothesis proposed as the conduit's proximal end never come within 4.5 Å of TPP1, and five never come within 7 Å. A side-chain-level conduit through 760-765 is therefore not supported by the coordinates; only 766 is at the interface, and 766 is where the strand-to-turn transition sits (see §5).

## 4. Which domains hold the modelled primer?

| domain | primer-contacting residues (<= 6.0 Å) | count |
|---|---|---|
| TEN (1-196) | [] | 0 |
| TRBD (197-325) | [] | 0 |
| RT (326-739) | [329, 500, 501, 502, 570, 571, 631, 681, 682, 683, 712, 717] | 12 |
| IFD-TRAP (740-800) | [756, 757] | 2 |
| CTE (801-1132) | [839, 866, 867, 868, 869, 930, 931, 932, 944, 945, 948, 949] ... (23 total) | 23 |

The TEN domain contributes **0 primer-contacting residues**; Tyr176, the residue at the centre of the original bipartite-clamp hypothesis, is **11.1 Å** from the nearest primer atom. The primer is held by the RT domain, the IFD-TRAP, and the CTE. This is the measurement that already invalidated the original two-hand TEN/TPP1 clamp model, and it is unchanged by the present re-analysis.

## 5. What secondary structure connects the anchor residues to the TPP1 footprint?

Backbone dihedrals are computed from the coordinates and classified with plain Ramachandran boxes (alpha-helix: -140 < phi < -30 and -80 < psi < -5; beta-strand: -180 < phi < -40 and 90 < psi < 180; anything else is a loop). A single residue's dihedrals are noisy, so the reported assignment is the majority of each residue and its two neighbours. The naive CA(i)->CA(i+4) distance test is tabulated beside it because it is the shortcut that produced the wrong answer in this repository's earlier write-up.

| residue | name | phi / deg | psi / deg | raw | smoothed | CA(i->i+4) / Å |
|---|---|---|---|---|---|---|
| 750 | ALA |  -89.3 |  -58.3 | H | L | 5.7 |
| 751 | ALA | -162.9 | -168.9 | L | H | 6.3 |
| 752 | HIS |  -66.8 |  -47.0 | H | L | 10.9 |
| 753 | GLY |  -88.8 |    3.9 | L | L | 13.4 |
| 754 | HIS | -151.6 |  152.9 | E | E | 12.9 |
| 755 | VAL | -129.7 |  141.5 | E | E | 13.2 |
| 756 | ARG |  -97.8 |  133.4 | E | E | 12.1 |
| 757 | LYS | -124.2 |  109.8 | E | E | 12.9 |
| 758 | ALA | -104.5 |  126.4 | E | E | 12.4 |
| 759 | PHE |  -94.9 |   85.9 | L | L | 13.5 |
| 760 | LYS |  -85.6 | -163.1 | L | L | 13.7 |
| 761 | SER | -162.7 |  125.6 | E | E | 13.8 |
| 762 | HIS | -105.4 |  142.9 | E | E | 13.3 |
| 763 | VAL | -132.8 |  142.5 | E | E | 12.0 |
| 764 | SER | -151.0 |  175.9 | E | E | 6.6 |
| 765 | THR | -122.1 |  169.0 | E | E | 5.7 |
| 766 | LEU |  -60.0 |  -26.7 | H | H | 9.1 |
| 767 | THR |  -67.6 |  -47.8 | H | H | 11.3 |
| 768 | ASP | -123.5 |   48.0 | L | L | 12.0 |
| 769 | LEU |  -74.2 |  168.0 | E | L | 10.7 |
| 770 | GLN | -127.1 |   79.8 | L | L | 8.9 |
| 771 | PRO |  -61.6 |  -54.3 | H | L | 8.1 |
| 772 | TYR |  -68.0 |  157.2 | E | H | 6.3 |
| 773 | MET |  -91.7 |  -30.4 | H | H | 6.1 |
| 774 | ARG |  -71.7 |  -40.8 | H | H | 6.1 |
| 775 | GLN |  -64.3 |  -35.7 | H | H | 6.2 |
| 776 | PHE |  -61.2 |  -46.7 | H | H | 6.1 |
| 777 | VAL |  -72.1 |  -29.3 | H | H | 6.1 |
| 778 | ALA |  -63.3 |  -45.5 | H | H | 6.1 |
| 779 | HIS |  -57.8 |  -44.1 | H | H | 6.0 |
| 780 | LEU |  -65.2 |  -37.7 | H | H | 5.7 |

* **beta-strand** (745-805 only): 745-749, 754-758, 761-765, 785-800
* **alpha-helix** (745-805 only): 751, 766-767, 772-783, 802-805
* **loop / other** (745-805 only): 750-753, 759-760, 768-771, 784, 799-801

**This corrects a claim made earlier in this repository.** `results/allostery/ALLOSTERIC_COUPLING_ANALYSIS.md` and `results/verification/corrected_breakthrough_claim.md` section 8 describe TPP1 as binding *the same IFD-TRAP helix* that carries the DNA anchor and describe 756-775 as one continuous helix. That description is wrong, and the CA(i)->CA(i+4) distances in the last column show why a distance-only check would have caught it: over 752-765 they sit at 12-14 Å, which is extended strand geometry, not helical (5.0-6.5 Å). What the coordinates actually contain:

* the anchor residues **K757 / F759 lie on an extended beta-strand** (strand-like dihedrals across 754-758 and 761-765) - beta17 of the IFD sheet;
* **766-767 is a single helical turn**; 766 is the one residue in the 760-766 window that touches TPP1 (3.6 Å);
* **768-771 is a loop, and 771 is a proline** at the N-cap of the following helix;
* **772-780 is a continuous alpha-helix** (CA(i->i+4) 6.0-6.3 Å) - helix alpha25, which carries TPP1's closest contact, R774 at 2.3 Å.

So the DNA anchor and the TPP1 footprint are joined by a **strand -> turn -> loop -> helix junction**, not by one helix. The distance gradient (26 Å at the strand to 2.3 Å at alpha25) survives this correction and remains the strongest structural fact in the repository. The mechanism wording built on top of it does not: a conduit reading of the gradient cannot be *helix propagation*, because there is no single helix to propagate along. Any surviving version of that idea has to be a junction/hinge mechanism that crosses a secondary-structure boundary and a proline - a weaker and far less specific claim, and the claim that work package 3 has to test.

## 6. What differs between 7TRE (TPP1-bound) and 7TRD (no TPP1)?

* Modelled TERT C-alpha atoms in both models: **934**, same author chain (A) in both, so the comparison is like-for-like.
* Global C-alpha RMSD after superposition: **0.3557 Å**
* Chain composition: {'7TRE': {'B': 217, 'A': 934, 'D': 7, 'C': 152}, '7TRD': {'B': 217, 'A': 934, 'N': 7}}
* The short nucleic-acid chain used as the primer is **D** in 7TRE but **N** in 7TRD; anything that hard-codes chain D will silently analyse the wrong chain in 7TRD.

The two deposits differ by the presence of TPP1 (chain C in 7TRE) and nothing else that the coordinate error can resolve: TERT is the same conformation in both. 7TRD is therefore **not** an independent observation of a TPP1-free state, and treating the pair as a bound/unbound couple for a mechanical model would build a trajectory on a difference of a third of an ångström. This is the reason the mechanical analysis cannot be given a bound/unbound axis from these structures alone.

## 7. Can TPP1 reach the nucleic acid directly?

* TPP1 to DNA primer: **28.73 Å** (minimum heavy-atom distance)
* TPP1 to hTR: **14.24 Å** (minimum heavy-atom distance)

Both are far outside contact range, which reproduces the published statement of Sekne et al. (fig. S10) that TPP1 does not contact the DNA. TPP1 therefore cannot stabilise the primer or hold the hybrid directly; any effect of TPP1 on processivity has to be carried through TERT or through a change in TERT dynamics.

---

Generated by `scripts/structural_validation.py`; machine-readable summary in `structural_validation_data.csv`.
