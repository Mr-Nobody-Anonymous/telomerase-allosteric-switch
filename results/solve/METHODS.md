<!-- generated, do not edit by hand: python scripts/solve/wp49_master_record.py -->

# Methods — how every number was produced

    python scripts/solve/wp49_master_record.py

Each section states the question, the input, the method, the command, the
output, and the published source the method follows.

## 0. Environment

```bash
bash scripts/setup_environment.sh
```

The venv lives at `/home/user/.venv-t` — **outside the repository root**, so it
is not persisted between sessions and must be rebuilt. The branch ref can also
be reset to the base commit; the same script restores it from the remote.

| package | role |
|---|---|
| `gemmi` | mmCIF parsing; all coordinate measurement |
| `numpy`/`scipy` | geometry, distances, statistics |
| `biopython` | sequence descriptors (GRAVY, pI, charge, instability) |
| `rdkit` | SMILES, Lipinski/PAINS/BRENK, QED, SA score, SMARTS |
| `sympy` | symbolic checks in the mechanism analysis |
| `requests`/`pyyaml` | the `src/telotrade` pipeline |

**Unavailable, and the reason several questions stay open:** OpenMM, any docking
engine, Rosetta, MODELLER, a PDB-wide mirror.

## 1. Structure handling — the three rules that were learned the hard way

### Rule 1: identify chains by sequence, never by chain letter

7TRE chain C is TPP1 and chain A is TERT. Both contain residues numbered 168
and 215, so numbering alone is ambiguous between the two proteins. Every script
finds TPP1 by testing that 166–172 reads `DWEEKEF`.

### Rule 2: identify TERT by residue TYPE, not by numbering

A real bug, caught in WP48: matching "the chain containing 78 and 132" selected
**the nucleotide chain**, which is also numbered past 132. The check must be
`LYS78` and `ARG132`.

### Rule 3: state the completeness of the deposit

| structure | chain | residues | heavy atoms/residue | with a side chain |
|---|---|---|---|---|
| 7TRE | TPP1 | 152 | 7.7 | 75% |
| 7TRE | TERT | 934 | 8.0 | 78% |
| 2I46 | TPP1 | 203 | 6.0 | 75% |
| 2I46 | TERT | — | — | absent |

These are **coarse models**. Any clash clearance computed against modelled atoms
alone is an **upper bound**, because missing atoms cannot clash. WP48 added
truncation-aware padding to correct this.

## 2. Interface measurement

**Question:** is the TPP1 TEL patch one surface or two?

**Method:** salt-bridge distances between TPP1 carboxylate oxygens and TERT
cationic nitrogens; then the shortest through-space path between the two
anchors, with a straight-line obstruction test.

**Command:** `python scripts/solve/wp42_designs.py`

| result | value |
|---|---|
| K78–R132 | 4.35 Å |
| straight anchor-to-anchor | 20.26 Å |
| bent path | 25.54 Å, 119° bend, tortuosity 1.26 |
| E215 to K78 / R132 | 4.1 / 4.07 Å |

A **salt bridge is judged on the nearer oxygen.** Quoting only OE1 overstates
the E215–R132 distance by 1.6 Å (5.69 vs 4.07).

## 3. Span pricing — why a small molecule cannot do it

**Method:** convert the geometric span into a ring count at 3.2 Å per phenyl
(para-substituted), then to MW and logP with standard increments.

**Result:** ~7.8 phenyls → MW ≈ 683, logP ≈ 11 against Lipinski 500/5.
**This is a class exclusion.**

## 4. Helix-register analysis — and why it is now deprecated for the triad

**Method:** ideal α-helix geometry, Cα at `(R cos θ, R sin θ, rise·i)` with
R = 2.3 Å, rise = 1.5 Å, twist 100°/residue.

**Valid where the segment is actually helical.** For K78/R132 it is: i+3 =
5.05 Å against a measured 4.35 Å.

**Invalid for triad residues 167–171.** ψ(E169) = +24.4°, which breaks a helix,
and CA168–CA171 = 7.68 Å in 7TRE and 7.70 Å in 2I46. Do not fit helix registers
to the triad (retraction X2).

## 5. Backbone grafting (the method that solved S3)

**Source:** Azoitei et al., *Science* 2011, **334**:373 — computation-guided
backbone grafting of a discontinuous motif. `Multigraft` searches 4-insert
positions for backbone continuity without clash; `Multigraft Design` then
models the connecting segments.

**Why it applies:** the TEL patch epitope is **discontinuous** — the triad and
the K78/R132 grip are 20 Å apart. The paper states epitope grafting had been
limited to *continuous* epitopes "even though most antibody epitopes are
discontinuous". That is exactly this problem.

**Command:** `python scripts/solve/wp45_graft.py`

**Method:** keep both native segments verbatim; measure the junction span;
require linker length ≥ span/3.5 Å; build the construct and recompute every
descriptor from first principles.

## 6. Crosslink screening (four filters)

**Command:** `python scripts/solve/wp47_close_open.py`

Reach ranges, from the stapling literature:

| chemistry | Cβ–Cβ range | brace |
|---|---|---|
| disulfide | 3.5–5.0 Å | native |
| lactam | 4.5–6.0 Å | 7–8 atom |
| hydrocarbon (RCM) | 6.0–7.5 Å | 8 atom, i,i+4 |
| 1,3-diyne (Glaser) | 6.5–7.5 Å | 7 Å diynyl |
| thioether | 9.0–11.5 Å | 10–12 atom |
| biaryl | 10.0–12.5 Å | i,i+7 |

The **direction filter** is the non-obvious one: a staple grows outward along
the native side-chain vector, so a residue whose side chain points away cannot
form one without being rebuilt. Screened as cos > 0 against the gap vector.

## 7. Filter calibration and its ceiling

**Command:** `python scripts/solve/wp48_new_designs.py`

Two things are measured separately and must not be conflated:

1. **The electrophile layer** — SMARTS patterns for Michael acceptors, alkyl
   halides, epoxides, aldehydes, isocyanates, quinones. Scored alone.
2. **The gate stack** — MW, logP, HBD, HBA, SA. Scored end-to-end.

Evidence comes from a **third panel** used by neither WP43 nor WP47. Validation
is defined on that panel only; authoring-panel scores are labelled
demonstrations.

**A rule that was added and removed:** a `MW ≥ 60` floor, which only ever
compensated for missing SMARTS. Once the SMARTS caught acrylonitrile and
formaldehyde it changed no decision, so it was deleted.

## 8. Ledgers, registers, and their enforcement

| artefact | generator | enforcement |
|---|---|---|
| `CLOSURE.md` | `wp46_closure.py` | every status re-checked against the owning JSON |
| `RETRACTIONS.md` | `wp38_archive.py` | a removed marker or unregistered claim fails the suite |
| `STATUS_LEDGER.md` | `wp41_status_ledger.py` | regenerated and compared |
| `MASTER_RECORD.md` | `wp49_master_record.py` | assembled from the above |

## 9. Every self-test in the tree

`wp49_master_record` is excluded here: it runs inside this build and cannot
test a document it is still writing. Its check lives in `tests/`.

| script | result |
|---|---|
| `wp28_operating_point` | **pass** — 8 checks, 0 failed |
| `wp29_maintenance` | **pass** — 11 checks, 0 failed |
| `wp30_ligand_completion` | **pass** — 25 checks, 0 failed |
| `wp31_epitope` | **pass** — 22 checks, 0 failed |
| `wp32_physics_scoring` | **pass** — 26 checks, 0 failed |
| `wp33_valid_scoring` | **pass** — 36 checks, 0 failed |
| `wp34_status` | **pass** — 18 checks, 0 failed |
| `wp35_p2_identifiability` | **pass** — WP35 self-test: 54 passed, 0 failed |
| `wp36_circularity` | **pass** — WP36 self-test: 34 passed, 0 failed |
| `wp37_missing` | **pass** — WP37 self-test: 28 passed, 0 failed |
| `wp38_archive` | **pass** — WP38 self-test: 61 passed, 0 failed |
| `wp39_executable` | **pass** — WP39 self-test: 39 passed, 0 failed |
| `wp40_design_to_fifty` | **pass** — WP40 self-test: 45 passed, 0 failed |
| `wp41_status_ledger` | **pass** — WP41 self-test: 42 passed, 0 failed |
| `wp42_designs` | **pass** — WP42 self-test: 48 checks, 0 failed |
| `wp43_test_battery` | **pass** — WP43 self-test: 48 checks, 0 failed |
| `wp44_modular_designs` | **pass** — WP44 self-test: 39 checks, 0 failed |
| `wp45_graft` | **pass** — WP45 self-test: 30 checks, 0 failed |
| `wp46_closure` | **pass** — WP46 self-test: 26 checks, 0 failed |
| `wp47_close_open` | **pass** — WP47 self-test: 33 checks, 0 failed |
| `wp48_new_designs` | **pass** — WP48 self-test: 25 checks, 0 failed |

## 10. Commit history of this branch

| commit | date | subject |
|---|---|---|
| `c63548d` | 2026-11-24 | WP48: a new design per failure, four isolated agents, and a critic |
| `9916cb2` | 2026-11-23 | WP47: P6 gets a real crosslink, P10 is proven unfixable as framed |
| `7e4af73` | 2026-11-23 | WP46: closure ledger - no, and here is exactly what is left |
| `6cd6832` | 2026-11-23 | WP45: solve that site, and correct the premise that made it look unsolvable |
| `79c3a88` | 2026-11-23 | WP44: one module per site, merged, and the merge stress-tested |
| `7651ca7` | 2026-11-23 | WP43: test all ten designs, and separate the three success rates |
| `55e4078` | 2026-11-23 | WP42: design new molecules and new proteins for the TEL patch |
| `9fdab63` | 2026-11-23 | Merge remote-tracking branch 'origin/arena/7b6ab982-telomerase-allosteric-switch' into arena/7b6ab982-telomerase-allosteric-switch |
| `1f976b6` | 2026-11-23 | Resolve the contradicting test and the uncertain claims |
| `ae92678` | 2026-11-21 | Ledger: add CLM-025 to CLM-029 from the no-lab audit |
| `61c2432` | 2026-11-21 | Merge remote-tracking branch 'origin/arena/7b6ab982-telomerase-allosteric-switch' into arena/7b6ab982-telomerase-allosteric-switch Co-authored-by: arena-agent <297053741+arena-agent@users.noreply.github.com> |
| `0483471` | 2026-11-21 | Ledger: add CLM-025 to CLM-029 from the no-lab audit |
| `1440d2e` | 2026-11-21 | No-lab audit: K757 geometry, panel re-analysis, power, mixture recovery |
| `29c3d79` | 2026-11-19 | Merge pull request #12 from Mr-Nobody-Anonymous/arena/86583e69-telomerase-allosteric-switch |
