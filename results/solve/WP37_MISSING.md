# WP37 — Adding the four things that were missing

    python scripts/solve/wp37_missing.py
    -> results/solve/wp37_missing.csv

    python scripts/solve/wp37_missing.py --self-test
    -> 28 checks, 0 failed

WP36 ended by naming what the programme lacks. This adds each one, and is
explicit about which are **closed** and which remain **impossible here**.

| # | missing | after WP37 |
|---|---|---|
| 1 | a measurement of `n` | **CLOSED as a design** — one chase cannot; two chases can, exactly |
| 2 | `δ*` was defined twice | **CLOSED** — the physical quantity needs no exponent |
| 3 | a second enzyme | **PARTLY** — a second independent R½ pair exists; still no measured fold |
| 4 | a TEL-patch affinity | **NOT POSSIBLE HERE** — route + falsifiable docking proposal added |

**Quarantine unchanged.** No therapeutic claim.

---

## 1. A measurement of `n` — the two-chase inversion

WP35's P2 measured **one** chase and tried to infer τ *and* the fold. WP36 showed
that is impossible: `n` and `τ` are confounded, and one fold cannot fix either.

Two chases on the same samples cancel the confound **exactly**:

```
untrapped  m = κ(n,τ) · (1 + τ)      stalling keeps complexes on the primer
trapped    c = κ(n,τ) · (1 − τ)      a trap removes the retries

m/c = (1+τ)/(1−τ)    ⇒    τ = (m/c − 1)/(m/c + 1)
```

`κ` and `n` **cancel**, so τ needs no model assumption at all. Then
`κ = c/(1−τ)` and `n` follows from κ at that τ.

| n_true | τ_true | m | c | τ_rec | κ_rec | n_rec |
|---|---|---|---|---|---|---|
| 2 | 0.1000 | 3.0714 | 2.5130 | 0.1000 | 2.7922 | **2.000** |
| 2 | 0.3500 | 3.7187 | 1.7905 | 0.3500 | 2.7546 | **2.000** |
| 6 | 0.1833 | 3.1954 | 2.2054 | 0.1833 | 2.7004 | **6.000** |
| 24 | 0.1833 | 3.1605 | 2.1813 | 0.1833 | 2.6709 | **24.000** |
| 24 | 0.3500 | 3.5403 | 1.7046 | 0.3500 | 2.6225 | **24.000** |

Exact across a 25-point grid. **This is the design P2 should have been.**

> **Do not shortcut this to `κ = √(m·c)`.** `m·c = κ²(1−τ²)`, so that shortcut is
> biased low by exactly the stall share. The division by `(1−τ)` is the whole
> correction — the same "apply it in the right direction" trap the first draft of
> WP35 fell into. Pinned by a test.

**A free qualitative prediction:** stalling makes the *untrapped* reaction look
**more** processive, so `m > c` always, and at τ = 0 the two coincide exactly
(both 2.733486). **A trap that produced `c > m` would refute the stall model
outright, before any fitting.** That is a stronger pre-registration than WP35's.

---

## 2. `δ*` was two quantities under one symbol

`N_INCORP = 6` has been used for two different things:

| quantity | value | status |
|---|---|---|
| nt per telomeric repeat | 6 | **physical, measured, fixed** |
| escape opportunities per repeat | 6 | **modelling choice, never measured** |

| R½ | 6·R½ (physical) | 1/log₂(1+a), n=3 | n=6 | n=24 |
|---|---|---|---|---|
| 2.963 | **17.78** | 8.89 | 17.78 | 71.11 |
| 8.000 | **48.00** | 24.00 | 48.00 | 192.00 |

The model expression equals the physical length **only at n = 6**. So the gate
must be stated exponent-free:

```
δ* = 6 · R½ = 48 bp/PD at the +PT anchor     →  below the 50 floor → CLOSED
a 90 bp/PD cell needs R½ ≥ 15 repeats        →  have 8
```

**This is not a new assumption — it is WP29's own objective**, `6·c·R½ − δ`.
WP29 already had it right. WP35's "OPEN at 57.3" came from substituting a model
expression into an objective that was already exponent-free.

### The one escape hatch, stated as an assumption

`gain = 6·c·R½`, so a cell at R½ = 8 can hold `δ ≤ 48·c`:

| δ (bp/PD) | c needed |
|---|---|
| 50 | **1.04** |
| 60 | 1.25 |
| 90 | **1.88** |

A 50 bp/PD cell needs 1.04 extensions per division; a HUVEC-like 90 bp/PD cell
needs two independent overhangs extended per cycle. **That is the only route past
the ceiling, and nothing in this repository establishes it.** It is recorded as an
open assumption, not claimed as a mechanism.

---

## 3. A second system exists — and M5 predicts a *different* fold

A second, independent R½ pair is reported in *Nature Communications* **11**:5442
(2020), under the same POT1–TPP1 manipulation:

| system | R½ pair | M5 predicted fold |
|---|---|---|
| L&C 2010 (HEK super-telomerase) | 2.963 → 8.0 | **2.7335** |
| *Nat Commun* 2020 | 2.5 → 4.5 | **1.8187** |

M5 is falsifiable **per system** — it predicts 2.73× here and 1.82× there, not a
universal 2.7×. A 2020-system fold near 2.7× would refute it.

**C6 still stands.** Neither is a *measured* off-rate fold, so the identity's only
measured test remains the L&C datum it was built on. What changed is that there
are now **two** independent targets rather than one.

---

## 4. The TEL-patch affinity — what is added, and what is not

**No affinity is claimed.** TLN-amide, recomputed rather than quoted:

```
O=C1C(=S)N(c2ccc(F)cc2CCC(=O)N)C(=O)N(c2ccccc2)C1CC(C)C
C23H24FN3O3S   MW 441.53   logP 4.00   TPSA 83.71
HBD 1  HBA 4  rotB 7  heavy 31  charge 0   Lipinski PASS
```

This also **confirms** WP34's 441.53 (WP37's first draft asserted 429.51 and was
wrong; the computed value agrees with WP34).

### Retrosynthesis — 7 steps

| # | transformation | material |
|---|---|---|
| 1 | Fmoc-Phe(4-F)-OH + methylamine, HATU/DIPEA | Fmoc-Phe(4-F)-NHMe |
| 2 | piperidine deprotection | 4-(2-aminoethyl)fluorobenzene — the ambipolar arm |
| 3 | + thiocarbonyl diimidazole | Ar–CH₂CH₂–N=C=S |
| 4 | + isobutylamine | thiourea |
| 5 | + phenyl isocyanate | unsymmetrical urea |
| 6 | oxalyl chloride / cyclisation | 2-thioxo-imidazolidinone core |
| 7 | + 3-aminopropionamide | **TLN-amide** (racemic — needs chiral SFC) |

**Provenance:** the linker and one-carboxylate anchor are the lifitegrast lesson
(WP30). WP33 showed that analogy transferred the **linker** correctly and the
**anchor** incorrectly — hence the neutral amide here.

### A falsifiable docking proposal, replacing "no compound binds the patch"

- **predicts:** the neutral amide **loses** its interface preference when the
  ambipolar patch is neutralised
- **mutant:** TERT R92Q/K170Q/R180Q (and the reciprocal D166N/E168Q/E169Q/E215Q)
  — patch net charge −1 → 0
- **control:** a purely hydrophobic analogue (no amide) should be **insensitive**
  to the charge-neutralising mutant
- **refuted if:** the amide keeps its interface preference (Δ < −0.15 per heavy
  atom) in the neutralised patch

That is a prediction with a negative control, which is more than a score.

---

## 5. What is still missing, after all of this

| # | still missing | why it cannot be closed here |
|---|---|---|
| 1 | **`n`, measured** | the design exists; the experiment does not |
| 2 | **a measured second-enzyme fold** | no published k_off fold for any second system |
| 3 | **P2, as a two-chase** | needs wet lab; the single-chase version is superseded |
| 4 | **any affinity at the TEL patch** | no FEP/MM-GBSA; rigid-body surrogate only |
| 5 | **`c > 1` extensions/division** | would open the gate; nothing establishes it |
| 6 | **Gate 3 per-band densitometry** | not deposited; provably unreachable from what exists |

The programme now has: one falsifiable mechanism (M5, P ≈ 35 %), one closed gate
with a named escape hatch, one design that would measure the last free parameter,
and one candidate with a route and a negative control. It has **no drug**.

---

## 6. Tests

28 checks, 0 failed. They pin:

- the two-chase inversion recovers `τ`, `κ` **and** `n` exactly across a 25-point grid;
- `κ = √(m·c)` is **not** the fold, and the test asserts the size of that bias;
- at `τ = 0` the two chases coincide, and `m > c` for every `τ > 0`;
- one chase's fold is shared by ≥3 different `(n, τ)` pairs — the confound;
- `δ*`'s two definitions agree **iff** `n = 6`, and `6·R½ = 48` at the anchor;
- the two systems give materially different folds;
- TLN-amide is `C23H24FN3O3S`, MW 441.53, neutral, Lipinski PASS;
- the route has 7 steps and ends at the amide;
- the wild-type patch is net −1 and the proposed mutant net 0.
