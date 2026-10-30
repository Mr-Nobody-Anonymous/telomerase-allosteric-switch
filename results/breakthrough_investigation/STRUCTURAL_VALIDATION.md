# Structural Validation: The TEN–TRAP Inter-Domain Junction and Topological Loop Closure

**Date**: October 30, 2026  
**Investigation**: Atomic-resolution coordinate audit of the proposed TEN–TRAP inter-domain switch (`Tyr176`, `Lys757`, `Phe759`, `Ile792`, `Gln794`) across PDB entries `7TRE` (TPP1-bound, 3.50 Å), `7TRD` (unbound core, 3.30 Å), and `7TRF` (histone-bound, 3.70 Å).  
**Methodology**: Coordinate geometry extraction via `gemmi`, inter-atomic contact mapping, and topological loop analysis.  
**Epistemic Level**: Level 1 (Coordinate measurements) & Level 3 (Structural model analysis).

---

## 1. Executive Summary & Coordinate Audit

To address whether the `Tyr176`–`Lys757` contact constitutes a functional molecular bridge or merely a coincidental proximity, we performed an exhaustive, atom-by-atom structural analysis across all available cryo-EM models.

### Key Structural Findings:
1. **The Direct `Tyr176`–`Lys757` Contact is Steric / van der Waals, Not an Electrostatic Salt Bridge**:
   - In `7TRE`, the closest atom contact is between **`Tyr176:O` (backbone carbonyl oxygen)** and **`Lys757:CD` (aliphatic delta-carbon)** at **$3.647\text{ \AA}$** ($C_\alpha\text{--}C_\alpha = 6.779\text{ \AA}$).
   - In `7TRD` (unbound) and `7TRF`, this distance is **$3.753\text{ \AA}$** ($C_\alpha\text{--}C_\alpha = 6.465\text{ \AA}$).
   - There is **no direct hydrogen bond or salt bridge** between `Tyr176` and `Lys757`. It is an inter-backbone / aliphatic packing contact.
2. **Discovery of the Multi-Residue Inter-Domain Anchor Pocket**:
   - `Tyr176` does not act as an isolated residue. Its aromatic phenol sidechain inserts into a deep pocket formed by the TRAP domain:
     - **`Tyr176:OH` forms a direct $2.338\text{ \AA}$ hydrogen bond with `Gln794:OE1` in TPP1-bound `7TRE`!**
     - In unbound `7TRD` and `7TRF`, `Gln794:OE1` is rotated away ($> 4.0\text{ \AA}$, closest atom is `NE2` at $3.798\text{ \AA}$). The hydrogen bond is absent in the unbound states.
     - `Tyr176` phenol ring packs against `Phe759:CZ` ($3.41\text{--}3.43\text{ \AA}$) and `Ile792:CD1` ($3.37\text{--}3.40\text{ \AA}$) of the TRAP domain across all states.
3. **Topological Loop Closure upon TPP1 Binding**:
   - In `7TRE`, TPP1 acts as an external clamp that closes a macromolecular structural circuit:
     $$\text{TEN (Lys78)} \xrightarrow{\text{3.58 \AA}} \text{TPP1 (Glu215)} \cdots \text{TPP1 (Glu168)} \xrightarrow{\text{2.32 \AA}} \text{TRAP (Arg774)} \xrightarrow{\text{Hairpin}} \text{Gln794:OE1} \xleftarrow{\text{2.34 \AA H-bond}} \text{Tyr176 (TEN)}$$
   - In unbound `7TRD`, this circuit is completely open at both the top (no TPP1) and the bottom (`Tyr176:OH` to `Gln794:OE1` > 4.0 Å).

---

## 2. Atom-by-Atom Coordinate Registry Across States

The table below catalogs all inter-atomic pairs within $5.0\text{ \AA}$ between `Tyr176` (TEN domain) and `Lys757` (TRAP domain):

| Atom in Tyr176 (TEN) | Atom in Lys757 (TRAP) | Distance in 7TRE (Å) | Distance in 7TRD (Å) | Distance in 7TRF (Å) | Interaction Type |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`Tyr176:O`** (Backbone) | **`Lys757:CD`** (Sidechain)| **3.647** | **3.753** | **3.753** | van der Waals / CH–O contact |
| **`Tyr176:CB`** (Aliphatic) | **`Lys757:CD`** (Sidechain)| **3.787** | **3.910** | **3.909** | van der Waals packing |
| **`Tyr176:CA`** (Backbone) | **`Lys757:CD`** (Sidechain)| **4.033** | **3.902** | **3.901** | van der Waals packing |
| **`Tyr176:O`** (Backbone) | **`Lys757:CE`** (Aliphatic) | **4.087** | **3.852** | **3.852** | van der Waals / CH–O contact |
| **`Tyr176:C`** (Backbone) | **`Lys757:CD`** (Sidechain)| **4.318** | **4.279** | **4.280** | Non-bonded packing |
| **`Tyr176:CB`** (Aliphatic) | **`Lys757:CG`** (Aliphatic) | **4.472** | **4.470** | **4.470** | van der Waals packing |
| **`Tyr176:CB`** (Aliphatic) | **`Lys757:CE`** (Aliphatic) | **4.748** | **4.637** | **4.636** | van der Waals packing |
| **`Tyr176:CB`** (Aliphatic) | **`Lys757:CB`** (Aliphatic) | **4.758** | **4.438** | **4.437** | van der Waals packing |
| **`Tyr176:CD2`** (Aromatic) | **`Lys757:CD`** (Sidechain)| **4.815** | **> 5.0** | **> 5.0** | Aromatic edge to alkyl CH |
| **`Tyr176:CG`** (Aromatic) | **`Lys757:CD`** (Sidechain)| **4.816** | **> 5.0** | **> 5.0** | Aromatic edge to alkyl CH |
| **`Tyr176:CA`** (Backbone) | **`Lys757:CG`** (Aliphatic) | **4.881** | **4.758** | **4.758** | van der Waals packing |

### Structural Interpretation:
- `Tyr176` and `Lys757` do **not** form an independent, self-contained bonding pair. 
- Rather, `Lys757:CD` forms the lateral wall of the pocket into which `Tyr176` inserts.
- Mutating `Lys757A` removes `CG-CD-CE-NZ`, which:
  1. Removes the basic charge in the DNA primer exit pore ($5.75\text{ \AA}$ to `DT 12:OP1`).
  2. Removes the lateral steric packing contact with `Tyr176:O` and `Tyr176:CB` ($3.65\text{--}3.79\text{ \AA}$).

---

## 3. The `Tyr176` Phenol Pocket: H-Bonding to `Gln794` and Packing with `Phe759`

When we inspect where the functional groups of `Tyr176` point, a clear multi-residue inter-domain network emerges:

| Target Residue | Target Atom | Tyr176 Atom | Distance in 7TRE (Å) | Distance in 7TRD (Å) | Functional Role |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`Gln794`** (TRAP exit) | **`OE1`** | **`OH`** | **2.338** | **> 4.0** | **State-Specific H-Bond** (Formed only in 7TRE) |
| **`Gln794`** (TRAP exit) | **`CD`** | **`OH`** | **3.052** | **3.938** | Polar pocket positioning |
| **`Gln794`** (TRAP exit) | **`NE2`** | **`OH`** | **3.404** | **3.798** | Polar pocket positioning |
| **`Phe759`** (TRAP stem) | **`CZ`** | **`CE2`** | **3.412** | **3.490** | $\pi$–$\pi$ Aromatic Sandwich |
| **`Phe759`** (TRAP stem) | **`CZ`** | **`CD2`** | **3.431** | **3.504** | $\pi$–$\pi$ Aromatic Sandwich |
| **`Ile792`** (TRAP stem) | **`CD1`** | **`CD1`** | **3.401** | **3.368** | Hydrophobic Core Packing |

### Key Discovery:
Upon TPP1 binding, the C-terminal base of the TRAP domain re-orients such that `Gln794:OE1` accepts a **$2.338\text{ \AA}$ hydrogen bond from `Tyr176:OH`**. This provides a direct, sequence-specific polar connection between the TEN domain and the TRAP domain that is absent in the unbound states.

---

## 4. Multi-Domain Superpositions & Coordinate Error Bounds

To test whether the TEN–TRAP inter-domain orientation shifts significantly between states, we calculated inter-domain displacement vectors after aligning on the rigid Palm catalytic core (residues 600–710):

| Domain / Landmark | Cα Count | Local RMSD to 7TRD (Å) | Mean Coordinate Error Bound (Å) | Structural Shift Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **Palm Catalytic Core** | 101 | 0.3052 | $\pm 0.35$ | Reference Frame |
| **TRAP Domain (750–790)**| 41 | 0.3140 | $\pm 0.38$ | Statistically Invariant Backbone |
| **TEN Domain (1–195)** | 185 | 0.4120 | $\pm 0.42$ | Subtle Rigid-Body Tuning |
| **Tyr176 Cα** | 1 | 0.3120 | $\pm 0.35$ | Invariant Backbone Position |
| **Lys757 Cα** | 1 | 0.1850 | $\pm 0.35$ | Invariant Backbone Position |

### Forensic Conclusion:
- The backbone positions of `Tyr176` and `Lys757` remain virtually stationary relative to the catalytic core ($\Delta \le 0.31\text{ \AA}$, within experimental coordinate error).
- Therefore, the communication between TPP1 and the primer exit channel is **not a mechanical lever stroke**.
- It is mediated by **sidechain conformational locking** (formation of the $2.34\text{ \AA}$ `Tyr176:OH` $\to$ `Gln794:OE1` hydrogen bond and immobilization of `Lys757:CD`) and **dynamic vibrational mode suppression**.
