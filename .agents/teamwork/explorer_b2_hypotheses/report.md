# Investigation Report: Research Question & Falsifiable Hypotheses Formulation

**Author**: `explorer_b2_hypotheses` (Research Question & Falsifiable Hypotheses Specialist)  
**Date**: October 2026  
**Milestone**: Benchmark Milestone B2 (Phase 2, Requirement R3)  
**Deliverables Produced**: `results/research/question_and_hypotheses.md`  
**Working Directory**: `C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_b2_hypotheses`  

---

## 1. Executive Summary

This report establishes the foundational research question and competing scientific hypotheses governing the benchmark evaluation of the human telomerase TERT–TPP1 recruitment interface. Human telomerase reverse transcriptase (TERT) is recruited to telomeres via the oligosaccharide/oligonucleotide-binding (OB)-fold domain of the shelterin component TPP1. 

Structural investigations across cryo-EM models (PDB: `7TRE`, 3.5 Å; `7TRD`, 3.3 Å; `7TRF`, 3.7 Å) and X-ray crystallography (PDB: `2I46`, 2.7 Å) reveal apparent tensions between rigid structural conservation, extensive localized conformational remodeling, and sharp computational sensitivity to threshold choices. 

To transform narrative interpretations into falsifiable science, we have formulated the **Primary Research Question** and defined **Three Competing, Falsifiable Scientific Hypotheses**:
- **H1 (Evolutionary Conservation / Rigid Architectural Joint / Experimental Artifact)**: The interface is a pre-organized, rigid joint. Coordinate differences reflect experimental resolution limits, crystallographic lattice contacts in `2I46`, and construct tags rather than biological plasticity.
- **H2 (Genuine State-Dependent Variability / Induced Fit)**: The interface undergoes genuine conformational remodeling upon recruitment, specifically localized flexibility in the TPP1 DE-loop (peaking at Lys170 with 4.293 Å Cα shift) and bipartite TERT clamping (TEN domain excursion 0.434 Å, TRAP/IFD excursion 0.427 Å), mediating functional regulation and dynamic allostery.
- **H3 (Computational & Parameter Sensitivity)**: Apparent structural differences, contact networks, and buried surface area are predominantly driven by computational parameter choices (distance cutoffs 3.5 vs 4.0 vs 4.5 Å, atomic vdW radius sets, and alignment reference selections).

Each hypothesis is equipped with explicit, mathematically operationalized **support criteria** and **falsification thresholds**. The complete deliverable has been written to `results/research/question_and_hypotheses.md`.

---

## 2. Context & Background Synthesis

### 2.1 The Structural Data Envelope
The comparative structural dataset comprises:
- `7TRE`: Human telomerase catalytic core in complex with TPP1 OB-fold domain (3.5 Å cryo-EM). Contains TERT chain A (residues 7–1132 modeled, 934 Cα) and TPP1 chain C (residues 90–241 modeled, 152 Cα).
- `7TRD`: Free human telomerase catalytic core (3.3 Å cryo-EM). TERT chain A (residues 7–1132 modeled, 934 Cα).
- `7TRF`: Telomerase catalytic core RNP bound to histone H2A–H2B (3.7 Å cryo-EM). TERT chain A (residues 7–1132 modeled, 934 Cα).
- `2I46`: Unbound human TPP1 OB-fold homodimer (2.7 Å X-ray). Chain A has 152 residues modeled (90–241) with DE-loop resolved in crystal contacts; Chain B has 134 residues modeled with DE-loop (168–171) completely disordered.

### 2.2 Reconciling Prior-Art and Evidence Ledger Audits
The prior-art audit (`results/literature/novelty_gap_report.md`) and evidence ledger (`results/audit/evidence_ledger.csv`) established crucial physical and methodological boundaries:
1. **Positive Controls**: The canonical salt bridges `Arg774–Glu168` (2.32 Å) and `Lys78–Glu215` (3.58 Å), the hydrophobic anchors `Leu769–Leu183` (3.67 Å) and `Pro771–Leu212` (3.52 Å), and the buried surface area (856.82 Å²) are confirmatory positive controls derived from published structures, not computational discoveries.
2. **Deconstruction of the "Zero Allostery" Deduction (CLM-022)**: The 0.338 Å catalytic core Cα RMSD between 7TRE and 7TRD is a measure of static coordinate identity in frozen vitreous ice. It falls within experimental coordinate error (~0.35–0.50 Å) and cannot biophysically prove an absence of dynamic allostery, conformational entropy modulation, or altered processivity (RAP).
3. **Deconstruction of the Sequential Induced-Fit Model (CLM-021)**: Inferring a chronological sequence ("Leu183 docks first, then DE-loop collapses") from static crystal and cryo-EM endpoints is scientifically invalid without time-resolved kinetic measurements.
4. **Quarantine of TELIN-1 (CLM-023)**: TELIN-1 is a purely theoretical in silico docking proposal with zero literature citations, zero synthesis, and zero binding affinity data.
5. **Multi-Cutoff Volatility (CLM-002, 003, 004)**: Expanding distance thresholds from 3.5 Å to 4.5 Å increases contact counts from 9 to 123 (>13-fold), demonstrating that qualitative descriptions of interface size are heavily dependent on arbitrary computational parameters.

---

## 3. Detailed Hypothesis Formulations & Falsification Protocols

### 3.1 Primary Research Question
The primary inquiry of this project is:
> **Does human telomerase recruitment via the TERT–TPP1 interface operate as an evolutionarily conserved, rigid structural joint whose observed conformational variations are experimental and crystallographic artifacts (H1); or does it function as a genuine state-dependent dynamic interface driven by induced-fit DE-loop remodeling and bipartite domain clamping that mediates functional regulation (H2); or are the reported structural differences and contact architectures predominantly artifacts of computational cutoff and parameter sensitivities (H3)?**

---

### 3.2 Hypothesis 1: Evolutionary Conservation & Rigid Architectural Joint
- **Premise**: The TERT–TPP1 recruitment interface is an evolutionarily invariant, pre-formed structural dock. Telomerase recruitment occurs via rigid key-in-lock recognition of pre-organized hydrophobic knobs (Leu183, Leu212) and complementary electrostatic surfaces. All structural variations across 7TRE, 7TRD, and 2I46 are artifacts of coordinate uncertainty, crystal packing, and construct tags.
- **Support Criteria**:
  1. Catalytic core Cα RMSD between 7TRE and 7TRD $\le 0.40\ \text{Å}$ (observed: 0.338 Å), falling strictly within Cruickshank DPI coordinate error ($\sigma_{\text{coord}} \approx 0.35\text{--}0.50\ \text{Å}$).
  2. Anchor residues Leu183 and Leu212 show pre-organized rigidity between 2I46 and 7TRE ($\Delta r_{\text{C}\alpha}(\text{Leu183}) = 0.512\ \text{Å} \le 0.60\ \text{Å}$; $\Delta r_{\text{C}\alpha}(\text{Leu212}) = 1.040\ \text{Å} \le 1.10\ \text{Å}$).
  3. Cα displacements across 2I46 Chain A correlate strongly with crystallographic B-factors (Pearson $r \ge 0.55, p < 0.001$) and proximity to crystal lattice neighbors ($< 4.0\ \text{Å}$ for residues 168–171).
  4. MD simulations of apo-TPP1 in solution sample the 7TRE-bound DE-loop conformation without crossing a significant thermodynamic barrier ($\Delta G^{\ddagger} \le 1.0\ \text{kcal/mol}$).
- **Falsification Criteria**:
  1. Cryo-EM map directional coordinate uncertainty analysis demonstrates that TEN (0.434 Å) or TRAP (0.427 Å) domain movements exceed isotropic noise with $p < 0.01$.
  2. Solution NMR or cryo-EM of isolated apo-TPP1 shows the DE-loop stably adopts an open conformation deviating from 7TRE by Cα RMSD $> 2.50\ \text{Å}$.
  3. Transitioning the DE-loop from apo to bound state requires crossing a free energy barrier $\Delta G^{\ddagger} > 3.0\ \text{kcal/mol}$.

---

### 3.3 Hypothesis 2: Genuine State-Dependent Remodeling & Induced-Fit Clamping
- **Premise**: Telomerase recruitment is a dynamic, multi-step induced-fit process. The intrinsically mobile TPP1 DE-loop undergoes substantial conformational remodeling to clamp into the basic cradle, coordinated with a subtle rigid-body hinge clamping of the TEN domain (0.434 Å) and TRAP hairpin (0.427 Å). This transition mediates functional allosteric regulation across the catalytic platform.
- **Support Criteria**:
  1. DE-loop Cα displacement exceeds combined coordinate error ($> 4\sigma$): Lys170 shift = $4.293\ \text{Å} > 6.6\sigma$; mean DE-loop shift = $3.169\ \text{Å} > 4.8\sigma$.
  2. Rigid-body domain displacements of TEN and TRAP/IFD exhibit significant vectorial alignment toward the docked TPP1 center-of-mass ($\cos \theta > 0.70, p < 0.05$).
  3. Normal mode analysis or atomistic MD demonstrates significant remodeling of low-frequency vibrational eigenvalues ($\Delta \omega_1 / \omega_1 > 15\%$) and template-cleft RMSF ($> 0.5\ \text{Å}$) upon TPP1 docking, confirming dynamic allostery.
  4. Real-time binding kinetics (stopped-flow / SPR) fit a two-step induced-fit model ($A + B \rightleftharpoons AB^* \rightleftharpoons AB$) significantly better than a one-step model ($\Delta\text{AIC} > 10.0, k_{\text{isomerization}} > 0$).
- **Falsification Criteria**:
  1. Solution HDX-MS shows zero change in backbone amide exchange ($\Delta\text{HDX} < 2\%, p > 0.20$) across the TERT catalytic core upon TPP1 binding.
  2. Solution NMR of apo-TPP1 reveals the DE-loop is pre-organized in the 7TRE conformation (RMSD $\le 1.0\ \text{Å}$).
  3. Association kinetics are strictly linear with ligand concentration ($R^2 \ge 0.995$), ruling out an isomerization/induced-fit step.
  4. Disulfide-locked DE-loop mutants display binding affinity and processivity identical to wild-type ($\Delta K_d \le 1.2$-fold, $\Delta\text{RAP} \le 10\%$).

---

### 3.4 Hypothesis 3: Computational Parameter Sensitivity & Methodological Artifact
- **Premise**: Reported structural distinctions, contact networks, buried surface area, and domain excursions are computational artifacts arising from threshold choices, radius conventions, and reference alignment selections.
- **Support Criteria**:
  1. Contact count increases $> 10$-fold between 3.5 Å and 4.5 Å ($9 \rightarrow 46 \rightarrow 123$, factor of 13.67x), with $> 60\%$ of contacts categorized as Threshold Sensitive (62.6% appear only at 4.5 Å).
  2. Buried surface area varies by $> 15\%$ or $> 200\ \text{Å}^2$ across standard vdW radius sets (Bondi vs ProtOr) and probe radii (1.2–1.6 Å).
  3. Domain excursion metrics vary by $> 25\%$ between Global TERT (934 Cα) and Catalytic Core (770 Cα) superpositions.
  4. Independent software implementations (Gemmi, Biopython, FreeSASA) diverge by $> 100\ \text{Å}^2$ in BSA when using default settings.
- **Falsification Criteria**:
  1. The 5 Core Robust contact pairs account for $> 80\%$ of electrostatic binding free energy ($\Delta G_{\text{elec}}$), rendering threshold-sensitive contacts functionally negligible.
  2. Under harmonized radius definitions, BSA across independent algorithms converges within $\le 2.0\%$ ($\le 17\ \text{Å}^2$).
  3. Domain mobility ranking ($\text{TEN} > \text{TRAP} > \text{CTE} > \text{RT} > \text{TRBD}$) is perfectly invariant (Spearman $\rho = 1.0, p < 0.001$) across all alignment methods.

---

## 4. Synthesis & Translational Relevance (Axis A vs Axis B)

The formulated hypotheses directly govern the translational feasibility of targeting the TERT–TPP1 interface for cancer therapeutics:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        TRANSLATIONAL IMPACT OF HYPOTHESES                              │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                        │
│   Under H1 (Rigid Conservation):                                                       │
│   • The interface is an invariant evolutionary joint.                                  │
│   • Disruption by small molecules will inevitably phenocopy germline DC mutations      │
│     (Lys170del, Glu168del), causing bone marrow aplasia and pulmonary fibrosis.        │
│   • Therapeutic window depends entirely on kinetic cell-turnover differentials.        │
│                                                                                        │
│   Under H2 (State-Dependent Remodeling):                                               │
│   • Conformational intermediates (e.g. open DE-loop, unlatched TEN) exist in solution. │
│   • Conformation-specific allosteric inhibitors could selectively modulate binding     │
│     kinetics or processivity (RAP) without completely abolishing recruitment.          │
│   • Opens potential for a viable therapeutic window sparing healthy stem cells.        │
│                                                                                        │
│   Under H3 (Parameter Sensitivity):                                                    │
│   • Nominal docked poses (e.g. TELIN-1) are computational artifacts.                   │
│   • Target engagement and affinity must be re-evaluated using rigorous biophysical     │
│     force fields and wet-lab validation before drug discovery commitments.             │
│                                                                                        │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 5. Verification & Deliverable Status

The authoritative deliverable `results/research/question_and_hypotheses.md` has been successfully created in the repository root. All 3 competing hypotheses have been fully articulated with quantitative metrics, testable predictions, and unambiguous falsification thresholds, completing Milestone B2 Requirement R3.
