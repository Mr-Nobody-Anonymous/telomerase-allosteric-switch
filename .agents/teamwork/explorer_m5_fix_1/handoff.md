# Handoff Report: Dynamic Remediation Plan for Figure 2 (`src/telotrade/figures.py`)

**Agent**: `explorer_m5_fix_1`  
**Milestone**: Milestone 5 Remediation  
**Target Module**: `src/telotrade/figures.py:generate_figure2_cross_state_rmsd()`  
**Recipient**: Parent Orchestrator (`eb46fecb-7562-4ce7-aa98-998216cee2c4`) / Remediation Worker  
**Date**: 2026-10-16  

---

## 1. Observation

### 1.1 Verbatim Code Inspection of `src/telotrade/figures.py`
In `src/telotrade/figures.py`, lines 218–322:

1. **Unassigned Table Load on Line 239**:
   ```python
   235:     domain_path = tables_dir / "domain_rmsd.csv"
   236:     cross_path = tables_dir / "cross_state_rmsd.csv"
   237: 
   238:     df_domain = pd.read_csv(domain_path) if domain_path.exists() else pd.DataFrame()
   239:     pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
   ```
   Line 239 evaluates `pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()` as a standalone expression statement. The returned DataFrame is discarded and never assigned to any variable.

2. **Hardcoded Fallback Metric Arrays in Panel A on Lines 262–263**:
   ```python
   254:     if not df_domain.empty:
   255:         df_7trd = df_domain[df_domain["reference_structure"] == "7TRD"].set_index("domain_name")
   256:         df_7trf = df_domain[df_domain["reference_structure"] == "7TRF"].set_index("domain_name")
   257: 
   258:         domain_keys = ["TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"]
   259:         rmsds_7trd = [df_7trd.loc[k, "independent_rmsd_angstrom"] if k in df_7trd.index else 0.0 for k in domain_keys]
   260:         rmsds_7trf = [df_7trf.loc[k, "independent_rmsd_angstrom"] if k in df_7trf.index else 0.0 for k in domain_keys]
   261:     else:
   262:         rmsds_7trd = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]
   263:         rmsds_7trf = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]
   ```
   Lines 262–263 embed hardcoded numeric floats representing pre-computed RMSD constants.

3. **Hardcoded Comparison RMSD Array in Panel B on Lines 284–294**:
   ```python
   284:     comparisons = [
   285:         "7TRE vs 7TRD\n(Global TERT, 934 Cα)",
   286:         "7TRE vs 7TRD\n(Catalytic Core, 770 Cα)",
   287:         "7TRE vs 7TRD\n(TEN on Core Fit, 160 Cα)",
   288:         "7TRE vs 7TRF\n(Catalytic Core, 770 Cα)",
   289:         "7TRE vs 2I46 Chain A\n(TPP1 OB-fold, 152 Cα)",
   290:         "7TRE vs 2I46 Chain A\n(TEL-patch Core, 6 Cα)",
   291:     ]
   292:     comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]
   293:     comp_colors = ["#2A9D8F", "#264653", "#E76F51", "#264653", "#E63946", "#F4A261"]
   ```
   `comp_rmsds` is hardcoded as a fixed list literal `[0.356, 0.338, 0.434, 0.338, 1.966, 1.651]`.

### 1.2 Inspection of Underlying Data Source `results/tables/cross_state_rmsd.csv`
`results/tables/cross_state_rmsd.csv` contains genuine, dynamically calculated RMSD and atom count metrics:
```csv
target_structure,reference_structure,alignment_selection,num_aligned_ca_atoms,superposition_method,rmsd_angstrom,mean_ca_displacement_angstrom,max_ca_displacement_angstrom,max_displacement_residue
7TRE,7TRD,Global,934,Kabsch SVD,0.356,0.315,1.226,ALA612
7TRE,7TRD,Catalytic Core,770,Kabsch SVD,0.338,0.299,1.232,ALA612
7TRE,7TRD,TEN Domain,160,Kabsch SVD,0.411,0.376,1.0,CYS7
7TRE,7TRD,TEN on Core Fit,160,Kabsch SVD,0.434,0.394,0.947,CYS7
7TRE,7TRF,Global,934,Kabsch SVD,0.356,0.315,1.226,ALA612
7TRE,7TRF,Catalytic Core,770,Kabsch SVD,0.338,0.299,1.232,ALA612
7TRE,7TRF,TEN Domain,160,Kabsch SVD,0.411,0.376,1.0,CYS7
7TRE,7TRF,TEN on Core Fit,160,Kabsch SVD,0.434,0.394,0.947,CYS7
7TRE,2I46,TPP1 OB-fold,152,Kabsch SVD,1.966,1.616,7.001,GLU107
7TRE,2I46,TEL-patch Core,6,Kabsch SVD,1.651,1.531,2.207,GLU171
7TRE,2I46,TEL-patch Extended,9,Kabsch SVD,1.746,1.57,2.814,GLU171
```
Every single value in the hardcoded list `[0.356, 0.338, 0.434, 0.338, 1.966, 1.651]` and its aligned atom count `[934, 770, 160, 770, 152, 6]` exists with 100% precision in `cross_state_rmsd.csv`.

---

## 2. Logic Chain

1. **Root Cause Analysis (Observations 1.1.1, 1.1.3)**:
   - Line 239 attempted to read `cross_path` but neglected to assign the result to `df_cross`.
   - Because `df_cross` was unavailable in local scope, the original author inserted a static list literal `comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]` on line 292.
2. **Identification of Redundant Hardcoded Fallbacks (Observation 1.1.2)**:
   - Lines 262–263 contain hardcoded fallback arrays `[0.411, 0.270, 0.295, 0.262, 0.238, 0.338]` for `rmsds_7trd` and `rmsds_7trf`.
   - If `df_domain.empty` occurs, fallbacks should be initialized to zeros (`[0.0] * len(domain_keys)`), not pre-computed scientific constants.
3. **Dynamic Query Formulation (Observations 1.1, 1.2)**:
   - Assign `df_cross = pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()`.
   - Also allow fallback to `domain_rmsds.csv` if `domain_rmsd.csv` is requested.
   - For Panel B, iterate over structured tuples:
     `(reference_structure, alignment_selection, label_template, hex_color)`
   - Filter `df_cross[(df_cross["reference_structure"] == ref) & (df_cross["alignment_selection"] == sel)]`.
   - Dynamically extract `rmsd_angstrom` and `num_aligned_ca_atoms` for each comparison.
   - Dynamically format labels with `{n}` Cα atoms.
   - Dynamically compute x-axis limits (`ax2.set_xlim(0, max(max_rmsd_b + 0.35, 2.3))`) and y-axis limits for Panel A (`ax1.set_ylim(0, max(max_rmsd_a * 1.25, 0.55))`).
4. **Resolution of Integrity Violation**:
   - Zero hardcoded scientific numbers remain in `generate_figure2_cross_state_rmsd()`.
   - The plot is 100% sensitive to underlying changes in `cross_state_rmsd.csv` and `domain_rmsd.csv`.

---

## 3. Caveats

1. **Scope Boundary**: This report focuses strictly on `generate_figure2_cross_state_rmsd()` in `src/telotrade/figures.py`. A separate issue in Figure 4 (`atom_counts = [9, 46, 123]`) is assigned to Explorer 2 (`explorer_m5_fix_2`).
2. **Missing Row Fallback**: In the rare event that `df_cross` exists but is missing one of the 6 comparison combinations, the logic gracefully appends `0.0` with `n=0` to ensure array lengths match across bars, labels, and colors without crashing.
3. **Execution Safety**: Per constraints, this explorer has NOT modified `src/telotrade/figures.py`. The patch is provided in `proposed_figure2_patch.patch` and detailed below for the remediation worker to apply.

---

## 4. Conclusion & Actionable Recommendation

The remediation plan eliminates all hardcoded constants and unassigned expressions in `generate_figure2_cross_state_rmsd()`.

### Exact Replacement Instructions for Worker

In `src/telotrade/figures.py`:

#### 1. Lines 235–240: Replace unassigned read with `df_cross` assignment
**Before**:
```python
    domain_path = tables_dir / "domain_rmsd.csv"
    cross_path = tables_dir / "cross_state_rmsd.csv"

    df_domain = pd.read_csv(domain_path) if domain_path.exists() else pd.DataFrame()
    pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
```

**After**:
```python
    domain_path = tables_dir / "domain_rmsd.csv"
    if not domain_path.exists():
        domain_path = tables_dir / "domain_rmsds.csv"
    cross_path = tables_dir / "cross_state_rmsd.csv"

    df_domain = pd.read_csv(domain_path) if domain_path.exists() else pd.DataFrame()
    df_cross = pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
```

#### 2. Lines 254–264: Eliminate hardcoded fallback arrays
**Before**:
```python
    if not df_domain.empty:
        df_7trd = df_domain[df_domain["reference_structure"] == "7TRD"].set_index("domain_name")
        df_7trf = df_domain[df_domain["reference_structure"] == "7TRF"].set_index("domain_name")

        domain_keys = ["TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"]
        rmsds_7trd = [df_7trd.loc[k, "independent_rmsd_angstrom"] if k in df_7trd.index else 0.0 for k in domain_keys]
        rmsds_7trf = [df_7trf.loc[k, "independent_rmsd_angstrom"] if k in df_7trf.index else 0.0 for k in domain_keys]
    else:
        rmsds_7trd = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]
        rmsds_7trf = [0.411, 0.270, 0.295, 0.262, 0.238, 0.338]
```

**After**:
```python
    domain_keys = ["TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"]
    if not df_domain.empty:
        df_7trd = df_domain[df_domain["reference_structure"] == "7TRD"].set_index("domain_name")
        df_7trf = df_domain[df_domain["reference_structure"] == "7TRF"].set_index("domain_name")

        rmsds_7trd = [float(df_7trd.loc[k, "independent_rmsd_angstrom"]) if k in df_7trd.index else 0.0 for k in domain_keys]
        rmsds_7trf = [float(df_7trf.loc[k, "independent_rmsd_angstrom"]) if k in df_7trf.index else 0.0 for k in domain_keys]
    else:
        rmsds_7trd = [0.0] * len(domain_keys)
        rmsds_7trf = [0.0] * len(domain_keys)
```

#### 3. Line 273: Dynamic Y-Limit for Panel A
**Before**:
```python
    ax1.set_ylim(0, 0.55)
```

**After**:
```python
    max_rmsd_a = max(rmsds_7trd + rmsds_7trf) if (rmsds_7trd and rmsds_7trf) else 0.4
    ax1.set_ylim(0, max(max_rmsd_a * 1.25, 0.55))
```

#### 4. Lines 284–294: Replace hardcoded `comp_rmsds` with dynamic query
**Before**:
```python
    comparisons = [
        "7TRE vs 7TRD\n(Global TERT, 934 Cα)",
        "7TRE vs 7TRD\n(Catalytic Core, 770 Cα)",
        "7TRE vs 7TRD\n(TEN on Core Fit, 160 Cα)",
        "7TRE vs 7TRF\n(Catalytic Core, 770 Cα)",
        "7TRE vs 2I46 Chain A\n(TPP1 OB-fold, 152 Cα)",
        "7TRE vs 2I46 Chain A\n(TEL-patch Core, 6 Cα)",
    ]
    comp_rmsds = [0.356, 0.338, 0.434, 0.338, 1.966, 1.651]
    comp_colors = ["#2A9D8F", "#264653", "#E76F51", "#264653", "#E63946", "#F4A261"]
```

**After**:
```python
    comparison_targets = [
        ("7TRD", "Global", "7TRE vs 7TRD\n(Global TERT, {n} Cα)", "#2A9D8F"),
        ("7TRD", "Catalytic Core", "7TRE vs 7TRD\n(Catalytic Core, {n} Cα)", "#264653"),
        ("7TRD", "TEN on Core Fit", "7TRE vs 7TRD\n(TEN on Core Fit, {n} Cα)", "#E76F51"),
        ("7TRF", "Catalytic Core", "7TRE vs 7TRF\n(Catalytic Core, {n} Cα)", "#264653"),
        ("2I46", "TPP1 OB-fold", "7TRE vs 2I46 Chain A\n(TPP1 OB-fold, {n} Cα)", "#E63946"),
        ("2I46", "TEL-patch Core", "7TRE vs 2I46 Chain A\n(TEL-patch Core, {n} Cα)", "#F4A261"),
    ]

    comparisons: list[str] = []
    comp_rmsds: list[float] = []
    comp_colors: list[str] = []

    if not df_cross.empty:
        for ref, sel, label_tpl, col in comparison_targets:
            match = df_cross[
                (df_cross["reference_structure"] == ref)
                & (df_cross["alignment_selection"] == sel)
            ]
            if not match.empty:
                row = match.iloc[0]
                n_atoms = int(row.get("num_aligned_ca_atoms", 0))
                rmsd_val = float(row["rmsd_angstrom"])
                comparisons.append(label_tpl.format(n=n_atoms))
                comp_rmsds.append(rmsd_val)
                comp_colors.append(col)
            else:
                comparisons.append(label_tpl.format(n=0))
                comp_rmsds.append(0.0)
                comp_colors.append(col)
    else:
        comparisons = [tpl.format(n=0) for _, _, tpl, _ in comparison_targets]
        comp_rmsds = [0.0] * len(comparison_targets)
        comp_colors = [col for _, _, _, col in comparison_targets]
```

#### 5. Line 315: Dynamic X-Limit for Panel B
**Before**:
```python
    ax2.set_xlim(0, 2.3)
```

**After**:
```python
    max_rmsd_b = max(comp_rmsds) if comp_rmsds else 2.0
    ax2.set_xlim(0, max(max_rmsd_b + 0.35, 2.3))
```

### Complete Unified Patch File
The complete unified patch is saved at:  
`C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\explorer_m5_fix_1\proposed_figure2_patch.patch`

---

## 5. Verification Method

To independently verify this plan before and after applying the patch:

1. **Verify No Hardcoded Arrays Remain**:
   ```bash
   grep -n "comp_rmsds =" src/telotrade/figures.py
   # Expected: comp_rmsds: list[float] = [] (no literal float list)

   grep -n "rmsds_7trd =" src/telotrade/figures.py
   # Expected: rmsds_7trd = [float(...)] and rmsds_7trd = [0.0] * len(domain_keys) (no literal floats)

   grep -n "df_cross =" src/telotrade/figures.py
   # Expected: df_cross = pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
   ```

2. **Verify Execution with Real Pipeline Data**:
   ```bash
   python -c "from telotrade.figures import generate_figure2_cross_state_rmsd; generate_figure2_cross_state_rmsd('results/tables', 'results/figures/figure2_cross_state_rmsd.png')"
   ```
   Check that `results/figures/figure2_cross_state_rmsd.png` is generated (> 1024 bytes) without warnings or exceptions.

3. **Verify Runtime Sensitivity**:
   Run an inline test that creates a temporary copy of `cross_state_rmsd.csv` with perturbed values (e.g., setting `rmsd_angstrom` to `0.999` and `num_aligned_ca_atoms` to `9999` for Global 7TRD), and confirm that the extracted comparison text contains `9999` and the bar value is `0.999`.

4. **Verify Linter and Tests**:
   ```bash
   ruff check src/telotrade/figures.py
   pytest tests/e2e/test_tier4_scenarios.py -k test_t4_figures_directory_and_files
   ```
