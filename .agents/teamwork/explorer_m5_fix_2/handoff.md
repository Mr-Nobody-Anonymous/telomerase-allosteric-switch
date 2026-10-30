# Milestone 5 Remediation Explorer 2 Handoff Report: Figure 4 Dynamic Architecture Plan

**Role**: Milestone 5 Remediation Explorer 2 (`explorer_m5_fix_2`)  
**Target File**: `src/telotrade/figures.py` (`generate_figure4_cutoff_robustness`, lines 458–580)  
**Parent Orchestrator**: `eb46fecb-7562-4ce7-aa98-998216cee2c4`  
**Date**: 2026-10-16  

---

## 1. Observation

Direct examination of `src/telotrade/figures.py` (lines 458–580) and `results/tables/` revealed the following specific hardcoding and facade issues:

1. **Panel A Contact Counts Hardcoding (lines 491–493)**:
   ```python
   # Panel A: Monotonic Escalation
   cutoffs = [3.5, 4.0, 4.5]
   atom_counts = [9, 46, 123]
   res_pair_counts = [5, 21, 32]
   ```
   `atom_counts` and `res_pair_counts` are static integer lists rather than dynamically computed from the contact tables.
   Furthermore, line 515 sets a static y-axis upper limit:
   ```python
   ax1.set_ylim(0, 145)
   ```

2. **Panel B Hardcoded Defaults and Fallback Array (lines 518–533)**:
   ```python
   if not df_rob.empty and "robustness_class" in df_rob.columns:
       rob_counts = df_rob["robustness_class"].value_counts()
       labels = [
           f"Core Robust (3.5 Å)\n[n={rob_counts.get('Core Robust', 9)}]",
           f"Intermediate (4.0 Å)\n[n={rob_counts.get('Intermediate', 37)}]",
           f"Sensitive (4.5 Å)\n[n={rob_counts.get('Threshold Sensitive', 77)}]",
       ]
       values = [
           rob_counts.get("Core Robust", 9),
           rob_counts.get("Intermediate", 37),
           rob_counts.get("Threshold Sensitive", 77),
       ]
   else:
       labels = ["Core Robust (3.5 Å)\n[n=9]", "Intermediate (4.0 Å)\n[n=37]", "Sensitive (4.5 Å)\n[n=77]"]
       values = [9, 37, 77]
   ```
   `rob_counts.get(..., default)` uses hardcoded integers (`9`, `37`, `77`), and the `else:` branch contains literal fallback lists `values = [9, 37, 77]` and hardcoded labels.

3. **Panel C Synthetic Coordinate Series Fallbacks (lines 550–556)**:
   ```python
   if not df_rob.empty and "min_atom_distance_angstrom" in df_rob.columns:
       ten_dists = df_rob[df_rob["tert_domain"] == "TEN"]["min_atom_distance_angstrom"].dropna()
       trap_dists = df_rob[df_rob["tert_domain"] == "TRAP_IFD"]["min_atom_distance_angstrom"].dropna()
   else:
       ten_dists = pd.Series([3.35, 3.45, 3.51, 3.58, 3.62, 4.12, 4.28])
       trap_dists = pd.Series([2.32, 3.24, 3.36, 3.40, 3.48, 3.52, 3.56, 3.67])
   ```
   The `else:` branch defines synthetic mock distance lists for TEN and TRAP_IFD domains instead of enforcing data availability or raising an informative error.

4. **Resource Management and File Flushing (lines 479, 577–579)**:
   ```python
   fig = plt.figure(figsize=(14, 6))
   ...
   plt.tight_layout(rect=[0, 0, 1, 0.95])
   plt.savefig(output_path)
   plt.close()
   return output_path
   ```
   - `plt.close()` is not enclosed within a `try...finally` block, risking unreleased memory and file handles if plotting fails.
   - Calling global `plt.close()` rather than `plt.close(fig)` can close arbitrary figures.
   - No post-save verification ensures `output_path.stat().st_size > 1024`, contributing to intermittent 0-byte file issues observed in E2E tests (`test_t4_figures_directory_and_files`).

5. **Empirical Verification of Available Data Sources**:
   Direct query of `results/tables/` confirmed:
   - `interface_contacts_7tre_3.5A.csv`: 9 rows (atomic contacts), 5 unique residue pairs `(tert_resnum, tpp1_resnum)`
   - `interface_contacts_7tre_4.0A.csv`: 46 rows (atomic contacts), 21 unique residue pairs `(tert_resnum, tpp1_resnum)`
   - `interface_contacts_7tre_4.5A.csv`: 123 rows (atomic contacts), 32 unique residue pairs `(tert_resnum, tpp1_resnum)`
   - `interface_contacts_robustness.csv`: 123 rows; `robustness_class` counts: `Core Robust` = 9, `Intermediate` = 37, `Threshold Sensitive` = 77; `tert_domain` counts: `TEN` = 49 (distances 3.35–4.50 Å), `TRAP_IFD` = 74 (distances 2.32–4.50 Å).
   - Distances in `interface_contacts_7tre_4.5A.csv["distance_angstrom"]` and `interface_contacts_robustness.csv["min_atom_distance_angstrom"]` are 100% identical.

---

## 2. Logic Chain

1. **Root Cause**: The developer originally used hardcoded literals in `generate_figure4_cutoff_robustness()` to expedite figure drafting or guard against potential table absences.
2. **Standard Violation**: Hardcoding numbers that mimic calculation outputs violates zero-tolerance integrity requirements and renders Figure 4 completely insensitive to changes in underlying structural contacts.
3. **Data Parity & Feasibility**:
   - For Panel A: `len(df_cutoff)` gives exact atomic contact counts (`9`, `46`, `123`). Grouping by or dropping duplicates on `['tert_resnum', 'tpp1_resnum']` yields exact residue pair counts (`5`, `21`, `32`).
   - If individual cutoff files are absent, querying `interface_contacts_robustness.csv` via boolean columns `present_at_3.5A`, `present_at_4.0A`, `present_at_4.5A` and counting unique residue pairs yields identical dynamic results.
   - For Panel B: `df_rob["robustness_class"].value_counts()` yields exact counts (`Core Robust` = 9, `Intermediate` = 37, `Threshold Sensitive` = 77). Setting fallback default to `0` via `rob_counts.get("Core Robust", 0)` ensures accurate reporting without synthetic inflation.
   - For Panel C: `interface_contacts_7tre_4.5A.csv` contains all 123 contacts partitioned by `tert_domain` (`TEN` vs `TRAP_IFD`). Querying either `df_c45["distance_angstrom"]` or `df_rob["min_atom_distance_angstrom"]` produces the exact continuous distribution without any synthetic series.
4. **Error Handling & Resource Safety**:
   - If required CSV files are absent, the function must raise `FileNotFoundError` with a clear explanation rather than fabricating fake data.
   - Wrapping figure generation in `try...finally: plt.close(fig)` guarantees resource release.
   - Asserting `output_path.exists() and output_path.stat().st_size > 1024` guarantees non-zero atomic image generation.

---

## 3. Caveats

1. **Scope Boundary**: This report focuses strictly on Figure 4 (`generate_figure4_cutoff_robustness`). Figure 2 remediation is assigned to `explorer_m5_fix_1`, and general resource wrapping across Figures 1 and 3 is covered by `explorer_m5_fix_3`.
2. **Column Naming Conventions**: In `interface_contacts_7tre_*.csv`, residue numbers are stored in `tert_resnum` and `tpp1_resnum`. In `interface_contacts_robustness.csv`, residues are represented as strings in `tert_residue` and `tpp1_residue`. The dynamic extraction logic supports both schemas seamlessly.
3. **Zero Read-Only Violation**: No source code was modified by this explorer agent. All code proposals are provided as recommendations for the Worker agent.

---

## 4. Conclusion & Concrete Implementation Plan

### Recommended Drop-in Replacement for Worker

Replace lines 458–580 of `src/telotrade/figures.py` with the following clean, robust, and zero-hardcoded implementation:

```python
def generate_figure4_cutoff_robustness(
    tables_dir: Path | str,
    output_path: Path | str,
) -> Path:
    """
    Generate Figure 4: Interface contact count and persistence distribution
    across distance cutoffs (3.5 Å, 4.0 Å, 4.5 Å).

    Panels:
    - (A) Contact Escalation across Cutoffs (Atomic Contacts & Residue Pairs).
    - (B) Persistence Classification Breakdown (Core Robust vs Intermediate vs Sensitive).
    - (C) Minimum Atom Distance Distribution across TERT Domains.
    """
    set_publication_style()
    tables_dir = Path(tables_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rob_path = tables_dir / "interface_contacts_robustness.csv"
    df_rob = pd.read_csv(rob_path) if rob_path.exists() else pd.DataFrame()

    fig = plt.figure(figsize=(14, 6))
    try:
        fig.suptitle(
            "Figure 4: Computational Robustness and Distance Sensitivity of TERT–TPP1 Contacts",
            fontweight="bold",
            y=0.98,
        )
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.2])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])

        # Panel A: Contact Escalation across Cutoffs (Dynamic Table Ingestion)
        cutoffs = [3.5, 4.0, 4.5]
        atom_counts = []
        res_pair_counts = []

        for c in cutoffs:
            c_path = tables_dir / f"interface_contacts_7tre_{c:.1f}A.csv"
            if not c_path.exists():
                c_path = tables_dir / f"interface_contacts_{c:.1f}A.csv"

            if c_path.exists():
                df_c = pd.read_csv(c_path)
                atom_counts.append(len(df_c))
                if "tert_resnum" in df_c.columns and "tpp1_resnum" in df_c.columns:
                    n_pairs = len(df_c[["tert_resnum", "tpp1_resnum"]].drop_duplicates())
                elif "tert_residue" in df_c.columns and "tpp1_residue" in df_c.columns:
                    n_pairs = len(df_c[["tert_residue", "tpp1_residue"]].drop_duplicates())
                else:
                    n_pairs = len(df_c)
                res_pair_counts.append(n_pairs)
            elif not df_rob.empty and f"present_at_{c:.1f}A" in df_rob.columns:
                sub = df_rob[df_rob[f"present_at_{c:.1f}A"] == True]
                atom_counts.append(len(sub))
                if "tert_residue" in sub.columns and "tpp1_residue" in sub.columns:
                    n_pairs = len(sub[["tert_residue", "tpp1_residue"]].drop_duplicates())
                else:
                    n_pairs = len(sub)
                res_pair_counts.append(n_pairs)
            else:
                raise FileNotFoundError(
                    f"Contact table for cutoff {c:.1f} Å not found (expected interface_contacts_7tre_{c:.1f}A.csv "
                    f"or interface_contacts_robustness.csv in {tables_dir})"
                )

        x = np.arange(len(cutoffs))
        width = 0.35

        b1 = ax1.bar(x - width / 2, atom_counts, width, label="Atomic Contacts", color="#1D3557", edgecolor="black")
        b2 = ax1.bar(x + width / 2, res_pair_counts, width, label="Residue Pairs", color="#E63946", edgecolor="black")

        ax1.set_xlabel("Distance Cutoff (Å)")
        ax1.set_ylabel("Count")
        ax1.set_title("A. Contact Count Escalation", fontweight="bold")
        ax1.set_xticks(x)
        ax1.set_xticklabels([f"{c:.1f} Å" for c in cutoffs])
        ax1.legend(loc="upper left", framealpha=0.9)

        for bar in b1:
            h = bar.get_height()
            ax1.annotate(f"{int(h)}", xy=(bar.get_x() + bar.get_width() / 2, h), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontweight="bold", fontsize=8.5)
        for bar in b2:
            h = bar.get_height()
            ax1.annotate(f"{int(h)}", xy=(bar.get_x() + bar.get_width() / 2, h), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontweight="bold", fontsize=8.5)

        max_count = max(atom_counts + res_pair_counts) if (atom_counts and res_pair_counts) else 100
        ax1.set_ylim(0, max_count * 1.18)

        # Panel B: Persistence Classification Breakdown (Dynamic Value Counts)
        if df_rob.empty or "robustness_class" not in df_rob.columns:
            raise FileNotFoundError(
                f"interface_contacts_robustness.csv with 'robustness_class' column not found in {tables_dir}"
            )

        rob_counts = df_rob["robustness_class"].value_counts()
        n_core = int(rob_counts.get("Core Robust", 0))
        n_intermediate = int(rob_counts.get("Intermediate", 0))
        n_sensitive = int(rob_counts.get("Threshold Sensitive", 0))

        labels = [
            f"Core Robust (3.5 Å)\n[n={n_core}]",
            f"Intermediate (4.0 Å)\n[n={n_intermediate}]",
            f"Sensitive (4.5 Å)\n[n={n_sensitive}]",
        ]
        values = [n_core, n_intermediate, n_sensitive]
        pie_colors = ["#D90429", "#F4A261", "#2A9D8F"]

        if sum(values) > 0:
            _wedges, _texts, autotexts = ax2.pie(
                values,
                labels=labels,
                autopct="%1.1f%%",
                startangle=140,
                colors=pie_colors,
                wedgeprops={"width": 0.45, "edgecolor": "black", "linewidth": 1},
                textprops={"fontsize": 8},
            )
            for autotext in autotexts:
                autotext.set_fontsize(8.5)
                autotext.set_fontweight("bold")
        else:
            ax2.text(0.5, 0.5, "No robustness contacts", ha="center", va="center")
        ax2.set_title("B. Persistence Partition", fontweight="bold")

        # Panel C: Minimum Atom Distance Distribution across TERT Domains (Dynamic Query)
        c45_path = tables_dir / "interface_contacts_7tre_4.5A.csv"
        if not c45_path.exists():
            c45_path = tables_dir / "interface_contacts_4.5A.csv"

        if c45_path.exists():
            df_c45 = pd.read_csv(c45_path)
            ten_dists = df_c45[df_c45["tert_domain"] == "TEN"]["distance_angstrom"].dropna()
            trap_dists = df_c45[df_c45["tert_domain"] == "TRAP_IFD"]["distance_angstrom"].dropna()
        elif not df_rob.empty and "min_atom_distance_angstrom" in df_rob.columns:
            ten_dists = df_rob[df_rob["tert_domain"] == "TEN"]["min_atom_distance_angstrom"].dropna()
            trap_dists = df_rob[df_rob["tert_domain"] == "TRAP_IFD"]["min_atom_distance_angstrom"].dropna()
        else:
            raise FileNotFoundError(
                f"Distance data not found in {tables_dir}: "
                f"neither interface_contacts_7tre_4.5A.csv nor interface_contacts_robustness.csv is available."
            )

        if len(trap_dists) > 0 or len(ten_dists) > 0:
            all_dists = pd.concat([trap_dists, ten_dists])
            d_min = max(0.0, float(all_dists.min()) - 0.2)
            d_max = float(all_dists.max()) + 0.2
            bins = np.linspace(min(2.0, d_min), max(4.6, d_max), 14)

            ax3.hist(
                [trap_dists, ten_dists],
                bins=bins,
                stacked=True,
                label=["TRAP/IFD Domain Contacts", "TEN Domain Contacts"],
                color=["#1D3557", "#E76F51"],
                edgecolor="black",
                alpha=0.85,
            )
            ax3.axvline(x=3.5, color="#D90429", linestyle="--", lw=1.5, label="3.5 Å (Core Gate)")
            ax3.axvline(x=4.0, color="#F4A261", linestyle="--", lw=1.5, label="4.0 Å (Intermediate Gate)")
            ax3.axvline(x=4.5, color="#2A9D8F", linestyle="--", lw=1.5, label="4.5 Å (Outer Shell)")
            ax3.legend(loc="upper left", framealpha=0.9, fontsize=8)
        else:
            ax3.text(0.5, 0.5, "No contact distances available", ha="center", va="center")

        ax3.set_xlabel("Minimum Interatomic Distance (Å)")
        ax3.set_ylabel("Number of Atomic Contacts")
        ax3.set_title("C. Interatomic Distance Distribution by TERT Domain", fontweight="bold")

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    finally:
        plt.close(fig)

    if not output_path.exists() or output_path.stat().st_size <= 1024:
        raise RuntimeError(f"Failed to generate valid figure at {output_path}")

    return output_path
```

---

## 5. Verification Method

To independently verify the proposed remediation:

1. **Verify Complete Elimination of Hardcoded Lists**:
   After applying the patch, verify that grep returns 0 hits for hardcoded lists in Figure 4:
   ```bash
   grep -n "atom_counts = \[" src/telotrade/figures.py
   # Expected: 0 hits (replaced with dynamic list comprehension/append)

   grep -n "res_pair_counts = \[" src/telotrade/figures.py
   # Expected: 0 hits

   grep -n "ten_dists = pd.Series" src/telotrade/figures.py
   # Expected: 0 hits (all synthetic series removed)
   ```

2. **Verify Dynamic Sensitivity to Table Modification**:
   Execute a test script that appends a synthetic contact to `interface_contacts_7tre_3.5A.csv`:
   - Verify that the height of bar 1 at 3.5 Å increases from 9 to 10.
   - Verify that Panel B and Panel C reflect altered counts and distributions.

3. **Verify Resource Cleanup and Artifact Integrity**:
   Execute:
   ```bash
   python -c "from telotrade.figures import generate_figure4_cutoff_robustness; p = generate_figure4_cutoff_robustness('results/tables', 'results/figures/figure4_cutoff_robustness.png'); print('Size:', p.stat().st_size)"
   ```
   Confirm file size > 300,000 bytes and zero memory/file-handle leaks.
