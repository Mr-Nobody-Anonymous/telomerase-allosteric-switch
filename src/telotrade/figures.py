"""
Publication Figures Generator for Telotrade.

Generates 4 publication-quality figures (PNG, 300 DPI, >= 1024 bytes) in results/figures/:
1. figure1_interface_contacts.png:
   TERT-TPP1 interface contact map and polar/hydrophobic interaction network.
2. figure2_cross_state_rmsd.png:
   Domain-by-domain RMSD comparison and catalytic core rigidity across 7TRE, 7TRD, 7TRF, and 2I46.
3. figure3_tpp1_conformational_shifts.png:
   TPP1 Cα displacement profile (DE-loop remodeling vs Leu183 rigidity).
4. figure4_cutoff_robustness.png:
   Interface contact count and persistence distribution across distance cutoffs.
"""

import os
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def set_publication_style() -> None:
    """Set clean, publication-grade matplotlib styling."""
    plt.style.use("default")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 14,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": "--",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def save_figure_safely(fig: plt.Figure, output_path: Path | str, dpi: int = 300) -> Path:
    """
    Atomically save a matplotlib figure and guarantee non-zero file size (> 1024 bytes).

    Guarantees:
    1. Writes to a temporary file on the same filesystem/directory.
    2. Flushes Python stream buffers and forces OS disk cache commit via os.fsync.
    3. Validates file size > 1024 bytes before commitment.
    4. Atomically replaces target path via os.replace (preventing 0-byte visibility).
    5. Cleans up temporary file if any error occurs.
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}_{time.time_ns()}.tmp")

    try:
        with open(tmp_path, "wb") as fh:
            fig.savefig(fh, format="png", dpi=dpi, bbox_inches="tight")
            fh.flush()
            os.fsync(fh.fileno())

        size = tmp_path.stat().st_size
        if size <= 1024:
            raise OSError(
                f"Generated figure {output_path.name} is unexpectedly small ({size} bytes, expected > 1024 bytes)"
            )

        os.replace(tmp_path, output_path)

        committed_size = output_path.stat().st_size
        if committed_size <= 1024:
            raise OSError(
                f"Committed figure {output_path.name} has invalid size ({committed_size} bytes)"
            )
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass

    return output_path



def generate_figure1_interface_contacts(
    tables_dir: Path | str,
    output_path: Path | str,
) -> Path:
    """
    Generate Figure 1: TERT-TPP1 interface contact map and polar/hydrophobic interaction network.

    Panels:
    - (A) 2D Contact Map of 7TRE TERT vs TPP1 residues color-coded by interaction type.
    - (B) Polar and Salt-Bridge Interaction Network with key interatomic distances.
    """
    set_publication_style()
    tables_dir = Path(tables_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    contacts_path = tables_dir / "interface_contacts_7tre_4.5A.csv"
    polar_path = tables_dir / "polar_interactions_7tre.csv"

    if not contacts_path.exists():
        contacts_path = tables_dir / "interface_contacts_4.5A.csv"
    if not polar_path.exists():
        polar_path = tables_dir / "polar_interactions.csv"

    df_contacts = pd.read_csv(contacts_path)
    df_polar = pd.read_csv(polar_path)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5), gridspec_kw={"width_ratios": [1.2, 1]})
    try:
        fig.suptitle(
            "Figure 1: High-Resolution Anatomy of Human TERT–TPP1 Recruitment Interface (PDB 7TRE)",
            fontweight="bold",
            y=0.98,
        )

        # Panel A: Contact Map
        type_colors = {
            "Salt Bridge": "#D90429",       # Deep Red
            "Polar/H-Bond": "#0077B6",      # Strong Blue
            "Hydrophobic": "#2B9348",       # Forest Green
            "Van der Waals": "#ADB5BD",     # Slate Gray
        }
        type_markers = {
            "Salt Bridge": "s",
            "Polar/H-Bond": "^",
            "Hydrophobic": "D",
            "Van der Waals": "o",
        }

        # Group by residue pairs to get minimum distance and interaction type
        pair_summary = (
            df_contacts.groupby(["tert_resnum", "tert_resname", "tert_domain", "tpp1_resnum", "tpp1_resname", "interaction_type"])
            .agg(min_dist=("distance_angstrom", "min"), is_tel=("is_tel_patch", "any"))
            .reset_index()
        )

        # Sort so prominent interactions draw on top
        draw_order = ["Van der Waals", "Hydrophobic", "Polar/H-Bond", "Salt Bridge"]
        pair_summary["order"] = pair_summary["interaction_type"].map(lambda x: draw_order.index(x) if x in draw_order else 0)
        pair_summary = pair_summary.sort_values("order")

        for itype in draw_order:
            subset = pair_summary[pair_summary["interaction_type"] == itype]
            if not subset.empty:
                ax1.scatter(
                    subset["tpp1_resnum"],
                    subset["tert_resnum"],
                    c=type_colors.get(itype, "#6C757D"),
                    marker=type_markers.get(itype, "o"),
                    s=70 if itype != "Van der Waals" else 35,
                    alpha=0.9 if itype != "Van der Waals" else 0.5,
                    edgecolors="black" if itype != "Van der Waals" else "none",
                    linewidths=0.8,
                    label=f"{itype} (n={len(subset)})",
                )

        # Annotate Key TEL-Patch residues on TPP1 (x-axis)
        tel_patch_x = [168, 169, 170, 171, 183, 212]
        for x in tel_patch_x:
            ax1.axvline(x=x, color="#EF233C", linestyle=":", alpha=0.3, zorder=0)

        def _get_pair_dist(t_res: int, p_res: int) -> float | None:
            m = pair_summary[(pair_summary["tert_resnum"] == t_res) & (pair_summary["tpp1_resnum"] == p_res)]
            return float(m["min_dist"].min()) if not m.empty else None

        # Highlight Key Salt Bridge
        sb_dist = _get_pair_dist(774, 168)
        sb_text = f"Arg774 -- Glu168\n({sb_dist:.2f} Å Salt Bridge)" if sb_dist is not None else "Arg774 -- Glu168\n(Salt Bridge)"
        ax1.annotate(
            sb_text,
            xy=(168, 774),
            xytext=(145, 830),
            arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
        )

        # Highlight Lys78 -- Glu215
        latch_dist = _get_pair_dist(78, 215)
        latch_text = f"Lys78 -- Glu215\n({latch_dist:.2f} Å Latch)" if latch_dist is not None else "Lys78 -- Glu215\n(Latch)"
        ax1.annotate(
            latch_text,
            xy=(215, 78),
            xytext=(185, 130),
            arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
        )

        # Highlight Leu769 -- Leu183
        anchor_dist = _get_pair_dist(769, 183)
        anchor_text = (
            f"Leu769 -- Leu183\n({anchor_dist:.2f} Å Hydrophobic Anchor)"
            if anchor_dist is not None
            else "Leu769 -- Leu183\n(Hydrophobic Anchor)"
        )
        ax1.annotate(
            anchor_text,
            xy=(183, 769),
            xytext=(195, 710),
            arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 1},
        )

        ax1.set_xlabel("TPP1 OB-fold Residue Index (Canonical UniProt Q96AP0)")
        ax1.set_ylabel("TERT Residue Index (Canonical UniProt O14746)")
        ax1.set_title("A. Bipartite Interface Contact Map (4.5 Å Enclosing Shell)", fontweight="bold")
        ax1.legend(loc="upper left", framealpha=0.9)
        ax1.set_xlim(85, 245)
        ax1.set_ylim(20, 860)

        # Panel B: Polar & Salt Bridge Distance Hierarchy
        df_polar_unique = (
            df_polar.groupby(["tert_resnum", "tert_resname", "tert_atom", "tpp1_resnum", "tpp1_resname", "tpp1_atom", "interaction_type"])
            .agg(dist=("distance_angstrom", "min"), is_tel=("is_tel_patch", "first"), domain=("tert_domain", "first"))
            .reset_index()
            .sort_values("dist", ascending=True)
        )

        labels = [
            f"{r['tert_resname']}{r['tert_resnum']}:{r['tert_atom']} ↔ {r['tpp1_resname']}{r['tpp1_resnum']}:{r['tpp1_atom']}"
            for _, r in df_polar_unique.iterrows()
        ]
        distances = df_polar_unique["dist"].values
        bar_colors = [
            "#D90429" if "Salt Bridge" in t else "#0077B6"
            for t in df_polar_unique["interaction_type"]
        ]

        y_pos = np.arange(len(labels))
        bars = ax2.barh(y_pos, distances, color=bar_colors, edgecolor="black", height=0.65, alpha=0.85)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(labels, fontsize=8)
        ax2.invert_yaxis()  # Shortest at top
        ax2.set_xlabel("Interatomic Distance (Å)")
        ax2.set_title("B. Polar Network & Salt-Bridge Distance Hierarchy", fontweight="bold")
        ax2.axvline(x=3.5, color="gray", linestyle="--", alpha=0.7, label="Standard Polar Cutoff (3.5 Å)")
        ax2.axvline(x=4.0, color="orange", linestyle=":", alpha=0.7, label="Cryo-EM Distance Gate (4.0 Å)")

        for bar, d, itype in zip(bars, distances, df_polar_unique["interaction_type"]):
            ax2.text(
                d + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{d:.2f} Å ({'SB' if 'Salt' in itype else 'HB'})",
                va="center",
                fontsize=7.5,
                fontweight="bold",
            )

        ax2.set_xlim(0, 4.6)
        ax2.legend(loc="lower right", framealpha=0.9, fontsize=8)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path


def generate_figure2_cross_state_rmsd(
    tables_dir: Path | str,
    output_path: Path | str,
) -> Path:
    """
    Generate Figure 2: Domain-by-domain RMSD comparison and catalytic core rigidity
    across 7TRE, 7TRD, 7TRF, and 2I46.

    Panels:
    - (A) Grouped Bar Chart of Domain RMSD across structural states.
    - (B) Catalytic Core Rigidity vs Peripheral Domain Mobility.
    """
    set_publication_style()
    tables_dir = Path(tables_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    domain_path = tables_dir / "domain_rmsd.csv"
    if not domain_path.exists():
        domain_path = tables_dir / "domain_rmsds.csv"
    cross_path = tables_dir / "cross_state_rmsd.csv"

    df_domain = pd.read_csv(domain_path) if domain_path.exists() else pd.DataFrame()
    df_cross = pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    try:
        fig.suptitle(
            "Figure 2: Cross-State Structural Conservation and Catalytic Core Invariance",
            fontweight="bold",
            y=0.98,
        )

        # Panel A: Domain RMSD Comparison (7TRE vs 7TRD and 7TRE vs 7TRF)
        # Target standard domains: TEN, TRBD, RT, IFD_TRAP, CTE, Catalytic Core
        domains = ["TEN\n(7-178)", "TRBD\n(326-540)", "RT\n(541-935)", "IFD_TRAP\n(755-785)", "CTE\n(936-1132)", "Catalytic Core\n(326-1132)"]
        x = np.arange(len(domains))
        width = 0.35

        domain_keys = ["TEN", "TRBD", "RT", "IFD_TRAP", "CTE", "Catalytic Core"]
        if not df_domain.empty:
            df_7trd = df_domain[df_domain["reference_structure"] == "7TRD"].set_index("domain_name")
            df_7trf = df_domain[df_domain["reference_structure"] == "7TRF"].set_index("domain_name")

            rmsds_7trd = [float(df_7trd.loc[k, "independent_rmsd_angstrom"]) if k in df_7trd.index else 0.0 for k in domain_keys]
            rmsds_7trf = [float(df_7trf.loc[k, "independent_rmsd_angstrom"]) if k in df_7trf.index else 0.0 for k in domain_keys]
        else:
            rmsds_7trd = [0.0] * len(domain_keys)
            rmsds_7trf = [0.0] * len(domain_keys)

        rects1 = ax1.bar(x - width / 2, rmsds_7trd, width, label="7TRE vs 7TRD (Unbound Core)", color="#1D3557", edgecolor="black")
        rects2 = ax1.bar(x + width / 2, rmsds_7trf, width, label="7TRE vs 7TRF (+ H2A/H2B)", color="#457B9D", edgecolor="black")

        ax1.set_ylabel("Cα RMSD (Å)")
        ax1.set_title("A. Domain-Specific Cα RMSD (Kabsch SVD)", fontweight="bold")
        ax1.set_xticks(x)
        ax1.set_xticklabels(domains)
        ax1.axhline(y=0.4, color="#E63946", linestyle="--", alpha=0.7, label="Cryo-EM Coordinate Uncertainty (~0.4 Å)")
        max_rmsd_a = max(rmsds_7trd + rmsds_7trf) if (rmsds_7trd and rmsds_7trf) else 0.4
        ax1.set_ylim(0, max(max_rmsd_a * 1.25, 0.55))
        ax1.legend(loc="upper right", framealpha=0.9)

        for rect in rects1:
            h = rect.get_height()
            ax1.annotate(f"{h:.3f}", xy=(rect.get_x() + rect.get_width() / 2, h), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7.5)
        for rect in rects2:
            h = rect.get_height()
            ax1.annotate(f"{h:.3f}", xy=(rect.get_x() + rect.get_width() / 2, h), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=7.5)

        # Panel B: Macromolecular Comparison across all complexes (including 2I46 TPP1)
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

        y_pos = np.arange(len(comparisons))
        bars2 = ax2.barh(y_pos, comp_rmsds, color=comp_colors, edgecolor="black", height=0.6, alpha=0.9)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(comparisons, fontsize=8.5)
        ax2.invert_yaxis()
        ax2.set_xlabel("Cα RMSD (Å)")
        ax2.set_title("B. Catalytic Core Invariance vs TPP1 OB-fold Remodeling", fontweight="bold")
        ax2.axvline(x=0.5, color="gray", linestyle=":", alpha=0.7, label="Rigid Threshold (<0.5 Å)")
        ax2.axvline(x=1.5, color="red", linestyle=":", alpha=0.7, label="Extensive Plasticity (>1.5 Å)")

        for bar, val in zip(bars2, comp_rmsds):
            ax2.text(
                val + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.3f} Å",
                va="center",
                fontsize=8.5,
                fontweight="bold",
            )

        max_rmsd_b = max(comp_rmsds) if comp_rmsds else 2.0
        ax2.set_xlim(0, max(max_rmsd_b + 0.35, 2.3))
        ax2.legend(loc="lower right", framealpha=0.9)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path


def generate_figure3_tpp1_conformational_shifts(
    tables_dir: Path | str,
    output_path: Path | str,
) -> Path:
    """
    Generate Figure 3: TPP1 Cα displacement profile highlighting DE-loop plasticity
    (Glu168–Glu171 peaking at Lys170 at 4.29 Å) vs Leu183 rigidity (0.51 Å).
    """
    set_publication_style()
    tables_dir = Path(tables_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tpp1_path = tables_dir / "tpp1_bound_vs_unbound_rmsd.csv"
    if not tpp1_path.exists():
        tpp1_path = tables_dir / "tpp1_conformational_shifts.csv"

    df_tpp1 = pd.read_csv(tpp1_path)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    try:
        fig.suptitle(
            "Figure 3: Induced-Fit Plasticity of the TPP1 TEL Patch: DE-Loop Remodeling vs Invariant Leu183",
            fontweight="bold",
            y=0.98,
        )

        # Panel A: Cα Displacement along sequence
        resnums = df_tpp1["auth_seq_id"].values
        displacements = df_tpp1["ca_distance_7tre_vs_2i46_angstrom"].values
        is_tel = df_tpp1["is_tel_patch"].values

        # Highlight background zones
        ax1.axhspan(0, 1.0, color="#E8F5E9", alpha=0.5, label="Rigid Anchor (< 1.0 Å)")
        ax1.axhspan(1.0, 2.5, color="#FFF3E0", alpha=0.5, label="Moderate Shift (1.0 – 2.5 Å)")
        ax1.axhspan(2.5, 7.5, color="#FFEBEE", alpha=0.5, label="Extensive Remodeling (> 2.5 Å)")

        # Plot continuous line and points
        ax1.plot(resnums, displacements, color="#1D3557", lw=1.5, zorder=2)
        ax1.scatter(
            resnums[~is_tel],
            displacements[~is_tel],
            color="#457B9D",
            s=25,
            alpha=0.7,
            zorder=3,
            label="Non-TEL Patch Residues",
        )
        ax1.scatter(
            resnums[is_tel],
            displacements[is_tel],
            color="#D90429",
            s=80,
            edgecolors="black",
            lw=1.2,
            zorder=4,
            label="TEL-Patch Residues (168-171, 180, 183, 212)",
        )

        # Shaded band for DE-loop (168-171)
        ax1.axvspan(167.5, 171.5, color="#EF233C", alpha=0.2, label="DE-Loop Window (168-171)")

        disp_map = (
            df_tpp1.set_index("auth_seq_id")["ca_distance_7tre_vs_2i46_angstrom"].to_dict()
            if "auth_seq_id" in df_tpp1.columns and "ca_distance_7tre_vs_2i46_angstrom" in df_tpp1.columns
            else {}
        )

        # Lys170 (Peak)
        d_170 = disp_map.get(170)
        if d_170 is not None:
            ax1.annotate(
                f"Lys170\n(Peak: {d_170:.2f} Å)",
                xy=(170, d_170),
                xytext=(160, max(d_170 + 1.5, 5.8)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
                fontsize=8.5,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
            )

        # Glu169
        d_169 = disp_map.get(169)
        if d_169 is not None:
            ax1.annotate(
                f"Glu169\n({d_169:.2f} Å)",
                xy=(169, d_169),
                xytext=(145, max(d_169 + 0.7, 4.5)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#FFF0F3", "ec": "#D90429", "lw": 0.8},
            )

        # Glu168
        d_168 = disp_map.get(168)
        if d_168 is not None:
            ax1.annotate(
                f"Glu168\n({d_168:.2f} Å)",
                xy=(168, d_168),
                xytext=(145, max(d_168 + 0.9, 2.2)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#FFF0F3", "ec": "#D90429", "lw": 0.8},
            )

        # Leu183 (Invariant Anchor)
        d_183 = disp_map.get(183)
        if d_183 is not None:
            ax1.annotate(
                f"Leu183 Rigid Anchor\n({d_183:.2f} Å Displacement)",
                xy=(183, d_183),
                xytext=(195, 2.8),
                arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.5, "headwidth": 6},
                fontsize=8.5,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.3", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 1},
            )

        # Leu212
        d_212 = disp_map.get(212)
        if d_212 is not None:
            ax1.annotate(
                f"Leu212\n({d_212:.2f} Å)",
                xy=(212, d_212),
                xytext=(218, 2.2),
                arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 0.8},
            )

        ax1.set_ylabel("Cα Displacement (Å)\n[7TRE Bound vs 2I46 Chain A Unbound]")
        ax1.set_title("A. Per-Residue Cα Displacement across the TPP1 OB-fold Domain", fontweight="bold")
        ax1.set_ylim(-0.2, 7.5)
        ax1.legend(loc="upper right", framealpha=0.9, fontsize=8)

        # Panel B: B-Factor Comparison (Cryo-EM 7TRE vs Crystal 2I46 Chain A)
        b_7tre = df_tpp1["local_b_factor_7tre"].values
        b_2i46 = df_tpp1["local_b_factor_2i46"].values

        ax2.plot(resnums, b_7tre, color="#0077B6", lw=1.2, label="7TRE Cryo-EM (TPP1-Bound)")
        ax2.plot(resnums, b_2i46, color="#E76F51", lw=1.2, linestyle="--", label="2I46 X-ray Crystal (Unbound)")
        ax2.axvspan(167.5, 171.5, color="#EF233C", alpha=0.2)
        ax2.set_xlabel("TPP1 Author Residue Number (Canonical UniProt Q96AP0)")
        ax2.set_ylabel("B-Factor (Å²)")
        ax2.set_title("B. Local Coordinate Temperature Factor Profiles (Cα B-iso)", fontweight="bold")
        ax2.legend(loc="upper right", framealpha=0.9, fontsize=8)
        if len(resnums) > 0:
            ax2.set_xlim(resnums.min() - 2, resnums.max() + 2)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path


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
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path


def generate_all_figures(
    tables_dir: Path | str,
    figures_dir: Path | str,
) -> dict[str, Path]:
    """
    Generate all 4 publication figures in figures_dir.

    Returns:
        dict[str, Path]: Mapping of figure filename to output path.
    """
    tables_dir = Path(tables_dir)
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig1_path = figures_dir / "figure1_interface_contacts.png"
    fig2_path = figures_dir / "figure2_cross_state_rmsd.png"
    fig3_path = figures_dir / "figure3_tpp1_conformational_shifts.png"
    fig4_path = figures_dir / "figure4_cutoff_robustness.png"

    generate_figure1_interface_contacts(tables_dir, fig1_path)
    generate_figure2_cross_state_rmsd(tables_dir, fig2_path)
    generate_figure3_tpp1_conformational_shifts(tables_dir, fig3_path)
    generate_figure4_cutoff_robustness(tables_dir, fig4_path)

    figures = {
        "figure1_interface_contacts.png": fig1_path,
        "figure2_cross_state_rmsd.png": fig2_path,
        "figure3_tpp1_conformational_shifts.png": fig3_path,
        "figure4_cutoff_robustness.png": fig4_path,
    }

    # Verify all generated artifacts are present and valid (> 1024 bytes)
    for fname, fpath in figures.items():
        if not fpath.exists():
            raise FileNotFoundError(f"Expected publication figure {fname} was not found at {fpath}")
        size = fpath.stat().st_size
        if size <= 1024:
            raise OSError(f"Figure {fname} is corrupted or empty ({size} bytes, expected > 1024 bytes)")

    return figures
