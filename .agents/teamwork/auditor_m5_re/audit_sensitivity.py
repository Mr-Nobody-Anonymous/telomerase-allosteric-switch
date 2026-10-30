"""
Forensic Dynamic Sensitivity Test Suite for src/telotrade/figures.py.
Milestone 5 Remediation Re-Audit.

Tests that every figure generator dynamically ingests input CSV tables
and that modifying underlying data produces immediate, measurable changes
in the resulting graphical objects.
"""

import hashlib
import os
import shutil
import tempfile
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from telotrade.figures import (
    generate_figure1_interface_contacts,
    generate_figure2_cross_state_rmsd,
    generate_figure3_tpp1_conformational_shifts,
    generate_figure4_cutoff_robustness,
)

def test_dynamic_sensitivity():
    real_tables_dir = Path("results/tables").resolve()
    assert real_tables_dir.exists(), f"Missing tables dir: {real_tables_dir}"

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        tables_copy = tmp_path / "tables"
        shutil.copytree(real_tables_dir, tables_copy)
        out_dir = tmp_path / "figures"
        out_dir.mkdir()

        results = {}

        # -------------------------------------------------------------
        # Test 1: Figure 2 Sensitivity (cross_state_rmsd.csv)
        # -------------------------------------------------------------
        base_fig2 = out_dir / "fig2_base.png"
        generate_figure2_cross_state_rmsd(tables_copy, base_fig2)
        base_fig2_hash = hashlib.sha256(base_fig2.read_bytes()).hexdigest()

        # Mutate cross_state_rmsd.csv
        cross_csv = tables_copy / "cross_state_rmsd.csv"
        df_cross = pd.read_csv(cross_csv)
        orig_val = df_cross.loc[0, "rmsd_angstrom"]
        df_cross.loc[0, "rmsd_angstrom"] = 9.999  # Massive perturbation
        df_cross.to_csv(cross_csv, index=False)

        mut_fig2 = out_dir / "fig2_mutated.png"
        generate_figure2_cross_state_rmsd(tables_copy, mut_fig2)
        mut_fig2_hash = hashlib.sha256(mut_fig2.read_bytes()).hexdigest()

        fig2_sensitive = (base_fig2_hash != mut_fig2_hash)
        results["Figure 2 (cross_state_rmsd)"] = {
            "sensitive": fig2_sensitive,
            "details": f"Original RMSD: {orig_val:.3f} -> Mutated: 9.999, Base hash != Mutated hash: {fig2_sensitive}",
        }
        # Restore
        shutil.copyfile(real_tables_dir / "cross_state_rmsd.csv", cross_csv)

        # -------------------------------------------------------------
        # Test 2: Figure 4 Sensitivity (interface_contacts_7tre_3.5A.csv)
        # -------------------------------------------------------------
        base_fig4 = out_dir / "fig4_base.png"
        generate_figure4_cutoff_robustness(tables_copy, base_fig4)
        base_fig4_hash = hashlib.sha256(base_fig4.read_bytes()).hexdigest()

        # Mutate 3.5A contacts by duplicating rows (simulating 50 new contacts)
        c35_csv = tables_copy / "interface_contacts_7tre_3.5A.csv"
        df_c35 = pd.read_csv(c35_csv)
        orig_len = len(df_c35)
        df_c35_mut = pd.concat([df_c35, df_c35, df_c35], ignore_index=True)
        df_c35_mut.to_csv(c35_csv, index=False)

        mut_fig4 = out_dir / "fig4_mutated.png"
        generate_figure4_cutoff_robustness(tables_copy, mut_fig4)
        mut_fig4_hash = hashlib.sha256(mut_fig4.read_bytes()).hexdigest()

        fig4_sensitive = (base_fig4_hash != mut_fig4_hash)
        results["Figure 4 (contact cutoffs)"] = {
            "sensitive": fig4_sensitive,
            "details": f"Original 3.5A contacts: {orig_len} -> Mutated: {len(df_c35_mut)}, Base hash != Mutated hash: {fig4_sensitive}",
        }
        # Restore
        shutil.copyfile(real_tables_dir / "interface_contacts_7tre_3.5A.csv", c35_csv)

        # -------------------------------------------------------------
        # Test 3: Figure 3 Sensitivity (tpp1_bound_vs_unbound_rmsd.csv)
        # -------------------------------------------------------------
        base_fig3 = out_dir / "fig3_base.png"
        generate_figure3_tpp1_conformational_shifts(tables_copy, base_fig3)
        base_fig3_hash = hashlib.sha256(base_fig3.read_bytes()).hexdigest()

        tpp1_csv = tables_copy / "tpp1_bound_vs_unbound_rmsd.csv"
        df_tpp1 = pd.read_csv(tpp1_csv)
        # Mutate Lys170 displacement
        idx_170 = df_tpp1[df_tpp1["auth_seq_id"] == 170].index[0]
        orig_disp = df_tpp1.loc[idx_170, "ca_distance_7tre_vs_2i46_angstrom"]
        df_tpp1.loc[idx_170, "ca_distance_7tre_vs_2i46_angstrom"] = 12.345
        df_tpp1.to_csv(tpp1_csv, index=False)

        mut_fig3 = out_dir / "fig3_mutated.png"
        generate_figure3_tpp1_conformational_shifts(tables_copy, mut_fig3)
        mut_fig3_hash = hashlib.sha256(mut_fig3.read_bytes()).hexdigest()

        fig3_sensitive = (base_fig3_hash != mut_fig3_hash)
        results["Figure 3 (tpp1 displacement)"] = {
            "sensitive": fig3_sensitive,
            "details": f"Original Lys170 displacement: {orig_disp:.3f} -> Mutated: 12.345, Base hash != Mutated hash: {fig3_sensitive}",
        }
        # Restore
        shutil.copyfile(real_tables_dir / "tpp1_bound_vs_unbound_rmsd.csv", tpp1_csv)

        # -------------------------------------------------------------
        # Test 4: Figure 1 Sensitivity (polar_interactions_7tre.csv)
        # -------------------------------------------------------------
        base_fig1 = out_dir / "fig1_base.png"
        generate_figure1_interface_contacts(tables_copy, base_fig1)
        base_fig1_hash = hashlib.sha256(base_fig1.read_bytes()).hexdigest()

        polar_csv = tables_copy / "polar_interactions_7tre.csv"
        df_polar = pd.read_csv(polar_csv)
        orig_polar_len = len(df_polar)
        # Drop half the interactions
        df_polar.iloc[: len(df_polar) // 2].to_csv(polar_csv, index=False)

        mut_fig1 = out_dir / "fig1_mutated.png"
        generate_figure1_interface_contacts(tables_copy, mut_fig1)
        mut_fig1_hash = hashlib.sha256(mut_fig1.read_bytes()).hexdigest()

        fig1_sensitive = (base_fig1_hash != mut_fig1_hash)
        results["Figure 1 (polar interactions)"] = {
            "sensitive": fig1_sensitive,
            "details": f"Original polar count: {orig_polar_len} -> Mutated: {orig_polar_len // 2}, Base hash != Mutated hash: {fig1_sensitive}",
        }

        print("=== DYNAMIC TABLE SENSITIVITY RESULTS ===")
        all_passed = True
        for name, data in results.items():
            status = "PASS" if data["sensitive"] else "FAIL"
            if not data["sensitive"]:
                all_passed = False
            print(f"[{status}] {name}: {data['details']}")

        return all_passed

if __name__ == "__main__":
    success = test_dynamic_sensitivity()
    import sys
    sys.exit(0 if success else 1)
