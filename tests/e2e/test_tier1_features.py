"""Tier 1: Feature Coverage E2E Tests.

Verifies all core features (Features 1 to 18) declared in PROJECT.md:
- Raw data acquisition, SHA-256 integrity, data manifest
- Sequence-aware chain mapping and expression tag offsets
- 2I46 dimer asymmetry handling
- Interface contact networks at 3.5, 4.0, 4.5 A
- Buried surface area (BSA)
- Polar interactions and salt bridges
- Hydrophobic packing and TEL-patch mapping
- Cross-state SVD Kabsch superposition and RMSDs
- Standardized CSV tables and report generation
"""

import hashlib

import numpy as np
import pandas as pd
import pytest

from tests.e2e.conftest import (
    extract_ca_atoms,
    extract_heavy_atoms,
    find_contacts,
    kabsch_fit,
)


def compute_distance(coord1, coord2):
    return float(np.linalg.norm(np.array(coord1) - np.array(coord2)))


@pytest.mark.tier1
class TestTier1RawDataAndManifest:
    """Features 1, 2, 3: Raw data acquisition, SHA-256 verification, manifest."""

    def test_f1_all_cif_files_exist(self, resolve_cif):
        """Feature 1: Verify all 5 required mmCIF structures exist and are non-empty."""
        for pid in ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]:
            cif_path = resolve_cif(pid)
            assert cif_path.exists(), f"CIF file for {pid} does not exist at {cif_path}"
            assert cif_path.stat().st_size > 100_000, f"CIF file {pid} is unexpectedly small"

    def test_f1_sha256_checksums(self, resolve_cif, ref_data):
        """Feature 1: Verify SHA-256 hashes of all 5 files match authoritative values."""
        for pid, meta in ref_data["structures"].items():
            cif_path = resolve_cif(pid)
            h = hashlib.sha256()
            with open(cif_path, "rb") as f:
                while chunk := f.read(65536):
                    h.update(chunk)
            computed = h.hexdigest()
            expected = meta["sha256"]
            assert computed == expected, f"SHA-256 mismatch for {pid}: got {computed}, expected {expected}"

    def test_f2_data_manifest_schema(self, project_root):
        """Feature 2: Verify data_manifest.csv exists and contains required schema."""
        manifest_path = project_root / "data_manifest.csv"
        alt_path = project_root / "data" / "manifest" / "data_manifest.csv"
        target_path = manifest_path if manifest_path.exists() else alt_path
        if not target_path.exists():
            pytest.skip("data_manifest.csv not yet generated (M1 in progress)")

        df = pd.read_csv(target_path)
        required_cols = [
            "pdb_id", "experimental_method", "resolution", "sha256",
            "source_url", "chain_mappings"
        ]
        for col in required_cols:
            assert col in df.columns, f"Missing required column '{col}' in manifest"
        assert len(df) >= 5, "Manifest must contain all 5 target structures"

    def test_f2_manifest_resolutions(self, project_root, ref_data):
        """Feature 2: Verify manifest contains accurate experimental resolutions."""
        manifest_path = project_root / "data_manifest.csv"
        alt_path = project_root / "data" / "manifest" / "data_manifest.csv"
        target_path = manifest_path if manifest_path.exists() else alt_path
        if not target_path.exists():
            pytest.skip("data_manifest.csv not yet generated (M1 in progress)")

        df = pd.read_csv(target_path)
        df.set_index("pdb_id", inplace=True)
        for pid in ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]:
            assert pid in df.index, f"{pid} not found in manifest"
            res = float(df.loc[pid, "resolution"])
            expected_res = ref_data["structures"][pid]["resolution"]
            assert abs(res - expected_res) < 0.20, f"Resolution mismatch for {pid}: got {res}, expected {expected_res}"

    def test_f3_validate_module_contract(self):
        """Feature 3: Verify validate module interface contracts."""
        try:
            from telotrade import validate
            assert hasattr(validate, "verify_checksums") or hasattr(validate, "download_structures") or hasattr(validate, "validate_dataset"), \
                "validate module must expose verification functions"
        except ImportError:
            pytest.skip("telotrade.validate not yet implemented (M1 in progress)")


@pytest.mark.tier1
class TestTier1SequenceAndChainMapping:
    """Features 4, 5, 6: Sequence-aware chain mapping, tag offsets, 2I46 asymmetry."""

    def test_f4_author_vs_label_chain_identifiers(self, gemmi_structures):
        """Feature 4: Verify 7TRE macromolecular chains have correct author vs label IDs."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        model = st[0]
        # In 7TRE: TERT is author chain 'A', label chain 'B'
        tert_chain = model["A"]
        assert tert_chain is not None, "TERT author chain 'A' must exist in 7TRE"
        assert tert_chain[0].subchain == "B", "TERT CIF label chain must be 'B'"

        # TPP1 is author chain 'C', label chain 'D'
        tpp1_chain = model["C"]
        assert tpp1_chain is not None, "TPP1 author chain 'C' must exist in 7TRE"
        assert tpp1_chain[0].subchain == "D", "TPP1 CIF label chain must be 'D'"

    def test_f4_homologous_tert_chain_lengths(self, gemmi_structures):
        """Feature 4: Verify 7TRE, 7TRD, 7TRF each contain exactly 934 modeled TERT residues."""
        for pid in ["7TRE", "7TRD", "7TRF"]:
            st = gemmi_structures.get(pid)
            assert st is not None
            ca_dict = extract_ca_atoms(st[0]["A"])
            assert len(ca_dict) == 934, f"{pid} TERT must contain exactly 934 modeled C-alpha atoms"

    def test_f5_tert_expression_tag_offset(self, gemmi_structures):
        """Feature 5: Verify TERT expression tag offset +35 between label_seq and auth_seq."""
        st = gemmi_structures.get("7TRE")
        model = st[0]
        tert_chain = model["A"]
        first_res = tert_chain[0]
        # First modeled residue is author 7, label 42 (offset +35)
        assert first_res.seqid.num == 7, "First modeled TERT author residue must be 7"
        assert first_res.label_seq == 42, "First modeled TERT label_seq must be 42 (+35 Twin-Strep tag offset)"
        assert first_res.label_seq - first_res.seqid.num == 35

    def test_f5_tpp1_expression_tag_offset(self, gemmi_structures):
        """Feature 5: Verify 7TRE TPP1 expression tag offset +86."""
        st = gemmi_structures.get("7TRE")
        model = st[0]
        tpp1_chain = model["C"]
        first_res = tpp1_chain[0]
        # First modeled TPP1 residue in 7TRE is author 90, label 4 (offset = 90 - 4 = 86 due to GAG prefix)
        assert first_res.seqid.num == 90, "First modeled TPP1 author residue in 7TRE must be 90"
        assert first_res.label_seq == 4, "First modeled TPP1 label_seq in 7TRE must be 4"
        assert first_res.seqid.num - first_res.label_seq == 86

    def test_f6_2i46_dimer_asymmetry_chain_a_resolved(self, gemmi_structures):
        """Feature 6: Verify 2I46 Chain A has TEL-patch residues 168-171 fully resolved."""
        st = gemmi_structures.get("2I46")
        assert st is not None
        ca_dict = extract_ca_atoms(st[0]["A"])
        for res_num in [168, 169, 170, 171]:
            assert res_num in ca_dict, f"2I46 Chain A must resolve TEL-patch residue {res_num}"

    def test_f6_2i46_dimer_asymmetry_chain_b_disordered(self, gemmi_structures):
        """Feature 6: Verify 2I46 Chain B has TEL-patch DE-loop (168-171) disordered/unmodeled."""
        st = gemmi_structures.get("2I46")
        assert st is not None
        ca_dict = extract_ca_atoms(st[0]["B"])
        for res_num in [168, 169, 170, 171]:
            assert res_num not in ca_dict, f"2I46 Chain B must NOT resolve disordered residue {res_num}"


@pytest.mark.tier1
class TestTier1InterfaceGeometry:
    """Features 7, 8, 9, 10: Multi-cutoff contacts, BSA, polar networks, TEL patch."""

    def test_f7_contact_networks_cutoff_progression(self, gemmi_structures):
        """Feature 7: Verify contact network size strictly increases across 3.5, 4.0, 4.5 A."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        contacts_35 = find_contacts(tert_heavy, tpp1_heavy, 3.5)
        contacts_40 = find_contacts(tert_heavy, tpp1_heavy, 4.0)
        contacts_45 = find_contacts(tert_heavy, tpp1_heavy, 4.5)

        assert len(contacts_35) == 9, f"Expected 9 atomic contacts at 3.5 A, got {len(contacts_35)}"
        assert len(contacts_40) == 46, f"Expected 46 atomic contacts at 4.0 A, got {len(contacts_40)}"
        assert len(contacts_45) == 123, f"Expected 123 atomic contacts at 4.5 A, got {len(contacts_45)}"
        assert len(contacts_35) < len(contacts_40) < len(contacts_45)

    def test_f8_buried_surface_area_formula_and_value(self, ref_data):
        """Feature 8: Verify buried surface area formula BSA = (SASA_A + SASA_B - SASA_AB) / 2."""
        bsa_ref = ref_data["bsa"]
        sasa_a = bsa_ref["sasa_tert"]
        sasa_b = bsa_ref["sasa_tpp1"]
        sasa_ab = bsa_ref["sasa_complex"]
        delta_sasa = sasa_a + sasa_b - sasa_ab
        bsa = delta_sasa / 2.0

        assert abs(delta_sasa - bsa_ref["delta_sasa"]) < 0.1, "Delta SASA calculation mismatch"
        assert abs(bsa - bsa_ref["bsa_interface"]) < 0.1, "Interface BSA calculation mismatch"
        assert 800.0 <= bsa <= 1200.0, f"BSA ({bsa} A^2) must lie in canonical PPI range"

    def test_f9_primary_salt_bridge_arg774_glu168(self, gemmi_structures):
        """Feature 9: Verify primary recruitment salt bridge TERT Arg774 -- TPP1 Glu168."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        arg774 = next((r for r in st[0]["A"] if r.seqid.num == 774 and r.name == "ARG"), None)
        glu168 = next((r for r in st[0]["C"] if r.seqid.num == 168 and r.name == "GLU"), None)
        assert arg774 is not None, "TERT Arg774 must exist"
        assert glu168 is not None, "TPP1 Glu168 must exist"

        min_dist = float("inf")
        for a1 in arg774:
            if a1.name in ["NH1", "NH2", "NE"]:
                for a2 in glu168:
                    if a2.name in ["OE1", "OE2"]:
                        d = compute_distance([a1.pos.x, a1.pos.y, a1.pos.z], [a2.pos.x, a2.pos.y, a2.pos.z])
                        min_dist = min(min_dist, d)
        assert min_dist < 2.45, f"TERT Arg774 -- TPP1 Glu168 salt bridge distance {min_dist} A > 2.45 A"
        assert abs(min_dist - 2.32) < 0.1, f"Expected salt bridge distance ~2.32 A, got {min_dist} A"

    def test_f9_secondary_salt_bridge_lys78_glu215(self, gemmi_structures):
        """Feature 9: Verify secondary electrostatic clamp TERT Lys78 -- TPP1 Glu215."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        lys78 = next((r for r in st[0]["A"] if r.seqid.num == 78 and r.name == "LYS"), None)
        glu215 = next((r for r in st[0]["C"] if r.seqid.num == 215 and r.name == "GLU"), None)
        assert lys78 is not None and glu215 is not None

        nz = next((a for a in lys78 if a.name == "NZ"), None)
        oe = [a for a in glu215 if a.name in ["OE1", "OE2"]]
        assert nz is not None and len(oe) > 0
        min_dist = min(compute_distance([nz.pos.x, nz.pos.y, nz.pos.z], [a.pos.x, a.pos.y, a.pos.z]) for a in oe)
        assert min_dist <= 3.8, f"Lys78 -- Glu215 distance {min_dist} A > 3.8 A"
        assert abs(min_dist - 3.58) < 0.1, f"Expected Lys78 -- Glu215 distance ~3.58 A, got {min_dist} A"

    def test_f10_hydrophobic_packing_groove(self, gemmi_structures):
        """Feature 10: Verify hydrophobic packing between TERT (Leu769, Pro771) and TPP1 (Leu183, Leu212)."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        leu769 = next((r for r in st[0]["A"] if r.seqid.num == 769), None)
        pro771 = next((r for r in st[0]["A"] if r.seqid.num == 771), None)
        leu183 = next((r for r in st[0]["C"] if r.seqid.num == 183), None)
        leu212 = next((r for r in st[0]["C"] if r.seqid.num == 212), None)
        assert all([leu769, pro771, leu183, leu212]), "Hydrophobic core residues must all exist"

        d1 = min(compute_distance([a1.pos.x, a1.pos.y, a1.pos.z], [a2.pos.x, a2.pos.y, a2.pos.z]) for a1 in leu769 for a2 in leu183)
        d2 = min(compute_distance([a1.pos.x, a1.pos.y, a1.pos.z], [a2.pos.x, a2.pos.y, a2.pos.z]) for a1 in pro771 for a2 in leu212)
        assert d1 <= 4.0, f"Leu769 -- Leu183 distance {d1} A > 4.0 A"
        assert d2 <= 4.0, f"Pro771 -- Leu212 distance {d2} A > 4.0 A"


@pytest.mark.tier1
class TestTier1CrossStateStructuralAnalysis:
    """Features 11, 12, 13, 14: Cross-state Kabsch RMSD, TEN-TRAP rearrangements, TPP1 remodeling."""

    def test_f11_kabsch_svd_proper_rotation_determinant(self):
        """Feature 11: Verify Kabsch SVD guarantees proper rotation (det(R) == +1.0)."""
        coords_a = np.array([
            [1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0], [2.0, 0.0, 4.0]
        ], dtype=float)
        theta = np.pi / 4.0
        rot_z = np.array([
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0]
        ])
        coords_b = np.dot(coords_a, rot_z.T) + np.array([10.0, -5.0, 3.0])

        rmsd, r, _t, _ = kabsch_fit(coords_b, coords_a)
        assert rmsd < 1e-10, f"Kabsch RMSD should be ~0 for rigid transform, got {rmsd}"
        assert abs(np.linalg.det(r) - 1.0) < 1e-5, "Rotation matrix determinant must be +1.0 strictly"

    def test_f11_catalytic_core_rigidity_7tre_vs_7trd(self, gemmi_structures):
        """Feature 11: Verify catalytic core (residues >= 326) RMSD between 7TRE and 7TRD is ~0.338 A."""
        st_tre = gemmi_structures.get("7TRE")
        st_trd = gemmi_structures.get("7TRD")
        assert st_tre is not None and st_trd is not None

        ca_tre = extract_ca_atoms(st_tre[0]["A"])
        ca_trd = extract_ca_atoms(st_trd[0]["A"])

        core_res = sorted([r for r in ca_tre if r >= 326 and r in ca_trd])
        assert len(core_res) == 770

        p = np.array([ca_tre[r][1] for r in core_res])
        q = np.array([ca_trd[r][1] for r in core_res])

        rmsd, _, _, _ = kabsch_fit(p, q)
        assert abs(rmsd - 0.338) < 0.05, f"Expected core RMSD ~0.338 A, got {rmsd} A"
        assert rmsd < 0.45, "Catalytic core must remain rigidly locked (<0.45 A)"

    def test_f12_ten_domain_displacement_relative_to_core(self, gemmi_structures):
        """Feature 12: Verify TEN domain displacement on Core-fitted coordinates."""
        st_tre = gemmi_structures.get("7TRE")
        st_trd = gemmi_structures.get("7TRD")
        assert st_tre is not None and st_trd is not None

        ca_tre = extract_ca_atoms(st_tre[0]["A"])
        ca_trd = extract_ca_atoms(st_trd[0]["A"])
        common_res = sorted(set(ca_tre.keys()) & set(ca_trd.keys()))

        # Align on catalytic core (res >= 326)
        core_res = [r for r in common_res if r >= 326]
        p_core = np.array([ca_tre[r][1] for r in core_res])
        q_core = np.array([ca_trd[r][1] for r in core_res])
        _, r_core, t_core, _ = kabsch_fit(p_core, q_core)

        # Evaluate TEN domain (res < 200) displacement when fitted by core
        ten_res = [r for r in common_res if r < 200]
        assert len(ten_res) == 160
        p_ten = np.array([ca_tre[r][1] for r in ten_res])
        q_ten = np.array([ca_trd[r][1] for r in ten_res])

        p_ten_fitted = np.dot(p_ten, r_core.T) + t_core
        ten_disp = float(np.sqrt(np.mean(np.sum((p_ten_fitted - q_ten) ** 2, axis=1))))
        assert abs(ten_disp - 0.434) < 0.05, f"Expected TEN domain displacement ~0.434 A, got {ten_disp} A"

    def test_f14_tpp1_bound_vs_unbound_tel_patch_remodeling(self, gemmi_structures, ref_data):
        """Feature 14: Verify bound TPP1 (7TRE) vs unbound TPP1 (2I46 Chain A) global RMSD ~1.966 A and loop remodeling."""
        st_tre = gemmi_structures.get("7TRE")
        st_2i46 = gemmi_structures.get("2I46")
        assert st_tre is not None and st_2i46 is not None

        ca_tre = extract_ca_atoms(st_tre[0]["C"])
        ca_2i46 = extract_ca_atoms(st_2i46[0]["A"])

        common_res = sorted(set(ca_tre.keys()) & set(ca_2i46.keys()))
        assert len(common_res) == 152, f"Expected 152 common C-alpha residues, got {len(common_res)}"

        p = np.array([ca_tre[r][1] for r in common_res])
        q = np.array([ca_2i46[r][1] for r in common_res])
        rmsd, _r, _t, p_rot = kabsch_fit(p, q)

        assert abs(rmsd - 1.966) < 0.05, f"Expected global TPP1 RMSD ~1.966 A, got {rmsd} A"

        # Check TEL patch per-residue displacements
        res_disp = {
            common_res[i]: float(np.linalg.norm(p_rot[i] - q[i]))
            for i in range(len(common_res))
        }

        # Leu183 is invariant rigid hydrophobic anchor (<0.8 A displacement)
        assert res_disp[183] < 0.8, f"Leu183 displacement {res_disp[183]} A > 0.8 A"
        assert abs(res_disp[183] - 0.512) < 0.1
        # Glu169 and Lys170 undergo extensive remodeling (>3.0 A)
        assert res_disp[169] > 3.0, f"Glu169 displacement {res_disp[169]} A <= 3.0 A"
        assert res_disp[170] > 3.5, f"Lys170 displacement {res_disp[170]} A <= 3.5 A"


@pytest.mark.tier1
class TestTier1PipelineAndOutputArtifacts:
    """Features 15, 16, 17, 18: Standardized CSV tables, report, figures, pipeline runner."""

    def test_f15_expected_tables_list(self, ref_data):
        """Feature 15: Verify all 11 required CSV table names are tracked."""
        expected = ref_data["expected_tables"]
        assert len(expected) == 11
        assert "interface_contacts_7tre_4.0A.csv" in expected
        assert "buried_surface_area.csv" in expected
        assert "cross_state_rmsd.csv" in expected

    def test_f16_synthesis_report_sections(self, project_root):
        """Feature 16: Verify report.md contains required sections if present."""
        report_path = project_root / "report.md"
        if not report_path.exists():
            pytest.skip("report.md not yet generated (M5 in progress)")

        content = report_path.read_text(encoding="utf-8")
        assert "AXIS A" in content or "Axis A" in content, "Report must feature Axis A (Molecular Mechanism)"
        assert "AXIS B" in content or "Axis B" in content, "Report must feature Axis B (Therapeutic Selectivity)"
        assert "7TRE" in content and "7TRD" in content and "2I46" in content

    def test_f18_run_pipeline_syntax_or_api(self, project_root):
        """Feature 18: Verify run_pipeline.py exists or pipeline module is callable."""
        runner_path = project_root / "run_pipeline.py"
        if not runner_path.exists():
            pytest.skip("run_pipeline.py not yet implemented (M5 in progress)")
        assert runner_path.stat().st_size > 0
