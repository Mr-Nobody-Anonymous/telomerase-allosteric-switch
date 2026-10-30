"""
Unit and integration test suite for telotrade interface geometry (Milestone 3).

Covers:
- Contact networks across discrete distance thresholds (3.5 Å, 4.0 Å, 4.5 Å).
- Chain ID disambiguation (author A/C vs label B/D) for interface calculations.
- Buried Surface Area (BSA) via Shrake-Rupley / FreeSASA and mathematical identity.
- Polar interaction profiling: primary salt bridge Arg774-Glu168 (2.32 Å) and secondary
  clamp Lys78-Glu215 (3.58 Å), plus backbone-sidechain hydrogen bonds.
- TEL-patch spatial mapping: hydrophobic anchors Leu769-Leu183 (3.67 Å) and Pro771-Leu212 (3.52 Å).
- Distance cutoff robustness scoring: Core Robust (9), Intermediate (37), Threshold Sensitive (77).
- Standardized CSV table generation matching data schemas in results/tables/.
"""

from pathlib import Path

import pandas as pd
import pytest

from telotrade import (
    HeavyAtomInfo,
    classify_interaction_type,
    compute_bsa,
    compute_contact_robustness,
    compute_interface_contacts,
    compute_polar_interactions,
    generate_interface_tables,
    get_interface_tert_domain,
    map_tel_patch_contacts,
)
from telotrade.interface_geometry import (
    TEL_PATCH_CORE,
    TEL_PATCH_EXTENDED,
    find_project_root,
)


@pytest.fixture(scope="module")
def project_root() -> Path:
    return find_project_root()


@pytest.fixture(scope="module")
def cif_7tre(project_root: Path) -> Path:
    p = project_root / "data" / "raw" / "7TRE.cif"
    if not p.is_file():
        p = project_root / ".agents" / "teamwork" / "survey_pdb" / "cifs" / "7TRE.cif"
    assert p.is_file(), f"7TRE mmCIF file not found at {p}"
    return p


# =============================================================================
# 1. Contact Network Analysis
# =============================================================================


class TestInterfaceContacts:
    """Tests contact network computation across discrete distance cutoffs."""

    def test_classify_interaction_type_unit_cases(self):
        """Direct unit tests for interaction classification logic."""
        # Salt bridge
        assert classify_interaction_type("ARG", "NH1", "N", "GLU", "OE1", "O", 2.32) == "Salt Bridge"
        assert classify_interaction_type("LYS", "NZ", "N", "GLU", "OE2", "O", 3.58) == "Salt Bridge"
        assert classify_interaction_type("ASP", "OD1", "O", "ARG", "NH2", "N", 3.0) == "Salt Bridge"
        # Salt bridge beyond 4.0 A falls back to Polar/H-Bond or VdW
        assert classify_interaction_type("ARG", "NH1", "N", "GLU", "OE1", "O", 4.05) == "Van der Waals"
        # Polar/H-bond
        assert classify_interaction_type("ALA", "N", "N", "GLU", "OE1", "O", 3.45) == "Polar/H-Bond"
        # Hydrophobic
        assert classify_interaction_type("LEU", "CD1", "C", "LEU", "CD2", "C", 3.67) == "Hydrophobic"
        assert classify_interaction_type("PRO", "CG", "C", "LEU", "CD2", "C", 3.52) == "Hydrophobic"
        # Van der Waals
        assert classify_interaction_type("ALA", "CB", "C", "GLU", "OE1", "O", 3.8) == "Van der Waals"

    def test_get_interface_tert_domain_boundaries(self):
        """Domain assignment boundaries for interface residues."""
        assert get_interface_tert_domain(10) == "TEN"
        assert get_interface_tert_domain(78) == "TEN"
        assert get_interface_tert_domain(196) == "TEN"
        assert get_interface_tert_domain(755) == "TRAP_IFD"
        assert get_interface_tert_domain(774) == "TRAP_IFD"
        assert get_interface_tert_domain(800) == "TRAP_IFD"
        assert get_interface_tert_domain(600) == "RT"
        assert get_interface_tert_domain(400) == "TRBD"
        assert get_interface_tert_domain(1000) == "CTE"
        assert get_interface_tert_domain(250) == "LINKER"
        assert get_interface_tert_domain(9999) == "UNKNOWN"

    def test_heavy_atom_info_dataclass(self):
        """Instantiating HeavyAtomInfo preserves fields."""
        import numpy as np

        atom = HeavyAtomInfo(
            auth_seq_id=774,
            label_seq_id=809,
            auth_asym_id="A",
            label_asym_id="B",
            resname="ARG",
            atom_name="NH1",
            element="N",
            pos=np.array([1.0, 2.0, 3.0], dtype=np.float32),
            b_iso=45.0,
        )
        assert atom.auth_seq_id == 774
        assert atom.resname == "ARG"
        assert atom.atom_name == "NH1"

    def test_tel_patch_core_and_extended_constants(self):
        """TEL_PATCH_CORE is a strict subset of TEL_PATCH_EXTENDED."""
        assert set(TEL_PATCH_CORE).issubset(set(TEL_PATCH_EXTENDED))
        assert {168, 169, 170, 171, 183, 212} == set(TEL_PATCH_CORE)

    def test_contacts_7tre_counts_across_cutoffs(self, cif_7tre: Path):
        """Feature 7 / R3: Exactly 9, 46, and 123 atomic contacts at 3.5, 4.0, and 4.5 A."""
        contacts = compute_interface_contacts(cif_7tre, cutoffs=[3.5, 4.0, 4.5])
        assert isinstance(contacts, dict)

        df_35 = contacts[3.5]
        df_40 = contacts[4.0]
        df_45 = contacts[4.5]

        assert len(df_35) == 9, f"Expected 9 contacts at 3.5 A, got {len(df_35)}"
        assert len(df_40) == 46, f"Expected 46 contacts at 4.0 A, got {len(df_40)}"
        assert len(df_45) == 123, f"Expected 123 contacts at 4.5 A, got {len(df_45)}"

    def test_contacts_residue_pair_counts(self, cif_7tre: Path):
        """Feature 7: Exactly 5, 21, and 32 residue-residue pairs at 3.5, 4.0, and 4.5 A."""
        contacts = compute_interface_contacts(cif_7tre, cutoffs=[3.5, 4.0, 4.5])
        assert isinstance(contacts, dict)

        res_pairs_35 = {(r["tert_resnum"], r["tpp1_resnum"]) for _, r in contacts[3.5].iterrows()}
        res_pairs_40 = {(r["tert_resnum"], r["tpp1_resnum"]) for _, r in contacts[4.0].iterrows()}
        res_pairs_45 = {(r["tert_resnum"], r["tpp1_resnum"]) for _, r in contacts[4.5].iterrows()}

        assert len(res_pairs_35) == 5, f"Expected 5 residue pairs at 3.5 A, got {len(res_pairs_35)}"
        assert len(res_pairs_40) == 21, f"Expected 21 residue pairs at 4.0 A, got {len(res_pairs_40)}"
        assert len(res_pairs_45) == 32, f"Expected 32 residue pairs at 4.5 A, got {len(res_pairs_45)}"

    def test_contacts_strict_subset_monotonicity(self, cif_7tre: Path):
        """Feature 7: C_3.5 is a strict subset of C_4.0, which is a strict subset of C_4.5."""
        contacts = compute_interface_contacts(cif_7tre, cutoffs=[3.5, 4.0, 4.5])
        assert isinstance(contacts, dict)

        pairs_35 = {
            (r["tert_resnum"], r["tert_atom"], r["tpp1_resnum"], r["tpp1_atom"])
            for _, r in contacts[3.5].iterrows()
        }
        pairs_40 = {
            (r["tert_resnum"], r["tert_atom"], r["tpp1_resnum"], r["tpp1_atom"])
            for _, r in contacts[4.0].iterrows()
        }
        pairs_45 = {
            (r["tert_resnum"], r["tert_atom"], r["tpp1_resnum"], r["tpp1_atom"])
            for _, r in contacts[4.5].iterrows()
        }

        assert pairs_35.issubset(pairs_40), "Every 3.5 A contact must be in 4.0 A network"
        assert pairs_40.issubset(pairs_45), "Every 4.0 A contact must be in 4.5 A network"
        assert len(pairs_35) < len(pairs_40) < len(pairs_45)

    def test_contacts_chain_disambiguation_author_vs_label(self, cif_7tre: Path):
        """Feature 11: Author chains (A/C) and label chains (B/D) produce identical contacts."""
        contacts_auth = compute_interface_contacts(
            cif_7tre, chain_tert="A", chain_tpp1="C", cutoffs=4.0, use_auth=True
        )
        contacts_label = compute_interface_contacts(
            cif_7tre, chain_tert="B", chain_tpp1="D", cutoffs=4.0, use_auth=False
        )

        assert isinstance(contacts_auth, pd.DataFrame)
        assert isinstance(contacts_label, pd.DataFrame)
        assert len(contacts_auth) == len(contacts_label) == 46

        # Check distances and residue pairs match identically
        assert list(contacts_auth["distance_angstrom"]) == list(contacts_label["distance_angstrom"])
        assert list(contacts_auth["tert_resnum"]) == list(contacts_label["tert_resnum"])
        assert list(contacts_auth["tpp1_resnum"]) == list(contacts_label["tpp1_resnum"])

    def test_contacts_single_cutoff_argument_returns_dataframe(self, cif_7tre: Path):
        """Single numeric cutoff should return a single DataFrame."""
        df_40 = compute_interface_contacts(cif_7tre, cutoffs=4.0)
        assert isinstance(df_40, pd.DataFrame)
        assert len(df_40) == 46

    def test_contacts_zero_and_negative_cutoffs_return_empty(self, cif_7tre: Path):
        """Zero or negative distance thresholds yield empty DataFrames without errors."""
        df_zero = compute_interface_contacts(cif_7tre, cutoffs=0.0)
        assert isinstance(df_zero, pd.DataFrame)
        assert len(df_zero) == 0
        assert "distance_angstrom" in df_zero.columns

        df_neg = compute_interface_contacts(cif_7tre, cutoffs=-1.5)
        assert isinstance(df_neg, pd.DataFrame)
        assert len(df_neg) == 0

    def test_contacts_missing_file_raises_filenotfound(self):
        """Nonexistent mmCIF path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            compute_interface_contacts("nonexistent_path/foo.cif")

    def test_contacts_invalid_chain_raises_keyerror(self, cif_7tre: Path):
        """Invalid chain ID raises KeyError."""
        with pytest.raises(KeyError, match="Chain 'Z' not found"):
            compute_interface_contacts(cif_7tre, chain_tert="Z", chain_tpp1="C")

    def test_contacts_schema_and_types(self, cif_7tre: Path):
        """Verify output table adheres strictly to Section 4.5 data schema."""
        df = compute_interface_contacts(cif_7tre, cutoffs=4.0)
        assert isinstance(df, pd.DataFrame)
        expected_cols = [
            "tert_auth_asym",
            "tert_label_asym",
            "tert_resnum",
            "tert_resname",
            "tert_atom",
            "tert_domain",
            "tpp1_auth_asym",
            "tpp1_label_asym",
            "tpp1_resnum",
            "tpp1_resname",
            "tpp1_atom",
            "is_tel_patch",
            "distance_angstrom",
            "interaction_type",
        ]
        assert list(df.columns) == expected_cols
        assert (df["tert_auth_asym"] == "A").all()
        assert (df["tert_label_asym"] == "B").all()
        assert (df["tpp1_auth_asym"] == "C").all()
        assert (df["tpp1_label_asym"] == "D").all()
        assert (df["distance_angstrom"] <= 4.0).all()

    def test_interaction_types_present(self, cif_7tre: Path):
        """Verify presence of all 4 biophysical interaction categories."""
        df = compute_interface_contacts(cif_7tre, cutoffs=4.5)
        assert isinstance(df, pd.DataFrame)
        types = set(df["interaction_type"])
        assert "Salt Bridge" in types
        assert "Polar/H-Bond" in types
        assert "Hydrophobic" in types
        assert "Van der Waals" in types


# =============================================================================
# 2. Buried Surface Area (BSA)
# =============================================================================


class TestBuriedSurfaceArea:
    """Tests Buried Surface Area (BSA) calculation."""

    def test_compute_bsa_7tre_empirical_benchmark(self, cif_7tre: Path):
        """Feature 8 / R3: BSA = 856.82 +/- 5.0 A^2; Delta SASA = 1713.65 +/- 5.0 A^2."""
        res = compute_bsa(cif_7tre)
        assert isinstance(res, dict)

        bsa_val = res["bsa_interface_half_angstrom2"]
        delta_val = res["delta_sasa_total_angstrom2"]
        sasa_a = res["sasa_partner_a_alone_angstrom2"]
        sasa_b = res["sasa_partner_b_alone_angstrom2"]
        sasa_ab = res["sasa_complex_angstrom2"]

        assert abs(bsa_val - 856.82) <= 5.0, f"Expected BSA ~856.82 A^2, got {bsa_val}"
        assert abs(delta_val - 1713.65) <= 5.0, f"Expected Delta SASA ~1713.65 A^2, got {delta_val}"
        assert abs(sasa_a - 49776.0) <= 20.0, f"TERT SASA {sasa_a} differs from 49776"
        assert abs(sasa_b - 10651.36) <= 20.0, f"TPP1 SASA {sasa_b} differs from 10651.36"
        assert abs(sasa_ab - 58713.7) <= 20.0, f"Complex SASA {sasa_ab} differs from 58713.7"
        assert 800.0 <= bsa_val <= 1200.0, "BSA must fall squarely in physiological PPI range"

    def test_compute_bsa_mathematical_identity(self, cif_7tre: Path):
        """Mathematical identity: BSA = (SASA_A + SASA_B - SASA_AB) / 2."""
        res = compute_bsa(cif_7tre)
        sasa_a = res["sasa_partner_a_alone_angstrom2"]
        sasa_b = res["sasa_partner_b_alone_angstrom2"]
        sasa_ab = res["sasa_complex_angstrom2"]
        delta = sasa_a + sasa_b - sasa_ab
        expected_bsa = delta / 2.0

        assert abs(res["delta_sasa_total_angstrom2"] - delta) < 0.05
        assert abs(res["bsa_interface_half_angstrom2"] - expected_bsa) < 0.05

    def test_compute_bsa_label_chains(self, cif_7tre: Path):
        """BSA calculation succeeds with CIF label chain IDs B and D."""
        res = compute_bsa(cif_7tre, chain_tert="B", chain_tpp1="D", use_auth=False)
        assert abs(res["bsa_interface_half_angstrom2"] - 856.82) <= 5.0

    def test_compute_bsa_missing_file_raises(self):
        """Nonexistent mmCIF raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            compute_bsa("nonexistent_path/7TRE.cif")


# =============================================================================
# 3. Polar Interaction Profiling (Salt Bridges & H-Bonds)
# =============================================================================


class TestPolarInteractions:
    """Tests polar network profiling, salt bridges, and hydrogen bonds."""

    def test_polar_interactions_count_and_schema(self, cif_7tre: Path):
        """Exactly 11 polar atom contacts <= 4.0 A with matching Section 4.8 schema."""
        df = compute_polar_interactions(cif_7tre)
        assert len(df) == 11

        expected_cols = [
            "interaction_type",
            "tert_resnum",
            "tert_resname",
            "tert_atom",
            "tert_domain",
            "tpp1_resnum",
            "tpp1_resname",
            "tpp1_atom",
            "distance_angstrom",
            "is_tel_patch",
        ]
        assert list(df.columns) == expected_cols

    def test_primary_salt_bridge_arg774_glu168(self, cif_7tre: Path):
        """Feature 9: Primary electrostatic anchor TERT Arg774 to TPP1 Glu168 at ~2.32 A."""
        df = compute_polar_interactions(cif_7tre)
        sb1 = df[(df["tert_resnum"] == 774) & (df["tpp1_resnum"] == 168)]
        assert len(sb1) >= 1, "Primary salt bridge Arg774-Glu168 not found"

        min_dist = sb1["distance_angstrom"].min()
        assert abs(min_dist - 2.32) < 0.1, f"Expected salt bridge distance ~2.32 A, got {min_dist}"
        # Validate bidentate interaction atoms
        nh1_oe1 = sb1[(sb1["tert_atom"] == "NH1") & (sb1["tpp1_atom"] == "OE1")]
        assert len(nh1_oe1) == 1
        assert nh1_oe1.iloc[0]["interaction_type"] == "Salt Bridge"

    def test_secondary_clamp_lys78_glu215(self, cif_7tre: Path):
        """Feature 9: Secondary electrostatic clamp TERT Lys78 to TPP1 Glu215 at ~3.58 A."""
        df = compute_polar_interactions(cif_7tre)
        sb2 = df[(df["tert_resnum"] == 78) & (df["tpp1_resnum"] == 215)]
        assert len(sb2) == 1, "Secondary salt bridge Lys78-Glu215 not found"
        row = sb2.iloc[0]
        assert row["tert_atom"] == "NZ"
        assert row["tpp1_atom"] == "OE2"
        assert abs(row["distance_angstrom"] - 3.58) < 0.1
        assert row["interaction_type"] == "Salt Bridge"
        assert row["tert_domain"] == "TEN"
        assert bool(row["is_tel_patch"]) is True

    def test_polar_hbond_backbone_sidechain(self, cif_7tre: Path):
        """Polar interactions include Ala45-Glu171 and Tyr772-Leu211."""
        df = compute_polar_interactions(cif_7tre)
        # Ala45 (TEN backbone N) to Glu171 (OE1)
        ala45 = df[(df["tert_resnum"] == 45) & (df["tpp1_resnum"] == 171)]
        assert len(ala45) >= 1
        assert abs(ala45.iloc[0]["distance_angstrom"] - 3.45) < 0.1
        assert ala45.iloc[0]["interaction_type"] == "H-Bond Backbone-Sidechain"

        # Tyr772 (OH) to Leu211 (O backbone)
        tyr772 = df[(df["tert_resnum"] == 772) & (df["tpp1_resnum"] == 211)]
        assert len(tyr772) >= 1
        assert abs(tyr772.iloc[0]["distance_angstrom"] - 3.36) < 0.1
        assert tyr772.iloc[0]["interaction_type"] == "H-Bond Backbone-Sidechain"


# =============================================================================
# 4. TEL-Patch Epitope Mapping
# =============================================================================


class TestTELPatchMapping:
    """Tests dedicated mapping of TPP1 TEL-patch contacts and hydrophobic anchors."""

    def test_map_tel_patch_contacts_filters_correctly(self, cif_7tre: Path):
        """Extracts contacts for canonical and extended TEL-patch residues."""
        df_45 = compute_interface_contacts(cif_7tre, cutoffs=4.5)
        assert isinstance(df_45, pd.DataFrame)
        df_tel = map_tel_patch_contacts(df_45)

        assert not df_tel.empty
        assert set(df_tel["tpp1_resnum"]).issubset(set(TEL_PATCH_EXTENDED))
        assert "tel_patch_subregion" in df_tel.columns

    def test_hydrophobic_anchors_leu769_and_pro771(self, cif_7tre: Path):
        """Feature 10: Leu769-Leu183 (3.67 A) and Pro771-Leu212 (3.52 A)."""
        df_45 = compute_interface_contacts(cif_7tre, cutoffs=4.5)
        assert isinstance(df_45, pd.DataFrame)
        df_tel = map_tel_patch_contacts(df_45)

        # Anchor 1: Leu769 - Leu183
        h1 = df_tel[(df_tel["tert_resnum"] == 769) & (df_tel["tpp1_resnum"] == 183)]
        assert len(h1) >= 1, "Hydrophobic anchor Leu769-Leu183 missing"
        min_d1 = h1["distance_angstrom"].min()
        assert abs(min_d1 - 3.67) < 0.1, f"Expected Leu769-Leu183 distance ~3.67 A, got {min_d1}"
        assert (h1["interaction_type"] == "Hydrophobic").any()

        # Anchor 2: Pro771 - Leu212
        h2 = df_tel[(df_tel["tert_resnum"] == 771) & (df_tel["tpp1_resnum"] == 212)]
        assert len(h2) >= 1, "Hydrophobic anchor Pro771-Leu212 missing"
        min_d2 = h2["distance_angstrom"].min()
        assert abs(min_d2 - 3.52) < 0.1, f"Expected Pro771-Leu212 distance ~3.52 A, got {min_d2}"
        assert (h2["interaction_type"] == "Hydrophobic").any()


# =============================================================================
# 5. Distance Cutoff Robustness Analysis
# =============================================================================


class TestContactRobustness:
    """Tests classification of contact persistence across distance thresholds."""

    def test_robustness_partition_sum(self, cif_7tre: Path):
        """Feature 13: Core Robust (9) + Intermediate (37) + Threshold Sensitive (77) == 123."""
        df_rob = compute_contact_robustness(cif_7tre)
        assert len(df_rob) == 123

        core = df_rob[df_rob["robustness_class"] == "Core Robust"]
        intermediate = df_rob[df_rob["robustness_class"] == "Intermediate"]
        sensitive = df_rob[df_rob["robustness_class"] == "Threshold Sensitive"]

        assert len(core) == 9, f"Expected 9 Core Robust, got {len(core)}"
        assert len(intermediate) == 37, f"Expected 37 Intermediate, got {len(intermediate)}"
        assert len(sensitive) == 77, f"Expected 77 Threshold Sensitive, got {len(sensitive)}"
        assert len(core) + len(intermediate) + len(sensitive) == 123

    def test_robustness_schema(self, cif_7tre: Path):
        """Verifies Section 4.6 data schema compliance."""
        df_rob = compute_contact_robustness(cif_7tre)
        expected_cols = [
            "contact_pair_id",
            "tert_residue",
            "tert_domain",
            "tpp1_residue",
            "is_tel_patch",
            "min_atom_distance_angstrom",
            "present_at_3.5A",
            "present_at_4.0A",
            "present_at_4.5A",
            "persistence_count",
            "robustness_class",
        ]
        assert list(df_rob.columns) == expected_cols

    def test_primary_salt_bridge_is_core_robust(self, cif_7tre: Path):
        """Arg774-Glu168 salt bridge is Core Robust with persistence score 3."""
        df_rob = compute_contact_robustness(cif_7tre)
        row = df_rob[df_rob["contact_pair_id"].str.contains("ARG774:NH1--TPP1_GLU168:OE1")]
        assert len(row) == 1
        assert row.iloc[0]["persistence_count"] == 3
        assert row.iloc[0]["robustness_class"] == "Core Robust"
        assert bool(row.iloc[0]["present_at_3.5A"]) is True
        assert bool(row.iloc[0]["present_at_4.0A"]) is True
        assert bool(row.iloc[0]["present_at_4.5A"]) is True


# =============================================================================
# 6. Standardized Table Generation
# =============================================================================


class TestStandardizedTableGenerators:
    """Tests file writing and export for standardized results tables."""

    def test_generate_interface_tables_creates_all_files(self, cif_7tre: Path, tmp_path: Path):
        """Verifies all required Milestone 3 CSV tables are generated and non-empty."""
        out_dir = tmp_path / "tables"
        generated = generate_interface_tables(cif_path=cif_7tre, output_dir=out_dir)

        expected_files = [
            "interface_contacts_7tre_3.5A.csv",
            "interface_contacts_3.5A.csv",
            "interface_contacts_7tre_4.0A.csv",
            "interface_contacts_4.0A.csv",
            "interface_contacts_7tre_4.5A.csv",
            "interface_contacts_4.5A.csv",
            "interface_contacts_robustness.csv",
            "buried_surface_area.csv",
            "polar_interactions_7tre.csv",
            "polar_interactions.csv",
            "tel_patch_contacts.csv",
        ]

        for fname in expected_files:
            assert fname in generated, f"Table {fname} not reported in generated dictionary"
            p = generated[fname]
            assert p.is_file(), f"Table file {p} does not exist"
            assert p.stat().st_size > 0, f"Table file {p} is empty"

        # Verify row counts of generated CSVs
        df_35 = pd.read_csv(generated["interface_contacts_7tre_3.5A.csv"])
        df_40 = pd.read_csv(generated["interface_contacts_7tre_4.0A.csv"])
        df_45 = pd.read_csv(generated["interface_contacts_7tre_4.5A.csv"])
        df_rob = pd.read_csv(generated["interface_contacts_robustness.csv"])
        df_bsa = pd.read_csv(generated["buried_surface_area.csv"])
        df_polar = pd.read_csv(generated["polar_interactions_7tre.csv"])
        df_tel = pd.read_csv(generated["tel_patch_contacts.csv"])

        assert len(df_35) == 9
        assert len(df_40) == 46
        assert len(df_45) == 123
        assert len(df_rob) == 123
        assert len(df_bsa) == 1
        assert len(df_polar) == 11
        assert len(df_tel) >= 40


# =============================================================================
# 7. Adversarial & Edge Case Robustness
# =============================================================================


class TestAdversarialInterfaceGeometry:
    """Stress tests and boundary condition validation for interface geometry."""

    def test_unsorted_cutoffs_in_robustness(self, cif_7tre: Path):
        """Unsorted cutoffs are sorted properly and partition persists."""
        df_rob = compute_contact_robustness(cif_7tre, cutoffs=[4.5, 3.5, 4.0])
        assert len(df_rob) == 123
        assert len(df_rob[df_rob["robustness_class"] == "Core Robust"]) == 9
        assert len(df_rob[df_rob["robustness_class"] == "Intermediate"]) == 37
        assert len(df_rob[df_rob["robustness_class"] == "Threshold Sensitive"]) == 77

    def test_non_protein_chain_returns_empty_contacts(self, cif_7tre: Path):
        """Passing RNA chain (author 'B') returns 0 protein-protein contacts without crashing."""
        df = compute_interface_contacts(cif_7tre, chain_tert="B", chain_tpp1="C", cutoffs=4.0)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_custom_probe_radius_bsa_scaling(self, cif_7tre: Path):
        """Smaller probe radius yields larger SASA, larger probe yields smaller SASA."""
        bsa_12 = compute_bsa(cif_7tre, probe_radius=1.2, n_sphere_points=480)
        bsa_16 = compute_bsa(cif_7tre, probe_radius=1.6, n_sphere_points=480)
        assert bsa_12["sasa_tert"] > bsa_16["sasa_tert"]
        assert bsa_12["sasa_tpp1"] > bsa_16["sasa_tpp1"]

    def test_tel_patch_mapping_empty_and_custom_inputs(self):
        """map_tel_patch_contacts handles empty DataFrame or custom residues."""
        empty_df = pd.DataFrame(columns=["tpp1_resnum", "distance_angstrom"])
        res_empty = map_tel_patch_contacts(empty_df)
        assert res_empty.empty

        mock_df = pd.DataFrame(
            [
                {"tpp1_resnum": 168, "distance_angstrom": 2.32},
                {"tpp1_resnum": 999, "distance_angstrom": 3.50},
            ]
        )
        res_mock = map_tel_patch_contacts(mock_df, tel_patch_residues=(168,))
        assert len(res_mock) == 1
        assert res_mock.iloc[0]["tpp1_resnum"] == 168
        assert res_mock.iloc[0]["tel_patch_subregion"] == "Acidic DE-loop"

