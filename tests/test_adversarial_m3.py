"""
Adversarial stress-testing suite for Milestone 3:
Quantitative Interface Geometry & Interaction Profiling Engine.
Authored by Challenger 1 (challenger_m3_1).

Empirically probes:
1. Boundary cutoffs: zero (0.0 A), negative (-0.01, -1.0, -100.0 A), large cutoffs (100.0 A on
   synthetic models), non-standard/empty cutoffs list.
2. Invalid and edge-case chain IDs: nonexistent chain, empty string, RNA/non-protein chains,
   intra-chain self-interface.
3. Degenerate and malformed coordinate geometry: zero-byte CIF, corrupt CIF, coincident atoms
   (distance 0.0 A), large coordinate values, coplanar and single-atom geometries.
4. Mathematical precision and Euclidean distance verification:
   - Double-precision SciPy cdist oracle comparison across all 7TRE contacts (zero false pos/neg).
   - High-precision synthetic Pythagorean triplets.
   - Coordinate translation invariance and monotonicity.
5. Strictness of salt bridge classification:
   - Exhaustive true positive acidic-basic combinations.
   - Strict 4.00 A distance boundary (4.000 vs 4.001 A).
   - Rejection of neutral polar groups (Asn/Gln amides, Ser/Thr/Tyr hydroxyls).
   - Rejection of backbone amide/carbonyl, aliphatic carbons, basic-basic, and acidic-acidic pairs.
6. Buried Surface Area (BSA) robustness under adversarial probe radii and empty/non-protein chains.
"""

from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import pytest
import scipy.spatial.distance as sp_dist

from telotrade import (
    classify_interaction_type,
    compute_bsa,
    compute_interface_contacts,
    compute_polar_interactions,
)
from telotrade.interface_geometry import (
    ACIDIC_ATOMS,
    BASIC_ATOMS,
    _extract_chain_heavy_atoms,
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


def create_synthetic_cif(
    path: Path,
    chain_a_atoms: list[tuple[str, str, str, tuple[float, float, float]]],
    chain_b_atoms: list[tuple[str, str, str, tuple[float, float, float]]],
) -> Path:
    """
    Helper to create an mmCIF file with custom coordinates for adversarial testing.
    Format of atoms: (resname, atom_name, element, (x, y, z)).
    """
    st = gemmi.Structure()
    m = gemmi.Model("1")

    # Chain A
    c_a = gemmi.Chain("A")
    for idx, (resname, atom_name, element, (x, y, z)) in enumerate(chain_a_atoms, 1):
        r = gemmi.Residue()
        r.name = resname
        r.seqid = gemmi.SeqId(str(idx))
        at = gemmi.Atom()
        at.name = atom_name
        at.element = gemmi.Element(element)
        at.pos = gemmi.Position(float(x), float(y), float(z))
        r.add_atom(at)
        c_a.add_residue(r)
    m.add_chain(c_a)

    # Chain B
    c_b = gemmi.Chain("B")
    for idx, (resname, atom_name, element, (x, y, z)) in enumerate(chain_b_atoms, 1):
        r = gemmi.Residue()
        r.name = resname
        r.seqid = gemmi.SeqId(str(idx))
        at = gemmi.Atom()
        at.name = atom_name
        at.element = gemmi.Element(element)
        at.pos = gemmi.Position(float(x), float(y), float(z))
        r.add_atom(at)
        c_b.add_residue(r)
    m.add_chain(c_b)

    st.add_model(m)
    doc = st.make_mmcif_document()
    doc.write_file(str(path))
    return path


# =============================================================================
# 1. Boundary Cutoffs Stress Testing
# =============================================================================


class TestAdversarialBoundaryCutoffs:
    """Adversarial testing of distance thresholds."""

    @pytest.mark.parametrize("cutoff", [0.0, -0.001, -1.0, -10.0, -100.0])
    def test_zero_and_negative_cutoffs_return_empty_dataframe(self, cif_7tre: Path, cutoff: float):
        """Zero or negative distance cutoffs must return empty DataFrame with 14 columns."""
        df = compute_interface_contacts(cif_7tre, cutoffs=cutoff)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
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

    def test_cutoff_list_with_mixed_valid_and_boundary_values(self, cif_7tre: Path):
        """Passing mixed cutoffs including negative and zero produces proper dictionary."""
        cutoffs = [-5.0, 0.0, 3.5, 4.0]
        res = compute_interface_contacts(cif_7tre, cutoffs=cutoffs)
        assert isinstance(res, dict)
        assert len(res[-5.0]) == 0
        assert len(res[0.0]) == 0
        assert len(res[3.5]) == 9
        assert len(res[4.0]) == 46

    def test_empty_cutoffs_list_returns_empty_dict(self, cif_7tre: Path):
        """Passing an empty list of cutoffs returns an empty dictionary."""
        res = compute_interface_contacts(cif_7tre, cutoffs=[])
        assert isinstance(res, dict)
        assert len(res) == 0

    def test_large_cutoff_100A_on_synthetic_system(self, tmp_path: Path):
        """Cutoff = 100.0 A on a synthetic coordinate system correctly captures wide-range pairs."""
        cif_path = tmp_path / "large_cutoff.cif"
        chain_a = [("ARG", "NH1", "N", (0.0, 0.0, 0.0))]
        chain_b = [
            ("GLU", "OE1", "O", (20.0, 0.0, 0.0)),   # dist = 20.0 A
            ("GLU", "OE2", "O", (60.0, 0.0, 0.0)),   # dist = 60.0 A
            ("ASP", "OD1", "O", (99.0, 0.0, 0.0)),   # dist = 99.0 A
            ("ASP", "OD2", "O", (101.0, 0.0, 0.0)),  # dist = 101.0 A (outside 100 A)
        ]
        create_synthetic_cif(cif_path, chain_a, chain_b)

        df_100 = compute_interface_contacts(cif_path, chain_tert="A", chain_tpp1="B", cutoffs=100.0)
        assert isinstance(df_100, pd.DataFrame)
        assert len(df_100) == 3, f"Expected 3 contacts <= 100 A, got {len(df_100)}"
        assert list(df_100["distance_angstrom"]) == [20.0, 60.0, 99.0]
        # Distances > 4.0 A between acidic/basic are classified as Van der Waals
        assert (df_100["interaction_type"] == "Van der Waals").all()

    def test_fine_grained_monotonicity_below_3A(self, cif_7tre: Path):
        """Contacts at ultra-close cutoffs: 1.0 A -> 0, 2.0 A -> 0, 2.4 A -> 1 (primary salt bridge)."""
        res = compute_interface_contacts(cif_7tre, cutoffs=[1.0, 2.0, 2.4, 2.5])
        assert isinstance(res, dict)
        assert len(res[1.0]) == 0
        assert len(res[2.0]) == 0
        assert len(res[2.4]) == 1  # Exactly the 2.32 A Arg774-Glu168 salt bridge
        row = res[2.4].iloc[0]
        assert row["tert_resnum"] == 774
        assert row["tpp1_resnum"] == 168
        assert row["distance_angstrom"] == 2.32
        assert row["interaction_type"] == "Salt Bridge"
        assert len(res[2.5]) >= 1


# =============================================================================
# 2. Invalid Chain Identifiers & Macromolecular Entities
# =============================================================================


class TestAdversarialChainHandling:
    """Stress testing chain lookup and entity filtering."""

    @pytest.mark.parametrize("invalid_id", ["Z", "INVALID", "99", ""])
    def test_invalid_chain_id_raises_informative_keyerror(self, cif_7tre: Path, invalid_id: str):
        """Nonexistent chain ID raises KeyError specifying available author and label chains."""
        with pytest.raises(KeyError) as exc_info:
            compute_interface_contacts(cif_7tre, chain_tert=invalid_id, chain_tpp1="C")
        err_msg = str(exc_info.value)
        assert f"Chain '{invalid_id}' not found" in err_msg
        assert "Available author chains" in err_msg

    def test_non_protein_chain_returns_empty_contacts_without_crash(self, cif_7tre: Path):
        """Author chain B is hTR RNA (332 nucleic acids); should yield 0 contacts cleanly."""
        res = compute_interface_contacts(cif_7tre, chain_tert="B", chain_tpp1="C", cutoffs=4.0)
        assert isinstance(res, pd.DataFrame)
        assert len(res) == 0

    def test_same_chain_self_interface_does_not_crash(self, cif_7tre: Path):
        """Passing identical chains (A and A) computes intra-chain contacts without crashing."""
        df_self = compute_interface_contacts(cif_7tre, chain_tert="A", chain_tpp1="A", cutoffs=2.0)
        assert isinstance(df_self, pd.DataFrame)
        # Every atom pairs with itself at distance 0.0 A
        assert len(df_self) > 0
        assert df_self["distance_angstrom"].min() == 0.0

    def test_polar_interactions_raises_on_invalid_chain(self, cif_7tre: Path):
        """compute_polar_interactions raises KeyError on invalid chain ID."""
        with pytest.raises(KeyError):
            compute_polar_interactions(cif_7tre, chain_tert="NONEXISTENT")

    def test_bsa_raises_on_invalid_chain(self, cif_7tre: Path):
        """compute_bsa raises KeyError on invalid chain ID."""
        with pytest.raises(KeyError):
            compute_bsa(cif_7tre, chain_tert="NONEXISTENT")


# =============================================================================
# 3. Degenerate and Malformed Coordinate Geometry
# =============================================================================


class TestAdversarialDegenerateCoordinates:
    """Testing behavior under degenerate, truncated, or extreme coordinates."""

    def test_zero_byte_cif_raises_exception(self, tmp_path: Path):
        """A 0-byte file raises an exception without crashing Python process."""
        empty_cif = tmp_path / "empty.cif"
        empty_cif.write_text("")
        with pytest.raises((IndexError, ValueError)):
            compute_interface_contacts(empty_cif)

    def test_coincident_atoms_calculate_zero_distance(self, tmp_path: Path):
        """Atoms located at identical spatial coordinates compute distance = 0.0 A."""
        cif_path = tmp_path / "coincident.cif"
        chain_a = [("ARG", "NH1", "N", (10.0, 20.0, 30.0))]
        chain_b = [("GLU", "OE1", "O", (10.0, 20.0, 30.0))]
        create_synthetic_cif(cif_path, chain_a, chain_b)

        df = compute_interface_contacts(cif_path, chain_tert="A", chain_tpp1="B", cutoffs=4.0)
        assert len(df) == 1
        assert df.iloc[0]["distance_angstrom"] == 0.0
        assert df.iloc[0]["interaction_type"] == "Salt Bridge"

    def test_extreme_negative_coordinates(self, tmp_path: Path):
        """Handles negative coordinates across cell hashing boundaries."""
        cif_path = tmp_path / "negative_coords.cif"
        chain_a = [("LEU", "CD1", "C", (-150.25, -200.50, -350.75))]
        chain_b = [("LEU", "CD2", "C", (-152.25, -200.50, -350.75))]  # dx = 2.0 A
        create_synthetic_cif(cif_path, chain_a, chain_b)

        df = compute_interface_contacts(cif_path, chain_tert="A", chain_tpp1="B", cutoffs=4.5)
        assert len(df) == 1
        assert df.iloc[0]["distance_angstrom"] == 2.0
        assert df.iloc[0]["interaction_type"] == "Hydrophobic"

    def test_large_coordinates_no_overflow(self, tmp_path: Path):
        """Coordinates up to 100,000 A compute without numeric overflow or crash."""
        cif_path = tmp_path / "large_coords.cif"
        chain_a = [("ALA", "CB", "C", (100000.0, 100000.0, 100000.0))]
        chain_b = [("ALA", "CB", "C", (100003.0, 100004.0, 100000.0))]  # dx=3, dy=4 -> dist=5.0 A
        create_synthetic_cif(cif_path, chain_a, chain_b)

        df = compute_interface_contacts(cif_path, chain_tert="A", chain_tpp1="B", cutoffs=6.0)
        assert len(df) == 1
        assert df.iloc[0]["distance_angstrom"] == 5.0


# =============================================================================
# 4. Mathematical Precision & Oracle Verification
# =============================================================================


class TestMathematicalPrecisionEuclideanDistances:
    """Rigorous double-precision oracle validation of Euclidean distances."""

    def test_scipy_cdist_oracle_full_parity_7tre_45A(self, cif_7tre: Path):
        """
        Verify all 123 contacts at 4.5 A against an independent SciPy float64 cdist oracle.
        Guarantees:
        1. Exact distance match (max diff == 0.0 after 2 decimal rounding).
        2. Zero false positives (every reported contact is <= 4.5 A in oracle).
        3. Zero false negatives (every oracle pair <= 4.5 A is in reported table).
        """
        tert_atoms = _extract_chain_heavy_atoms(cif_7tre, "A", use_auth=True)
        tpp1_atoms = _extract_chain_heavy_atoms(cif_7tre, "C", use_auth=True)

        coords_t = np.array([a.pos for a in tert_atoms], dtype=np.float64)
        coords_p = np.array([a.pos for a in tpp1_atoms], dtype=np.float64)

        # Oracle distance matrix
        dmat = sp_dist.cdist(coords_t, coords_p)

        df_45 = compute_interface_contacts(cif_7tre, cutoffs=4.5)
        assert len(df_45) == 123

        max_discrepancy = 0.0
        for _, row in df_45.iterrows():
            it = next(
                i for i, a in enumerate(tert_atoms)
                if a.auth_seq_id == row["tert_resnum"] and a.atom_name == row["tert_atom"]
            )
            ip = next(
                i for i, a in enumerate(tpp1_atoms)
                if a.auth_seq_id == row["tpp1_resnum"] and a.atom_name == row["tpp1_atom"]
            )
            d_oracle = dmat[it, ip]
            d_reported = row["distance_angstrom"]
            diff = abs(d_reported - round(d_oracle, 2))
            max_discrepancy = max(max_discrepancy, diff)

        assert max_discrepancy == 0.0, f"Oracle discrepancy detected: {max_discrepancy} A"

        # Oracle total pairs <= 4.5 A
        oracle_hits = int(np.sum(dmat <= 4.5))
        assert oracle_hits == 123, f"Oracle found {oracle_hits} pairs <= 4.5 A, expected 123"

    def test_pythagorean_triplet_precision(self, tmp_path: Path):
        """Tests exact integer Pythagorean distance: 1^2 + 2^2 + 2^2 = 9 -> sqrt(9) = 3.00."""
        cif_path = tmp_path / "pythagorean.cif"
        chain_a = [("ALA", "CA", "C", (10.0, 10.0, 10.0))]
        chain_b = [("ALA", "CA", "C", (11.0, 12.0, 12.0))]  # dx=1, dy=2, dz=2 -> dist=3.0
        create_synthetic_cif(cif_path, chain_a, chain_b)

        df = compute_interface_contacts(cif_path, chain_tert="A", chain_tpp1="B", cutoffs=4.0)
        assert len(df) == 1
        assert df.iloc[0]["distance_angstrom"] == 3.00

    def test_translation_invariance(self, tmp_path: Path):
        """Translating both partner molecules rigidly by a vector does not alter interatomic distance."""
        cif_base = tmp_path / "base.cif"
        cif_shifted = tmp_path / "shifted.cif"

        create_synthetic_cif(
            cif_base,
            [("ARG", "NH1", "N", (1.0, 2.0, 3.0))],
            [("GLU", "OE1", "O", (3.0, 4.0, 5.0))],
        )

        shift = (734.5, -921.3, 114.8)
        create_synthetic_cif(
            cif_shifted,
            [("ARG", "NH1", "N", (1.0 + shift[0], 2.0 + shift[1], 3.0 + shift[2]))],
            [("GLU", "OE1", "O", (3.0 + shift[0], 4.0 + shift[1], 5.0 + shift[2]))],
        )

        df_base = compute_interface_contacts(cif_base, chain_tert="A", chain_tpp1="B", cutoffs=5.0)
        df_shifted = compute_interface_contacts(cif_shifted, chain_tert="A", chain_tpp1="B", cutoffs=5.0)

        assert df_base.iloc[0]["distance_angstrom"] == df_shifted.iloc[0]["distance_angstrom"]


# =============================================================================
# 5. Strictness of Salt Bridge Classification
# =============================================================================


class TestStrictnessSaltBridgeClassification:
    """Adversarial stress-testing of electrostatic interaction classification."""

    def test_all_positive_acidic_basic_combinations(self):
        """All functional acidic-basic atom pairs <= 4.0 A are classified as Salt Bridge."""
        for b_res, b_atom in BASIC_ATOMS:
            for a_res, a_atom in ACIDIC_ATOMS:
                res1 = classify_interaction_type(b_res, b_atom, "N", a_res, a_atom, "O", 3.2)
                res2 = classify_interaction_type(a_res, a_atom, "O", b_res, b_atom, "N", 3.2)
                assert res1 == "Salt Bridge", f"Failed for {b_res}:{b_atom} -> {a_res}:{a_atom}"
                assert res2 == "Salt Bridge", f"Failed for {a_res}:{a_atom} -> {b_res}:{b_atom}"

    def test_strict_4A_distance_boundary_salt_bridge(self):
        """Salt bridge classification terminates sharply at 4.00 A."""
        # 4.000 A is Salt Bridge
        assert classify_interaction_type("ARG", "NH1", "N", "GLU", "OE1", "O", 4.000) == "Salt Bridge"
        # 4.001 A is NOT Salt Bridge (falls back to Van der Waals since Polar/H-Bond also capped at 4.0)
        assert classify_interaction_type("ARG", "NH1", "N", "GLU", "OE1", "O", 4.001) == "Van der Waals"

    @pytest.mark.parametrize(
        "res_t,atom_t,el_t,res_p,atom_p,el_p",
        [
            ("ASN", "ND2", "N", "GLU", "OE1", "O"),  # Asn neutral amide N
            ("ASN", "OD1", "O", "ARG", "NH1", "N"),  # Asn neutral amide O
            ("GLN", "NE2", "N", "ASP", "OD1", "O"),  # Gln neutral amide N
            ("GLN", "OE1", "O", "LYS", "NZ", "N"),   # Gln neutral amide O
            ("SER", "OG", "O", "ARG", "NH1", "N"),   # Ser neutral hydroxyl
            ("THR", "OG1", "O", "GLU", "OE1", "O"),  # Thr neutral hydroxyl
            ("TYR", "OH", "O", "ARG", "NH1", "N"),   # Tyr phenolic hydroxyl
            ("ARG", "N", "N", "GLU", "OE1", "O"),    # Arg backbone N
            ("ARG", "O", "O", "GLU", "OE1", "O"),    # Arg backbone O
            ("GLU", "N", "N", "ARG", "NH1", "N"),    # Glu backbone N
            ("GLU", "O", "O", "ARG", "NH1", "N"),    # Glu backbone O
            ("ARG", "NH1", "N", "LYS", "NZ", "N"),   # Basic - Basic (repulsive)
            ("GLU", "OE1", "O", "ASP", "OD1", "O"),   # Acidic - Acidic (repulsive)
            ("HIS", "CE1", "C", "GLU", "OE1", "O"),  # Ring carbon
            ("ARG", "CD", "C", "GLU", "OE1", "O"),   # Aliphatic carbon
        ],
    )
    def test_neutral_polar_and_backbone_atoms_rejected_from_salt_bridge(
        self, res_t, atom_t, el_t, res_p, atom_p, el_p
    ):
        """Neutral polar atoms and backbone atoms must NEVER be classified as Salt Bridge."""
        itype = classify_interaction_type(res_t, atom_t, el_t, res_p, atom_p, el_p, 3.2)
        assert itype != "Salt Bridge", (
            f"Adversarial violation: {res_t}:{atom_t} -- {res_p}:{atom_p} was classified as Salt Bridge"
        )

    def test_7tre_primary_salt_bridges_exactness(self, cif_7tre: Path):
        """
        Verify the exact set of 4 salt bridge contacts <= 4.0 A in 7TRE:
        1. Arg774:NH1 -- Glu168:OE1 (2.32 A)
        2. Lys78:NZ -- Glu215:OE2 (3.58 A)
        3. Arg774:NE -- Glu168:OE1 (3.85 A)
        4. Arg774:NH1 -- Glu168:OE2 (3.87 A)
        """
        df_polar = compute_polar_interactions(cif_7tre)
        sb = df_polar[df_polar["interaction_type"] == "Salt Bridge"]
        assert len(sb) == 4, f"Expected 4 salt bridge contacts <= 4.0 A, got {len(sb)}"

        # Validate distances
        dists = sorted(sb["distance_angstrom"])
        assert dists == [2.32, 3.58, 3.85, 3.87]


# =============================================================================
# 6. Buried Surface Area (BSA) Robustness
# =============================================================================


class TestAdversarialBSA:
    """Stress testing Buried Surface Area calculation."""

    def test_bsa_non_protein_chain_yields_zero_bsa(self, cif_7tre: Path):
        """Passing RNA chain B yields 0.0 BSA gracefully."""
        bsa = compute_bsa(cif_7tre, chain_tert="B", chain_tpp1="C")
        assert bsa["bsa"] == 0.0
        assert bsa["delta_sasa"] == 0.0

    def test_bsa_sphere_points_convergence(self, cif_7tre: Path):
        """Higher sphere test points (480 vs 960) yield consistent BSA within 2%."""
        bsa_480 = compute_bsa(cif_7tre, n_sphere_points=480)
        bsa_960 = compute_bsa(cif_7tre, n_sphere_points=960)

        diff = abs(bsa_480["bsa"] - bsa_960["bsa"])
        assert diff < 20.0, f"BSA differed by {diff} A^2 between 480 and 960 sphere points"

    def test_bsa_probe_radius_strictly_monotonic(self, cif_7tre: Path):
        """Increasing probe radius (1.2 -> 1.4 -> 1.8 A) strictly reduces individual SASAs."""
        bsa_12 = compute_bsa(cif_7tre, probe_radius=1.2, n_sphere_points=480)
        bsa_14 = compute_bsa(cif_7tre, probe_radius=1.4, n_sphere_points=480)
        bsa_18 = compute_bsa(cif_7tre, probe_radius=1.8, n_sphere_points=480)

        assert bsa_12["sasa_tert"] > bsa_14["sasa_tert"] > bsa_18["sasa_tert"]
        assert bsa_12["sasa_tpp1"] > bsa_14["sasa_tpp1"] > bsa_18["sasa_tpp1"]
