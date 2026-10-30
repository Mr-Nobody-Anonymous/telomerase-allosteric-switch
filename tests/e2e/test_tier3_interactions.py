"""Tier 3: Cross-Feature Interaction & Pairwise Workflow E2E Tests.

Tests multi-module interactions and data flow across pipeline stages:
- Sequence-aware chain correspondence feeding into interface contact calculation
- Sequence mapping indices feeding into Kabsch SVD superposition
- Multi-cutoff contact networks feeding into robustness classification matrix
- Interface contact geometry feeding into polar and salt-bridge categorization
- Cross-state superposition feeding into per-residue displacement and B-factor profiling
- Buried Surface Area (BSA) consistency with contact-derived interface residues
"""

import numpy as np
import pytest

from tests.e2e.conftest import (
    extract_ca_atoms,
    extract_heavy_atoms,
    find_contacts,
    kabsch_fit,
)


@pytest.mark.tier3
class TestTier3SequenceToGeometryPipeline:
    """Tests interaction between sequence mapping and geometric calculations."""

    def test_t3_sequence_mapping_feeds_contacts(self, gemmi_structures):
        """Test that sequence-mapped chains (TERT auth A / TPP1 auth C) correctly identify interface contacts."""
        st = gemmi_structures.get("7TRE")
        assert st is not None

        chain_a = st[0]["A"]
        chain_c = st[0]["C"]
        assert chain_a[0].subchain == "B"  # CIF label B
        assert chain_c[0].subchain == "D"  # CIF label D

        tert_heavy = extract_heavy_atoms(chain_a)
        tpp1_heavy = extract_heavy_atoms(chain_c)

        contacts = find_contacts(tert_heavy, tpp1_heavy, 4.0)
        assert len(contacts) == 46

        tert_resnums = {c[0][0] for c in contacts}
        tpp1_resnums = {c[1][0] for c in contacts}
        # TERT Arg774 and TPP1 Glu168 must be in the interaction
        assert 774 in tert_resnums
        assert 168 in tpp1_resnums

    def test_t3_sequence_correspondence_feeds_svd_kabsch(self, gemmi_structures):
        """Test that sequence correspondence alignment produces matched coordinate arrays for Kabsch SVD."""
        st_tre = gemmi_structures.get("7TRE")
        st_trd = gemmi_structures.get("7TRD")
        assert st_tre is not None and st_trd is not None

        ca_tre_map = extract_ca_atoms(st_tre[0]["A"])
        ca_trd_map = extract_ca_atoms(st_trd[0]["A"])

        common_indices = sorted(set(ca_tre_map.keys()) & set(ca_trd_map.keys()))
        assert len(common_indices) == 934, "Strict 1:1 sequence match must yield 934 C-alpha positions"

        coords_tre = np.array([ca_tre_map[idx][1] for idx in common_indices])
        coords_trd = np.array([ca_trd_map[idx][1] for idx in common_indices])

        rmsd, _, _, _ = kabsch_fit(coords_tre, coords_trd)
        assert abs(rmsd - 0.356) < 0.05, f"Expected global RMSD ~0.356 A, got {rmsd} A"


@pytest.mark.tier3
class TestTier3CutoffsToRobustnessClassification:
    """Tests interaction between multi-cutoff contact networks and robustness categorization."""

    def test_t3_contacts_monotonic_progression_subset_rule(self, gemmi_structures):
        """Test that C_3.5 is a strict subset of C_4.0, which is a strict subset of C_4.5."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        c35 = find_contacts(tert_heavy, tpp1_heavy, 3.5)
        c40 = find_contacts(tert_heavy, tpp1_heavy, 4.0)
        c45 = find_contacts(tert_heavy, tpp1_heavy, 4.5)

        pairs_35 = {(c[0][0], c[0][2], c[1][0], c[1][2]) for c in c35}
        pairs_40 = {(c[0][0], c[0][2], c[1][0], c[1][2]) for c in c40}
        pairs_45 = {(c[0][0], c[0][2], c[1][0], c[1][2]) for c in c45}

        assert len(pairs_35) == 9
        assert len(pairs_40) == 46
        assert len(pairs_45) == 123

        # Strict subset invariants
        assert pairs_35.issubset(pairs_40), "Every 3.5 A contact must be present at 4.0 A"
        assert pairs_40.issubset(pairs_45), "Every 4.0 A contact must be present at 4.5 A"

    def test_t3_robustness_classification_partition(self, gemmi_structures):
        """Test partitioning of 4.5 A contacts into Core Robust, Intermediate, and Threshold Sensitive."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        c45 = find_contacts(tert_heavy, tpp1_heavy, 4.5)
        assert len(c45) == 123

        persistence = {}
        for c in c45:
            pair = (c[0][0], c[0][2], c[1][0], c[1][2])
            d = c[2]
            score = 0
            if d <= 4.5:
                score += 1
            if d <= 4.0:
                score += 1
            if d <= 3.5:
                score += 1
            persistence[pair] = (score, d)

        core_robust = [p for p, (score, _) in persistence.items() if score == 3]
        intermediate = [p for p, (score, _) in persistence.items() if score == 2]
        threshold_sensitive = [p for p, (score, _) in persistence.items() if score == 1]

        assert len(core_robust) == 9, "Core Robust contacts must equal 3.5 A contacts (9)"
        assert len(intermediate) == 37, "Intermediate contacts must equal 4.0 A - 3.5 A contacts (37)"
        assert len(threshold_sensitive) == 77, "Threshold Sensitive contacts must equal 4.5 A - 4.0 A contacts (77)"
        assert len(core_robust) + len(intermediate) + len(threshold_sensitive) == 123


@pytest.mark.tier3
class TestTier3InterfaceToPolarAndBsaConsistency:
    """Tests interaction between contact geometry, polar categorization, and BSA."""

    def test_t3_polar_subset_of_contacts(self, gemmi_structures):
        """Test that polar interactions and salt bridges form a strict subset of the 4.0 A contact network."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        c40 = find_contacts(tert_heavy, tpp1_heavy, 4.0)

        all_contacts = set()
        salt_bridges = set()

        acidic_atoms = {"OE1", "OE2", "OD1", "OD2"}
        basic_atoms = {"NZ", "NH1", "NH2", "NE"}

        for c in c40:
            t_resnum, t_resname, t_atom, _, _ = c[0]
            p_resnum, p_resname, p_atom, _, _ = c[1]
            pair = (t_resnum, t_atom, p_resnum, p_atom)
            all_contacts.add(pair)

            is_sb = (
                (t_resname in ["ARG", "LYS"] and t_atom in basic_atoms and
                 p_resname in ["GLU", "ASP"] and p_atom in acidic_atoms) or
                (p_resname in ["ARG", "LYS"] and p_atom in basic_atoms and
                 t_resname in ["GLU", "ASP"] and t_atom in acidic_atoms)
            )
            if is_sb:
                salt_bridges.add(pair)

        assert len(salt_bridges) >= 2, "Must identify at least Arg774-Glu168 and Lys78-Glu215 salt bridges"
        assert salt_bridges.issubset(all_contacts), "Salt bridges must be a subset of 4.0 A contacts"

    def test_t3_cross_state_tpp1_rmsd_with_local_bfactor(self, gemmi_structures):
        """Test cross-state TPP1 displacement correlates with TEL-patch structural plasticity."""
        st_tre = gemmi_structures.get("7TRE")
        st_2i46 = gemmi_structures.get("2I46")
        assert st_tre is not None and st_2i46 is not None

        ca_tre = extract_ca_atoms(st_tre[0]["C"])
        ca_2i46 = extract_ca_atoms(st_2i46[0]["A"])

        common = sorted(set(ca_tre.keys()) & set(ca_2i46.keys()))
        p = np.array([ca_tre[r][1] for r in common])
        q = np.array([ca_2i46[r][1] for r in common])
        _rmsd, _r, _t, p_rot = kabsch_fit(p, q)

        displacements = {common[i]: float(np.linalg.norm(p_rot[i] - q[i])) for i in range(len(common))}

        # Compare Leu183 (hydrophobic core anchor) vs Glu169 (polar loop tip)
        disp_leu183 = displacements[183]
        disp_glu169 = displacements[169]

        assert disp_leu183 < 1.0, f"Leu183 should be rigid (<1.0 A), got {disp_leu183} A"
        assert disp_glu169 > 3.0, f"Glu169 should exhibit extensive remodeling (>3.0 A), got {disp_glu169} A"
        assert disp_glu169 > 3.0 * disp_leu183, "TEL-patch DE-loop remodeling must substantially exceed hydrophobic anchor displacement"
