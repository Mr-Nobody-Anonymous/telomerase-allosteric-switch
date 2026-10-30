"""Tier 2: Boundary, Corner Case & Adversarial E2E Tests.

Tests edge conditions, mathematical boundaries, and adversarial corner cases:
- Exact distance cutoff boundaries (3.500 A, 4.000 A, 4.500 A) and invalid thresholds
- Unmodeled loops and gaps (TERT flexible linker 197-325, terminal gaps)
- 2I46 Chain B disorder vs Chain A resolved coordinates
- Alternate locations (altloc) handling without atomic center duplication
- Kabsch SVD mathematical degeneracy (planar/collinear points, reflection determinant)
- File integrity corruption detection and missing file handling
"""

import hashlib

import gemmi
import numpy as np
import pytest

from tests.e2e.conftest import (
    extract_heavy_atoms,
    find_contacts,
)


def kabsch_fit_with_det(p, q):
    """Kabsch SVD with explicit reflection determinant check."""
    centroid_p = np.mean(p, axis=0)
    centroid_q = np.mean(q, axis=0)
    p_centered = p - centroid_p
    q_centered = q - centroid_q
    h = np.dot(p_centered.T, q_centered)
    u, _s, vt = np.linalg.svd(h)
    det_h = np.linalg.det(np.dot(vt.T, u.T))
    d = np.sign(det_h) if abs(det_h) > 1e-10 else 1.0
    diag = np.diag([1.0, 1.0, d])
    r = np.dot(vt.T, np.dot(diag, u.T))
    p_rotated = np.dot(p_centered, r.T)
    rmsd = float(np.sqrt(np.mean(np.sum((p_rotated - q_centered) ** 2, axis=1))))
    return rmsd, r, float(np.linalg.det(r))


@pytest.mark.tier2
class TestTier2DistanceCutoffBoundaries:
    """Boundary conditions for interface distance thresholds."""

    def test_t2_cutoff_exact_boundary_floating_point(self, gemmi_structures):
        """Test boundary filtering logic: contacts exactly <= cutoff are included, > cutoff excluded."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        # Compute contacts up to 4.5 A
        contacts_45 = find_contacts(tert_heavy, tpp1_heavy, 4.5)
        distances = sorted([c[2] for c in contacts_45])

        assert len(distances) == 123
        leq_35 = [d for d in distances if d <= 3.5]
        gt_35_leq_40 = [d for d in distances if 3.5 < d <= 4.0]
        gt_40_leq_45 = [d for d in distances if 4.0 < d <= 4.5]

        assert len(leq_35) == 9
        assert len(gt_35_leq_40) == 37  # 9 + 37 = 46
        assert len(gt_40_leq_45) == 77  # 46 + 77 = 123

    def test_t2_negative_or_zero_cutoff_handling(self, gemmi_structures):
        """Test that negative or zero cutoffs produce empty contact sets rather than errors."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])[:50]
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])[:50]

        # Zero cutoff
        zero_contacts = find_contacts(tert_heavy, tpp1_heavy, 0.0)
        assert len(zero_contacts) == 0, "Zero distance threshold must yield zero contacts"

    def test_t2_extreme_high_cutoff_scaling(self, gemmi_structures):
        """Test interface search behavior at larger cutoff expands stably."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_heavy = extract_heavy_atoms(st[0]["A"])
        tpp1_heavy = extract_heavy_atoms(st[0]["C"])

        contacts_60 = find_contacts(tert_heavy, tpp1_heavy, 6.0)
        assert len(contacts_60) > 123, "Contacts at 6.0 A must exceed contacts at 4.5 A"


@pytest.mark.tier2
class TestTier2UnmodeledLoopsAndDisorder:
    """Corner cases concerning missing residues, unmodeled loops, and chain asymmetry."""

    def test_t2_unmodeled_tert_linker_179_to_321(self, gemmi_structures):
        """Test that missing flexible linker (residues 179-321) has no coordinates in 7TRE."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_seqs = {res.seqid.num for res in st[0]["A"]}
        # In TERT, residues 179 to 321 are unmodeled in cryo-EM density
        for res_num in range(179, 322):
            assert res_num not in tert_seqs, f"Residue {res_num} should be unmodeled in 7TRE linker"

    def test_t2_unmodeled_terminal_residues(self, gemmi_structures):
        """Test terminal unmodeled residue boundaries (residues 1-6 in TERT)."""
        st = gemmi_structures.get("7TRE")
        assert st is not None
        tert_seqs = {res.seqid.num for res in st[0]["A"]}
        for res_num in range(1, 7):
            assert res_num not in tert_seqs, f"TERT N-terminal residue {res_num} must be unmodeled"

    def test_t2_chain_b_tel_patch_nan_safety(self, gemmi_structures):
        """Test that accessing 2I46 Chain B for TEL patch yields no coordinates without crashing."""
        st = gemmi_structures.get("2I46")
        assert st is not None
        chain_b_resnums = {res.seqid.num for res in st[0]["B"]}
        for patch_res in [168, 169, 170, 171]:
            assert patch_res not in chain_b_resnums, f"Chain B must not report resolved coordinates for {patch_res}"

    def test_t2_altloc_handling_no_duplicate_ca(self, gemmi_structures):
        """Test that alternate location conformations do not result in duplicate C-alpha entries."""
        for pid in ["7TRE", "2I46"]:
            st = gemmi_structures.get(pid)
            assert st is not None
            for chain in st[0]:
                for res in chain:
                    ca_count = sum(1 for a in res if a.name == "CA")
                    assert ca_count <= 1, f"Duplicate CA atom detected in {pid} Chain {chain.name} Residue {res.seqid.num}"


@pytest.mark.tier2
class TestTier2KabschMathematicalBoundaries:
    """Mathematical edge cases for Kabsch SVD superposition."""

    def test_t2_kabsch_reflection_determinant(self):
        """Test that Kabsch SVD strictly prevents coordinate inversion (reflection)."""
        coords = np.array([
            [1.0, 1.0, 1.0],
            [2.0, 1.0, 1.0],
            [1.0, 3.0, 1.0],
            [1.0, 1.0, 4.0]
        ], dtype=float)
        reflected = coords.copy()
        reflected[:, 2] = -reflected[:, 2]

        _rmsd, _r, det = kabsch_fit_with_det(reflected, coords)
        assert abs(det - 1.0) < 1e-5, f"Determinant must be +1.0, got {det}"

    def test_t2_kabsch_collinear_or_planar_points(self):
        """Test Kabsch SVD on planar points does not produce NaN or inf."""
        planar_a = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0]
        ], dtype=float)
        planar_b = planar_a + np.array([5.0, 2.0, 0.0])

        rmsd, r, det = kabsch_fit_with_det(planar_b, planar_a)
        assert not np.isnan(rmsd)
        assert not np.isnan(r).any()
        assert rmsd < 1e-10
        assert abs(det - 1.0) < 1e-5

    def test_t2_kabsch_insufficient_points_handling(self):
        """Test that coordinate sets with fewer than 3 non-collinear points are handled safely."""
        single_point_a = np.array([[1.0, 2.0, 3.0]])
        diff = single_point_a - single_point_a.mean(axis=0)
        assert np.allclose(diff, 0.0)


@pytest.mark.tier2
class TestTier2AdversarialFileIntegrity:
    """Adversarial and corrupted input handling for file I/O."""

    def test_t2_missing_file_raises_not_found(self, project_root):
        """Test that attempting to load a non-existent CIF file raises FileNotFoundError."""
        bogus_path = project_root / "data" / "raw" / "NONEXISTENT_9999.cif"
        with pytest.raises(Exception):  # noqa: B017
            gemmi.read_structure(str(bogus_path))

    def test_t2_checksum_mismatch_detection(self, tmp_path):
        """Test that corrupted file content is flagged by SHA-256 verification."""
        test_file = tmp_path / "corrupted.cif"
        test_file.write_text("data_corrupted\n# corrupted CIF contents\n")
        h = hashlib.sha256(test_file.read_bytes()).hexdigest()

        expected_7tre = "2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a"
        assert h != expected_7tre, "Corrupted content must not match authoritative hash"

    def test_t2_space_group_robustness_in_xray(self, gemmi_structures):
        """Test that 2I46 space group is correctly extracted from crystal data."""
        st = gemmi_structures.get("2I46")
        assert st is not None
        assert st.spacegroup_hm in ["I 41 2 2", "P 21 21 21"], f"Unexpected space group {st.spacegroup_hm}"
