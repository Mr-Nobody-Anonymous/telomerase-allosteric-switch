"""Shared pytest fixtures, paths, and reference ground truth values for E2E tests."""

from pathlib import Path

import gemmi
import numpy as np
import pytest
from Bio.PDB import MMCIFParser

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw"
DATA_MANIFEST_PATH = PROJECT_ROOT / "data_manifest.csv"
DATA_MANIFEST_ALT_PATH = PROJECT_ROOT / "data" / "manifest" / "data_manifest.csv"
RESULTS_TABLES_DIR = PROJECT_ROOT / "results" / "tables"
RESULTS_FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
REPORT_PATH = PROJECT_ROOT / "report.md"

# Fallback survey cache directory (used if data/raw has not yet been populated)
SURVEY_CIFS_DIR = PROJECT_ROOT / ".agents" / "teamwork" / "survey_pdb" / "cifs"

# Authoritative reference metadata
REFERENCE_DATA = {
    "structures": {
        "7TRE": {
            "sha256": "2c681322533afea566f7acd014b6b98c7dbfb51748e6207712dc4e65819a996a",
            "method": "ELECTRON MICROSCOPY",
            "resolution": 3.50,
            "tert_auth_chain": "A",
            "tert_label_chain": "B",
            "tpp1_auth_chain": "C",
            "tpp1_label_chain": "D",
            "tert_uniprot": "O14746",
            "tpp1_uniprot": "Q96AP0",
            "tert_modeled_count": 934,
            "tpp1_modeled_count": 152,
            "tert_modeled_range": (7, 1132),
            "tpp1_modeled_range": (90, 241),
        },
        "7TRD": {
            "sha256": "faf6d1aa6401d3dd0f66d97b818e95918ffbb5f9d05694751b8f107aafee67e4",
            "method": "ELECTRON MICROSCOPY",
            "resolution": 3.30,
            "tert_auth_chain": "A",
            "tert_label_chain": "B",
            "tert_modeled_count": 934,
            "tert_modeled_range": (7, 1132),
        },
        "7TRF": {
            "sha256": "e62419d5c2aac19d9f6c01ce983d4443ebbda8c8db55edb9ac93890f6defaadd",
            "method": "ELECTRON MICROSCOPY",
            "resolution": 3.70,
            "tert_auth_chain": "A",
            "tert_label_chain": "D",
            "tert_modeled_count": 934,
            "tert_modeled_range": (7, 1132),
        },
        "7TRC": {
            "sha256": "99b61793d50b3c79e3d65acc628a0764056c0c4f17682920ef2630e6ed282aea",
            "method": "ELECTRON MICROSCOPY",
            "resolution": 3.30,
        },
        "2I46": {
            "sha256": "4f4402ce3e9aadfa51456e795b8c71765d1420107a311677acb4488fddcaf779",
            "method": "X-RAY DIFFRACTION",
            "resolution": 2.70,
            "space_group": "I 41 2 2",
            "chain_a_tel_patch_resolved": True,
            "chain_b_tel_patch_disordered": True,
            "chain_a_modeled_count": 152,
        },
    },
    "contacts": {
        3.5: {"atom_contacts": 9, "residue_pairs": 5},
        4.0: {"atom_contacts": 46, "residue_pairs": 21},
        4.5: {"atom_contacts": 123, "residue_pairs": 32},
    },
    "bsa": {
        "sasa_tert": 49776.00,
        "sasa_tpp1": 10651.36,
        "sasa_complex": 58713.70,
        "delta_sasa": 1713.65,
        "bsa_interface": 856.82,
        "tolerance": 5.0,
    },
    "key_salt_bridges": [
        {"tert": ("ARG", 774, "NH1"), "tpp1": ("GLU", 168, "OE1"), "dist_approx": 2.32},
        {"tert": ("LYS", 78, "NZ"), "tpp1": ("GLU", 215, "OE2"), "dist_approx": 3.58},
    ],
    "tel_patch_residues": [168, 169, 170, 171, 183, 212],
    "rmsd_benchmarks": {
        "7tre_vs_7trd_global": 0.356,
        "7tre_vs_7trd_core": 0.338,
        "7tre_vs_2i46_tpp1": 1.966,
        "displacements_tpp1": {
            168: 1.253,
            169: 3.809,
            170: 4.293,
            171: 3.320,
            183: 0.512,  # Invariant rigid hydrophobic anchor
            212: 1.040,
        },
    },
    "expected_tables": [
        "data_manifest.csv",
        "chain_entity_mapping.csv",
        "residue_correspondence_tert.csv",
        "residue_correspondence_tpp1.csv",
        "interface_contacts_7tre_3.5A.csv",
        "interface_contacts_7tre_4.0A.csv",
        "interface_contacts_7tre_4.5A.csv",
        "interface_contacts_robustness.csv",
        "buried_surface_area.csv",
        "polar_interactions_7tre.csv",
        "cross_state_rmsd.csv",
    ],
}


def extract_heavy_atoms(chain):
    """Extract non-hydrogen atoms from a Gemmi chain."""
    atoms = []
    for res in chain:
        for atom in res:
            if atom.element.name != "H":
                atoms.append((
                    res.seqid.num,
                    res.name,
                    atom.name,
                    [atom.pos.x, atom.pos.y, atom.pos.z],
                    atom.b_iso
                ))
    return atoms


def extract_ca_atoms(chain):
    """Extract C-alpha atoms from a Gemmi chain as dict: auth_seq_num -> (resname, pos, b_iso)."""
    ca_dict = {}
    for res in chain:
        ca = res.find_atom("CA", "*")
        if ca:
            ca_dict[res.seqid.num] = (
                res.name,
                np.array([ca.pos.x, ca.pos.y, ca.pos.z], dtype=float),
                ca.b_iso
            )
    return ca_dict


def find_contacts(atoms_a, atoms_b, cutoff):
    """Pure-NumPy spatial cell hashing neighbor search between heavy atoms."""
    if cutoff <= 0.0:
        return []
    coords_a = np.array([a[3] for a in atoms_a], dtype=np.float32)
    coords_b = np.array([a[3] for a in atoms_b], dtype=np.float32)
    cell_size = cutoff
    grid = {}
    for idx_b, pt in enumerate(coords_b):
        cell = (int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size))
        grid.setdefault(cell, []).append(idx_b)

    contacts = []
    cutoff_sq = cutoff * cutoff
    for idx_a, pt in enumerate(coords_a):
        cx, cy, cz = int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor_cell = (cx + dx, cy + dy, cz + dz)
                    if neighbor_cell in grid:
                        for idx_b in grid[neighbor_cell]:
                            d_sq = np.sum((pt - coords_b[idx_b]) ** 2)
                            if d_sq <= cutoff_sq:
                                contacts.append((atoms_a[idx_a], atoms_b[idx_b], float(np.sqrt(d_sq))))
    return contacts


def kabsch_fit(p, q):
    """Kabsch algorithm for rigid superposition of mobile p onto target q."""
    centroid_p = np.mean(p, axis=0)
    centroid_q = np.mean(q, axis=0)
    p_centered = p - centroid_p
    q_centered = q - centroid_q
    h = np.dot(p_centered.T, q_centered)
    u, _s, vt = np.linalg.svd(h)
    det_val = np.linalg.det(np.dot(vt.T, u.T))
    d = np.sign(det_val) if abs(det_val) > 1e-10 else 1.0
    diag = np.diag([1.0, 1.0, d])
    r = np.dot(vt.T, np.dot(diag, u.T))
    t = centroid_q - np.dot(r, centroid_p)
    p_rotated = np.dot(p, r.T) + t
    rmsd = float(np.sqrt(np.mean(np.sum((p_rotated - q) ** 2, axis=1))))
    return rmsd, r, t, p_rotated


@pytest.fixture(scope="session")
def project_root() -> Path:
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def ref_data() -> dict:
    return REFERENCE_DATA


@pytest.fixture(scope="session")
def resolve_cif():
    """Returns a function that resolves the path to an mmCIF file."""
    def _resolver(entry_id: str) -> Path:
        raw_path = DATA_RAW_DIR / f"{entry_id}.cif"
        if raw_path.exists():
            return raw_path
        survey_path = SURVEY_CIFS_DIR / f"{entry_id}.cif"
        if survey_path.exists():
            return survey_path
        raise FileNotFoundError(
            f"CIF file for {entry_id} not found in {DATA_RAW_DIR} or {SURVEY_CIFS_DIR}"
        )
    return _resolver


@pytest.fixture(scope="session")
def gemmi_structures(resolve_cif):
    """Loads all 5 structures via Gemmi once per session."""
    structures = {}
    for pid in ["7TRE", "7TRD", "7TRF", "7TRC", "2I46"]:
        try:
            path = resolve_cif(pid)
            structures[pid] = gemmi.read_structure(str(path))
        except FileNotFoundError:
            structures[pid] = None
    return structures


@pytest.fixture(scope="session")
def biopython_parser():
    return MMCIFParser(QUIET=True)


@pytest.fixture
def manifest_dataframe():
    """Loads data_manifest.csv if available."""
    import pandas as pd
    if DATA_MANIFEST_PATH.exists():
        return pd.read_csv(DATA_MANIFEST_PATH)
    if DATA_MANIFEST_ALT_PATH.exists():
        return pd.read_csv(DATA_MANIFEST_ALT_PATH)
    return None
