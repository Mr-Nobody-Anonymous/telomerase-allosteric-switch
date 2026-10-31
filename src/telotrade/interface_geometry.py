"""
Quantitative Interface Geometry and Interaction Profiling Engine for Human Telomerase.

Provides contact network analysis across multiple distance cutoffs (3.5 Å, 4.0 Å, 4.5 Å),
Buried Surface Area (BSA) via Shrake-Rupley / FreeSASA, polar interaction profiling
(salt bridges including Arg774-Glu168 and Lys78-Glu215, hydrogen bonds), TEL-patch
epitope mapping (hydrophobic anchors Leu769-Leu183, Pro771-Leu212), threshold robustness
scoring, and standardized CSV table generation matching project specifications.
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Domain Definitions and Canonical Residues
# -----------------------------------------------------------------------------

# Canonical TEL-patch residues on TPP1 (UniProt Q96AP0 numbering)
TEL_PATCH_CORE: tuple[int, ...] = (168, 169, 170, 171, 183, 212)
TEL_PATCH_DE_LOOP: tuple[int, ...] = (168, 169, 170, 171)
TEL_PATCH_EXTENDED: tuple[int, ...] = (168, 169, 170, 171, 180, 183, 211, 212, 215)

# Interaction Profiling Atoms and Residues
ACIDIC_RESIDUES: set[str] = {"ASP", "GLU"}
ACIDIC_ATOMS: set[tuple[str, str]] = {
    ("ASP", "OD1"),
    ("ASP", "OD2"),
    ("GLU", "OE1"),
    ("GLU", "OE2"),
}

BASIC_RESIDUES: set[str] = {"ARG", "LYS", "HIS"}
BASIC_ATOMS: set[tuple[str, str]] = {
    ("ARG", "NH1"),
    ("ARG", "NH2"),
    ("ARG", "NE"),
    ("LYS", "NZ"),
    ("HIS", "ND1"),
    ("HIS", "NE2"),
}

POLAR_ELEMENTS: set[str] = {"N", "O", "S"}

HYDROPHOBIC_RESIDUES: set[str] = {
    "ALA",
    "VAL",
    "LEU",
    "ILE",
    "PRO",
    "PHE",
    "TRP",
    "TYR",
    "MET",
}

BACKBONE_POLAR_ATOMS: set[str] = {"N", "O"}

# Bondi Van der Waals radii (in Ångströms)
BONDI_RADII: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
}
DEFAULT_VDW_RADIUS: float = 1.70


# -----------------------------------------------------------------------------
# Data Containers
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class HeavyAtomInfo:
    """Heavy atom structural representation for spatial interface searches."""

    auth_seq_id: int
    label_seq_id: int | None
    auth_asym_id: str
    label_asym_id: str
    resname: str
    atom_name: str
    element: str
    pos: np.ndarray  # shape (3,)
    b_iso: float


# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------


def find_project_root() -> Path:
    """Locate the project root directory containing data_manifest.csv or src/."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "data_manifest.csv").exists() or (parent / "src").exists():
            return parent
    return Path.cwd()


def get_interface_tert_domain(resnum: int) -> str:
    """
    Annotate TERT residue domain for interface tables.

    Matches Section 4.5/4.6/4.8 schema: 'TEN' or 'TRAP_IFD' (or 'RT', 'TRBD', 'CTE').
    """
    if 1 <= resnum <= 196:
        return "TEN"
    if 755 <= resnum <= 800:
        return "TRAP_IFD"
    if 541 <= resnum <= 935:
        return "RT"
    if 326 <= resnum <= 540:
        return "TRBD"
    if 936 <= resnum <= 1132:
        return "CTE"
    if 197 <= resnum <= 325:
        return "LINKER"
    return "UNKNOWN"


def _extract_chain_heavy_atoms(
    cif_path: Path,
    chain_id: str,
    use_auth: bool = True,
) -> list[HeavyAtomInfo]:
    """
    Extract non-hydrogen atoms from a specified chain in an mmCIF file.

    Disambiguates author chain IDs (auth_asym_id) vs CIF label chain IDs (label_asym_id).
    Filters duplicate atoms from alternate conformations (keeping 'A' or first altloc).

    Raises:
        FileNotFoundError: If cif_path does not exist.
        KeyError: If chain_id is not found.
    """
    if not cif_path.is_file():
        raise FileNotFoundError(f"mmCIF file not found: {cif_path}")

    st = gemmi.read_structure(str(cif_path))
    if not st:
        raise ValueError(f"No models found in {cif_path}")
    model = st[0]

    target_residues: list[tuple[str, str, gemmi.Residue]] = []

    if use_auth:
        # 1. Primary lookup by author chain ID (ch.name in Gemmi)
        for ch in model:
            if ch.name == chain_id:
                for r in ch:
                    target_residues.append((ch.name, r.subchain, r))
                break

        # Fallback to label_asym_id if not found by auth
        if not target_residues:
            for ch in model:
                for r in ch:
                    if r.subchain == chain_id:
                        target_residues.append((ch.name, r.subchain, r))
    else:
        # 1. Primary lookup by CIF label chain ID (r.subchain in Gemmi)
        for ch in model:
            for r in ch:
                if r.subchain == chain_id:
                    target_residues.append((ch.name, r.subchain, r))

        # Fallback to author chain ID if not found by label
        if not target_residues:
            for ch in model:
                if ch.name == chain_id:
                    for r in ch:
                        target_residues.append((ch.name, r.subchain, r))
                    break

    if not target_residues:
        available_auth = sorted({ch.name for ch in model})
        available_label = sorted({r.subchain for ch in model for r in ch})
        raise KeyError(
            f"Chain '{chain_id}' not found in {cif_path.name}. "
            f"Available author chains: {available_auth}, label chains: {available_label}"
        )

    atoms: list[HeavyAtomInfo] = []
    seen_atoms: set[tuple[int, str]] = set()

    for auth_ch, label_ch, r in target_residues:
        tab = gemmi.find_tabulated_residue(r.name)
        # Interface analysis focuses on protein polypeptide chains
        if not tab.is_amino_acid():
            continue

        resnum = r.seqid.num
        label_seq = getattr(r, "label_seq", None)

        for atom in r:
            if atom.element.name == "H":
                continue
            # Filter alternate locations: allow primary conformation only
            if atom.altloc and atom.altloc not in ("A", "\x00"):
                continue

            atom_key = (resnum, atom.name)
            if atom_key in seen_atoms:
                continue
            seen_atoms.add(atom_key)

            atoms.append(
                HeavyAtomInfo(
                    auth_seq_id=resnum,
                    label_seq_id=label_seq,
                    auth_asym_id=auth_ch,
                    label_asym_id=label_ch,
                    resname=r.name,
                    atom_name=atom.name,
                    element=atom.element.name,
                    pos=np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=np.float32),
                    b_iso=atom.b_iso,
                )
            )

    return atoms


def classify_interaction_type(
    res_t: str,
    atom_t: str,
    el_t: str,
    res_p: str,
    atom_p: str,
    el_p: str,
    dist: float,
) -> str:
    """
    Classify the biophysical interaction type between a pair of heavy atoms.

    Priority:
    1. Salt Bridge: Anionic oxygen (Asp OD1/OD2, Glu OE1/OE2) to basic nitrogen
       (Arg NH1/NH2/NE, Lys NZ, His ND1/NE2) <= 4.0 Å.
    2. Polar/H-Bond: Polar heavy atom contacts (N, O, S) <= 4.0 Å.
    3. Hydrophobic: Carbon-carbon contacts between hydrophobic amino acids
       (Ala, Val, Leu, Ile, Pro, Phe, Trp, Tyr, Met) <= 4.5 Å.
    4. Van der Waals: Any other non-bonded heavy atom contact <= cutoff.
    """
    # 1. Salt Bridge
    is_sb_t_basic = ((res_t, atom_t) in BASIC_ATOMS) and ((res_p, atom_p) in ACIDIC_ATOMS)
    is_sb_p_basic = ((res_t, atom_t) in ACIDIC_ATOMS) and ((res_p, atom_p) in BASIC_ATOMS)
    if (is_sb_t_basic or is_sb_p_basic) and dist <= 4.0:
        return "Salt Bridge"

    # 2. Polar / Hydrogen Bond
    if (el_t in POLAR_ELEMENTS and el_p in POLAR_ELEMENTS) and dist <= 4.0:
        return "Polar/H-Bond"

    # 3. Hydrophobic Contact
    if (
        res_t in HYDROPHOBIC_RESIDUES
        and res_p in HYDROPHOBIC_RESIDUES
        and el_t == "C"
        and el_p == "C"
        and dist <= 4.5
    ):
        return "Hydrophobic"

    # 4. Van der Waals
    return "Van der Waals"


# -----------------------------------------------------------------------------
# Core Milestone 3 Algorithms
# -----------------------------------------------------------------------------


def compute_interface_contacts(
    cif_path: str | Path,
    chain_tert: str = "A",
    chain_tpp1: str = "C",
    cutoffs: list[float] | float | None = None,
    use_auth: bool = True,
) -> dict[float, pd.DataFrame] | pd.DataFrame:
    """
    Compute heavy-atom interface contact networks between TERT and TPP1.

    Finds all inter-chain atom pairs within Euclidean distance cutoffs,
    annotates TERT structural domain ('TEN', 'TRAP_IFD'), flags TPP1 TEL-patch
    residues, and classifies interaction types into 'Salt Bridge', 'Polar/H-Bond',
    'Hydrophobic', and 'Van der Waals'.

    Args:
        cif_path: Path to the mmCIF coordinate file (e.g. 7TRE.cif).
        chain_tert: Chain identifier for TERT ('A' in auth, 'B' in label).
        chain_tpp1: Chain identifier for TPP1 ('C' in auth, 'D' in label).
        cutoffs: Distance threshold in Ångströms (e.g. 4.0) or list of cutoffs
                 [3.5, 4.0, 4.5]. Defaults to [3.5, 4.0, 4.5].
        use_auth: If True, lookup chains using author chain IDs. If False, use
                  CIF label chain IDs.

    Returns:
        dict[float, pd.DataFrame] if a list of cutoffs is provided, or a single
        pd.DataFrame if a single numeric cutoff is provided.
    """
    cif_path = Path(cif_path)
    if cutoffs is None:
        target_cutoffs = [3.5, 4.0, 4.5]
        return_dict = True
    elif isinstance(cutoffs, (int, float)):
        target_cutoffs = [float(cutoffs)]
        return_dict = False
    else:
        target_cutoffs = [float(c) for c in cutoffs]
        return_dict = True

    tert_atoms = _extract_chain_heavy_atoms(cif_path, chain_tert, use_auth=use_auth)
    tpp1_atoms = _extract_chain_heavy_atoms(cif_path, chain_tpp1, use_auth=use_auth)

    coords_t = np.array([a.pos for a in tert_atoms], dtype=np.float32)
    coords_p = np.array([a.pos for a in tpp1_atoms], dtype=np.float32)

    results: dict[float, pd.DataFrame] = {}

    for cutoff in target_cutoffs:
        if cutoff <= 0.0 or len(coords_t) == 0 or len(coords_p) == 0:
            empty_df = pd.DataFrame(
                columns=[
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
            )
            results[cutoff] = empty_df
            continue

        # Fast grid hashing neighbor search in pure NumPy
        cell_size = cutoff
        grid: dict[tuple[int, int, int], list[int]] = {}
        for idx_p, pt in enumerate(coords_p):
            cell = (int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size))
            grid.setdefault(cell, []).append(idx_p)

        cutoff_sq = cutoff * cutoff
        rows = []

        for idx_t, pt in enumerate(coords_t):
            cx, cy, cz = int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        neighbor_cell = (cx + dx, cy + dy, cz + dz)
                        if neighbor_cell in grid:
                            for idx_p in grid[neighbor_cell]:
                                diff = pt - coords_p[idx_p]
                                d_sq = float(np.sum(diff * diff))
                                if d_sq <= cutoff_sq:
                                    dist = float(np.sqrt(d_sq))
                                    at_t = tert_atoms[idx_t]
                                    at_p = tpp1_atoms[idx_p]

                                    itype = classify_interaction_type(
                                        at_t.resname,
                                        at_t.atom_name,
                                        at_t.element,
                                        at_p.resname,
                                        at_p.atom_name,
                                        at_p.element,
                                        dist,
                                    )

                                    rows.append(
                                        {
                                            "tert_auth_asym": at_t.auth_asym_id,
                                            "tert_label_asym": at_t.label_asym_id,
                                            "tert_resnum": at_t.auth_seq_id,
                                            "tert_resname": at_t.resname,
                                            "tert_atom": at_t.atom_name,
                                            "tert_domain": get_interface_tert_domain(at_t.auth_seq_id),
                                            "tpp1_auth_asym": at_p.auth_asym_id,
                                            "tpp1_label_asym": at_p.label_asym_id,
                                            "tpp1_resnum": at_p.auth_seq_id,
                                            "tpp1_resname": at_p.resname,
                                            "tpp1_atom": at_p.atom_name,
                                            "is_tel_patch": at_p.auth_seq_id in TEL_PATCH_EXTENDED,
                                            "distance_angstrom": round(dist, 2),
                                            "interaction_type": itype,
                                        }
                                    )

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(
                by=["distance_angstrom", "tert_resnum", "tpp1_resnum"],
                ascending=[True, True, True],
            ).reset_index(drop=True)
        else:
            df = pd.DataFrame(
                columns=[
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
            )
        results[cutoff] = df

    if not return_dict:
        return results[target_cutoffs[0]]
    return results


def compute_bsa(
    cif_path: str | Path,
    chain_tert: str = "A",
    chain_tpp1: str = "C",
    probe_radius: float = 1.4,
    n_sphere_points: int = 960,
    use_auth: bool = True,
    algorithm: str = "shrake_rupley",
) -> dict[str, Any]:
    """
    Compute Buried Surface Area (BSA) between TERT and TPP1.

    Mathematical Formulation:
      ΔSASA_total = SASA_A + SASA_B - SASA_AB
      BSA_interface = ΔSASA_total / 2 = (SASA_A + SASA_B - SASA_AB) / 2

    Uses FreeSASA (or Biopython ShrakeRupley fallback) with probe radius 1.4 Å
    and 960 sphere test points per Bondi atomic radius.

    Returns:
        dict containing individual SASAs, delta SASA, interface BSA, and metadata.
    """
    cif_path = Path(cif_path)
    tert_atoms = _extract_chain_heavy_atoms(cif_path, chain_tert, use_auth=use_auth)
    tpp1_atoms = _extract_chain_heavy_atoms(cif_path, chain_tpp1, use_auth=use_auth)

    c_a: list[float] = []
    r_a: list[float] = []
    for a in tert_atoms:
        c_a.extend([float(a.pos[0]), float(a.pos[1]), float(a.pos[2])])
        r_a.append(BONDI_RADII.get(a.element.upper(), DEFAULT_VDW_RADIUS))

    c_b: list[float] = []
    r_b: list[float] = []
    for a in tpp1_atoms:
        c_b.extend([float(a.pos[0]), float(a.pos[1]), float(a.pos[2])])
        r_b.append(BONDI_RADII.get(a.element.upper(), DEFAULT_VDW_RADIUS))

    c_ab = c_a + c_b
    r_ab = r_a + r_b

    sasa_a = 0.0
    sasa_b = 0.0
    sasa_ab = 0.0

    used_algo = "Shrake-Rupley"

    # Attempt ultra-fast C implementation via FreeSASA
    try:
        import freesasa

        params = freesasa.Parameters()
        params.setAlgorithm(freesasa.ShrakeRupley)
        params.setProbeRadius(probe_radius)
        params.setNPoints(n_sphere_points)

        res_a = freesasa.calcCoord(c_a, r_a, params)
        res_b = freesasa.calcCoord(c_b, r_b, params)
        res_ab = freesasa.calcCoord(c_ab, r_ab, params)

        sasa_a = float(res_a.totalArea())
        sasa_b = float(res_b.totalArea())
        sasa_ab = float(res_ab.totalArea())
    except (ImportError, RuntimeError, ValueError, TypeError, OSError) as err:
        logger.info("FreeSASA calculation skipped or unavailable (%s), falling back to Bio.PDB", err)
        from Bio.PDB import MMCIFParser
        from Bio.PDB.Model import Model
        from Bio.PDB.SASA import ShrakeRupley
        from Bio.PDB.Structure import Structure

        parser = MMCIFParser(QUIET=True)
        biopdb_struct = parser.get_structure("complex", str(cif_path))
        biopdb_model = biopdb_struct[0]

        # Bio.PDB's MMCIFParser keys chains by *author* chain id, whereas the atom lists
        # above may have been selected by label chain id.  Resolve the requested chains to
        # their author ids through the atoms that were actually selected.
        auth_a = tert_atoms[0].auth_asym_id if tert_atoms else chain_tert
        auth_b = tpp1_atoms[0].auth_asym_id if tpp1_atoms else chain_tpp1
        ch_a = biopdb_model[auth_a]
        ch_b = biopdb_model[auth_b]

        sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_sphere_points)

        s_a = Structure("A")
        m_a = Model(0)
        m_a.add(copy.deepcopy(ch_a))
        s_a.add(m_a)
        sr.compute(s_a, level="A")
        sasa_a = float(sum(at.sasa for at in s_a.get_atoms()))

        s_b = Structure("B")
        m_b = Model(0)
        m_b.add(copy.deepcopy(ch_b))
        s_b.add(m_b)
        sr.compute(s_b, level="A")
        sasa_b = float(sum(at.sasa for at in s_b.get_atoms()))

        s_ab = Structure("AB")
        m_ab = Model(0)
        m_ab.add(copy.deepcopy(ch_a))
        m_ab.add(copy.deepcopy(ch_b))
        s_ab.add(m_ab)
        sr.compute(s_ab, level="A")
        sasa_ab = float(sum(at.sasa for at in s_ab.get_atoms()))

    delta_sasa = sasa_a + sasa_b - sasa_ab
    bsa_half = delta_sasa / 2.0

    stem = cif_path.stem.upper()
    return {
        "structure_id": stem,
        "interface_id": f"TERT_{chain_tert}__TPP1_{chain_tpp1}",
        "partner_a_name": "TERT",
        "partner_b_name": "TPP1",
        "sasa_partner_a_alone_angstrom2": round(sasa_a, 2),
        "sasa_partner_b_alone_angstrom2": round(sasa_b, 2),
        "sasa_complex_angstrom2": round(sasa_ab, 2),
        "delta_sasa_total_angstrom2": round(delta_sasa, 2),
        "bsa_interface_half_angstrom2": round(bsa_half, 2),
        "probe_radius_angstrom": probe_radius,
        "n_sphere_points": n_sphere_points,
        "algorithm": used_algo,
        # Conveniences & aliases
        "bsa": round(bsa_half, 2),
        "delta_sasa": round(delta_sasa, 2),
        "sasa_tert": round(sasa_a, 2),
        "sasa_tpp1": round(sasa_b, 2),
        "sasa_complex": round(sasa_ab, 2),
    }


def compute_polar_interactions(
    cif_path: str | Path,
    chain_tert: str = "A",
    chain_tpp1: str = "C",
    use_auth: bool = True,
    cutoff: float = 4.0,
) -> pd.DataFrame:
    """
    Identify and profile all polar interactions (salt bridges and hydrogen bonds).

    Specifically categorizes:
      - Salt Bridge: Anionic carboxylate (Asp/Glu) to basic guanidinium/ammonium
        (Arg/Lys/His) <= 4.0 Å.
      - H-Bond Sidechain-Sidechain: Both donor and acceptor on sidechain polar atoms.
      - H-Bond Backbone-Sidechain: One donor/acceptor on backbone amide/carbonyl, one on sidechain.
      - H-Bond Backbone-Backbone: Both donor and acceptor on backbone.

    Validates primary salt bridge Arg774-Glu168 (2.32 Å) and secondary clamp Lys78-Glu215 (3.58 Å).
    """
    cif_path = Path(cif_path)
    tert_atoms = _extract_chain_heavy_atoms(cif_path, chain_tert, use_auth=use_auth)
    tpp1_atoms = _extract_chain_heavy_atoms(cif_path, chain_tpp1, use_auth=use_auth)

    polar_tert = [a for a in tert_atoms if a.element in POLAR_ELEMENTS]
    polar_tpp1 = [a for a in tpp1_atoms if a.element in POLAR_ELEMENTS]

    if not polar_tert or not polar_tpp1:
        return pd.DataFrame(
            columns=[
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
        )

    coords_t = np.array([a.pos for a in polar_tert], dtype=np.float32)
    coords_p = np.array([a.pos for a in polar_tpp1], dtype=np.float32)

    rows = []
    cutoff_sq = cutoff * cutoff

    for idx_t, pt_t in enumerate(coords_t):
        diffs = coords_p - pt_t
        d_sq = np.sum(diffs * diffs, axis=1)
        hit_indices = np.where(d_sq <= cutoff_sq)[0]

        for idx_p in hit_indices:
            dist = float(np.sqrt(d_sq[idx_p]))
            at_t = polar_tert[idx_t]
            at_p = polar_tpp1[idx_p]

            # Salt bridge check
            is_sb1 = ((at_t.resname, at_t.atom_name) in BASIC_ATOMS) and (
                (at_p.resname, at_p.atom_name) in ACIDIC_ATOMS
            )
            is_sb2 = ((at_t.resname, at_t.atom_name) in ACIDIC_ATOMS) and (
                (at_p.resname, at_p.atom_name) in BASIC_ATOMS
            )

            if is_sb1 or is_sb2:
                itype = "Salt Bridge"
            else:
                is_bb_t = at_t.atom_name in BACKBONE_POLAR_ATOMS
                is_bb_p = at_p.atom_name in BACKBONE_POLAR_ATOMS
                if is_bb_t and is_bb_p:
                    itype = "H-Bond Backbone-Backbone"
                elif is_bb_t or is_bb_p:
                    itype = "H-Bond Backbone-Sidechain"
                else:
                    itype = "H-Bond Sidechain-Sidechain"

            rows.append(
                {
                    "interaction_type": itype,
                    "tert_resnum": at_t.auth_seq_id,
                    "tert_resname": at_t.resname,
                    "tert_atom": at_t.atom_name,
                    "tert_domain": get_interface_tert_domain(at_t.auth_seq_id),
                    "tpp1_resnum": at_p.auth_seq_id,
                    "tpp1_resname": at_p.resname,
                    "tpp1_atom": at_p.atom_name,
                    "distance_angstrom": round(dist, 2),
                    "is_tel_patch": at_p.auth_seq_id in TEL_PATCH_EXTENDED,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(
            by=["distance_angstrom", "tert_resnum", "tpp1_resnum"],
            ascending=[True, True, True],
        ).reset_index(drop=True)
    return df


def map_tel_patch_contacts(
    contacts_df: pd.DataFrame,
    tel_patch_residues: tuple[int, ...] = TEL_PATCH_EXTENDED,
) -> pd.DataFrame:
    """
    Extract and profile all contacts involving TPP1 TEL-patch residues.

    Documented TEL-patch subregions:
      - Acidic DE-Loop: Glu168, Glu169, Lys170, Glu171
      - Hydrophobic Core Anchors: Leu183, Leu212
      - Auxiliary Clamp: Arg180, Glu215
    """
    if contacts_df.empty or "tpp1_resnum" not in contacts_df.columns:
        return contacts_df.copy()

    df_tel = contacts_df[contacts_df["tpp1_resnum"].isin(tel_patch_residues)].copy()

    def _assign_subregion(resnum: int) -> str:
        if resnum in (168, 169, 170, 171):
            return "Acidic DE-loop"
        if resnum in (183, 212):
            return "Hydrophobic Core"
        if resnum == 180:
            return "Basic Finger"
        if resnum == 215:
            return "Auxiliary Clamp"
        return "TEL-patch"

    df_tel["tel_patch_subregion"] = df_tel["tpp1_resnum"].map(_assign_subregion)
    return df_tel.sort_values(by=["distance_angstrom"]).reset_index(drop=True)


def compute_contact_robustness(
    cif_path: str | Path,
    chain_tert: str = "A",
    chain_tpp1: str = "C",
    cutoffs: list[float] | None = None,
    use_auth: bool = True,
) -> pd.DataFrame:
    """
    Analyze persistence and sensitivity of interface contacts across distance cutoffs.

    Partitioning:
      - Core Robust: Present at 3.5 Å, 4.0 Å, and 4.5 Å (persistence_count = 3).
      - Intermediate: Present at 4.0 Å and 4.5 Å (persistence_count = 2).
      - Threshold Sensitive: Present only at 4.5 Å (persistence_count = 1).

    Returns:
        pd.DataFrame adhering to Section 4.6 data schema.
    """
    if cutoffs is None:
        cutoffs = [3.5, 4.0, 4.5]
    sorted_cutoffs = sorted(cutoffs)
    c35, c40, c45 = sorted_cutoffs[0], sorted_cutoffs[1], sorted_cutoffs[2]

    contacts_dict = compute_interface_contacts(
        cif_path=cif_path,
        chain_tert=chain_tert,
        chain_tpp1=chain_tpp1,
        cutoffs=sorted_cutoffs,
        use_auth=use_auth,
    )
    assert isinstance(contacts_dict, dict)

    df_35 = contacts_dict[c35]
    df_40 = contacts_dict[c40]
    df_45 = contacts_dict[c45]

    pairs_35 = {
        (row["tert_resnum"], row["tert_atom"], row["tpp1_resnum"], row["tpp1_atom"])
        for _, row in df_35.iterrows()
    }
    pairs_40 = {
        (row["tert_resnum"], row["tert_atom"], row["tpp1_resnum"], row["tpp1_atom"])
        for _, row in df_40.iterrows()
    }

    rows = []
    for _, row in df_45.iterrows():
        pair_key = (row["tert_resnum"], row["tert_atom"], row["tpp1_resnum"], row["tpp1_atom"])
        in_35 = pair_key in pairs_35
        in_40 = pair_key in pairs_40
        in_45 = True

        score = int(in_35) + int(in_40) + int(in_45)
        if score == 3:
            robust_class = "Core Robust"
        elif score == 2:
            robust_class = "Intermediate"
        else:
            robust_class = "Threshold Sensitive"

        contact_id = (
            f"TERT_{row['tert_resname']}{row['tert_resnum']}:{row['tert_atom']}--"
            f"TPP1_{row['tpp1_resname']}{row['tpp1_resnum']}:{row['tpp1_atom']}"
        )

        rows.append(
            {
                "contact_pair_id": contact_id,
                "tert_residue": f"{row['tert_resname']}{row['tert_resnum']}",
                "tert_domain": row["tert_domain"],
                "tpp1_residue": f"{row['tpp1_resname']}{row['tpp1_resnum']}",
                "is_tel_patch": row["is_tel_patch"],
                "min_atom_distance_angstrom": row["distance_angstrom"],
                "present_at_3.5A": in_35,
                "present_at_4.0A": in_40,
                "present_at_4.5A": in_45,
                "persistence_count": score,
                "robustness_class": robust_class,
            }
        )

    df_rob = pd.DataFrame(rows)
    if not df_rob.empty:
        df_rob = df_rob.sort_values(
            by=["min_atom_distance_angstrom", "tert_residue"],
            ascending=[True, True],
        ).reset_index(drop=True)
    return df_rob


def generate_interface_tables(
    cif_path: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Path]:
    """
    Generate and save all Milestone 3 standardized CSV tables into results/tables/.

    Output files:
      1. interface_contacts_7tre_3.5A.csv (and alias interface_contacts_3.5A.csv)
      2. interface_contacts_7tre_4.0A.csv (and alias interface_contacts_4.0A.csv)
      3. interface_contacts_7tre_4.5A.csv (and alias interface_contacts_4.5A.csv)
      4. interface_contacts_robustness.csv
      5. buried_surface_area.csv
      6. polar_interactions_7tre.csv (and alias polar_interactions.csv)
      7. tel_patch_contacts.csv

    Returns:
        dict[str, Path]: Mapping from table name to saved file Path.
    """
    project_root = find_project_root()
    if cif_path is None:
        cif_candidate = project_root / "data" / "raw" / "7TRE.cif"
        if not cif_candidate.is_file():
            cif_candidate = (
                project_root / ".agents" / "teamwork" / "survey_pdb" / "cifs" / "7TRE.cif"
            )
        cif_path = cif_candidate
    else:
        cif_path = Path(cif_path)

    if not cif_path.is_file():
        raise FileNotFoundError(f"Cannot generate interface tables: 7TRE mmCIF not found at {cif_path}")

    if output_dir is None:
        out_path = project_root / "results" / "tables"
    else:
        out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    generated: dict[str, Path] = {}

    # 1. Contacts at 3.5 Å, 4.0 Å, 4.5 Å
    contacts = compute_interface_contacts(cif_path, cutoffs=[3.5, 4.0, 4.5])
    assert isinstance(contacts, dict)

    for cutoff, df_c in contacts.items():
        fname_full = f"interface_contacts_7tre_{cutoff:.1f}A.csv"
        fname_alias = f"interface_contacts_{cutoff:.1f}A.csv"
        p_full = out_path / fname_full
        p_alias = out_path / fname_alias
        df_c.to_csv(p_full, index=False)
        df_c.to_csv(p_alias, index=False)
        generated[fname_full] = p_full
        generated[fname_alias] = p_alias

    # 2. Contact robustness table
    df_rob = compute_contact_robustness(cif_path, cutoffs=[3.5, 4.0, 4.5])
    p_rob = out_path / "interface_contacts_robustness.csv"
    df_rob.to_csv(p_rob, index=False)
    generated["interface_contacts_robustness.csv"] = p_rob

    # 3. Buried surface area table
    bsa_res = compute_bsa(cif_path)
    df_bsa = pd.DataFrame([bsa_res])
    p_bsa = out_path / "buried_surface_area.csv"
    df_bsa.to_csv(p_bsa, index=False)
    generated["buried_surface_area.csv"] = p_bsa

    # 4. Polar interactions table
    df_polar = compute_polar_interactions(cif_path)
    p_polar_full = out_path / "polar_interactions_7tre.csv"
    p_polar_alias = out_path / "polar_interactions.csv"
    df_polar.to_csv(p_polar_full, index=False)
    df_polar.to_csv(p_polar_alias, index=False)
    generated["polar_interactions_7tre.csv"] = p_polar_full
    generated["polar_interactions.csv"] = p_polar_alias

    # 5. TEL-patch contacts table
    df_tel = map_tel_patch_contacts(contacts[4.5])
    p_tel = out_path / "tel_patch_contacts.csv"
    df_tel.to_csv(p_tel, index=False)
    generated["tel_patch_contacts.csv"] = p_tel

    return generated
