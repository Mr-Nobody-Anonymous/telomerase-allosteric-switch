"""
Cross-State Comparative Structural Analysis & Robustness Testing Engine.

Provides rigid-body Singular Value Decomposition (SVD) Kabsch superposition,
cross-state C-alpha RMSD comparisons (7TRE vs 7TRD, 7TRF, and 2I46 Chain A),
domain-specific RMSD partitions (TEN, TRBD, RT, IFD_TRAP, CTE, Catalytic Core),
per-residue core-fit displacement profiling, TPP1 TEL-patch conformational shifts
(quantifying DE-loop plasticity up to 4.3 Å vs Leu183 rigidity at 0.51 Å),
interface cutoff robustness persistence scoring, and CSV table generation matching
project specifications.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, NamedTuple

import gemmi
import numpy as np
import pandas as pd

from .interface_geometry import compute_contact_robustness
from .sequence_mapping import (
    TEL_PATCH_EXTENDED,
    TEL_PATCH_RESIDUES,
    get_tert_domain,
    select_valid_tpp1_chain,
)

logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Canonical Domain Definitions and Constants
# -----------------------------------------------------------------------------

# TERT structural domains (author residue numbers matching UniProt O14746)
CANONICAL_DOMAIN_RANGES: dict[str, tuple[int, int]] = {
    "TEN": (1, 196),
    "TRBD": (326, 540),
    "RT": (541, 935),
    "IFD_TRAP": (755, 785),
    "CTE": (936, 1132),
    "Catalytic Core": (326, 1132),
}

# Empirical benchmark expectations
BENCHMARK_RMSDS: dict[str, float] = {
    "7tre_vs_7trd_global": 0.356,
    "7tre_vs_7trd_core": 0.338,
    "7tre_vs_7trd_ten_alone": 0.411,
    "7tre_vs_7trd_ten_core_fit": 0.434,
    "7tre_vs_7trf_global": 0.356,
    "7tre_vs_7trf_core": 0.338,
    "7tre_vs_2i46_tpp1_global": 1.966,
    "tpp1_leu183_displacement": 0.512,
    "tpp1_lys170_displacement": 4.293,
}


class KabschResult(NamedTuple):
    """
    Result of SVD Kabsch rigid-body superposition.

    Attributes:
        aligned_coords: Superimposed coordinates of mobile structure (N x 3).
        rotation: 3x3 orthogonal rotation matrix in SO(3).
        translation: 3D translation vector.
        rmsd: Root-mean-square deviation in Ångströms.
    """

    aligned_coords: np.ndarray
    rotation: np.ndarray
    translation: np.ndarray
    rmsd: float


# -----------------------------------------------------------------------------
# Path and Structure Resolution Helpers
# -----------------------------------------------------------------------------


def find_project_root() -> Path:
    """Locate the project root directory containing data_manifest.csv or src/."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "data_manifest.csv").exists() or (parent / "src").exists():
            return parent
    return Path.cwd()


def resolve_cif_path(
    entry_id_or_path: str | Path,
    raw_dir: str | Path | None = None,
) -> Path:
    """
    Resolve an mmCIF file path from either an entry ID (e.g. '7TRE') or explicit path.

    Checks:
    1. Direct file path if valid.
    2. raw_dir / {entry_id}.cif.
    3. project_root / 'data' / 'raw' / {entry_id}.cif.
    4. project_root / '.agents' / 'teamwork' / 'survey_pdb' / 'cifs' / {entry_id}.cif.

    Raises:
        FileNotFoundError: If the mmCIF file cannot be located.
    """
    p = Path(entry_id_or_path)
    if p.is_file():
        return p

    pid = p.stem.upper()
    candidates: list[Path] = []

    if raw_dir is not None:
        rdir = Path(raw_dir)
        candidates.append(rdir / f"{pid}.cif")
        candidates.append(rdir / f"{pid.lower()}.cif")

    root = find_project_root()
    candidates.extend(
        [
            root / "data" / "raw" / f"{pid}.cif",
            root / "data" / "raw" / f"{pid.lower()}.cif",
            root / ".agents" / "teamwork" / "survey_pdb" / "cifs" / f"{pid}.cif",
            root / ".agents" / "teamwork" / "survey_pdb" / "cifs" / f"{pid.lower()}.cif",
        ]
    )

    for cand in candidates:
        if cand.is_file():
            return cand

    raise FileNotFoundError(
        f"mmCIF file for '{entry_id_or_path}' not found. Searched locations: {candidates}"
    )


def extract_ca_coordinates(
    cif_path: str | Path,
    chain_id: str,
    use_auth: bool = True,
) -> dict[int, tuple[str, np.ndarray, float]]:
    """
    Extract C-alpha atomic coordinates from an mmCIF file for a specified chain.

    Disambiguates author chain IDs (auth_asym_id) vs CIF label chain IDs (label_asym_id).
    Filters duplicate atoms from alternate conformations (retaining primary 'A' or '\\x00').

    Args:
        cif_path: Path or identifier for mmCIF file.
        chain_id: Chain identifier (e.g. 'A', 'C', etc.).
        use_auth: If True, lookup chain by author ID and use author residue numbers.
                  If False, lookup by CIF label ID and use label residue numbers.

    Returns:
        dict[int, tuple[str, np.ndarray, float]]: Mapping of residue number ->
            (3-letter resname, coordinate array shape (3,), isotropic B-factor).

    Raises:
        FileNotFoundError: If cif_path does not exist.
        KeyError: If chain_id is not found in the structure.
    """
    resolved_path = resolve_cif_path(cif_path)
    st = gemmi.read_structure(str(resolved_path))
    if not st:
        raise ValueError(f"No models found in {resolved_path}")
    model = st[0]

    target_residues: list[gemmi.Residue] = []

    if use_auth:
        for ch in model:
            if ch.name == chain_id:
                for r in ch:
                    tab = gemmi.find_tabulated_residue(r.name)
                    if tab.is_amino_acid():
                        # Exclude solvent heteroatoms sharing author chain ID in 2I46
                        if "2I46" in resolved_path.name and r.subchain not in ("A", "B"):
                            continue
                        target_residues.append(r)
                break

        if not target_residues:
            for ch in model:
                for r in ch:
                    if r.subchain == chain_id:
                        tab = gemmi.find_tabulated_residue(r.name)
                        if tab.is_amino_acid():
                            target_residues.append(r)
    else:
        for ch in model:
            for r in ch:
                if r.subchain == chain_id:
                    tab = gemmi.find_tabulated_residue(r.name)
                    if tab.is_amino_acid():
                        target_residues.append(r)

        if not target_residues:
            for ch in model:
                if ch.name == chain_id:
                    for r in ch:
                        tab = gemmi.find_tabulated_residue(r.name)
                        if tab.is_amino_acid():
                            target_residues.append(r)
                    break

    if not target_residues:
        available_auth = sorted({ch.name for ch in model})
        available_label = sorted({r.subchain for ch in model for r in ch})
        raise KeyError(
            f"Chain '{chain_id}' not found in {resolved_path.name}. "
            f"Available author chains: {available_auth}, label chains: {available_label}"
        )

    ca_dict: dict[int, tuple[str, np.ndarray, float]] = {}

    for r in target_residues:
        resnum = r.seqid.num if use_auth else (r.label_seq if r.label_seq is not None else r.seqid.num)
        if resnum in ca_dict:
            continue

        ca_atom = None
        for atom in r:
            if atom.name == "CA":
                if atom.altloc and atom.altloc not in ("A", "\x00"):
                    continue
                ca_atom = atom
                break

        if ca_atom is not None:
            pos = np.array([ca_atom.pos.x, ca_atom.pos.y, ca_atom.pos.z], dtype=float)
            ca_dict[resnum] = (r.name, pos, float(ca_atom.b_iso))

    return ca_dict


# -----------------------------------------------------------------------------
# Rigid-Body SVD Kabsch Superposition Engine
# -----------------------------------------------------------------------------


def svd_kabsch_align(
    coords_target: np.ndarray,
    coords_mobile: np.ndarray,
) -> KabschResult:
    """
    Perform optimal rigid-body superposition of mobile coordinates onto target coordinates
    using the Kabsch algorithm via Singular Value Decomposition (SVD).

    Formulation:
      1. Centroid centering:
         P_c = P - mean(P),  Q_c = Q - mean(Q)
      2. Cross-covariance matrix:
         H = P_c^T @ Q_c
      3. Singular Value Decomposition:
         H = U @ Sigma @ V^T
      4. Proper rotation enforcement in SO(3) (prevent inversion/reflection):
         d = sign(det(V @ U^T))
         R = V @ diag(1, 1, d) @ U^T
      5. Translation vector:
         t = mean(Q) - R @ mean(P)
      6. Transformed mobile coordinates:
         P_aligned = P @ R^T + t
      7. Root-Mean-Square Deviation:
         RMSD = sqrt( (1/N) * sum(||P_aligned_i - Q_i||^2) )

    Args:
        coords_target: (N, 3) reference/target coordinates.
        coords_mobile: (N, 3) mobile coordinates to superimpose onto target.

    Returns:
        KabschResult: NamedTuple containing:
            - aligned_coords: Superimposed mobile coordinates (N, 3).
            - rotation: (3, 3) rotation matrix strictly satisfying det(R) == +1.0.
            - translation: (3,) translation vector.
            - rmsd: Root-mean-square deviation in Ångströms.

    Raises:
        ValueError: If coordinate arrays are empty, non-2D, do not have shape (N, 3),
                    or have mismatched point counts.
    """
    p = np.asarray(coords_mobile, dtype=float)
    q = np.asarray(coords_target, dtype=float)

    if p.size == 0 or q.size == 0:
        raise ValueError("Cannot perform Kabsch alignment on empty coordinate arrays.")

    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"coords_mobile must have shape (N, 3), got {p.shape}")

    if q.ndim != 2 or q.shape[1] != 3:
        raise ValueError(f"coords_target must have shape (N, 3), got {q.shape}")

    if p.shape[0] != q.shape[0]:
        raise ValueError(
            f"Coordinate point count mismatch: {p.shape[0]} mobile points vs "
            f"{q.shape[0]} target points."
        )

    n_points = p.shape[0]

    # Handle 1-point trivial case
    if n_points == 1:
        t = q[0] - p[0]
        r = np.eye(3, dtype=float)
        aligned = p + t
        return KabschResult(aligned_coords=aligned, rotation=r, translation=t, rmsd=0.0)

    # 1. Compute centroids and center point clouds
    centroid_p = np.mean(p, axis=0)
    centroid_q = np.mean(q, axis=0)

    p_centered = p - centroid_p
    q_centered = q - centroid_q

    # 2. Cross-covariance matrix H = P_c^T @ Q_c
    h = np.dot(p_centered.T, q_centered)

    # 3. SVD of cross-covariance matrix
    # np.linalg.svd returns u, s, vt such that h = u @ diag(s) @ vt
    u, _s, vt = np.linalg.svd(h)
    v = vt.T
    ut = u.T

    # 4. Enforce proper right-handed rotation (reflection check)
    det_vu = float(np.linalg.det(np.dot(v, ut)))
    d = np.sign(det_vu) if abs(det_vu) > 1e-10 else 1.0

    diag = np.diag([1.0, 1.0, d])
    r = np.dot(v, np.dot(diag, ut))

    # Guarantee numerical det(R) == +1.0
    if abs(np.linalg.det(r) - 1.0) > 1e-4:
        logger.warning("Rotation matrix determinant deviated from +1.0: %f", np.linalg.det(r))

    # 5. Translation vector
    t = centroid_q - np.dot(r, centroid_p)

    # 6. Transform mobile points
    aligned = np.dot(p, r.T) + t

    # 7. Calculate RMSD
    diff = aligned - q
    rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))

    return KabschResult(aligned_coords=aligned, rotation=r, translation=t, rmsd=rmsd)


def apply_kabsch_transform(
    coords: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> np.ndarray:
    """
    Apply a rigid-body rotation and translation to an arbitrary set of 3D coordinates.

    Args:
        coords: (M, 3) coordinates.
        rotation: (3, 3) rotation matrix.
        translation: (3,) translation vector.

    Returns:
        np.ndarray: (M, 3) transformed coordinates.
    """
    c = np.asarray(coords, dtype=float)
    if c.size == 0:
        return c.copy()
    return np.dot(c, rotation.T) + translation


# -----------------------------------------------------------------------------
# Cross-State Superposition & RMSD Calculations
# -----------------------------------------------------------------------------


def compute_cross_state_rmsd(
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Compute comprehensive cross-state RMSD metrics across telomerase complexes.

    Evaluates:
      - 7TRE vs 7TRD:
        * Global C-alpha (n=934)
        * Catalytic Core (residues 326-1132, n=770, benchmark ~0.338 Å)
        * TEN Domain alone (residues 7-196, n=160, benchmark ~0.411 Å)
        * TEN Domain on Core Fit (displacement ~0.434 Å)
      - 7TRE vs 7TRF:
        * Global C-alpha (n=934)
        * Catalytic Core (n=770, benchmark ~0.338 Å)
        * TEN Domain alone (n=160, benchmark ~0.411 Å)
        * TEN Domain on Core Fit (displacement ~0.434 Å)
      - 7TRE vs 2I46 Chain A (strictly Chain A!):
        * TPP1 OB-fold / Global C-alpha (residues 90-241, n=152, benchmark ~1.966 Å)
        * TEL-patch Core (residues 168, 169, 170, 171, 183, 212, n=6)

    Returns:
        pd.DataFrame adhering to Section 4.9 cross_state_rmsd.csv schema:
            target_structure, reference_structure, alignment_selection,
            num_aligned_ca_atoms, superposition_method, rmsd_angstrom,
            mean_ca_displacement_angstrom, max_ca_displacement_angstrom,
            max_displacement_residue
    """
    cif_7tre = resolve_cif_path("7TRE", raw_dir=raw_dir)
    cif_7trd = resolve_cif_path("7TRD", raw_dir=raw_dir)
    cif_7trf = resolve_cif_path("7TRF", raw_dir=raw_dir)
    cif_2i46 = resolve_cif_path("2I46", raw_dir=raw_dir)

    # Extract TERT C-alphas
    ca_7tre_tert = extract_ca_coordinates(cif_7tre, "A", use_auth=True)
    ca_7trd_tert = extract_ca_coordinates(cif_7trd, "A", use_auth=True)
    ca_7trf_tert = extract_ca_coordinates(cif_7trf, "A", use_auth=True)

    # Extract TPP1 C-alphas (7TRE chain C vs 2I46 Chain A strictly)
    valid_2i46_chain = select_valid_tpp1_chain(cif_2i46, preferred_chain="A", require_tel_patch=True)
    ca_7tre_tpp1 = extract_ca_coordinates(cif_7tre, "C", use_auth=True)
    ca_2i46_tpp1 = extract_ca_coordinates(cif_2i46, valid_2i46_chain, use_auth=True)

    rows: list[dict[str, Any]] = []

    def _eval_alignment(
        target_name: str,
        ref_name: str,
        selection_name: str,
        ca_dict_t: dict[int, tuple[str, np.ndarray, float]],
        ca_dict_r: dict[int, tuple[str, np.ndarray, float]],
        res_subset: list[int] | None = None,
        use_core_fit_transform: tuple[np.ndarray, np.ndarray] | None = None,
    ) -> None:
        if res_subset is None:
            common_res = sorted(set(ca_dict_t.keys()) & set(ca_dict_r.keys()))
        else:
            common_res = sorted(set(res_subset) & set(ca_dict_t.keys()) & set(ca_dict_r.keys()))

        if not common_res:
            return

        p = np.array([ca_dict_t[r][1] for r in common_res])
        q = np.array([ca_dict_r[r][1] for r in common_res])

        if use_core_fit_transform is None:
            res = svd_kabsch_align(coords_target=q, coords_mobile=p)
            aligned_p = res.aligned_coords
            rmsd_val = res.rmsd
        else:
            r_mat, t_vec = use_core_fit_transform
            aligned_p = apply_kabsch_transform(p, r_mat, t_vec)
            diff = aligned_p - q
            rmsd_val = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))

        diffs = aligned_p - q
        displacements = np.sqrt(np.sum(diffs * diffs, axis=1))
        mean_disp = float(np.mean(displacements))
        max_idx = int(np.argmax(displacements))
        max_disp = float(displacements[max_idx])
        max_res_num = common_res[max_idx]
        max_res_name = ca_dict_t[max_res_num][0]
        max_res_str = f"{max_res_name}{max_res_num}"

        rows.append(
            {
                "target_structure": target_name,
                "reference_structure": ref_name,
                "alignment_selection": selection_name,
                "num_aligned_ca_atoms": len(common_res),
                "superposition_method": "Kabsch SVD",
                "rmsd_angstrom": round(rmsd_val, 3),
                "mean_ca_displacement_angstrom": round(mean_disp, 3),
                "max_ca_displacement_angstrom": round(max_disp, 3),
                "max_displacement_residue": max_res_str,
            }
        )

    # 1. 7TRE vs 7TRD comparisons
    # Core fit transformation for TEN displacement
    common_td_core = sorted(set(range(326, 1133)) & set(ca_7tre_tert.keys()) & set(ca_7trd_tert.keys()))
    p_core_td = np.array([ca_7tre_tert[r][1] for r in common_td_core])
    q_core_td = np.array([ca_7trd_tert[r][1] for r in common_td_core])
    res_core_td = svd_kabsch_align(coords_target=q_core_td, coords_mobile=p_core_td)
    td_core_transform = (res_core_td.rotation, res_core_td.translation)

    _eval_alignment("7TRE", "7TRD", "Global", ca_7tre_tert, ca_7trd_tert)
    _eval_alignment(
        "7TRE",
        "7TRD",
        "Catalytic Core",
        ca_7tre_tert,
        ca_7trd_tert,
        res_subset=list(range(326, 1133)),
    )
    _eval_alignment(
        "7TRE",
        "7TRD",
        "TEN Domain",
        ca_7tre_tert,
        ca_7trd_tert,
        res_subset=list(range(1, 197)),
    )
    _eval_alignment(
        "7TRE",
        "7TRD",
        "TEN on Core Fit",
        ca_7tre_tert,
        ca_7trd_tert,
        res_subset=list(range(1, 197)),
        use_core_fit_transform=td_core_transform,
    )

    # 2. 7TRE vs 7TRF comparisons
    common_tf_core = sorted(set(range(326, 1133)) & set(ca_7tre_tert.keys()) & set(ca_7trf_tert.keys()))
    p_core_tf = np.array([ca_7tre_tert[r][1] for r in common_tf_core])
    q_core_tf = np.array([ca_7trf_tert[r][1] for r in common_tf_core])
    res_core_tf = svd_kabsch_align(coords_target=q_core_tf, coords_mobile=p_core_tf)
    tf_core_transform = (res_core_tf.rotation, res_core_tf.translation)

    _eval_alignment("7TRE", "7TRF", "Global", ca_7tre_tert, ca_7trf_tert)
    _eval_alignment(
        "7TRE",
        "7TRF",
        "Catalytic Core",
        ca_7tre_tert,
        ca_7trf_tert,
        res_subset=list(range(326, 1133)),
    )
    _eval_alignment(
        "7TRE",
        "7TRF",
        "TEN Domain",
        ca_7tre_tert,
        ca_7trf_tert,
        res_subset=list(range(1, 197)),
    )
    _eval_alignment(
        "7TRE",
        "7TRF",
        "TEN on Core Fit",
        ca_7tre_tert,
        ca_7trf_tert,
        res_subset=list(range(1, 197)),
        use_core_fit_transform=tf_core_transform,
    )

    # 3. 7TRE vs 2I46 Chain A TPP1 comparisons
    _eval_alignment("7TRE", "2I46", "TPP1 OB-fold", ca_7tre_tpp1, ca_2i46_tpp1)
    _eval_alignment(
        "7TRE",
        "2I46",
        "TEL-patch Core",
        ca_7tre_tpp1,
        ca_2i46_tpp1,
        res_subset=list(TEL_PATCH_RESIDUES),
    )
    _eval_alignment(
        "7TRE",
        "2I46",
        "TEL-patch Extended",
        ca_7tre_tpp1,
        ca_2i46_tpp1,
        res_subset=list(TEL_PATCH_EXTENDED),
    )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Domain-Specific Structural Rearrangements
# -----------------------------------------------------------------------------


def compute_domain_rmsds(
    cif_ref_or_raw_dir: str | Path | None = None,
    cif_mobile: str | Path | None = None,
    domain_defs: dict[str, tuple[int, int]] | None = None,
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Compute domain-specific RMSDs independently and under catalytic core fit.

    Supports two calling conventions:
      1. Single pair: compute_domain_rmsds(cif_ref, cif_mobile, domain_defs)
      2. Directory: compute_domain_rmsds(raw_dir=...) or compute_domain_rmsds()
         (evaluates 7TRE vs 7TRD and 7TRE vs 7TRF across all canonical domains).

    Domain Definitions:
      - TEN: residues 1-196 (modeled 7-178, n=160)
      - TRBD: residues 326-540 (n=188)
      - RT: residues 541-935 (n=385)
      - IFD_TRAP: residues 755-785 (n=31)
      - CTE: residues 936-1132 (n=197)
      - Catalytic Core: residues 326-1132 (n=770)

    Returns:
        pd.DataFrame with columns:
            target_structure, reference_structure, domain_name, residue_range_auth,
            num_ca_atoms, independent_rmsd_angstrom, core_fit_displacement_angstrom,
            mean_displacement_angstrom, max_displacement_angstrom, max_displacement_residue
    """
    if domain_defs is None:
        domain_defs = CANONICAL_DOMAIN_RANGES

    # Determine whether a pair or directory was passed
    is_pair = False
    if cif_mobile is not None:
        is_pair = True
    elif cif_ref_or_raw_dir is not None:
        cand = Path(cif_ref_or_raw_dir)
        if cand.is_file() and cif_mobile is not None:
            is_pair = True

    pairs_to_evaluate: list[tuple[str, Path, str, Path]] = []

    if is_pair:
        assert cif_ref_or_raw_dir is not None
        assert cif_mobile is not None
        p_ref = resolve_cif_path(cif_ref_or_raw_dir)
        p_mob = resolve_cif_path(cif_mobile)
        pairs_to_evaluate.append((p_ref.stem.upper(), p_ref, p_mob.stem.upper(), p_mob))
    else:
        actual_raw = cif_ref_or_raw_dir if (cif_ref_or_raw_dir and Path(cif_ref_or_raw_dir).is_dir()) else raw_dir
        p_tre = resolve_cif_path("7TRE", raw_dir=actual_raw)
        p_trd = resolve_cif_path("7TRD", raw_dir=actual_raw)
        p_trf = resolve_cif_path("7TRF", raw_dir=actual_raw)
        pairs_to_evaluate.append(("7TRE", p_tre, "7TRD", p_trd))
        pairs_to_evaluate.append(("7TRE", p_tre, "7TRF", p_trf))

    rows: list[dict[str, Any]] = []

    for target_name, path_target, ref_name, path_ref in pairs_to_evaluate:
        ca_target = extract_ca_coordinates(path_target, "A", use_auth=True)
        ca_ref = extract_ca_coordinates(path_ref, "A", use_auth=True)

        common_all = sorted(set(ca_target.keys()) & set(ca_ref.keys()))

        # Compute catalytic core fit (residues 326-1132)
        core_res = [r for r in common_all if 326 <= r <= 1132]
        p_core = np.array([ca_target[r][1] for r in core_res])
        q_core = np.array([ca_ref[r][1] for r in core_res])
        res_core = svd_kabsch_align(coords_target=q_core, coords_mobile=p_core)

        for dname, (start, end) in domain_defs.items():
            dom_res = [r for r in common_all if start <= r <= end]
            if not dom_res:
                continue

            p_dom = np.array([ca_target[r][1] for r in dom_res])
            q_dom = np.array([ca_ref[r][1] for r in dom_res])

            # Independent Kabsch fit
            indep_res = svd_kabsch_align(coords_target=q_dom, coords_mobile=p_dom)

            # Core-fit displacement
            p_fitted = apply_kabsch_transform(p_dom, res_core.rotation, res_core.translation)
            diffs = p_fitted - q_dom
            disps = np.sqrt(np.sum(diffs * diffs, axis=1))

            core_fit_rmsd = float(np.sqrt(np.mean(disps * disps)))
            mean_disp = float(np.mean(disps))
            max_idx = int(np.argmax(disps))
            max_disp = float(disps[max_idx])
            max_res_num = dom_res[max_idx]
            max_res_str = f"{ca_target[max_res_num][0]}{max_res_num}"

            range_str = f"{min(dom_res)}-{max(dom_res)}"

            rows.append(
                {
                    "target_structure": target_name,
                    "reference_structure": ref_name,
                    "domain_name": dname,
                    "residue_range_auth": range_str,
                    "num_ca_atoms": len(dom_res),
                    "independent_rmsd_angstrom": round(indep_res.rmsd, 4),
                    "core_fit_displacement_angstrom": round(core_fit_rmsd, 4),
                    "mean_displacement_angstrom": round(mean_disp, 4),
                    "max_displacement_angstrom": round(max_disp, 4),
                    "max_displacement_residue": max_res_str,
                }
            )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Per-Residue C-alpha Displacement Profiling
# -----------------------------------------------------------------------------


def compute_per_residue_displacement(
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Compute per-residue C-alpha Euclidean displacements across 7TRE, 7TRD, and 7TRF.

    Calculates:
      - disp_7tre_vs_7trd_global_angstrom: displacement under global Kabsch fit.
      - disp_7tre_vs_7trd_core_fit_angstrom: displacement under catalytic core (326-1132) fit.
      - disp_7tre_vs_7trf_core_fit_angstrom: displacement under catalytic core (326-1132) fit.

    Demonstrates that the catalytic core (residues 326-1132) acts as a rigid monolith
    (<0.34 Å) while the TEN domain undergoes a slight rigid excursion (0.434 Å).

    Returns:
        pd.DataFrame adhering to Section 4.10 per_residue_displacement_tert.csv schema:
            auth_seq_id, resname, domain, disp_7tre_vs_7trd_global_angstrom,
            disp_7tre_vs_7trd_core_fit_angstrom, disp_7tre_vs_7trf_core_fit_angstrom
    """
    cif_7tre = resolve_cif_path("7TRE", raw_dir=raw_dir)
    cif_7trd = resolve_cif_path("7TRD", raw_dir=raw_dir)
    cif_7trf = resolve_cif_path("7TRF", raw_dir=raw_dir)

    ca_tre = extract_ca_coordinates(cif_7tre, "A", use_auth=True)
    ca_trd = extract_ca_coordinates(cif_7trd, "A", use_auth=True)
    ca_trf = extract_ca_coordinates(cif_7trf, "A", use_auth=True)

    common_all = sorted(set(ca_tre.keys()) & set(ca_trd.keys()) & set(ca_trf.keys()))
    core_res = [r for r in common_all if 326 <= r <= 1132]

    p_all = np.array([ca_tre[r][1] for r in common_all])
    q_trd_all = np.array([ca_trd[r][1] for r in common_all])
    q_trf_all = np.array([ca_trf[r][1] for r in common_all])

    # 1. Global fit 7TRE onto 7TRD
    res_glob_td = svd_kabsch_align(coords_target=q_trd_all, coords_mobile=p_all)
    disp_glob_td = np.sqrt(np.sum((res_glob_td.aligned_coords - q_trd_all) ** 2, axis=1))

    # 2. Catalytic Core fit 7TRE onto 7TRD
    p_core = np.array([ca_tre[r][1] for r in core_res])
    q_core_td = np.array([ca_trd[r][1] for r in core_res])
    res_core_td = svd_kabsch_align(coords_target=q_core_td, coords_mobile=p_core)
    p_all_core_td = apply_kabsch_transform(p_all, res_core_td.rotation, res_core_td.translation)
    disp_core_td = np.sqrt(np.sum((p_all_core_td - q_trd_all) ** 2, axis=1))

    # 3. Catalytic Core fit 7TRE onto 7TRF
    q_core_tf = np.array([ca_trf[r][1] for r in core_res])
    res_core_tf = svd_kabsch_align(coords_target=q_core_tf, coords_mobile=p_core)
    p_all_core_tf = apply_kabsch_transform(p_all, res_core_tf.rotation, res_core_tf.translation)
    disp_core_tf = np.sqrt(np.sum((p_all_core_tf - q_trf_all) ** 2, axis=1))

    rows: list[dict[str, Any]] = []
    for idx, r in enumerate(common_all):
        rows.append(
            {
                "auth_seq_id": r,
                "resname": ca_tre[r][0],
                "domain": get_tert_domain(r),
                "disp_7tre_vs_7trd_global_angstrom": round(float(disp_glob_td[idx]), 4),
                "disp_7tre_vs_7trd_core_fit_angstrom": round(float(disp_core_td[idx]), 4),
                "disp_7tre_vs_7trf_core_fit_angstrom": round(float(disp_core_tf[idx]), 4),
            }
        )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Bound TPP1 (7TRE) vs Unbound TPP1 (2I46 Chain A) Plasticity
# -----------------------------------------------------------------------------


def compute_tpp1_conformational_shifts(
    raw_dir: str | Path | None = None,
    cif_7tre: str | Path | None = None,
    cif_2i46: str | Path | None = None,
) -> pd.DataFrame:
    """
    Quantify conformational shifts between bound TPP1 (7TRE) and unbound TPP1 (2I46 Chain A).

    Crucial Requirements:
      - Strictly targets Chain A of 2I46 (Chain B has the TEL-patch DE-loop 168-171 disordered).
      - Performs global C-alpha Kabsch superposition of residues 90-241 (n=152, RMSD ~1.966 Å).
      - Quantifies per-residue C-alpha displacements.
      - Proves that Leu183 is a rigid pre-organized anchor (displacement 0.512 Å, 'Rigid Anchor').
      - Proves that the acidic DE-loop (Glu168-Glu171) undergoes plastic remodeling (up to 4.293 Å
        at Lys170, 'Extensive Remodeling').
      - Extracts local B-factors from both structures.
      - Classifies residues:
          * 'Rigid Anchor' (< 1.0 Å)
          * 'Moderate Shift' (1.0 - 2.5 Å)
          * 'Extensive Remodeling' (> 2.5 Å)

    Returns:
        pd.DataFrame adhering to Section 4.11 tpp1_bound_vs_unbound_rmsd.csv schema:
            auth_seq_id, resname, is_tel_patch, ca_distance_7tre_vs_2i46_angstrom,
            local_b_factor_7tre, local_b_factor_2i46, conformational_shift_flag
    """
    p_tre = resolve_cif_path(cif_7tre if cif_7tre else "7TRE", raw_dir=raw_dir)
    p_2i46 = resolve_cif_path(cif_2i46 if cif_2i46 else "2I46", raw_dir=raw_dir)

    # Strictly enforce Chain A on 2I46
    valid_chain_2i46 = select_valid_tpp1_chain(p_2i46, preferred_chain="A", require_tel_patch=True)

    ca_tre = extract_ca_coordinates(p_tre, "C", use_auth=True)
    ca_2i46 = extract_ca_coordinates(p_2i46, valid_chain_2i46, use_auth=True)

    common_res = sorted(set(ca_tre.keys()) & set(ca_2i46.keys()))
    if len(common_res) != 152:
        logger.warning("Expected 152 common C-alpha residues between 7TRE and 2I46, found %d", len(common_res))

    p = np.array([ca_tre[r][1] for r in common_res])
    q = np.array([ca_2i46[r][1] for r in common_res])

    # Global Kabsch superposition
    kabsch_res = svd_kabsch_align(coords_target=q, coords_mobile=p)
    diffs = kabsch_res.aligned_coords - q
    displacements = np.sqrt(np.sum(diffs * diffs, axis=1))

    rows: list[dict[str, Any]] = []

    for idx, r in enumerate(common_res):
        disp = float(displacements[idx])
        resname = ca_tre[r][0]
        b_tre = ca_tre[r][2]
        b_2i46 = ca_2i46[r][2]

        if disp < 1.0:
            flag = "Rigid Anchor"
        elif disp <= 2.5:
            flag = "Moderate Shift"
        else:
            flag = "Extensive Remodeling"

        is_tel = r in TEL_PATCH_EXTENDED

        rows.append(
            {
                "auth_seq_id": r,
                "resname": resname,
                "is_tel_patch": is_tel,
                "ca_distance_7tre_vs_2i46_angstrom": round(disp, 4),
                "local_b_factor_7tre": round(b_tre, 2),
                "local_b_factor_2i46": round(b_2i46, 2),
                "conformational_shift_flag": flag,
            }
        )

    return pd.DataFrame(rows)


# Alias for PROJECT.md interface contract compatibility
compute_tpp1_remodeling = compute_tpp1_conformational_shifts


# -----------------------------------------------------------------------------
# Interface Cutoff Robustness Analysis
# -----------------------------------------------------------------------------


def compute_cutoff_robustness(
    cif_path: str | Path | None = None,
    chain_tert: str = "A",
    chain_tpp1: str = "C",
    cutoffs: list[float] | None = None,
    use_auth: bool = True,
    raw_dir: str | Path | None = None,
) -> pd.DataFrame:
    """
    Evaluate persistence of interface contacts across distance thresholds (3.5, 4.0, 4.5 Å).

    Classifies each contact into:
      - Core Robust (persistence = 3): present at 3.5 Å, 4.0 Å, and 4.5 Å.
      - Intermediate (persistence = 2): present at 4.0 Å and 4.5 Å.
      - Threshold Sensitive (persistence = 1): present only at 4.5 Å.

    Wraps compute_contact_robustness from telotrade.interface_geometry.
    """
    if cif_path is None:
        cif_path = resolve_cif_path("7TRE", raw_dir=raw_dir)
    else:
        cif_path = resolve_cif_path(cif_path, raw_dir=raw_dir)

    return compute_contact_robustness(
        cif_path=cif_path,
        chain_tert=chain_tert,
        chain_tpp1=chain_tpp1,
        cutoffs=cutoffs,
        use_auth=use_auth,
    )


# -----------------------------------------------------------------------------
# Standard Results Table Generators
# -----------------------------------------------------------------------------


def generate_cross_state_tables(
    output_dir: str | Path | None = None,
    raw_dir: str | Path | None = None,
) -> dict[str, Path]:
    """
    Generate all standardized results CSV tables for Milestone 4 and save to output_dir.

    Outputs:
      - cross_state_rmsd.csv (Section 4.9)
      - domain_rmsd.csv and domain_rmsds.csv
      - per_residue_displacement_tert.csv (Section 4.10)
      - tpp1_bound_vs_unbound_rmsd.csv and tpp1_conformational_shifts.csv (Section 4.11)
      - interface_contacts_robustness.csv (Section 4.6)

    Returns:
        dict[str, Path]: Mapping from table name to saved file path.
    """
    if output_dir is None:
        output_dir = find_project_root() / "results" / "tables"
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    generated: dict[str, Path] = {}

    # 1. cross_state_rmsd.csv
    logger.info("Generating cross_state_rmsd.csv...")
    df_cross = compute_cross_state_rmsd(raw_dir=raw_dir)
    p_cross = out_path / "cross_state_rmsd.csv"
    df_cross.to_csv(p_cross, index=False)
    generated["cross_state_rmsd.csv"] = p_cross

    # 2. domain_rmsd.csv and domain_rmsds.csv
    logger.info("Generating domain_rmsd.csv...")
    df_domain = compute_domain_rmsds(raw_dir=raw_dir)
    p_domain = out_path / "domain_rmsd.csv"
    p_domain_alt = out_path / "domain_rmsds.csv"
    df_domain.to_csv(p_domain, index=False)
    df_domain.to_csv(p_domain_alt, index=False)
    generated["domain_rmsd.csv"] = p_domain
    generated["domain_rmsds.csv"] = p_domain_alt

    # 3. per_residue_displacement_tert.csv
    logger.info("Generating per_residue_displacement_tert.csv...")
    df_disp = compute_per_residue_displacement(raw_dir=raw_dir)
    p_disp = out_path / "per_residue_displacement_tert.csv"
    df_disp.to_csv(p_disp, index=False)
    generated["per_residue_displacement_tert.csv"] = p_disp

    # 4. tpp1_bound_vs_unbound_rmsd.csv and tpp1_conformational_shifts.csv
    logger.info("Generating tpp1_bound_vs_unbound_rmsd.csv...")
    df_tpp1 = compute_tpp1_conformational_shifts(raw_dir=raw_dir)
    p_tpp1 = out_path / "tpp1_bound_vs_unbound_rmsd.csv"
    p_tpp1_alt = out_path / "tpp1_conformational_shifts.csv"
    df_tpp1.to_csv(p_tpp1, index=False)
    df_tpp1.to_csv(p_tpp1_alt, index=False)
    generated["tpp1_bound_vs_unbound_rmsd.csv"] = p_tpp1
    generated["tpp1_conformational_shifts.csv"] = p_tpp1_alt

    # 5. interface_contacts_robustness.csv
    logger.info("Generating interface_contacts_robustness.csv...")
    df_rob = compute_cutoff_robustness(raw_dir=raw_dir)
    p_rob = out_path / "interface_contacts_robustness.csv"
    df_rob.to_csv(p_rob, index=False)
    generated["interface_contacts_robustness.csv"] = p_rob

    return generated
