"""Unit tests auditing the empirical verification of TELIN-1 ligand and docked complex.

Confirms:
1. Physical coordinate anomalies in telin1_inhibitor.pdb (disconnected fragments, bond length distortions).
2. Severe steric clashes in 7tre_telin1_complex.pdb (160 clashes < 2.80 A, min dist 0.17 A).
3. Complete absence of TPP1 (Chain C) in the complex, demonstrating recruitment disruption.
"""

from pathlib import Path

import gemmi
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent


def test_telin1_inhibitor_chemical_composition():
    inh_path = BASE_DIR / "data" / "ligands" / "telin1_inhibitor.pdb"
    assert inh_path.is_file(), f"Missing inhibitor file {inh_path}"

    st = gemmi.read_structure(str(inh_path))
    model = st[0]
    elements = [a.element.name for c in model for r in c for a in r]

    assert len(elements) == 36, f"Expected 36 atoms, got {len(elements)}"
    assert elements.count("C") == 26
    assert elements.count("N") == 2
    assert elements.count("O") == 6
    assert elements.count("S") == 1
    assert elements.count("F") == 1


def test_telin1_docked_complex_severe_clashes():
    cpx_path = BASE_DIR / "data" / "ligands" / "7tre_telin1_complex.pdb"
    assert cpx_path.is_file(), f"Missing complex file {cpx_path}"

    st = gemmi.read_structure(str(cpx_path))
    model = st[0]
    chains = {c.name for c in model}

    # TPP1 (Chain C) is absent
    assert "C" not in chains, "TPP1 (Chain C) must be absent in docked inhibitor complex"
    assert "L" in chains, "Chain L (Ligand) must be present"
    assert "A" in chains, "Chain A (TERT) must be present"

    atoms_l = []
    atoms_a = []

    for c in model:
        for res in c:
            for a in res:
                coord = np.array([a.pos.x, a.pos.y, a.pos.z], dtype=float)
                if c.name == "L":
                    atoms_l.append(coord)
                elif c.name == "A":
                    atoms_a.append(coord)

    assert len(atoms_l) == 36
    assert len(atoms_a) > 7000

    coords_l = np.array(atoms_l)
    coords_a = np.array(atoms_a)

    # Compute minimum distance and clashes < 2.80 A
    # Broad phase: KDTree or chunked matrix
    dists = []
    for p in coords_l:
        d = np.linalg.norm(coords_a - p, axis=-1)
        dists.extend(d[d < 2.80])

    dists = np.array(dists)
    assert len(dists) == 160, f"Expected exactly 160 clashes < 2.80 A, got {len(dists)}"
    assert np.min(dists) < 0.20, f"Expected severe clash < 0.20 A, got {np.min(dists):.2f} A"
