"""Adversarial challenger test suite 2 for Milestone B2.

Challenger: challenger_b2_2 (Milestone B2 Adversarial Challenger 2)
Focus:
- Verification of results/interface/sensitivity_analysis.csv:
  - 32 residue pairs partitioned into Core Persistent (5), Intermediate Boundary (16), and Peripheral/Sensitive (11).
  - Exact concordance with raw deposited coordinates in 7TRE (PDB: 7TRE).
  - Monotonicity, boolean presence flags, and heavy-atom contact breakdowns (9 -> 46 -> 123).
  - Positive controls and key interatomic distance benchmarks (Arg774-Glu168, Lys78-Glu215, Leu769-Leu183, Pro771-Leu212, Gln775-Asp166).
- Verification of results/research/question_and_hypotheses.md:
  - Competing hypotheses H1, H2, and H3.
  - Non-circular, falsifiable support and rejection criteria mapping directly to interface metrics.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import MMCIFParser
from scipy.spatial import cKDTree

RAW_DIR = Path("data/raw")
SENSITIVITY_CSV = Path("results/interface/sensitivity_analysis.csv")
HYPOTHESES_MD = Path("results/research/question_and_hypotheses.md")


@pytest.fixture(scope="module")
def sensitivity_df():
    assert SENSITIVITY_CSV.is_file(), f"Missing {SENSITIVITY_CSV}"
    df = pd.read_csv(SENSITIVITY_CSV)
    return df


@pytest.fixture(scope="module")
def raw_7tre_heavy_atoms():
    cif_path = RAW_DIR / "7tre.cif"
    assert cif_path.is_file(), f"Missing {cif_path}"
    parser = MMCIFParser(auth_chains=True, auth_residues=True, QUIET=True)
    structure = parser.get_structure("7TRE", str(cif_path))
    model = structure[0]

    tert_atoms = []  # (resnum, resname, atomname, coord)
    for res in model["A"].get_residues():
        hetflag, resseq, _ = res.get_id()
        if hetflag.strip():
            continue
        resname = res.get_resname().strip()
        for atom in res.get_atoms():
            el = atom.element.strip().upper() if atom.element else atom.name[0].upper()
            if el == "H":
                continue
            tert_atoms.append((resseq, resname, atom.name.strip(), atom.coord))

    tpp1_atoms = []
    for res in model["C"].get_residues():
        hetflag, resseq, _ = res.get_id()
        if hetflag.strip():
            continue
        resname = res.get_resname().strip()
        for atom in res.get_atoms():
            el = atom.element.strip().upper() if atom.element else atom.name[0].upper()
            if el == "H":
                continue
            tpp1_atoms.append((resseq, resname, atom.name.strip(), atom.coord))

    return tert_atoms, tpp1_atoms


def test_sensitivity_analysis_residue_pair_partitioning(sensitivity_df):
    """Verify that sensitivity_analysis.csv has 32 rows partitioned into 5 Core, 16 Intermediate, 11 Peripheral."""
    assert len(sensitivity_df) == 32, f"Expected 32 rows, got {len(sensitivity_df)}"

    core = sensitivity_df[sensitivity_df["persistence_class"] == "Core Persistent"]
    inter = sensitivity_df[sensitivity_df["persistence_class"] == "Intermediate Boundary"]
    periph = sensitivity_df[sensitivity_df["persistence_class"] == "Peripheral / Sensitive"]

    assert len(core) == 5, f"Expected 5 Core Persistent pairs, got {len(core)}"
    assert len(inter) == 16, f"Expected 16 Intermediate Boundary pairs, got {len(inter)}"
    assert len(periph) == 11, f"Expected 11 Peripheral / Sensitive pairs, got {len(periph)}"

    # Check Core properties
    for _, row in core.iterrows():
        d = float(row["min_distance_angstrom"])
        assert d <= 3.50, f"Core pair {row['contact_pair']} distance {d} > 3.50 A"
        assert row["present_at_3_5A"] is True or row["present_at_3_5A"] == "True"
        assert row["present_at_4_0A"] is True or row["present_at_4_0A"] == "True"
        assert row["present_at_4_5A"] is True or row["present_at_4_5A"] == "True"

    # Check Intermediate properties
    for _, row in inter.iterrows():
        d = float(row["min_distance_angstrom"])
        assert 3.50 < d <= 4.00, f"Intermediate pair {row['contact_pair']} distance {d} not in (3.5, 4.0]"
        assert row["present_at_3_5A"] is False or row["present_at_3_5A"] == "False"
        assert row["present_at_4_0A"] is True or row["present_at_4_0A"] == "True"
        assert row["present_at_4_5A"] is True or row["present_at_4_5A"] == "True"

    # Check Peripheral properties
    for _, row in periph.iterrows():
        d = float(row["min_distance_angstrom"])
        assert 4.00 < d <= 4.50, f"Peripheral pair {row['contact_pair']} distance {d} not in (4.0, 4.5]"
        assert row["present_at_3_5A"] is False or row["present_at_3_5A"] == "False"
        assert row["present_at_4_0A"] is False or row["present_at_4_0A"] == "False"
        assert row["present_at_4_5A"] is True or row["present_at_4_5A"] == "True"


def test_sensitivity_analysis_matches_raw_7tre_structure(sensitivity_df, raw_7tre_heavy_atoms):
    """Stress-test sensitivity table directly against raw 7tre.cif coordinates."""
    tert_atoms, tpp1_atoms = raw_7tre_heavy_atoms
    coords_a = np.array([x[3] for x in tert_atoms], dtype=np.float64)
    coords_c = np.array([x[3] for x in tpp1_atoms], dtype=np.float64)

    tree_a = cKDTree(coords_a)
    tree_c = cKDTree(coords_c)

    # Query all heavy atom contacts at 4.5 A
    pairs = tree_a.query_ball_tree(tree_c, r=4.5)

    emp_pairs = {}
    for i, js in enumerate(pairs):
        t_seq = tert_atoms[i][0]
        for j in js:
            p_seq = tpp1_atoms[j][0]
            dist = float(np.linalg.norm(coords_a[i] - coords_c[j]))
            key = (t_seq, p_seq)
            if key not in emp_pairs or dist < emp_pairs[key]:
                emp_pairs[key] = dist

    assert len(emp_pairs) == 32, f"Expected 32 unique empirical residue pairs, got {len(emp_pairs)}"

    # Check every row in sensitivity_df
    for _, row in sensitivity_df.iterrows():
        t_seq = int(row["resnum_tert"])
        p_seq = int(row["resnum_tpp1"])
        key = (t_seq, p_seq)
        assert key in emp_pairs, f"Pair {row['contact_pair']} not found in empirical 7TRE contact network"
        emp_dist = emp_pairs[key]
        csv_dist = float(row["min_distance_angstrom"])
        assert abs(emp_dist - csv_dist) <= 0.02, (
            f"Distance mismatch for {row['contact_pair']}: CSV={csv_dist}, Emp={emp_dist}"
        )


def test_heavy_atom_contacts_scaling_and_subsets(sensitivity_df, raw_7tre_heavy_atoms):
    """Verify that contact counts are exactly 9 (3.5 A), 46 (4.0 A), and 123 (4.5 A).

    Verify the distribution of contacts among persistence classes.
    """
    tert_atoms, tpp1_atoms = raw_7tre_heavy_atoms
    coords_a = np.array([x[3] for x in tert_atoms], dtype=np.float64)
    coords_c = np.array([x[3] for x in tpp1_atoms], dtype=np.float64)

    tree_a = cKDTree(coords_a)
    tree_c = cKDTree(coords_c)

    cnt_35 = tree_a.count_neighbors(tree_c, r=3.5)
    cnt_40 = tree_a.count_neighbors(tree_c, r=4.0)
    cnt_45 = tree_a.count_neighbors(tree_c, r=4.5)

    assert cnt_35 == 9, f"Expected 9 contacts at 3.5 A, got {cnt_35}"
    assert cnt_40 == 46, f"Expected 46 contacts at 4.0 A, got {cnt_40}"
    assert cnt_45 == 123, f"Expected 123 contacts at 4.5 A, got {cnt_45}"

    # Verify class contributions
    class_map = {
        (int(r["resnum_tert"]), int(r["resnum_tpp1"])): r["persistence_class"]
        for _, r in sensitivity_df.iterrows()
    }

    pairs_45 = tree_a.query_ball_tree(tree_c, r=4.5)
    counts_40 = {"Core Persistent": 0, "Intermediate Boundary": 0, "Peripheral / Sensitive": 0}
    counts_45 = {"Core Persistent": 0, "Intermediate Boundary": 0, "Peripheral / Sensitive": 0}

    for i, js in enumerate(pairs_45):
        t_seq = tert_atoms[i][0]
        for j in js:
            p_seq = tpp1_atoms[j][0]
            dist = np.linalg.norm(coords_a[i] - coords_c[j])
            p_cls = class_map[(t_seq, p_seq)]
            if dist <= 4.0:
                counts_40[p_cls] += 1
            if dist <= 4.5:
                counts_45[p_cls] += 1

    assert counts_40["Core Persistent"] == 18
    assert counts_40["Intermediate Boundary"] == 28
    assert counts_40["Peripheral / Sensitive"] == 0

    assert counts_45["Core Persistent"] == 30
    assert counts_45["Intermediate Boundary"] == 70
    assert counts_45["Peripheral / Sensitive"] == 23


def test_key_anchor_and_latch_interatomic_distances(raw_7tre_heavy_atoms):
    """Stress-test exact Euclidean distances of primary structural anchors in 7TRE."""
    tert_atoms, tpp1_atoms = raw_7tre_heavy_atoms

    def get_atom_coord(atoms, resnum, atomname):
        for seq, _, name, coord in atoms:
            if seq == resnum and name == atomname:
                return coord
        raise KeyError(f"Atom {atomname} of residue {resnum} not found")

    # 1. Primary Salt Bridge: Arg774:NH1 -- Glu168:OE1 == 2.32 A
    c_nh1 = get_atom_coord(tert_atoms, 774, "NH1")
    c_oe1 = get_atom_coord(tpp1_atoms, 168, "OE1")
    d_salt = float(np.linalg.norm(c_nh1 - c_oe1))
    assert round(d_salt, 2) == 2.32, f"Expected 2.32 A, got {d_salt}"

    # 2. Secondary Salt Latch: Lys78:NZ -- Glu215:OE2 == 3.58 A
    c_nz = get_atom_coord(tert_atoms, 78, "NZ")
    c_oe2 = get_atom_coord(tpp1_atoms, 215, "OE2")
    d_latch = float(np.linalg.norm(c_nz - c_oe2))
    assert round(d_latch, 2) == 3.58, f"Expected 3.58 A, got {d_latch}"

    # 3. Hydrophobic Anchor: Leu769:CD2 -- Leu183:CD2 == 3.67 A
    c_cd2_769 = get_atom_coord(tert_atoms, 769, "CD2")
    c_cd2_183 = get_atom_coord(tpp1_atoms, 183, "CD2")
    d_hydro = float(np.linalg.norm(c_cd2_769 - c_cd2_183))
    assert round(d_hydro, 2) == 3.67, f"Expected 3.67 A, got {d_hydro}"

    # 4. Hydrophobic Groove: Pro771:CG -- Leu212:CD2 == 3.52 A
    c_cg = get_atom_coord(tert_atoms, 771, "CG")
    c_cd2_212 = get_atom_coord(tpp1_atoms, 212, "CD2")
    d_groove = float(np.linalg.norm(c_cg - c_cd2_212))
    assert round(d_groove, 2) == 3.52, f"Expected 3.52 A, got {d_groove}"

    # 5. Peripheral Contact: Gln775:NE2 -- Asp166:O == 4.41 A
    c_ne2 = get_atom_coord(tert_atoms, 775, "NE2")
    c_o = get_atom_coord(tpp1_atoms, 166, "O")
    d_periph = float(np.linalg.norm(c_ne2 - c_o))
    assert round(d_periph, 2) == 4.41, f"Expected 4.41 A, got {d_periph}"


def test_question_and_hypotheses_falsifiability_and_non_circularity():
    """Verify results/research/question_and_hypotheses.md structure, non-circularity, and falsifiability."""
    assert HYPOTHESES_MD.is_file(), f"Missing {HYPOTHESES_MD}"
    content = HYPOTHESES_MD.read_text(encoding="utf-8")

    # Verify hypotheses exist
    assert "Hypothesis 1 (H1)" in content
    assert "Hypothesis 2 (H2)" in content
    assert "Hypothesis 3 (H3)" in content

    # Verify explicit support and falsification criteria headers
    assert "Quantitative Support Criteria" in content
    assert "Quantitative Falsification Criteria" in content

    # Verify mapped quantitative metrics
    metrics = [
        "0.338",  # Catalytic core RMSD
        "0.434",  # TEN excursion
        "0.427",  # TRAP excursion
        "4.293",  # Lys170 peak shift
        "0.512",  # Leu183 anchor
        "856.82",  # FreeSASA Bondi
        "ProtOr",  # FreeSASA ProtOr convention
        "Bondi",  # Bondi vdW convention
        "123",  # 4.5 A contacts
        "46",  # 4.0 A contacts
        "9",  # 3.5 A contacts
    ]
    for m in metrics:
        assert m in content, f"Metric/convention '{m}' not cited in {HYPOTHESES_MD}"

    # Verify key falsification modalities
    assert "HDX" in content or "HDX-MS" in content
    assert "NMR" in content
    assert "stopped-flow" in content or "Stopped-flow" in content
    assert "AIC" in content
    assert "Cruickshank" in content or "DPI" in content

    # Verify comparative decision matrix and roadmap
    assert "Multi-Hypothesis Comparative Decision Matrix" in content
    assert "Multi-Phase Empirical Discrimination Roadmap" in content
