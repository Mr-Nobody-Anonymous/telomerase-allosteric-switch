"""Regression guards for the 2026-10-30 claim-verification pass.

These tests do not re-run the full audit (see ``scripts/verify_all_claims.py``); they pin the
handful of *facts* the audit rests on, so that if the underlying data files ever change the
correction documents stop being silently wrong.

Run: ``PYTHONPATH=src pytest tests/test_claim_verification.py -q``
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import gemmi
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
VERIFY = ROOT / "results" / "verification"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _model(pdb_id: str) -> gemmi.Model:
    st = gemmi.read_structure(str(RAW / f"{pdb_id}.cif"))
    st.setup_entities()
    return st[0]


def _ca_map(model: gemmi.Model, chain: str = "A") -> dict[int, np.ndarray]:
    return {
        r.seqid.num: np.array([a.pos.x, a.pos.y, a.pos.z], float)
        for ch in model
        if ch.name == chain
        for r in ch
        for a in r
        if a.name == "CA"
    }


def test_audit_tables_are_present_and_complete():
    csv_path = VERIFY / "claim_recheck.csv"
    assert csv_path.is_file(), "run scripts/verify_all_claims.py to regenerate the audit"
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) >= 50
    verdicts = {r["verdict"] for r in rows}
    assert verdicts <= {"VERIFIED", "CONTRADICTED", "PARTIALLY_VERIFIED", "UNTESTED", "MINOR_DRIFT", "COMPUTED"}
    ids = [r["claim_id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate claim ids in the audit table"


def test_manifest_hashes_match_the_shipped_structures():
    """data_manifest.csv (not FINAL_RESEARCH_REPORT.md) carries the true hashes."""
    with open(ROOT / "data_manifest.csv") as fh:
        manifest = {r["pdb_id"]: r["sha256"] for r in csv.DictReader(fh)}
    for pdb_id, expected in manifest.items():
        path = RAW / f"{pdb_id}.cif"
        assert path.is_file(), f"missing {path}"
        assert _sha256(path) == expected, f"{pdb_id} hash drifted from the manifest"


def _kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    pc, qc = p - p.mean(0), q - q.mean(0)
    u, _, vt = np.linalg.svd(qc.T @ pc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    qt = qc @ r.T
    return float(np.sqrt(np.mean(np.sum((pc - qt) ** 2, axis=1))))


def test_seven_trd_and_seven_trf_are_the_same_tetr_model():
    """Documented in CLAIM_VERIFICATION_REPORT.md §3: the two 'unbound' entries share one model.

    The deposited coordinates are offset by a small rigid-body translation, so compare after
    superposition: the two TERT models are identical to numerical precision.
    """
    a, b = _ca_map(_model("7TRD")), _ca_map(_model("7TRF"))
    common = sorted(set(a) & set(b))
    assert len(common) == 934
    rmsd = _kabsch_rmsd(np.array([a[r] for r in common]), np.array([b[r] for r in common]))
    assert rmsd < 0.01, f"7TRD/7TRF TERT models are no longer superposable (RMSD {rmsd:.3f} Å)"


def test_lys757a_system_file_is_identical_to_wild_type():
    """Documented in §4.4: no K757A simulation system was ever built."""
    wt = ROOT / "data" / "junction_bound_fixed.pdb"
    mut = ROOT / "data" / "junction_bound_k757a_fixed.pdb"
    assert wt.is_file() and mut.is_file()
    assert _sha256(wt) == _sha256(mut)
    # The junction systems are renumbered 1..46; TERT Lys757 is the 8th residue of chain B.
    model = gemmi.read_structure(str(mut))[0]
    chain_b = next(ch for ch in model if ch.name == "B")
    assert list(chain_b)[7].name == "LYS", "expected an unmutated LYS at the 757-equivalent position"


def test_telin1_is_not_a_single_chemical_entity():
    """Documented in §6: two disconnected fragments, so it cannot be a lead compound."""
    path = ROOT / "data" / "ligands" / "telin1_inhibitor.pdb"
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if el == "H":
                continue
            atoms.append((el, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    assert len(atoms) == 36

    cov = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "F": 0.57}
    parent = list(range(len(atoms)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = float(np.linalg.norm(atoms[i][1] - atoms[j][1]))
            if d < cov.get(atoms[i][0], 0.8) + cov.get(atoms[j][0], 0.8) + 0.45:
                parent[find(i)] = find(j)
    assert len({find(i) for i in range(len(atoms))}) == 2


def test_interface_anchors_reproduce_from_raw_coordinates():
    """Spot-check of the numbers that survived verification (report §3)."""
    model = _model("7TRE")
    atoms = {
        (r.seqid.num, a.name): np.array([a.pos.x, a.pos.y, a.pos.z], float)
        for ch in model
        if ch.name in {"A", "C"}
        for r in ch
        for a in r
        if a.element.name != "H"
    }
    d_arg = np.linalg.norm(atoms[(774, "NH1")] - atoms[(168, "OE1")])
    d_lys = np.linalg.norm(atoms[(78, "NZ")] - atoms[(215, "OE2")])
    assert abs(d_arg - 2.320) < 0.01
    assert abs(d_lys - 3.584) < 0.01


@pytest.mark.parametrize("pdb_id,expected", [("2I46", 2.70), ("7TRE", 3.50), ("7TRD", 3.30), ("7TRF", 3.70)])
def test_deposited_resolutions(pdb_id, expected):
    st = gemmi.read_structure(str(RAW / f"{pdb_id}.cif"))
    assert abs(st.resolution - expected) < 0.01
