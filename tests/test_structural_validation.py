"""Guards for the work-package-1 structural validation (results/mechanism).

The report in results/mechanism/structural_validation.md is generated, not hand-written, but
generation alone does not stop the numbers from drifting when a structure file or the script
changes. These tests re-derive the load-bearing measurements independently of the report and
pin the ones the mechanism assessment depends on.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import gemmi
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "mechanism" / "structural_validation.md"
SCRIPT = ROOT / "scripts" / "structural_validation.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("structural_validation", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def sv():
    return _load_module()


@pytest.fixture(scope="module")
def model_7tre():
    st = gemmi.read_structure(str(ROOT / "data" / "raw" / "7TRE.cif"))
    st.setup_entities()
    return st[0]


def test_report_exists_and_states_the_correction():
    assert REPORT.exists()
    text = REPORT.read_text()
    for section in ("## 1. Chain identity", "## 2. Does the recomputed interface",
                    "## 3. What does TPP1 actually touch", "## 4. Which domains hold",
                    "## 5. What secondary structure", "## 6. What differs", "## 7. Can TPP1 reach"):
        assert section in text, section
    # the report must not silently drop the retraction of the earlier helix claim
    assert "strand -> turn -> loop -> helix junction" in text
    assert "one continuous helix" in text


def test_canonical_numbering(sv, model_7tre):
    names = sv.names(model_7tre, "A")
    assert names[762] == "HIS" and names[774] == "ARG"
    assert names[712] == "ASP" and names[868] == "ASP" and names[869] == "ASP"
    seq = sv.chain_sequence(model_7tre, "A")
    assert "".join(seq[i] for i in range(753, 762)) == "GHVRKAFKS"


def test_tpp1_does_not_contact_dna(sv, model_7tre):
    tpp1 = sv.flat(sv.heavy(model_7tre, "C"))
    dna = sv.flat(sv.heavy(model_7tre, "D"))
    rna = sv.flat(sv.heavy(model_7tre, "B"))
    assert sv.min_dist(tpp1, dna) > 20.0
    assert sv.min_dist(tpp1, rna) > 10.0


def test_ten_makes_no_primer_contact(sv, model_7tre):
    tert, dna = sv.heavy(model_7tre, "A"), sv.heavy(model_7tre, "D")
    dna_xyz = sv.flat(dna)
    ten = [r for r, pts in tert.items()
           if r <= 196 and sv.min_dist(dna_xyz, pts) <= sv.PRIMER_CUTOFF]
    assert ten == []
    assert sv.min_dist(dna_xyz, tert[176]) > 8.0


def test_only_766_of_the_conduit_window_is_at_the_interface(sv, model_7tre):
    tert, tpp1 = sv.heavy(model_7tre, "A"), sv.heavy(model_7tre, "C")
    tpp1_xyz = sv.flat(tpp1)
    touched = [r for r in range(760, 767) if sv.min_dist(tpp1_xyz, tert[r]) <= sv.CUTOFF]
    # the claim under test: the proximal end of the proposed conduit is not in contact
    assert touched == [766]
    assert sv.min_dist(tpp1_xyz, tert[762]) > 10.0


def test_anchor_to_footprint_gradient_is_monotone_in_the_core(sv, model_7tre):
    tert, tpp1 = sv.heavy(model_7tre, "A"), sv.heavy(model_7tre, "C")
    tpp1_xyz = sv.flat(tpp1)
    d = {r: sv.min_dist(tpp1_xyz, tert[r]) for r in (756, 759, 762, 765, 766, 771, 774)}
    assert d[756] > d[759] > d[762] > d[765] > d[766] > d[774]
    assert d[774] < 3.0 and d[756] > 20.0


def test_secondary_structure_contradicts_the_single_helix_claim(sv):
    # i -> i+4 distances in 752-763 are strand-like, not helical
    ca = sv.ca_map(sv.load("7TRE")[1], "A")
    d4 = [float(np.linalg.norm(ca[i + 4] - ca[i])) for i in range(752, 764)]
    assert min(d4) > 9.0, d4
    helix = [float(np.linalg.norm(ca[i + 4] - ca[i])) for i in range(772, 781)]
    assert max(helix) < 7.0, helix


def test_7trd_is_not_a_tpp1_free_conformation(sv):
    st_a, model_a = sv.load("7TRE")
    st_b, model_b = sv.load("7TRD")
    ca_a, ca_b = sv.ca_map(model_a, "A"), sv.ca_map(model_b, "A")
    shared = sorted(set(ca_a) & set(ca_b))
    rmsd = sv.kabsch_rmsd(np.array([ca_a[i] for i in shared]),
                          np.array([ca_b[i] for i in shared]))
    assert rmsd < 0.6, rmsd
    assert "N" in {ch.name for ch in model_b}
    assert "N" not in {ch.name for ch in model_a}
