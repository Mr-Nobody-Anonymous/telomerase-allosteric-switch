"""Guards for the mechanical-coupling analysis (results/allostery).

These pin the geometric facts the analysis report rests on, so the reframed claim cannot
silently drift if the structures or the script change.
"""

from __future__ import annotations

import json
from pathlib import Path

import gemmi
import numpy as np
import pytest
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "allostery"


def _model_7tre():
    st = gemmi.read_structure(str(ROOT / "data" / "raw" / "7TRE.cif"))
    st.setup_entities()
    return st[0]


def _coords(model, chain):
    out = {}
    for ch in model:
        if ch.name != chain:
            continue
        for r in ch:
            pts = [np.array([a.pos.x, a.pos.y, a.pos.z]) for a in r if a.element.name != "H"]
            if pts:
                out[r.seqid.num] = np.array(pts)
    return out


def test_distance_gradient_is_monotonic_along_the_anchor_helix():
    """TPP1 contact distance decays continuously from Lys757 to Arg774 (the 'no void' claim)."""
    model = _model_7tre()
    tpp1 = np.vstack(list(_coords(model, "C").values()))
    tert = _coords(model, "A")
    tree = cKDTree(tpp1)
    d = {rn: float(tree.query(tert[rn])[0].min()) for rn in (757, 759, 761, 765, 766, 774)}
    assert d[757] > 20.0                     # distal end of the helix, far from TPP1
    assert d[766] < 5.0                      # TPP1 footprint edge
    assert d[774] < 3.0                      # docking hotspot
    assert d[757] > d[759] > d[761] > d[765] > d[766] > d[774], d


def test_tpp1_never_touches_the_nucleic_acid():
    """The report's 28.7 A figure: TPP1 cannot act on DNA by direct contact in 7TRE."""
    model = _model_7tre()
    tpp1 = np.vstack(list(_coords(model, "C").values()))
    dna = np.vstack(list(_coords(model, "D").values()))
    rna = np.vstack(list(_coords(model, "B").values()))
    assert cKDTree(dna).query(tpp1)[0].min() > 25.0
    assert cKDTree(rna).query(tpp1)[0].min() > 12.0


def test_ten_domain_has_no_primer_contact():
    """The 'two-hand clamp' is absent: only the RT/IFD-TRAP side holds the modelled DNA."""
    model = _model_7tre()
    tert, dna = _coords(model, "A"), _coords(model, "D")
    tree = cKDTree(np.vstack(list(dna.values())))
    ten = [rn for rn in tert if rn <= 196]
    assert ten, "TEN domain residues missing"
    assert min(float(tree.query(tert[rn])[0].min()) for rn in ten) > 8.0


def test_analysis_outputs_exist_and_are_consistent():
    for name in ("anm_summary.json", "system_comparison.csv", "cutoff_robustness.csv",
                 "tpp1_distance_gradient.csv", "residue_coupling_ranked.csv"):
        assert (OUT / name).is_file(), f"missing {name}; run scripts/allosteric_coupling.py"
    summary = json.loads((OUT / "anm_summary.json").read_text())
    assert summary["H_switch_verdict"] in {"SUPPORTED", "NOT SUPPORTED"}
    assert summary["H_clamp_verdict"] in {"SUPPORTED", "NOT SUPPORTED", "MARGINAL"}
    # the positive control must be positive, otherwise the negative result is meaningless
    assert summary["positive_control_short_range"] > 0.01
    assert len(summary["table"]) == 3


def test_negative_result_is_not_overclaimed():
    """The report explicitly states the model cannot adjudicate; guard the wording."""
    text = (OUT / "ALLOSTERIC_COUPLING_ANALYSIS.md").read_text()
    assert "cannot adjudicate" in text
    assert "no coupling beyond" in text.lower() or "no coupling beyond proximity" in text.lower()
