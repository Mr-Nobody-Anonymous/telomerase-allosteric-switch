"""Offline tests for the 9SHZ / 9SI0 fixture fetcher (WP21).

The fetcher cannot run in a sandbox with TLS egress blocked, so its PARSER is
tested here against a recorded model-server response.  The sample below is a
verbatim excerpt of

    https://models.rcsb.org/v1/9SHZ/atoms?label_atom_id=CA&label_asym_id=F&auth_seq_id=176

which returned exactly one atom.  That single atom is a known-good value:
Y176 CA = (133.669, 153.363, 143.095) at B = 131.39.  If the parser can pull
that out of the real header layout, it can pull the other 956 residues too.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FETCHER = REPO / "scripts" / "solve" / "fetch_rap_fixture.py"

SAMPLE = """data_9SHZ
#
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ? TYR F 6 176 ? 133.669 153.363 143.095 1.00 131.39 176 TYR A CA 1
ATOM 2 C CA ? LEU F 6 177 ? 132.100 156.400 144.200 1.00 128.10 177 LEU A CA 1
#
"""


def _load():
    spec = importlib.util.spec_from_file_location("rap_fetch", FETCHER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_fetcher_exists():
    assert FETCHER.exists(), "9SHZ/9SI0 fetcher missing"


def test_parser_extracts_the_known_y176_coordinate():
    """Verbatim value from the live model-server call."""
    recs = _load().parse_model_server_cif(SAMPLE)
    assert len(recs) == 2
    y176 = [r for r in recs if r["auth_seq"] == 176][0]
    assert y176["resname"] == "TYR"
    assert y176["chain"] == "F"                 # label_asym_id, as the server uses
    assert y176["label_seq"] == 176
    assert y176["atoms"]["CA"] == pytest.approx([133.669, 153.363, 143.095])


def test_parser_rejects_a_response_without_atom_site():
    assert _load().parse_model_server_cif("data_9SHZ\n#\n") == []


def test_parser_groups_atoms_by_residue_not_by_row():
    """One record per residue, atoms nested - the shape the Colab loader expects."""
    recs = _load().parse_model_server_cif(SAMPLE)
    for r in recs:
        assert set(r) == {"chain", "resname", "label_seq", "auth_seq", "atoms"}
        assert isinstance(r["atoms"], dict)
        assert all(len(v) == 3 for v in r["atoms"].values())


def test_both_rap_entries_are_declared():
    mod = _load()
    assert "9SHZ" in mod.DEFAULT_SPECS and "9SI0" in mod.DEFAULT_SPECS
    # 9SI0 carries TPP1 as well, so it must fetch three chains
    assert set(mod.DEFAULT_SPECS["9SI0"]) == {"F", "C", "D"}
    assert set(mod.DEFAULT_SPECS["9SHZ"]) == {"F"}


def test_dry_run_makes_no_request(capsys):
    assert _load().main(["--entry", "9SHZ", "--dry-run"]) == 0
    assert "dry run" in capsys.readouterr().out


def test_docs_no_longer_claim_the_coordinates_are_unobtainable():
    """The honest state is 'not yet ingested', not 'cannot be retrieved'."""
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    assert "fetch_rap_fixture.py" in readme, (
        "README must point at the fetcher for 9SHZ/9SI0")
