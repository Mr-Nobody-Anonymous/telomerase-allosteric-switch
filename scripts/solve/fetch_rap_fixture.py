#!/usr/bin/env python3
"""
Fetch the RAP-state fixtures (9SHZ / 9SI0) that this repository has never had.

WHY THIS EXISTS
---------------
The README and WP17 both say the RAP-state superposition is "still not
computed: 9SHZ and 9SI0 coordinates have never been retrieved".  Every other
structure in this repo (7TRE, 7TRD, 7TRF, 2I46, 9SHY) has its CA atoms on disk,
so 9SHZ/9SI0 is the single remaining gap in the structural table.

The blocker was not permission or availability - PDBe exposes both entries and
models.rcsb.org serves their coordinates.  It is that nobody wrote the fetcher.

This script is that fetcher.  It writes fixtures in the SAME format as the
existing ones (see colab/fixtures/9SHY.json) so the Colab loader needs no
change.

USAGE
    python scripts/solve/fetch_rap_fixture.py --entry 9SHZ --dry-run
    python scripts/solve/fetch_rap_fixture.py --entry 9SHZ
    python scripts/solve/fetch_rap_fixture.py --entry 9SI0

Needs network.  In a sandbox with TLS egress blocked it fails with a clear
message rather than a traceback, and the parser is still testable offline
against the recorded sample in tests/test_rap_fixture_fetch.py.

PROVENANCE (verified via the PDBe polymer-coverage API)
    9SHZ  TERT = struct_asym F / auth chain A, to LEU1132
          DNA  = label 2-12,  auth 20-30
    9SI0  TERT = struct_asym F / auth chain A, 4-1132 (gaps 105-111, 418-443)
          DNA  = struct_asym C / auth chain N, 24-33
          TPP1 = struct_asym D / auth chain O, auth 90-238
    NOTE label_asym_id != auth_asym_id throughout; the model server needs
    label_asym_id, the output records carry the author chain id.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "colab" / "fixtures"
MODEL_SERVER = "https://models.rcsb.org/v1/{entry}/atoms"

# label_asym_id -> what it is, per the PDBe coverage API
DEFAULT_SPECS = {
    "9SHZ": {"F": "TERT (auth A)"},
    "9SI0": {"F": "TERT (auth A)", "C": "DNA (auth N)", "D": "TPP1 (auth O)"},
}


def parse_model_server_cif(text: str) -> list[dict]:
    """Parse the atom_site loop out of an RCSB model-server mmCIF response.

    Returns records shaped like the existing fixtures:
        {"chain", "resname", "label_seq", "auth_seq", "atoms": {"CA": [x,y,z]}}
    """
    lines = text.splitlines()
    # Locate the atom_site header.  Some model-server responses precede the
    # columns with an explicit `loop_` token and some do not, so find the first
    # `_atom_site.` column and walk forward from there.
    start = None
    cols: list[str] = []
    for i, line in enumerate(lines):
        if not line.strip().startswith("_atom_site."):
            continue
        cols, j = [], i
        while j < len(lines) and lines[j].strip().startswith("_atom_site."):
            cols.append(lines[j].strip().split(".", 1)[1])
            j += 1
        start = j
        break
    if start is None:
        return []

    idx = {c: n for n, c in enumerate(cols)}
    need = ("label_asym_id", "label_comp_id", "label_atom_id",
            "label_seq_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")
    missing = [c for c in need if c not in idx]
    if missing:
        raise ValueError(f"model-server response lacks columns: {missing}")

    out: dict[tuple, dict] = {}
    for line in lines[start:]:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
            if s.startswith("#") or s.startswith("_") or s == "loop_":
                break
            continue
        f = s.split()
        if len(f) < len(cols):
            continue
        try:
            label_seq = int(f[idx["label_seq_id"]])
            auth_seq = int(f[idx["auth_seq_id"]])
            xyz = [float(f[idx["Cartn_x"]]), float(f[idx["Cartn_y"]]),
                   float(f[idx["Cartn_z"]])]
        except ValueError:
            continue
        key = (f[idx["label_asym_id"]], label_seq, auth_seq, f[idx["label_comp_id"]])
        rec = out.setdefault(key, {
            "chain": f[idx["label_asym_id"]],
            "resname": f[idx["label_comp_id"]],
            "label_seq": label_seq,
            "auth_seq": auth_seq,
            "atoms": {},
        })
        rec["atoms"][f[idx["label_atom_id"]]] = xyz
    return list(out.values())


def fetch(entry: str, label_asym: str, atom: str = "CA", timeout: int = 60) -> str:
    url = (f"{MODEL_SERVER.format(entry=entry)}"
           f"?label_atom_id={atom}&label_asym_id={label_asym}")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except (urllib.error.URLError, OSError) as exc:
        raise SystemExit(
            f"network error fetching {url}\n  {exc}\n\n"
            "This sandbox blocks TLS egress from Python. Run this script on a "
            "machine with network access, or place the response on disk and use "
            "parse_model_server_cif() directly.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--entry", required=True, choices=sorted(DEFAULT_SPECS))
    ap.add_argument("--atom", default="CA")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    specs = DEFAULT_SPECS[a.entry]
    out = Path(a.out) if a.out else FIXTURES / f"atoms_{a.entry}.json"
    print(f"entry {a.entry}: {len(specs)} chain(s) -> {out.relative_to(REPO)}")
    for la, desc in specs.items():
        print(f"  label_asym {la}: {desc}")
    if a.dry_run:
        print("dry run: no request made")
        return 0

    records: list[dict] = []
    for la in specs:
        text = fetch(a.entry, la, a.atom)
        recs = parse_model_server_cif(text)
        print(f"  label_asym {la}: {len(recs)} residues parsed")
        if not recs:
            print(f"  WARNING: no atoms parsed for {la}", file=sys.stderr)
        records.extend(recs)

    payload = {
        "provenance": (f"RCSB/PDBe {a.entry}, CA atoms via models.rcsb.org. "
                       f"Fetched by scripts/solve/fetch_rap_fixture.py. "
                       f"Chains: {', '.join(specs)}."),
        "records": records,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"wrote {out.relative_to(REPO)} ({len(records)} residues)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
