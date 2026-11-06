#!/usr/bin/env python3
"""
TPP1 / telomerase structural test - a self-contained Colab runner.

WHY THIS FILE EXISTS
--------------------
Everything in scripts/solve/wp1-wp15 was computed from five mmCIF files that were already on
disk. The one measurement that could never be made is the one that matters most:

    does the TEN-domain anchor motif PLYQ (174-177) actually contact the DNA
    once the substrate is long enough to reach it?

7TRE / 7TRD / 7TRF all carry only a 7-nucleotide DNA, which stops two nucleotides short of
PLYQ. Extrapolating the 5' exit direction put Tyr176 at 6.8 A, but that was a straight line
through one nucleotide step - an inference, not an observation.

7QXS (30 nt) and 7QXA (30 nt) can answer it. They are unreachable from the sandbox where the
rest of this analysis was done, so this script is written to run anywhere with internet access.

WHAT IT DOES
------------
   1. downloads 7QXS 7QXA 9SHY 9SHZ 9SI0 7TRE 7TRD 7TRF
   2. identifies every chain by RESIDUE COMPOSITION, never by chain letter
   3. identifies TERT by contacting BOTH a DNA and an RNA chain
   4. identifies TPP1 BY CONTACT to TERT, never by chain length
      (7TRF carries histone H2A/H2B at 90/82 residues; a length rule mislabels them as TPP1)
   5. measures PLYQ / Y176 / TRAP-wedge / anchor strand / catalytic aspartate to DNA
   6. reports every TEN-domain residue that contacts DNA in each entry
   7. compares the three 2026 RAP states superposed on the RT domain
   8. applies a verdict rule fixed BEFORE the result is seen

USAGE
-----
    python colab/tpp1_colab.py              # full run (needs internet)
    python colab/tpp1_colab.py --selftest   # offline: verify against the 3 bundled CIFs
    python colab/tpp1_colab.py --offline    # skip downloads, use data/raw only

Outputs: colab_out/report.md, colab_out/results.json, colab_out/distance_table.csv

DESIGN NOTE
-----------
This script is allowed to fail, and says so. If a structure cannot be downloaded, the verdict
prints "the test did not run" - never "no contact". Three errors in this project came from
reading a missing measurement as a negative one.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CIF_DIR = REPO / "colab" / "cif"
OUT_DIR = REPO / "colab_out"

# ----------------------------------------------------------------------------- targets
STRUCTURES = [
    # pdb,    label,                                         expect_tpp1
    ("7QXS", "Sekne 2022, telomerase-DNA-TPP1-POT1", True),
    ("7QXA", "Sekne 2022, telomerase-DNA-TPP1", True),
    ("9SHY", "Balch 2026, initiation", True),
    ("9SHZ", "Balch 2026, elongation", True),
    ("9SI0", "Balch 2026, pre-termination", True),
    ("7TRE", "Liu 2022, TERT+hTR+TPP1+7nt DNA (local reference)", True),
    ("7TRD", "Liu 2022, apo (local reference, NO TPP1)", False),
    ("7TRF", "Liu 2022, apo + histones (local reference, NO TPP1)", False),
]

LOCAL_ONLY = {"7TRE", "7TRD", "7TRF"}

# hTERT author numbering as used in all these depositions
REGIONS = {
    "PLYQ_174_177": [174, 175, 176, 177],
    "Y176_only": [176],
    "TEN_1_196": list(range(1, 197)),
    "FGF_loop_99_101": [99, 100, 101],
    "G100": [100],
    "TRAP_wedge_752_754": [752, 753, 754],
    "anchor_strand_756_759": [756, 757, 758, 759],
    "IFD_766_775": list(range(766, 776)),
    "RT_channel_940_984": list(range(940, 985)),
    "catalytic_D712": [712],
    "catalytic_D868_869": [868, 869],
}

DOMAINS = {
    "TEN 1-196": range(1, 197),
    "TRBD 326-540": range(326, 541),
    "RT 541-935": range(541, 936),
    "IFD-TRAP 745-805": range(745, 806),
}

DNA_RES = {"DA", "DT", "DG", "DC", "DU"}
RNA_RES = {"A", "U", "G", "C", "I"}
AMINO = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
            "THR TRP TYR VAL".split())

CONTACT = 4.5          # A, heavy-atom contact cutoff
NEAR = 5.0             # A, "contact" for the anchor verdict
FAR = 10.0             # A, "no contact" for the anchor verdict

# Self-test expectations, measured independently in scripts/solve/wp1_anchor_geometry.py
SELFTEST = {
    "7TRE": {"PLYQ_174_177": 11.11, "anchor_strand_756_759": 5.02, "TPP1_present": True},
    "7TRD": {"PLYQ_174_177": 11.94, "anchor_strand_756_759": 4.35, "TPP1_present": False},
    "7TRF": {"PLYQ_174_177": 11.94, "anchor_strand_756_759": 4.34, "TPP1_present": False},
}


# ----------------------------------------------------------------------------- deps
def ensure_deps() -> None:
    missing = []
    for mod, pkg in (("numpy", "numpy"), ("gemmi", "gemmi"), ("scipy", "scipy")):
        try:
            __import__(mod)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"installing missing packages: {missing}")
        os.system(f"{sys.executable} -m pip -q install " + " ".join(missing))
        for mod in ("numpy", "gemmi", "scipy"):
            __import__(mod)


# ----------------------------------------------------------------------------- download
def fetch(pdb: str, offline: bool = False) -> tuple[Path | None, str]:
    """Return (path, status). Status is 'cached', 'downloaded' or a failure reason."""
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    dest = CIF_DIR / f"{pdb}.cif"

    if dest.exists() and dest.stat().st_size > 20_000:
        return dest, "cached"

    local = REPO / "data" / "raw" / f"{pdb}.cif"
    if local.exists() and local.stat().st_size > 20_000:
        dest.write_bytes(local.read_bytes())
        return dest, "copied from data/raw"

    if offline:
        return None, "offline mode, not in data/raw"

    urls = [
        f"https://files.rcsb.org/download/{pdb}.cif",
        f"https://models.rcsb.org/{pdb}.cif",
        f"https://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/"
        f"{pdb[1:3].lower()}/{pdb.lower()}.cif",
    ]
    last = "no URL tried"
    for url in urls:
        for attempt in range(2):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "telomerase-audit/1.0"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = r.read()
                if len(data) > 20_000:
                    dest.write_bytes(data)
                    return dest, "downloaded"
                last = f"{len(data)} bytes from {url}"
            except Exception as e:                              # noqa: BLE001
                last = f"{type(e).__name__}: {e}"
                time.sleep(1.5)
                continue
    return None, last


# ----------------------------------------------------------------------------- parsing
def load_model(path: Path):
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st[0]


def classify(chain) -> str:
    names = {r.name for r in chain}
    if names <= AMINO:
        return "PROTEIN"
    if names <= DNA_RES:
        return "DNA"
    if names <= RNA_RES:
        return "RNA"
    if names & DNA_RES:
        return "DNA"
    if names & RNA_RES:
        return "RNA"
    return "OTHER"


def heavy_atoms(chain):
    import numpy as np
    return [(r.seqid.num, r.name, a.name,
             np.array([a.pos.x, a.pos.y, a.pos.z]))
            for r in chain for a in r if a.element.name not in ("H", "D", "")]


def residue_xyz(chain) -> dict[int, "np.ndarray"]:
    import numpy as np
    out: dict[int, list] = {}
    for r in chain:
        for a in r:
            if a.element.name in ("H", "D", ""):
                continue
            out.setdefault(r.seqid.num, []).append([a.pos.x, a.pos.y, a.pos.z])
    return {k: np.array(v) for k, v in out.items()}


def ca_map(chain) -> dict[int, "np.ndarray"]:
    import numpy as np
    out = {}
    for r in chain:
        for a in r:
            if a.name == "CA":
                out[r.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out


def rigid_transform(P, Q):
    """Rotation+translation taking P onto Q (Kabsch)."""
    import numpy as np
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Q.mean(0) - R @ P.mean(0)


def domain_motion(ref_ca, tgt_ca, rt_range) -> float | None:
    """
    How far a domain moves RELATIVE to the TERT ring.

    The ring (RT) is superposed onto the target; the same transform is then applied to
    the domain, and its residual displacement is measured.  This is the quantity that
    matters: a domain that translates rigidly with the ring has NOT moved, while one
    that swings away from the ring has.

    Superposing each domain independently (plain per-domain RMSD) would absorb any rigid
    shift and report ~0 for exactly the motion we are looking for.  That bug was found
    by injecting a known 1.5 A shift and watching it vanish.
    """
    import numpy as np
    common_rt = [i for i in rt_range if i in ref_ca and i in tgt_ca]
    common_dm = [i for i in ref_ca if i in tgt_ca]
    if len(common_rt) < 20 or len(common_dm) < 10:
        return None
    A = np.array([ref_ca[i] for i in common_rt])
    B = np.array([tgt_ca[i] for i in common_rt])
    R, tr = rigid_transform(A, B)

    # apply the ring transform to the whole molecule, then look at the domain only
    dom = [i for i in common_dm if i in set(rt_range) or True]
    ref_all = np.array([ref_ca[i] for i in common_dm])
    tgt_all = np.array([tgt_ca[i] for i in common_dm])
    moved = (R @ ref_all.T).T + tr
    disp = np.linalg.norm(moved - tgt_all, axis=1)
    idx = {r: k for k, r in enumerate(common_dm)}
    return idx, disp, common_dm


# ----------------------------------------------------------------------------- analysis
def analyse_entry(pdb: str, path: Path, label: str, expect_tpp1: bool) -> dict:
    import numpy as np
    from scipy.spatial import cKDTree

    m = load_model(path)
    groups: dict[str, list] = {}
    for ch in m:
        groups.setdefault(classify(ch), []).append(ch)

    rec: dict = {"pdb": pdb, "label": label, "expect_tpp1": expect_tpp1,
                 "chains": {k: [(c.name, len(c)) for c in v] for k, v in groups.items()},
                 "warnings": []}

    prot = groups.get("PROTEIN", [])
    dna_chains = groups.get("DNA", [])
    rna_chains = groups.get("RNA", [])

    if not prot:
        rec["error"] = "no protein chain"
        return rec

    # --- TERT: the protein chain that contacts BOTH DNA and RNA, else the longest
    def contacts_any(ch, other_chains) -> int:
        if not other_chains:
            return 0
        P = np.array([x[3] for x in heavy_atoms(ch)])
        tree = cKDTree(np.vstack([np.array([x[3] for x in heavy_atoms(o)])
                                  for o in other_chains]))
        return int((tree.query(P)[0] <= CONTACT).sum())

    scores = [(ch, contacts_any(ch, dna_chains) + contacts_any(ch, rna_chains)) for ch in prot]
    tert, best = max(scores, key=lambda t: t[1])
    if best == 0:
        tert = max(prot, key=len)
        rec["warnings"].append("TERT identified by length; no chain contacts DNA+RNA")
    if len(tert) < 500:
        rec["warnings"].append(f"TERT candidate only {len(tert)} residues - may be incomplete")

    rec["tert_chain"] = f"{tert.name} ({len(tert)} res)"

    # --- DNA: longest DNA chain
    if not dna_chains:
        rec["error"] = "no DNA chain"
        return rec
    dna_c = max(dna_chains, key=len)
    dna_res = residue_xyz(dna_c)
    dna_ids = sorted(dna_res)
    dna_atoms = np.vstack([dna_res[i] for i in dna_ids])
    tree = cKDTree(dna_atoms)
    rec["n_dna_nt"] = len(dna_ids)
    rec["dna_range"] = f"{dna_ids[0]}-{dna_ids[-1]}"

    tert_res = residue_xyz(tert)

    # --- region distances
    for name, nums in REGIONS.items():
        pts = [tert_res[n] for n in nums if n in tert_res]
        rec[name] = (round(float(tree.query(np.vstack(pts))[0].min()), 2) if pts else None)

    # --- every TERT residue in contact with DNA
    contacting = []
    for resid, xyz in tert_res.items():
        if float(tree.query(xyz)[0].min()) <= CONTACT:
            contacting.append(resid)
    contacting.sort()
    rec["n_dna_contacting_residues"] = len(contacting)
    rec["dna_contacting_TEN_residues"] = [r for r in contacting if r <= 196]
    rec["dna_contacting_anchor_residues"] = [r for r in contacting if 170 <= r <= 180]

    # --- per-nucleotide profile: DNA residue -> min distance to PLYQ
    plyq = [tert_res[n] for n in REGIONS["PLYQ_174_177"] if n in tert_res]
    if plyq:
        plyq_tree = cKDTree(np.vstack(plyq))
        rec["plyq_per_nt"] = [
            {"nt": int(i), "d": round(float(plyq_tree.query(dna_res[i])[0].min()), 2)}
            for i in dna_ids]

    # --- TPP1 by CONTACT (never by length: 7TRF has histones at 90/82 res)
    tert_atoms = np.vstack(list(tert_res.values()))
    ttree = cKDTree(tert_atoms)
    cands = []
    for ch in prot:
        if ch is tert:
            continue
        P = np.array([x[3] for x in heavy_atoms(ch)])
        n_contacts = int((ttree.query(P)[0] <= CONTACT).sum())
        if n_contacts > 0:
            cands.append((ch, n_contacts))
    if cands:
        tpp1, ncon = max(cands, key=lambda t: t[1])
        rec["TPP1_chain"] = f"{tpp1.name} ({len(tpp1)} res, {ncon} TERT contacts)"
        rec["TPP1_present"] = True
        P = np.array([x[3] for x in heavy_atoms(tpp1)])
        rec["TPP1_to_DNA"] = round(float(tree.query(P)[0].min()), 2)
        rec["TPP1_to_PLYQ"] = (round(float(np.min(np.linalg.norm(
            P[:, None, :] - np.vstack(plyq)[None, :, :], axis=2))), 2) if plyq else None)
    else:
        rec["TPP1_chain"] = "ABSENT"
        rec["TPP1_present"] = False
        rec["TPP1_to_DNA"] = None
        rec["TPP1_to_PLYQ"] = None

    if rec["TPP1_present"] != expect_tpp1:
        rec["warnings"].append(
            f"TPP1 present={rec['TPP1_present']} but {expect_tpp1} was expected")

    rec["_ca"] = ca_map(tert)
    return rec


def compare_states(recs: list[dict]) -> dict:
    """
    Superpose the 2026 RAP states on the RT domain (the rigid TERT ring) and report how
    far each domain moves RELATIVE to that ring, in Angstrom.
    """
    import numpy as np
    states = [r for r in recs
              if r["pdb"] in ("9SHY", "9SHZ", "9SI0") and "_ca" in r]
    if len(states) < 2:
        return {"available": False,
                "reason": f"{len(states)} of 3 RAP states available"}

    ref = states[0]
    rt = list(DOMAINS["RT 541-935"])
    out = {"available": True, "reference": ref["pdb"], "domains": {},
           "tpp1_in_all_states": all(r["TPP1_present"] for r in states)}

    for dname, rng in DOMAINS.items():
        row = {}
        for s in states[1:]:
            common = [i for i in rng if i in ref["_ca"] and i in s["_ca"]]
            common_rt = [i for i in rt if i in ref["_ca"] and i in s["_ca"]]
            if len(common) < 10 or len(common_rt) < 20:
                row[s["pdb"]] = None
                continue
            A = np.array([ref["_ca"][i] for i in common_rt])
            B = np.array([s["_ca"][i] for i in common_rt])
            R, tr = rigid_transform(A, B)
            dom_ref = np.array([ref["_ca"][i] for i in common])
            dom_tgt = np.array([s["_ca"][i] for i in common])
            moved = (R @ dom_ref.T).T + tr          # ring-superposed position
            row[s["pdb"]] = round(float(np.linalg.norm(moved - dom_tgt, axis=1).mean()), 3)
        out["domains"][dname] = row
    return out


# ----------------------------------------------------------------------------- verdict
def verdict(recs: list[dict]) -> dict:
    """The rule is stated here, before the numbers are seen."""
    out: dict = {"status": "", "text": [], "long_substrate_available": False}

    usable = [r for r in recs if not r.get("error") and r.get("n_dna_nt")]
    long_dna = [r for r in usable if r["n_dna_nt"] >= 15]
    out["long_substrate_available"] = bool(long_dna)
    out["structures_used"] = [r["pdb"] for r in usable]

    if not long_dna:
        out["status"] = "TEST DID NOT RUN"
        out["text"] = [
            "No structure with >=15 nt of DNA was available, so the anchor question",
            "was not tested. This is NOT a negative result. Re-run with network access.",
        ]
        return out

    vals = [(r["pdb"], r["PLYQ_174_177"]) for r in long_dna if r.get("PLYQ_174_177") is not None]
    if not vals:
        out["status"] = "TEST DID NOT RUN"
        out["text"] = ["PLYQ was not modelled in any long-substrate structure."]
        return out

    best_pdb, best = min(vals, key=lambda t: t[1])
    out["closest_plyq"] = {"pdb": best_pdb, "distance_A": best}

    if best <= NEAR:
        out["status"] = "ANCHOR SITE ENGAGES DNA"
        out["text"] = [
            f"PLYQ comes within {best:.2f} A of the DNA in {best_pdb}.",
            "The TEN anchor site DOES contact DNA once the substrate is long enough.",
            "The 11.1 A figure measured a site the 7-mer never reached.",
        ]
    elif best >= FAR:
        out["status"] = "ANCHOR SITE DOES NOT ENGAGE DNA"
        out["text"] = [
            f"PLYQ stays >= {best:.2f} A from DNA even at {long_dna[0]['n_dna_nt']}+ nt.",
            "The anchor-site reading does NOT survive a long substrate.",
            "The repository's original conclusion was right on this point.",
        ]
    else:
        out["status"] = "INTERMEDIATE - APPROACH, NOT CONTACT"
        out["text"] = [
            f"PLYQ approaches to {best:.2f} A in {best_pdb} - between {NEAR} and {FAR} A.",
            "This is an approach, not a hydrogen-bonding contact.",
            "Report as such; do not upgrade it to a contact.",
        ]

    # supporting evidence, independent of the headline
    ten_hits = {r["pdb"]: r.get("dna_contacting_TEN_residues", []) for r in long_dna}
    out["TEN_residues_contacting_DNA"] = ten_hits
    out["any_TEN_contact"] = any(bool(v) for v in ten_hits.values())
    return out


# ----------------------------------------------------------------------------- output
def write_outputs(recs: list[dict], states: dict, verd: dict) -> None:
    import csv
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    strip = [{k: v for k, v in r.items() if not k.startswith("_")} for r in recs]
    (OUT_DIR / "results.json").write_text(
        json.dumps({"structures": strip, "states": states, "verdict": verd}, indent=2))

    cols = ["pdb", "n_dna_nt", "dna_range"] + list(REGIONS) + [
        "TPP1_present", "TPP1_to_DNA", "TPP1_to_PLYQ", "n_dna_contacting_residues"]
    with open(OUT_DIR / "distance_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(strip)

    L = []
    L.append("# TPP1 / telomerase structural test\n")
    L.append("Generated by `colab/tpp1_colab.py`.\n")
    L.append(f"## Verdict: **{verd['status']}**\n")
    L.append(" ".join(verd["text"]) + "\n")

    L.append("\n## Distance to DNA by region (A, heavy atom)\n")
    L.append("| PDB | DNA nt | PLYQ 174-177 | Y176 | anchor 756-759 | TRAP-wedge | IFD 766-775 | catalytic D868/9 | TPP1 present | TPP1-DNA |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in strip:
        if r.get("error"):
            L.append(f"| {r['pdb']} | - | ERROR: {r['error']} | | | | | | | |")
            continue
        def f(k):
            v = r.get(k)
            return f"{v:.2f}" if isinstance(v, (int, float)) else "-"
        L.append(f"| {r['pdb']} | {r.get('n_dna_nt','-')} | {f('PLYQ_174_177')} | {f('Y176_only')} | "
                 f"{f('anchor_strand_756_759')} | {f('TRAP_wedge_752_754')} | "
                 f"{f('IFD_766_775')} | {f('catalytic_D868_869')} | "
                 f"{'yes' if r.get('TPP1_present') else 'no'} | {f('TPP1_to_DNA')} |")

    L.append("\n## TEN-domain residues contacting DNA\n")
    for r in strip:
        if r.get("error") or r.get("n_dna_nt") is None:
            continue
        hits = r.get("dna_contacting_TEN_residues", [])
        L.append(f"- **{r['pdb']}** ({r['n_dna_nt']} nt): "
                 f"{len(hits)} TEN residues contact DNA" + (f" — {hits}" if hits else ""))

    if states.get("available"):
        L.append("\n## 2026 RAP states, superposed per domain\n")
        L.append(f"Reference: {states['reference']}. "
                 f"TPP1 present in all states: **{states['tpp1_in_all_states']}**\n")
        L.append("| domain | " + " | ".join(states["domains"][list(states["domains"])[0]]) + " |")
        L.append("|---" * (1 + len(states["domains"][list(states["domains"])[0]])) + "|")
        for d, row in states["domains"].items():
            cells = [f"{v:.2f} A" if v is not None else "-" for v in row.values()]
            L.append(f"| {d} | " + " | ".join(cells) + " |")
        L.append("\nIf TPP1 is present in all three states, these differences are CYCLE "
                 "motion, not TPP1-driven motion.")
    else:
        L.append(f"\n## 2026 RAP states\n\nNot compared: {states.get('reason')}\n")

    warn = [(r["pdb"], wmsg) for r in strip for wmsg in r.get("warnings", [])]
    if warn:
        L.append("\n## Warnings\n")
        for pdb, wm in warn:
            L.append(f"- {pdb}: {wm}")

    (OUT_DIR / "report.md").write_text("\n".join(L) + "\n")


# ----------------------------------------------------------------------------- selftest
def selftest() -> int:
    print("SELF-TEST against the three bundled CIFs (offline)")
    print("=" * 74)
    ok = True
    for pdb, exp in SELFTEST.items():
        path = REPO / "data" / "raw" / f"{pdb}.cif"
        if not path.exists():
            print(f"  {pdb}: MISSING {path}")
            ok = False
            continue
        r = analyse_entry(pdb, path, "selftest", exp["TPP1_present"])
        checks = []
        for key in ("PLYQ_174_177", "anchor_strand_756_759"):
            got, want = r.get(key), exp[key]
            good = got is not None and abs(got - want) <= 0.05
            checks.append((key, got, want, good))
        good_tpp1 = r.get("TPP1_present") == exp["TPP1_present"]
        print(f"  {pdb}: TPP1 present={r.get('TPP1_present')} "
              f"({'ok' if good_tpp1 else 'MISMATCH'})   {r.get('TPP1_chain')}")
        for key, got, want, good in checks:
            print(f"       {key:24s} got {got}  want {want}  {'ok' if good else 'MISMATCH'}")
            ok = ok and good
        ok = ok and good_tpp1
        if r.get("warnings"):
            print(f"       warnings: {r['warnings']}")
    print()
    print("VERDICT-BRANCH TESTS (fabricated records; every branch must be reachable)")
    def fake(pdb, ndna, plyq):
        return {"pdb": pdb, "n_dna_nt": ndna, "PLYQ_174_177": plyq,
                "dna_contacting_TEN_residues": [], "TPP1_present": True}

    branches = [
        ("contact",      [fake("X1", 30, 3.4)],   "ANCHOR SITE ENGAGES DNA"),
        ("no contact",   [fake("X2", 30, 12.1)],  "ANCHOR SITE DOES NOT ENGAGE DNA"),
        ("intermediate", [fake("X3", 30, 7.2)],   "INTERMEDIATE - APPROACH, NOT CONTACT"),
        ("no long DNA",  [fake("X4", 7, 11.1)],   "TEST DID NOT RUN"),
        ("missing PLYQ", [{"pdb": "X5", "n_dna_nt": 30, "PLYQ_174_177": None,
                           "dna_contacting_TEN_residues": [], "TPP1_present": True}],
         "TEST DID NOT RUN"),
    ]
    for name, recs, want in branches:
        got = verdict(recs)["status"]
        good = got == want
        ok = ok and good
        print(f"  {name:14s} -> {got:38s} {'ok' if good else 'MISMATCH'}")

    print()
    print("RAP-STATE COMPARISON TEST (inject a known motion; it must be recovered)")
    import numpy as np
    rec = analyse_entry("7TRE", REPO / "data" / "raw" / "7TRE.cif", "t", True)
    ca = rec["_ca"]

    def _perturb(ca, domain, shift):
        out = {k: v.copy() for k, v in ca.items()}
        for k in out:
            if k in domain:
                out[k] = out[k] + shift
        return out

    s1 = {"pdb": "9SHY", "_ca": ca, "TPP1_present": True}
    s2 = {"pdb": "9SHZ", "_ca": _perturb(ca, range(1, 197), np.array([1.5, 0.0, 0.0])),
          "TPP1_present": True}
    s3 = {"pdb": "9SI0", "_ca": _perturb(ca, range(1, 197), np.array([3.0, 0.0, 0.0])),
          "TPP1_present": True}
    cmp = compare_states([s1, s2, s3])
    ten = cmp["domains"]["TEN 1-196"]
    rt = cmp["domains"]["RT 541-935"]
    got = [ten["9SHZ"], ten["9SI0"], rt["9SHZ"], rt["9SI0"]]
    want = [1.5, 3.0, 0.0, 0.0]
    branch = all(abs(g - w) < 0.05 for g, w in zip(got, want))
    ok = ok and branch
    print(f"  TEN 9SHZ/9SI0 = {got[0]:.2f}/{got[1]:.2f} A (want 1.50/3.00)")
    print(f"  RT  9SHZ/9SI0 = {got[2]:.2f}/{got[3]:.2f} A (want 0.00/0.00)")
    print(f"  {'ok' if branch else 'MISMATCH - ring superposition is absorbing real motion'}")

    fallback = compare_states([s1])["available"] is False and compare_states([])["available"] is False
    ok = ok and fallback
    print(f"  1-state / 0-state fallbacks: {'ok' if fallback else 'MISMATCH'}")

    print()
    print("=" * 74)
    print("SELF-TEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


# ----------------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true", help="offline check against bundled CIFs")
    ap.add_argument("--offline", action="store_true", help="do not download")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    ensure_deps()

    print("=" * 96)
    print("TPP1 / telomerase structural test")
    print("=" * 96)

    recs, status = [], {}
    for pdb, label, expect in STRUCTURES:
        path, st = fetch(pdb, offline=args.offline)
        status[pdb] = st if path else f"FAILED ({st})"
        print(f"  {pdb:5s} {label[:52]:52s} {status[pdb]}")
        if path is None:
            recs.append({"pdb": pdb, "label": label, "error": f"not available: {st}"})
            continue
        try:
            recs.append(analyse_entry(pdb, path, label, expect))
        except Exception as e:                                   # noqa: BLE001
            recs.append({"pdb": pdb, "label": label, "error": f"{type(e).__name__}: {e}"})
            print(f"        analysis failed: {e}")

    print()
    hdr = (f"  {'PDB':5s} {'nt':>4s} {'PLYQ':>7s} {'Y176':>7s} {'anchor':>7s} "
           f"{'IFD':>7s} {'D868/9':>7s} {'TPP1':>5s} {'TPP1-DNA':>9s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in recs:
        if r.get("error"):
            print(f"  {r['pdb']:5s}  ERROR: {r['error'][:70]}")
            continue
        def g(k):
            v = r.get(k)
            return f"{v:7.2f}" if isinstance(v, (int, float)) else f"{'-':>7s}"
        print(f"  {r['pdb']:5s} {r.get('n_dna_nt', 0):4d} {g('PLYQ_174_177')} "
              f"{g('Y176_only')} {g('anchor_strand_756_759')} {g('IFD_766_775')} "
              f"{g('catalytic_D868_869')} "
              f"{'yes' if r.get('TPP1_present') else 'no':>5s} {g('TPP1_to_DNA')}")

    states = compare_states(recs)
    verd = verdict(recs)
    write_outputs(recs, states, verd)

    print()
    print("=" * 96)
    print(f"VERDICT: {verd['status']}")
    print("=" * 96)
    for line in verd["text"]:
        print("  " + line)

    if states.get("available"):
        print()
        print(f"2026 RAP states (reference {states['reference']}); "
              f"TPP1 in all: {states['tpp1_in_all_states']}")
        for d, row in states["domains"].items():
            cells = ", ".join(f"{k} {v:.2f}A" if v is not None else f"{k} -"
                              for k, v in row.items())
            print(f"  {d:18s} {cells}")

    print()
    print(f"wrote {OUT_DIR/'report.md'}")
    print(f"wrote {OUT_DIR/'results.json'}")
    print(f"wrote {OUT_DIR/'distance_table.csv'}")
    print()
    print("Paste colab_out/report.md (or just the VERDICT block) back into the conversation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
