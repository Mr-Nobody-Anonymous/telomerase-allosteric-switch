#!/usr/bin/env python3
"""
Allosteric / mechanical coupling analysis of human telomerase (TERT–TPP1).

Why this exists
---------------
`results/verification/CLAIM_VERIFICATION_REPORT.md` showed that the repository's mechanistic
claims rest on static coordinate differences, which cannot distinguish a mechanical link from
an independent-module arrangement. The open question is the bottleneck of the whole project:

    TPP1 docks 25-29 A away from the DNA primer and never touches it
    (TERT-side interface at 2.3 A, but TPP1->primer = 28.7 A, TPP1->hTR = 14.2 A).
    How does that docking event raise repeat-addition processivity 3-4 fold?

This script tests the two competing mechanical explanations with an anisotropic elastic network
model (ANM) of the deposited coordinates:

  H_switch  "signal propagation": the docking surface is mechanically coupled to the primer
            exit channel, so a ligand-induced change propagates 25 A along the protein.
  H_clamp   "two-hand retention": TPP1 does not send a signal anywhere; it staples the two
            DNA-binding modules of TERT (the TEN-domain DNA site and the IFD-TRAP DNA site)
            together, so that they move as one rigid unit and the primer is retained by two
            hands whose separation no longer fluctuates.

The two hypotheses make different predictions and are separately testable.

Systems (matched, same TERT coordinates to 0.34 A):
  BOUND  7TRE  TERT + hTR + DNA + TPP1
  APO    7TRD  TERT + hTR + DNA
  DECOY  BOUND with the TPP1 rigid body translated 250 A away -> same number of beads, same
         internal elasticity, no contact. Any BOUND-vs-DECOY difference is due to *specific*
         ligation rather than to adding a body to the network.

Observables
  * centroid-separation coupling index  K = 1 - Var(s_A - s_B) / (Var(s_A) + Var(s_B))
    K = 1 rigidly coupled, 0 independent, < 0 anti-correlated.
  * distance-normalised cross-correlation residual (proximity-corrected coupling)
  * site-to-site coupling with a distance-matched permutation null (z-score)
  * force-response (linear response) transfer with the same null
  * per-residue coupling ranking + strongest-coupling path (candidate transducer positions)
  * retrospective validation against the published loss-of-function positions

Limitations are stated in the emitted report. This is a harmonic single-structure model: it
reports *possible* mechanical communication, not energetics, and it cannot replace experiment.

Usage
    python scripts/allosteric_coupling.py [--cutoff 15] [--modes 60] [--out results/allostery]
"""

from __future__ import annotations

import argparse
import heapq
import json
import time
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import scipy.linalg
from scipy import sparse
from scipy.spatial import cKDTree

PROTEIN = "protein"
NUCLEIC = "nucleic"

# Published loss-of-function positions used ONLY for retrospective validation.
PUBLISHED_FUNCTIONAL = {
    176: "TEN PLYQ (Y176A/Q177A, Sekne 2022)",
    177: "TEN PLYQ (Y176A/Q177A, Sekne 2022)",
    757: "IFD-TRAP (K757A/F759A, Sekne 2022)",
    759: "IFD-TRAP (K757A/F759A, Sekne 2022)",
    794: "IFD-TRAP (Q794A, Sekne 2022)",
    749: "IFD-TRAP (K749E, Tesmer 2023)",
    756: "IFD-TRAP (R756E, Tesmer 2023)",
}
CATALYTIC = (712, 868, 869)

LOG: list[str] = []


def log(msg: str = "") -> None:
    LOG.append(msg)
    print(msg, flush=True)


# ---------------------------------------------------------------------------------------
# Structure loading
# ---------------------------------------------------------------------------------------


class Beads:
    def __init__(self, labels, xyz, kind):
        self.labels: list[tuple[str, int, str]] = labels
        self.xyz: np.ndarray = xyz
        self.kind: np.ndarray = np.array(kind)

    def __len__(self):
        return len(self.labels)

    def resnums(self, chain: str) -> dict[int, int]:
        return {r: i for i, (c, r, _n) in enumerate(self.labels) if c == chain}


def load_beads(cif_path: Path, chains: dict[str, str]) -> Beads:
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]
    labels, xyz, kinds = [], [], []
    for chain_id in chains.values():
        if chain_id is None:
            continue
        for ch in model:
            if ch.name != chain_id:
                continue
            for r in ch:
                tab = gemmi.find_tabulated_residue(r.name)
                is_aa, is_na = tab.is_amino_acid(), tab.is_nucleic_acid()
                if not (is_aa or is_na):
                    continue
                target = "CA" if is_aa else ("C1'" if any(a.name == "C1'" for a in r) else "P")
                pick = next((a for a in r if a.name == target), None)
                if pick is None:
                    continue
                labels.append((ch.name, r.seqid.num, r.name))
                xyz.append([pick.pos.x, pick.pos.y, pick.pos.z])
                kinds.append(PROTEIN if is_aa else NUCLEIC)
    return Beads(labels, np.array(xyz, float), kinds)


# ---------------------------------------------------------------------------------------
# Elastic network
# ---------------------------------------------------------------------------------------


def build_anm(xyz: np.ndarray, cutoff: float) -> sparse.csr_matrix:
    """Tirion ANM Hessian (unit springs) as a sparse 3N x 3N matrix."""
    n = len(xyz)
    pairs = cKDTree(xyz).query_pairs(cutoff, output_type="ndarray")
    rows, cols, vals = [], [], []
    diag = np.zeros((n, 3, 3))
    for i, j in pairs:
        d_vec = xyz[j] - xyz[i]
        dist = float(np.linalg.norm(d_vec))
        if dist < 1e-6:
            continue
        u = d_vec / dist
        block = np.outer(u, u)
        diag[i] += block
        diag[j] += block
        for a in range(3):
            for b in range(3):
                rows += [3 * i + a, 3 * j + a]
                cols += [3 * j + b, 3 * i + b]
                vals += [-block[a, b], -block[a, b]]
    for i in range(n):
        for a in range(3):
            for b in range(3):
                rows.append(3 * i + a)
                cols.append(3 * i + b)
                vals.append(diag[i, a, b])
    h = sparse.coo_matrix((vals, (rows, cols)), shape=(3 * n, 3 * n)).tocsr()
    return 0.5 * (h + h.T)


class Network:
    """ANM solution and residue-level observables."""

    def __init__(self, beads: Beads, cutoff: float, n_modes: int):
        self.beads = beads
        self.cutoff = cutoff
        t0 = time.time()
        h = build_anm(beads.xyz, cutoff)
        n_dof = h.shape[0]
        k = min(n_modes + 12, n_dof - 1)
        vals, vecs = scipy.linalg.eigh(h.toarray(), subset_by_index=[0, k])
        order = np.argsort(vals)
        vals, vecs = vals[order], vecs[:, order]
        keep = vals > 1e-8
        self.n_zero = int((~keep).sum())
        self.evals = vals[keep][:n_modes]
        self.evecs = vecs[:, keep][:, :n_modes]
        n_pairs = (h.nnz - 3 * len(beads)) // 18
        log(f"    ANM: {len(beads)} beads, cutoff {cutoff} A, {n_pairs} contacts "
            f"(mean degree {2 * n_pairs / len(beads):.1f}), {len(self.evals)} modes, "
            f"{self.n_zero} null modes discarded, {time.time() - t0:.1f} s")
        if self.n_zero > 8:
            log("      WARNING: more than the six rigid-body modes are null — the cutoff "
                "under-constrains this bead model; treat absolute values with care")

    # -- fluctuations ------------------------------------------------------------------
    def _mode_blocks(self) -> np.ndarray:
        return self.evecs.reshape(len(self.evals), -1, 3)

    def msf(self) -> np.ndarray:
        v = self._mode_blocks()
        return (1.0 / self.evals) @ (v ** 2).sum(axis=2)

    def centroid_msf(self, idx: np.ndarray) -> float:
        v = self._mode_blocks()[:, idx, :].mean(axis=1)          # (modes, 3)
        return float(np.sum((1.0 / self.evals) * (v ** 2).sum(axis=1)))

    def separation_variance(self, idx_a: np.ndarray, idx_b: np.ndarray) -> float:
        v = self._mode_blocks()
        va = v[:, idx_a, :].mean(axis=1) - v[:, idx_b, :].mean(axis=1)
        return float(np.sum((1.0 / self.evals) * (va ** 2).sum(axis=1)))

    def coupling_index(self, idx_a: np.ndarray, idx_b: np.ndarray) -> float:
        """
        K = 1 - Var(s_A - s_B) / (Var(s_A) + Var(s_B)).
        K = 1: the two modules move as one rigid unit (their separation is fixed).
        K = 0: their motions are independent.
        K < 0: their separation fluctuates more than independent motion would give.
        """
        var_ab = self.separation_variance(idx_a, idx_b)
        denom = self.centroid_msf(idx_a) + self.centroid_msf(idx_b)
        return float(1.0 - var_ab / denom) if denom > 0 else float("nan")

    # -- correlations ------------------------------------------------------------------
    def correlation_matrix(self, mask: np.ndarray) -> np.ndarray:
        v = self._mode_blocks()[:, mask, :]
        c = np.einsum("k,kia,kja->ij", 1.0 / self.evals, v, v, optimize=True)
        d = np.sqrt(np.diag(c))
        return c / np.outer(d, d)

    # -- response ----------------------------------------------------------------------
    def response(self, force: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
        proj = self.evecs.T @ force.reshape(-1)
        disp = (self.evecs * (proj / self.evals)).sum(axis=1)
        out = disp.reshape(-1, 3)
        return out[mask] if mask is not None else out


def coupling_residual(corr: np.ndarray, xyz: np.ndarray, n_bins: int = 40) -> np.ndarray:
    """|corr_ij| minus the median |corr| of pairs at the same bead-bead distance."""
    n = corr.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    d = np.linalg.norm(xyz[iu] - xyz[ju], axis=1)
    a = np.abs(corr[iu, ju])
    bins = np.linspace(np.percentile(d, 0.5), np.percentile(d, 99.5), n_bins + 1)
    idx = np.clip(np.digitize(d, bins) - 1, 0, n_bins - 1)
    expected = np.empty_like(a)
    for b in range(n_bins):
        m = idx == b
        expected[m] = np.median(a[m]) if m.sum() >= 5 else np.nan
    if not np.isfinite(expected).all():
        expected[~np.isfinite(expected)] = np.nanmedian(expected)
    resid = np.zeros((n, n))
    resid[iu, ju] = a - expected
    return resid + resid.T


def matched_null(resid, xyz, src_idx, tgt_idx, n_draws=1000, rng=None):
    """Null for a site-to-site score: random source residues at matched distance to the target."""
    rng = rng or np.random.default_rng(0)
    tgt_tree = cKDTree(xyz[tgt_idx])
    src_d = np.array([tgt_tree.query(xyz[i])[0] for i in src_idx])
    all_d = tgt_tree.query(xyz)[0]
    pool = np.setdiff1d(np.arange(len(xyz)), src_idx)
    draws = np.empty(n_draws)
    for k in range(n_draws):
        pick = []
        for dd in src_d:
            cand = pool[(all_d[pool] >= dd * 0.9) & (all_d[pool] <= dd * 1.1)]
            pick.append(rng.choice(cand) if len(cand) else rng.choice(pool))
        draws[k] = resid[np.ix_(np.array(pick), tgt_idx)].mean()
    return draws


def z_score(obs: float, null: np.ndarray) -> tuple[float, float, float]:
    return ((obs - null.mean()) / null.std() if null.std() > 0 else float("nan"),
            float(null.mean()), float(null.std()))


def strongest_path(resid, xyz, src, tgt, radius=22.0, penalty=0.02):
    """Dijkstra over the contact graph with edge cost 1/(coupling + penalty)."""
    n = resid.shape[0]
    contacts = cKDTree(xyz).query_ball_point(xyz, radius)
    w = 1.0 / (np.maximum(resid, 0.0) + penalty)
    adj = [[] for _ in range(n)]
    for i, nbrs in enumerate(contacts):
        for j in nbrs:
            if j > i:
                c = float(w[i, j])
                adj[i].append((j, c))
                adj[j].append((i, c))
    dist = np.full(n, np.inf)
    prev = np.full(n, -1)
    pq = []
    for s in np.atleast_1d(src).astype(int):
        dist[s] = 0.0
        heapq.heappush(pq, (0.0, int(s)))
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist[u]:
            continue
        for v_, c in adj[u]:
            nd = d + c
            if nd < dist[v_]:
                dist[v_], prev[v_] = nd, u
                heapq.heappush(pq, (nd, v_))
    ends = np.where(tgt)[0]
    if not len(ends):
        return []
    best = ends[np.argmin(dist[ends])]
    if not np.isfinite(dist[best]):
        return []
    path = [int(best)]
    while prev[path[-1]] != -1:
        path.append(int(prev[path[-1]]))
    return path[::-1]


# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------


def _subset(beads: Beads, keep: np.ndarray) -> Beads:
    """Sub-bean set with the same ordering convention (used to strip TPP1 for the APO system)."""
    labels = [lbl for lbl, k in zip(beads.labels, keep) if k]
    return Beads(labels, beads.xyz[keep], list(np.array(beads.kind)[keep]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--out", default="results/allostery")
    ap.add_argument("--cutoff", type=float, default=15.0)
    ap.add_argument("--modes", type=int, default=40)
    ap.add_argument("--cutoffs-robust", default="12,18")
    args = ap.parse_args()

    raw, out = Path(args.raw_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261030)

    log("=" * 96)
    log("TELOMERASE TERT MECHANICAL COUPLING — signal transmission, or a two-hand staple?")
    log("=" * 96)
    log("Design: all three systems share the IDENTICAL TERT/hTR/DNA coordinates (7TRE), so the")
    log("only variable is the presence and placement of the TPP1 OB fold. 7TRD and 7TRE TERT")
    log("superpose at 0.34 A (verified in results/verification), so 'APO' here means the same")
    log("conformation with the ligand removed, which isolates the ligand's mechanical effect.")

    full = load_beads(raw / "7TRE.cif", {"tert": "A", "rna": "B", "tpp1": "C", "dna": "D"})
    is_tpp1 = np.array([c == "C" for c, _r, _n in full.labels])
    apo = _subset(full, ~is_tpp1)
    log(f"\nBOUND {len(full)} beads ({int((full.kind == PROTEIN).sum())} protein, "
        f"{int((full.kind == NUCLEIC).sum())} nucleic acid)")
    log(f"APO   {len(apo)} beads (identical TERT/hTR/DNA, TPP1 deleted)")

    # ---- surfaces -----------------------------------------------------------------
    st = gemmi.read_structure(str(raw / "7TRE.cif"))
    heavy: dict[tuple[str, int], list[np.ndarray]] = {}
    for ch in st[0]:
        for r in ch:
            for a in r:
                if a.element.name != "H":
                    heavy.setdefault((ch.name, r.seqid.num), []).append(
                        np.array([a.pos.x, a.pos.y, a.pos.z]))
    tpp1_heavy = np.array([p for (c, _r), pts in heavy.items() if c == "C" for p in pts])
    dna_heavy = np.array([p for (c, _r), pts in heavy.items() if c == "D" for p in pts])
    tree_tpp1, tree_dna = cKDTree(tpp1_heavy), cKDTree(dna_heavy)

    dock_res = sorted({rn for (c, rn), pts in heavy.items()
                       if c == "A" and any(tree_tpp1.query(p)[0] <= 4.5 for p in pts)})
    primer_res = sorted({rn for (c, rn), pts in heavy.items()
                         if c == "A" and any(tree_dna.query(p)[0] <= 6.0 for p in pts)})
    ten_module = [r for r in range(7, 197)]
    trap_module = [r for r in range(740, 801)]
    catalytic = [712, 868, 869]
    # primer-channel residues: published anchor positions + everything within 8 A of the
    # modelled DNA 3' end (the primer that must be retained between repeats)
    dna3 = np.array([p for (c, rn), pts in heavy.items() if c == "D" and rn >= 14 for p in pts])
    tree3 = cKDTree(dna3)
    primer_channel = sorted({rn for (c, rn), pts in heavy.items()
                             if c == "A" and any(tree3.query(p)[0] <= 8.0 for p in pts)}
                            | {176, 177, 756, 757, 759, 794})
    # the docking surface minus anything that is itself part of the IFD-TRAP module, so that
    # "coupling from the dock to the primer channel" is not measuring self-coupling
    dock_only = [r for r in dock_res if not (740 <= r <= 800)]
    dock_trap = [r for r in dock_res if 740 <= r <= 800]

    log(f"\ndocking surface   (TERT <=4.5 A of TPP1, {len(dock_res)} res): {dock_res}")
    log(f"TEN module (DNA-binding, 7-196), mean min distance to DNA = "
        f"{np.mean([tree_dna.query(np.array([a.pos.x, a.pos.y, a.pos.z]))[0] for (c, rn), pts in heavy.items() if c == 'A' and rn in (176, 177) for a in []] or [0]):.1f}")
    log(f"  TEN DNA-anchor residues (<=6 A):  {[r for r in primer_res if r <= 196]}")
    log(f"  IFD-TRAP DNA-anchor residues:     {[r for r in primer_res if 740 <= r <= 800]}")
    log(f"  TPP1 -> DNA minimum heavy-atom distance: {tree_dna.query(tpp1_heavy)[0].min():.2f} A")
    log(f"  TPP1 -> TERT contacts at <=4.5 A: {len(dock_res)} residues = "
        f"{len(dock_only)} outside + {len(dock_trap)} inside the IFD-TRAP module "
        f"({dock_trap})")
    log(f"  primer-channel residues (<=8 A of the DNA 3' end + published anchors, "
        f"{len(primer_channel)}): {primer_channel}")

    # ---- the geometric premise: how far is TPP1 from every TERT residue? ----------
    tpp1_tree = cKDTree(tpp1_heavy)
    grad = []
    for (c, rn), pts in heavy.items():
        if c != "A":
            continue
        name = next((r.name for ch in st[0] if ch.name == "A" for r in ch if r.seqid.num == rn), "?")
        d = min(tpp1_tree.query(p)[0] for p in pts)
        grad.append({"resnum": rn, "resname": name, "min_dist_to_tpp1_A": round(float(d), 2)})
    grad_df = pd.DataFrame(grad).sort_values("resnum")
    grad_df.to_csv(out / "tpp1_distance_gradient.csv", index=False)

    def _ca_path_length(res_a: int, res_b: int) -> float:
        """Backbone path length between two TERT residues (successive modelled CA-CA)."""
        m = {rn: np.array([a.pos.x, a.pos.y, a.pos.z]) for ch in st[0] if ch.name == "A"
             for r in ch for a in r if a.name == "CA" for rn in [r.seqid.num]}
        lo, hi = sorted((res_a, res_b))
        pts = [m[r] for r in range(lo, hi + 1) if r in m]
        return float(sum(np.linalg.norm(pts[i + 1] - pts[i]) for i in range(len(pts) - 1)))

    log(f"\nTPP1 footprint reaches to within 4.5 A of TERT residues {min(dock_res)}-{max(dock_res)}")
    log("contact gradient along the IFD-TRAP helix (distance from TPP1):")
    for rn in range(756, 776):
        row = grad_df[grad_df.resnum == rn]
        if len(row):
            log(f"   {row.iloc[0]['resname']}{rn}: {row.iloc[0]['min_dist_to_tpp1_A']:5.2f} A")
    log(f"  backbone path length 757 -> 766 = {_ca_path_length(757, 766):.1f} A")

    def idx_of(beads: Beads, resnums, chain="A") -> np.ndarray:
        m = beads.resnums(chain)
        return np.array([m[r] for r in resnums if r in m])

    # ---- systems ------------------------------------------------------------------
    decoy_xyz = full.xyz.copy()
    decoy_xyz[is_tpp1] += np.array([250.0, 0.0, 0.0])
    decoy = Beads(list(full.labels), decoy_xyz, list(full.kind))
    systems = {"APO": apo, "BOUND": full, "DECOY": decoy}

    results, rows = {}, []
    for name, beads in systems.items():
        log(f"\n[{name}]")
        net = Network(beads, args.cutoff, args.modes)
        results[name] = net
        prot = net.beads.kind == PROTEIN
        prot_xyz = net.beads.xyz[prot]
        corr = net.correlation_matrix(prot)
        resid = coupling_residual(corr, prot_xyz)

        id_dock = idx_of(beads, dock_only)
        id_dock_all = idx_of(beads, dock_res)
        id_ten = idx_of(beads, ten_module)
        id_trap = idx_of(beads, trap_module)
        id_cat = idx_of(beads, catalytic)
        id_chan = idx_of(beads, primer_channel)
        # protein-mask coordinates for the residues of interest
        prot_labels = [lbl for lbl, m in zip(net.beads.labels, prot) if m]
        pos = {lbl: i for i, lbl in enumerate(prot_labels)}
        p_ten = np.array([pos[("A", r, _n)] for (c, r, _n) in beads.labels
                          if c == "A" and r in ten_module and ("A", r, _n) in pos])
        p_trap = np.array([pos[("A", r, _n)] for (c, r, _n) in beads.labels
                           if c == "A" and r in trap_module and ("A", r, _n) in pos])
        p_cat = np.array([pos[("A", r, _n)] for (c, r, _n) in beads.labels
                          if c == "A" and r in catalytic and ("A", r, _n) in pos])
        p_dock = np.array([pos[lbl] for lbl in (beads.labels[i] for i in id_dock)])
        p_chan = np.array([pos[lbl] for lbl in (beads.labels[i] for i in id_chan)
                           if lbl in pos])
        p_cat = np.array([pos[lbl] for lbl in (beads.labels[i] for i in id_cat)
                          if lbl in pos])
        p_tpp1 = [lbl for lbl in net.beads.labels if lbl[0] == "C"]
        p_tpp1 = np.array([pos[l] for l in p_tpp1 if l in pos])

        k_staple = net.coupling_index(p_ten, p_trap)
        k_control = net.coupling_index(p_cat, p_ten)
        k_dock = net.coupling_index(p_dock, p_ten)
        # positive control: two sub-sites of the same rigid fold must look coupled
        k_pos = net.coupling_index(p_cat, p_chan) if len(p_cat) and len(p_chan) else float("nan")
        # the switch metric, now measured from docking surface to primer channel
        k_switch = net.coupling_index(p_dock, p_chan) if len(p_dock) and len(p_chan) else float("nan")

        # site-to-site coupling with distance-matched null (dock -> primer channel)
        tgt_tree = cKDTree(prot_xyz[p_chan])
        all_d = tgt_tree.query(prot_xyz)[0]
        src_d = all_d[p_dock]
        pool = np.setdiff1d(np.arange(len(prot_xyz)), p_dock)
        obs = float(resid[np.ix_(p_dock, p_chan)].mean())
        draws = np.empty(400)
        for k in range(400):
            pick = []
            for dd in src_d:
                cand = pool[(all_d[pool] >= dd * 0.9) & (all_d[pool] <= dd * 1.1)]
                pick.append(rng.choice(cand) if len(cand) else rng.choice(pool))
            draws[k] = resid[np.ix_(np.array(pick), p_chan)].mean()
        z_switch = (obs - draws.mean()) / draws.std() if draws.std() > 0 else float("nan")

        # ---- positive controls: can this model detect coupling where it is known to exist?
        d_all = cKDTree(prot_xyz).query_pairs(8.0, output_type="ndarray")
        short_range = float(resid[d_all[:, 0], d_all[:, 1]].mean())
        # catalytic site -> DNA 3' end (a link that catalysis requires)
        dna3_idx = [i for i, (c, r, _n) in enumerate(beads.labels)
                    if c == "D" and r >= 16]
        full_corr = net.correlation_matrix(np.ones(len(beads), bool))
        cat_dna = float(np.abs(full_corr[np.ix_(id_cat, dna3_idx)]).mean()) if len(id_cat) and dna3_idx else float("nan")
        # ten-internal long-range pair (both in TEN, ~50 A apart)
        m_ten = beads.resnums("A")
        cat_dna_null = float(np.abs(full_corr[np.ix_(
            np.where(prot)[0][:len(id_cat)], dna3_idx)]).mean()) if dna3_idx else float("nan")

        row = {
            "system": name, "n_beads": len(beads), "n_modes": len(net.evals),
            "positive_control_short_range": short_range,
            "positive_control_catalytic_to_dna3": cat_dna,
            "null_modes": net.n_zero,
            "K_TEN_vs_IFDTRAP": k_staple,
            "K_catalytic_vs_TEN": k_control,
            "K_dock_vs_TEN": k_dock,
            "K_dock_vs_primer_channel": k_switch,
            "K_positive_control_cat_vs_channel": k_pos,
            "switch_coupling": obs, "switch_null_mean": float(draws.mean()),
            "switch_null_sd": float(draws.std()), "switch_z": z_switch,
        }
        rows.append(row)
        results[name]._aux = {"prot_labels": prot_labels, "pos": pos, "resid": resid,
                              "p_dock": p_dock, "p_ten": p_ten, "p_trap": p_trap,
                              "p_chan": p_chan, "id_dock": id_dock, "prot": prot}
        log(f"    K(TEN, IFD-TRAP)        = {k_staple:+.3f}   "
            f"[1 = rigid staple, 0 = independent modules]")
        log(f"    K(catalytic site, TEN)  = {k_control:+.3f}   [control: should not care about TPP1]")
        log(f"    K(docking surface, TEN) = {k_dock:+.3f}")
        log(f"    K(dock, primer channel) = {k_switch:+.3f}   "
            f"[H_switch would require a rise when TPP1 docks]")
        log(f"    positive control K(catalytic site, primer channel) = {k_pos:+.3f} "
            f"(overlapping sets — informational only)")
        log(f"    positive control: mean coupling residual for residue pairs <8 A apart = "
            f"{short_range:+.4f}  [must be strongly positive]")
        log(f"    positive control: |corr| catalytic site <-> DNA 3' end = {cat_dna:.3f}  "
            f"[must be well above the inter-chain background {cat_dna_null:.3f}]")
        log(f"    switch coupling dock->IFD-TRAP = {obs:+.4f} "
            f"(distance-matched null {draws.mean():+.4f} +- {draws.std():.4f}, z = {z_switch:+.2f})")

    table = pd.DataFrame(rows)
    table.to_csv(out / "system_comparison.csv", index=False)

    # ---- hypothesis test ----------------------------------------------------------
    b = table.set_index("system")
    log("\n" + "=" * 96)
    log("HYPOTHESIS TEST (all systems share the same TERT/hTR/DNA coordinates)")
    log("=" * 96)
    dK = b.loc["BOUND", "K_TEN_vs_IFDTRAP"] - b.loc["APO", "K_TEN_vs_IFDTRAP"]
    dK_decoy = b.loc["DECOY", "K_TEN_vs_IFDTRAP"] - b.loc["APO", "K_TEN_vs_IFDTRAP"]
    dC = b.loc["BOUND", "K_catalytic_vs_TEN"] - b.loc["APO", "K_catalytic_vs_TEN"]
    dC_decoy = b.loc["DECOY", "K_catalytic_vs_TEN"] - b.loc["APO", "K_catalytic_vs_TEN"]

    log(f"  H_clamp  (TPP1 staples the two DNA-binding modules together):")
    log(f"    K(TEN,IFD-TRAP):  APO {b.loc['APO','K_TEN_vs_IFDTRAP']:+.3f} -> "
        f"BOUND {b.loc['BOUND','K_TEN_vs_IFDTRAP']:+.3f}   Δ = {dK:+.3f}")
    log(f"    non-specific control (TPP1 present but 250 A away):       Δ = {dK_decoy:+.3f}")
    log(f"    specific ligation effect = {dK - dK_decoy:+.3f}")
    dKS = b.loc["BOUND", "K_dock_vs_primer_channel"] - b.loc["DECOY", "K_dock_vs_primer_channel"]
    log(f"    K(dock, primer channel): BOUND {b.loc['BOUND','K_dock_vs_primer_channel']:+.3f} "
        f"vs DECOY {b.loc['DECOY','K_dock_vs_primer_channel']:+.3f}   Δ = {dKS:+.3f}")
    log(f"  H_switch (signal reaches the primer channel from the dock):")
    log(f"    switch coupling Δ = "
        f"{b.loc['BOUND','switch_coupling'] - b.loc['DECOY','switch_coupling']:+.4f}; "
        f"z(BOUND) = {b.loc['BOUND','switch_z']:+.2f}, z(DECOY) = {b.loc['DECOY','switch_z']:+.2f}")
    log(f"  Positive control K(catalytic site, primer channel): BOUND "
        f"{b.loc['BOUND','K_positive_control_cat_vs_channel']:+.3f}, DECOY "
        f"{b.loc['DECOY','K_positive_control_cat_vs_channel']:+.3f}  "
        f"(the model must show strong coupling here or it lacks sensitivity)")
    log(f"  Control: catalytic-site/TEN coupling should be insensitive to TPP1:")
    log(f"    Δ = {dC:+.3f} (non-specific {dC_decoy:+.3f})")

    verdict_clamp = ("SUPPORTED" if (dK - dK_decoy) > 0.05 else
                     "MARGINAL" if (dK - dK_decoy) > 0.02 else "NOT SUPPORTED")
    dKS = b.loc["BOUND", "K_dock_vs_primer_channel"] - b.loc["DECOY", "K_dock_vs_primer_channel"]
    verdict_switch = ("SUPPORTED" if (dKS > 0.05 and b.loc["BOUND", "switch_z"] > 2.0)
                      else "NOT SUPPORTED")
    log(f"\n  => H_clamp : {verdict_clamp}")
    log(f"  => H_switch: {verdict_switch}")

    # ---- per-residue ranking (BOUND) ---------------------------------------------
    nb = results["BOUND"]
    aux = nb._aux
    resid, pos, p_dock, p_trap = aux["resid"], aux["pos"], aux["p_dock"], aux["p_trap"]
    msf_b = nb.msf()[aux["prot"]]
    na = results["APO"]
    msf_a = na.msf()[na._aux["prot"]]

    per_res = []
    for lbl, k in pos.items():
        if lbl[0] != "A":
            continue
        per_res.append({
            "resnum": lbl[1], "resname": lbl[2],
            "coupling_to_ifd_trap": float(resid[k, p_trap].mean()),
            "coupling_to_dock": float(resid[k, p_dock].mean()),
            "msf_bound": float(msf_b[k]),
            "msf_apo": float(msf_a[na._aux["pos"][lbl]]) if lbl in na._aux["pos"] else np.nan,
        })
    pr = pd.DataFrame(per_res)
    pr["rank_to_ifd_trap"] = pr["coupling_to_ifd_trap"].rank(ascending=False).astype(int)
    pr["published"] = pr["resnum"].map(lambda r: PUBLISHED_FUNCTIONAL.get(r, ""))
    pr["mobility_ratio_apo_bound"] = pr["msf_apo"] / pr["msf_bound"]
    pr.to_csv(out / "residue_coupling_ranked.csv", index=False)

    log("\n" + "=" * 96)
    log("RETROSPECTIVE VALIDATION — published loss-of-function positions")
    log("=" * 96)
    known = pr[pr["published"] != ""]
    rest = pr[pr["published"] == ""]
    log(known[["resnum", "resname", "rank_to_ifd_trap", "coupling_to_ifd_trap",
               "published"]].to_string(index=False))
    scores = np.concatenate([known["coupling_to_ifd_trap"].values,
                             rest["coupling_to_ifd_trap"].values])
    labels = np.concatenate([np.ones(len(known)), np.zeros(len(rest))])
    order = np.argsort(-scores)
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    n1, n0 = int(labels.sum()), int((1 - labels).sum())
    auc = 1 - ((ranks[labels == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    log(f"\n  AUC (published positions vs all other TERT residues) = {auc:.3f}")
    log("  Weak test only: 7 positions, several chosen because they contact DNA.")

    prot_xyz = nb.beads.xyz[aux["prot"]]
    path = strongest_path(resid, prot_xyz, p_dock, p_chan)
    log("\n" + "=" * 96)
    log("STRONGEST COUPLING PATH  docking surface -> IFD-TRAP")
    log("=" * 96)
    log("  " + " -> ".join(f"{aux['prot_labels'][i][2]}{aux['prot_labels'][i][1]}" for i in path))

    # ---- mobility change by domain ------------------------------------------------
    log("\n" + "=" * 96)
    log("MOBILITY CHANGE (ANM MSF ratio apo/bound) by domain — where TPP1 dampens motion")
    log("=" * 96)
    domains = {"TEN (7-196)": (7, 196), "TRBD (326-540)": (326, 540),
               "RT (541-739)": (541, 739), "IFD-TRAP (740-800)": (740, 800),
               "RT-after-IFD (801-935)": (801, 935), "CTE (936-1132)": (936, 1132)}
    dom = {}
    for dname, (lo, hi) in domains.items():
        sel = pr[(pr.resnum >= lo) & (pr.resnum <= hi)]
        dom[dname] = float(sel["mobility_ratio_apo_bound"].median()) if len(sel) else float("nan")
        log(f"  {dname:24s} n={len(sel):4d}  median MSF ratio = {dom[dname]:.3f}")

    # ---- robustness ---------------------------------------------------------------
    log("\n" + "=" * 96)
    log("ROBUSTNESS — repeat the staple metric at other contact cutoffs")
    log("=" * 96)
    robust = []
    for cut in [float(c) for c in args.cutoffs_robust.split(",")] + [args.cutoff]:
        rr, nn = {}, {}
        for name in ("APO", "BOUND", "DECOY"):
            beads = systems[name]
            net = Network(beads, cut, args.modes)
            prot_lbl = [lbl for lbl, m in zip(net.beads.labels, net.beads.kind == PROTEIN) if m]
            pm = {lbl: i for i, lbl in enumerate(prot_lbl)}
            p_ten = np.array([pm[("A", r, _n)] for (c, r, _n) in beads.labels
                              if c == "A" and r in ten_module and ("A", r, _n) in pm])
            p_trap = np.array([pm[("A", r, _n)] for (c, r, _n) in beads.labels
                               if c == "A" and r in trap_module and ("A", r, _n) in pm])
            rr[name] = net.coupling_index(p_ten, p_trap)
            nn[name] = net.n_zero
            del net
        rec = {"cutoff": cut, "K_apo": rr["APO"], "K_bound": rr["BOUND"],
               "K_decoy": rr["DECOY"], "null_modes_apo": nn["APO"],
               "null_modes_bound": nn["BOUND"], "null_modes_decoy": nn["DECOY"],
               "specific_effect": (rr["BOUND"] - rr["APO"]) - (rr["DECOY"] - rr["APO"])}
        robust.append(rec)
        log(f"  cutoff {cut:4.1f} A: K apo {rr['APO']:+.3f}, bound {rr['BOUND']:+.3f}, "
            f"decoy {rr['DECOY']:+.3f}  ->  specific effect {rec['specific_effect']:+.3f} "
            f"(null modes {nn['APO']}/{nn['BOUND']}/{nn['DECOY']})")
    pd.DataFrame(robust).to_csv(out / "cutoff_robustness.csv", index=False)

    summary = {
        "cutoff": args.cutoff, "modes": args.modes,
        "design": "APO/BOUND/DECOY share identical TERT-hTR-DNA coordinates (7TRE)",
        "dock_residues": dock_res, "primer_residues": primer_res,
        "ten_dna_anchor": [r for r in primer_res if r <= 196],
        "ifd_trap_dna_anchor": [r for r in primer_res if 740 <= r <= 800],
        "tpp1_min_distance_to_dna_A": float(tree_dna.query(tpp1_heavy)[0].min()),
        "table": table.to_dict(orient="records"),
        "H_clamp_specific_effect": float(dK - dK_decoy),
        "H_switch_specific_effect": float(b.loc["BOUND", "switch_coupling"]
                                          - b.loc["DECOY", "switch_coupling"]),
        "H_switch_K_effect": float(b.loc["BOUND", "K_dock_vs_primer_channel"]
                                   - b.loc["DECOY", "K_dock_vs_primer_channel"]),
        "positive_control_short_range": float(b.loc["BOUND", "positive_control_short_range"]),
        "positive_control_cat_to_dna3": float(b.loc["BOUND", "positive_control_catalytic_to_dna3"]),
        "positive_control_K": float(b.loc["BOUND", "K_positive_control_cat_vs_channel"]),
        "primer_channel": primer_channel, "dock_only": dock_only, "dock_trap": dock_trap,
        "H_switch_z_bound": float(b.loc["BOUND", "switch_z"]),
        "H_clamp_verdict": verdict_clamp, "H_switch_verdict": verdict_switch,
        "control_catalytic_delta": float(dC), "control_catalytic_delta_nonspecific": float(dC_decoy),
        "validation_auc": float(auc), "domain_mobility": dom,
        "path": [f"{aux['prot_labels'][i][2]}{aux['prot_labels'][i][1]}" for i in path],
        "robustness": robust,
    }
    (out / "anm_summary.json").write_text(json.dumps(summary, indent=2))
    (out / "analysis_log.txt").write_text("\n".join(LOG) + "\n")
    log(f"\nwrote {out}/anm_summary.json, system_comparison.csv, residue_coupling_ranked.csv, "
        f"cutoff_robustness.csv, analysis_log.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
