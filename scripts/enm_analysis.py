import gemmi
import numpy as np
from scipy import linalg

cif_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data/raw"

def get_ca_coords(structure, chain_id, min_res=None, max_res=None):
    model = structure[0]
    chain = model[chain_id]
    coords = []
    res_nums = []
    res_names = []
    b_factors = []
    for r in chain:
        if min_res and r.seqid.num < min_res:
            continue
        if max_res and r.seqid.num > max_res:
            continue
        ca = next((a for a in r if a.name == 'CA'), None)
        if ca:
            coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
            res_nums.append(r.seqid.num)
            res_names.append(r.name)
            b_factors.append(ca.b_iso)
    return np.array(coords), res_nums, res_names, np.array(b_factors)

print("Loading 7TRE and 7TRD...")
st_e = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/7TRE.cif").sole_block())
st_d = gemmi.make_structure_from_block(gemmi.cif.read_file(f"{cif_dir}/7TRD.cif").sole_block())

coords_e_tert, res_e, names_e, b_e = get_ca_coords(st_e, "A")
coords_d_tert, res_d, names_d, b_d = get_ca_coords(st_d, "A")
coords_e_tpp1, res_tpp1, names_tpp1, b_tpp1 = get_ca_coords(st_e, "C")

print(f"7TRE TERT Calpha: {len(coords_e_tert)}, TPP1 Calpha: {len(coords_e_tpp1)}")
print(f"7TRD TERT Calpha: {len(coords_d_tert)}")

# Ensure common residues for 1-to-1 comparison of TERT
common_res = sorted(set(res_e).intersection(set(res_d)))
idx_e = [res_e.index(r) for r in common_res]
idx_d = [res_d.index(r) for r in common_res]

N = len(common_res)
print(f"Common TERT residues: {N}")

coords_e_common = coords_e_tert[idx_e]
coords_d_common = coords_d_tert[idx_d]

def compute_gnm(coords, cutoff=10.0):
    N = coords.shape[0]
    # Pairwise distance matrix
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=-1))
    
    # Kirchhoff matrix
    gamma = np.zeros((N, N))
    contact = (dist <= cutoff) & (dist > 0.0)
    gamma[contact] = -1.0
    np.fill_diagonal(gamma, -np.sum(gamma, axis=1))
    
    # Eigendecomposition
    w, v = linalg.eigh(gamma)
    # Exclude 0th mode (rigid translation)
    # Covariance / cross-correlation
    inv_gamma = np.zeros((N, N))
    for k in range(1, N):
        inv_gamma += (1.0 / w[k]) * np.outer(v[:, k], v[:, k])
    
    # Normalize cross-correlation matrix
    diag = np.diag(inv_gamma)
    cross_corr = inv_gamma / np.sqrt(np.outer(diag, diag))
    return w, v, inv_gamma, cross_corr

print("\nComputing GNM for Unbound TERT (7TRD)...")
w_d, v_d, cov_d, cc_d = compute_gnm(coords_d_common, cutoff=10.0)

print("Computing GNM for Bound Complex (7TRE TERT + TPP1)...")
coords_complex = np.vstack([coords_e_common, coords_e_tpp1])
w_c, v_c, cov_c, cc_c = compute_gnm(coords_complex, cutoff=10.0)
cov_e = cov_c[:N, :N]
cc_e = cc_c[:N, :N]

# Identify indices of TEN domain (residues <= 195) and TRAP domain (residues 750-795)
ten_idx = [i for i, r in enumerate(common_res) if r <= 195]
trap_idx = [i for i, r in enumerate(common_res) if 750 <= r <= 795]

print(f"TEN domain residues: {len(ten_idx)} (Calpha index range {ten_idx[0]} to {ten_idx[-1]})")
print(f"TRAP domain residues: {len(trap_idx)} (Calpha index range {trap_idx[0]} to {trap_idx[-1]})")

# Cross-domain correlation between TEN and TRAP
mean_cc_unbound = np.mean(cc_d[np.ix_(ten_idx, trap_idx)])
mean_cc_bound = np.mean(cc_e[np.ix_(ten_idx, trap_idx)])

print("\n--- RESULTS OF ELASTIC NETWORK DYNAMICS ---")
print("Mean Dynamic Cross-Correlation between TEN and TRAP:")
print(f"  Unbound (7TRD): {mean_cc_unbound:.4f}")
print(f"  Bound (7TRE):   {mean_cc_bound:.4f}")
print(f"  Delta Coupling: +{(mean_cc_bound - mean_cc_unbound):.4f} ({(mean_cc_bound/mean_cc_unbound - 1)*100:+.1f}%)")

# Residue-specific correlations for Tyr176, Lys757, Gln794
i_y176 = common_res.index(176)
i_k757 = common_res.index(757)
i_q794 = common_res.index(794)

print("\nSpecific Pairwise Dynamic Cross-Correlations:")
print(f"  Tyr176 <---> Lys757: Unbound={cc_d[i_y176, i_k757]:.4f} | Bound={cc_e[i_y176, i_k757]:.4f} (Delta={cc_e[i_y176, i_k757]-cc_d[i_y176, i_k757]:+.4f})")
print(f"  Tyr176 <---> Gln794: Unbound={cc_d[i_y176, i_q794]:.4f} | Bound={cc_e[i_y176, i_q794]:.4f} (Delta={cc_e[i_y176, i_q794]-cc_d[i_y176, i_q794]:+.4f})")
print(f"  Lys757 <---> Gln794: Unbound={cc_d[i_k757, i_q794]:.4f} | Bound={cc_e[i_k757, i_q794]:.4f} (Delta={cc_e[i_k757, i_q794]-cc_d[i_k757, i_q794]:+.4f})")

# Mean Square Fluctuations (MSF) comparison
msf_d_trap = np.mean(np.diag(cov_d)[trap_idx])
msf_e_trap = np.mean(np.diag(cov_e)[trap_idx])
print("\nMean Square Fluctuations of TRAP Domain:")
print(f"  Unbound TRAP MSF: {msf_d_trap:.4f}")
print(f"  Bound TRAP MSF:   {msf_e_trap:.4f}")
print(f"  Suppression of TRAP breathing: {(1 - msf_e_trap/msf_d_trap)*100:.1f}% reduction in thermal variance!")
