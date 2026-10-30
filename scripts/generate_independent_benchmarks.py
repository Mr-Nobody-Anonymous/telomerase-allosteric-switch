"""Generate Independent Benchmark Datasets for SPR and RAP Automated Fitting Pipelines.

Outputs:
1. data/benchmarks/spr_independent_benchmark.csv:
   Multi-cycle SPR sensorgrams for WT, K757A, K757R, and E169A with experimental noise and baseline drift.
2. data/benchmarks/rap_independent_benchmark.csv:
   Multi-repeat direct extension gel band densitometries for WT, K757A, K757R, and D712A across minus/plus TPP1.
"""

from pathlib import Path

import numpy as np
import pandas as pd

BENCHMARK_DIR = Path("data/benchmarks")
BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. SPR Benchmark Dataset Generation
# -----------------------------------------------------------------------------
def generate_spr_benchmark():
    np.random.seed(42)
    time_assoc = np.linspace(0, 180, 45) # 45 pts during association
    time_dissoc = np.linspace(0, 600, 60) # 60 pts during dissociation
    
    # Constructs to benchmark
    constructs = [
        # (name, replicate, k_on, k_off, r_max)
        ("TERT_WT", 1, 4.3e4, 0.0062, 65.0), # Kd = 144 nM
        ("TERT_WT", 2, 4.1e4, 0.0065, 62.0), # Kd = 158 nM
        ("TERT_WT", 3, 4.4e4, 0.0061, 67.0), # Kd = 139 nM
        ("TERT_Lys757A", 1, 3.9e4, 0.0068, 64.0), # Kd = 174 nM
        ("TERT_Lys757A", 2, 3.8e4, 0.0071, 60.0), # Kd = 187 nM
        ("TERT_Lys757A", 3, 4.0e4, 0.0067, 63.0), # Kd = 168 nM
        ("TERT_Lys757R", 1, 4.2e4, 0.0064, 66.0), # Kd = 152 nM
        ("TERT_Lys757R", 2, 4.1e4, 0.0066, 64.0), # Kd = 161 nM
        ("TERT_Glu169A_neg", 1, 1.2e3, 0.0450, 25.0), # Kd = 37,500 nM (Disrupted)
    ]
    
    conc_nM = 240.0
    conc_M = conc_nM * 1e-9
    
    records = []
    for c_name, rep, k_on, k_off, r_max in constructs:
        # Association
        k_obs = k_on * conc_M + k_off
        r_eq = (r_max * conc_M) / (conc_M + (k_off / k_on))
        resp_a = r_eq * (1.0 - np.exp(-k_obs * time_assoc))
        resp_a += np.random.normal(0, 0.35, len(time_assoc)) # Noise
        resp_a += 0.002 * time_assoc # Baseline drift
        
        # Dissociation
        r0 = resp_a[-1]
        resp_d = r0 * np.exp(-k_off * time_dissoc)
        resp_d += np.random.normal(0, 0.35, len(time_dissoc))
        
        for t, r in zip(time_assoc, resp_a):
            records.append({
                "construct": c_name,
                "replicate": rep,
                "phase": "association",
                "time_sec": round(t, 2),
                "response_RU": round(r, 3),
                "conc_nM": conc_nM
            })
            
        for t, r in zip(time_dissoc, resp_d):
            records.append({
                "construct": c_name,
                "replicate": rep,
                "phase": "dissociation",
                "time_sec": round(t, 2),
                "response_RU": round(r, 3),
                "conc_nM": conc_nM
            })
            
    df_spr = pd.DataFrame(records)
    out_path = BENCHMARK_DIR / "spr_independent_benchmark.csv"
    df_spr.to_csv(out_path, index=False)
    print(f"Generated SPR benchmark: {out_path} ({len(df_spr)} points across {len(constructs)} runs)")

# -----------------------------------------------------------------------------
# 2. RAP Benchmark Dataset Generation
# -----------------------------------------------------------------------------
def generate_rap_benchmark():
    np.random.seed(123)
    repeats = 16 # 16 telomeric repeats
    rep_indices = np.arange(1, repeats + 1)
    
    # Profiles: (construct, rep, minus_k, plus_k, basal_scale)
    # k is decay rate: P_overall = exp(-k).
    # Basal (-TPP1): k ~ 1.05 -> P_overall ~ 0.35
    # Stimulated (+TPP1 for WT/Rescue): k ~ 0.16 -> P_overall ~ 0.85 -> chi_RAP ~ 2.4
    # Uncoupled (K757A): minus k ~ 1.05, plus k ~ 0.90 -> P_overall ~ 0.41 -> chi_RAP ~ 1.17 (Uncoupled!)
    # Catalytic Dead (D712A): basal_scale ~ 0.04
    configs = [
        ("TERT_WT", 1, 1.05, 0.16, 1.00),
        ("TERT_WT", 2, 1.08, 0.15, 0.98),
        ("TERT_WT", 3, 1.03, 0.17, 1.02),
        ("TERT_Lys757A", 1, 1.05, 0.90, 0.95),
        ("TERT_Lys757A", 2, 1.04, 0.92, 0.93),
        ("TERT_Lys757A", 3, 1.06, 0.88, 0.96),
        ("TERT_Lys757R", 1, 1.05, 0.16, 0.98),
        ("TERT_Lys757R", 2, 1.07, 0.17, 0.99),
        ("TERT_Asp712A_dead", 1, 3.50, 3.50, 0.04),
    ]

    records = []
    bg = 50.0  # phosphorimager background

    for c_name, rep, k_minus, k_plus, basal_scale in configs:
        # minus TPP1
        raw_minus = basal_scale * 1000.0 * np.exp(-k_minus * (rep_indices - 1)) * (3.0 * rep_indices) + bg
        raw_minus += np.random.normal(0, 1.0, repeats)
        raw_minus = np.maximum(bg, raw_minus)

        # plus TPP1
        raw_plus = basal_scale * 1000.0 * np.exp(-k_plus * (rep_indices - 1)) * (3.0 * rep_indices) + bg
        raw_plus += np.random.normal(0, 1.0, repeats)
        raw_plus = np.maximum(bg, raw_plus)
        
        for r_idx, v_m, v_p in zip(rep_indices, raw_minus, raw_plus):
            records.append({
                "construct": c_name,
                "replicate": rep,
                "repeat_number": r_idx,
                "raw_intensity_minus_tpp1": round(v_m, 1),
                "raw_intensity_plus_tpp1": round(v_p, 1),
                "background_intensity": bg,
                "internal_control_recovery": 1.0
            })
            
    df_rap = pd.DataFrame(records)
    out_path = BENCHMARK_DIR / "rap_independent_benchmark.csv"
    df_rap.to_csv(out_path, index=False)
    print(f"Generated RAP benchmark: {out_path} ({len(df_rap)} rows across {len(configs)} gel lanes)")

if __name__ == "__main__":
    generate_spr_benchmark()
    generate_rap_benchmark()
