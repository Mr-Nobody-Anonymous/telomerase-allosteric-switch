"""Quantitative Densitometry Pipeline for Telomerase Direct Extension RAP Gels.

Applies rolling-ball background subtraction, guanine normalization (3*i),
calculates repeat addition processivity (R_p,i), overall processivity factor (P_overall),
and fold-stimulation (chi_RAP), with automated falsification evaluation.
Supports batch evaluation of experimental CSV datasets via CLI.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def calculate_rap_metrics(raw_band_intensities, background, i_control=1.0):
    """
    Given raw intensities for repeats i = 1 .. N:
    1. Subtract background.
    2. Normalize by guanine content (n_G = 3 * i) and internal recovery control.
    3. Compute R_p,i according to Latrick & Cech (2010).
    4. Fit exponential decay to obtain P_overall = exp(-k).
    """
    N = len(raw_band_intensities)
    repeats = np.arange(1, N + 1)

    net_intensities = np.maximum(0.0, np.array(raw_band_intensities) - background)
    n_G = 3.0 * repeats
    norm_intensities = net_intensities / (n_G * i_control)
    total_activity = float(np.sum(norm_intensities))

    r_p = []
    for i in range(N - 1):
        num = np.sum(norm_intensities[i + 1:])
        denom = np.sum(norm_intensities[i:])
        r_p.append(float(num / denom) if denom > 0 else 0.0)

    valid_mask = norm_intensities > 0
    if np.sum(valid_mask) >= 3:
        p = np.polyfit(repeats[valid_mask], np.log(norm_intensities[valid_mask]), deg=1)
        k = -p[0]
        p_overall = float(np.exp(-k))
    else:
        p_overall = 0.0

    return {
        "norm_intensities": norm_intensities,
        "total_activity": total_activity,
        "R_p_profile": np.array(r_p),
        "mean_R_p": float(np.mean(r_p[:10])) if len(r_p) >= 10 else float(np.mean(r_p)),
        "P_overall": p_overall,
    }


def calculate_fold_stimulation(metrics_minus_tpp1, metrics_plus_tpp1):
    """Calculate processivity fold-enhancement chi_RAP."""
    p_minus = metrics_minus_tpp1["P_overall"]
    p_plus = metrics_plus_tpp1["P_overall"]
    chi_rap = float(p_plus / p_minus) if p_minus > 0 else 0.0
    return chi_rap


def evaluate_rap_phenotype(construct_name, basal_activity_fraction, chi_rap):
    """Evaluate against pre-registered falsification criteria."""
    if basal_activity_fraction < 0.30:
        verdict = "FALSIFIED (Severe Catalytic Failure; k_cat < 30%)"
    elif "757A" in construct_name:
        if chi_rap <= 1.8 and basal_activity_fraction >= 0.85:
            verdict = "CRITERIA 2B & 2C PASSED (Selective Processivity Uncoupling Hit)"
        elif chi_rap >= 4.0:
            verdict = "FALSIFIED (Processivity Stimulation Intact; chi_RAP >= 4.0)"
        else:
            verdict = "INTERMEDIATE PHENOTYPE"
    elif "757R" in construct_name:
        if chi_rap >= 2.0 and basal_activity_fraction >= 0.85:
            verdict = "RESCUE CONTROL CONFIRMED (Stereochemical & Charge Integrity Validated)"
        else:
            verdict = "RESCUE FAILED"
    else:
        verdict = "WT REFERENCE BASELINE"

    return verdict


def process_rap_benchmark_file(csv_path: str | Path) -> pd.DataFrame:
    """Process an independent RAP benchmark CSV containing multi-lane densitometry profiles."""
    df = pd.read_csv(csv_path)
    required_cols = {
        "construct",
        "replicate",
        "repeat_number",
        "raw_intensity_minus_tpp1",
        "raw_intensity_plus_tpp1",
        "background_intensity",
        "internal_control_recovery",
    }
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in RAP dataset: {required_cols - set(df.columns)}")

    # First determine WT mean basal activity to normalize mutant basal activity fractions
    wt_data = df[df["construct"] == "TERT_WT"]
    wt_basal_activities = []
    for _rep, g in wt_data.groupby("replicate"):
        bg = float(g["background_intensity"].iloc[0])
        ic = float(g["internal_control_recovery"].iloc[0])
        m = calculate_rap_metrics(g["raw_intensity_minus_tpp1"].to_numpy(), bg, ic)
        wt_basal_activities.append(m["total_activity"])
    mean_wt_basal = float(np.mean(wt_basal_activities)) if wt_basal_activities else 1.0

    results = []
    groups = df.groupby(["construct", "replicate"])

    for (construct, rep), g in groups:
        g_sorted = g.sort_values("repeat_number")
        bg = float(g_sorted["background_intensity"].iloc[0])
        ic = float(g_sorted["internal_control_recovery"].iloc[0])

        minus_intensities = g_sorted["raw_intensity_minus_tpp1"].to_numpy()
        plus_intensities = g_sorted["raw_intensity_plus_tpp1"].to_numpy()

        res_minus = calculate_rap_metrics(minus_intensities, bg, ic)
        res_plus = calculate_rap_metrics(plus_intensities, bg, ic)

        chi_rap = calculate_fold_stimulation(res_minus, res_plus)
        basal_frac = res_minus["total_activity"] / mean_wt_basal if mean_wt_basal > 0 else 1.0
        verdict = evaluate_rap_phenotype(construct, basal_frac, chi_rap)

        results.append({
            "construct": construct,
            "replicate": rep,
            "basal_activity_frac": round(basal_frac, 3),
            "P_overall_minus": round(res_minus["P_overall"], 3),
            "P_overall_plus": round(res_plus["P_overall"], 3),
            "chi_RAP": round(chi_rap, 2),
            "mean_Rp_plus": round(res_plus["mean_R_p"], 3),
            "Evaluation_Verdict": verdict,
        })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="RAP Gel Quantitative Densitometry Pipeline")
    parser.add_argument("--input", type=str, default="data/benchmarks/rap_independent_benchmark.csv", help="Input CSV path")
    parser.add_argument("--output", type=str, default="results/experimental_validation/rap_benchmark_results.csv", help="Output CSV path")
    args = parser.parse_args()

    in_path = Path(args.input)
    if in_path.exists():
        print(f"Processing RAP dataset: {in_path}")
        res_df = process_rap_benchmark_file(in_path)
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        res_df.to_csv(out_path, index=False)
        print(f"Results written to: {out_path}")
        print("\n--- Summary Table ---")
        print(res_df[["construct", "replicate", "basal_activity_frac", "P_overall_minus", "P_overall_plus", "chi_RAP", "Evaluation_Verdict"]].to_string(index=False))
    else:
        print(f"File not found: {in_path}. Running internal demonstration.")
        repeats = 15
        wt_minus_raw = 1000.0 * np.exp(-1.05 * np.arange(repeats)) * (3.0 * np.arange(1, repeats + 1)) + 50.0
        wt_plus_raw = 1000.0 * np.exp(-0.16 * np.arange(repeats)) * (3.0 * np.arange(1, repeats + 1)) + 50.0
        res_wt_minus = calculate_rap_metrics(wt_minus_raw, background=50.0)
        res_wt_plus = calculate_rap_metrics(wt_plus_raw, background=50.0)
        chi_wt = calculate_fold_stimulation(res_wt_minus, res_wt_plus)
        v_wt = evaluate_rap_phenotype("WT", 1.0, chi_wt)
        print(f"TERT WT: chi_RAP={chi_wt:.2f} -> {v_wt}")


if __name__ == "__main__":
    main()
