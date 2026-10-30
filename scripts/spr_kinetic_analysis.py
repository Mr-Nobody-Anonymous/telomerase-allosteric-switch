"""SPR Kinetic Analysis Pipeline for Telomerase Holoenzyme - TPP1 Binding.

Fits 1:1 Langmuir binding kinetics to multi-cycle or single-cycle SPR sensorgrams,
calculates association rate (k_on), dissociation rate (k_off), and Kd,
and evaluates pre-registered acceptance and falsification thresholds.
Supports batch evaluation of experimental CSV datasets via CLI.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def langmuir_1_to_1_association(t, k_on, k_off, r_max, conc):
    """Association phase sensorgram: R(t) = R_eq * (1 - exp(-(k_on*conc + k_off)*t))."""
    k_obs = k_on * conc + k_off
    r_eq = (r_max * conc) / (conc + (k_off / k_on))
    return r_eq * (1.0 - np.exp(-k_obs * t))


def langmuir_1_to_1_dissociation(t, r_0, k_off):
    """Dissociation phase sensorgram: R(t) = R_0 * exp(-k_off * (t - t_0))."""
    return r_0 * np.exp(-k_off * t)


def fit_spr_cycle(time_assoc, resp_assoc, time_dissoc, resp_dissoc, conc):
    """Fit a single SPR cycle to extract k_on, k_off, and calculate Kd."""
    r_0_guess = resp_dissoc[0] if len(resp_dissoc) > 0 else 50.0
    popt_dissoc, _ = curve_fit(
        langmuir_1_to_1_dissociation,
        time_dissoc,
        resp_dissoc,
        p0=[r_0_guess, 0.005],
        bounds=(0, [1000, 1.0]),
    )
    k_off_fit = popt_dissoc[1]

    def assoc_fixed_koff(t, k_on, r_max):
        return langmuir_1_to_1_association(t, k_on, k_off_fit, r_max, conc)

    popt_assoc, _ = curve_fit(
        assoc_fixed_koff,
        time_assoc,
        resp_assoc,
        p0=[5e4, 60.0],
        bounds=([1e2, 5.0], [1e7, 500.0]),
    )
    k_on_fit = popt_assoc[0]
    r_max_fit = popt_assoc[1]

    kd_molar = k_off_fit / k_on_fit
    kd_nm = kd_molar * 1e9

    pred_assoc = assoc_fixed_koff(time_assoc, k_on_fit, r_max_fit)
    chi2 = float(np.mean((resp_assoc - pred_assoc) ** 2))

    return {
        "k_on": float(k_on_fit),
        "k_off": float(k_off_fit),
        "Kd_nM": float(kd_nm),
        "R_max": float(r_max_fit),
        "Chi2": chi2,
    }


def evaluate_spr_candidate(construct_name, kd_nm, chi2, r_max):
    """Evaluate against pre-registered falsification criteria."""
    qc_pass = bool(chi2 <= (0.10 * r_max))

    if construct_name in ["WT", "TERT_WT"]:
        status = "BASELINE VALID" if (100.0 <= kd_nm <= 250.0) else "BASELINE ANOMALY"
    elif "757A" in construct_name:
        if kd_nm <= 300.0:
            status = "CRITERION 2A PASSED (Nanomolar Binding Preserved)"
        elif kd_nm > 1500.0:
            status = "FALSIFIED (Recruitment Disrupted)"
        else:
            status = "INTERMEDIATE ATTENUATION"
    elif "757R" in construct_name:
        status = "RESCUE CONTROL CONFIRMED" if kd_nm <= 300.0 else "UNEXPECTED BINDING LOSS"
    elif "169A" in construct_name:
        status = "RECRUITMENT DEFECT CONFIRMED" if kd_nm > 1500.0 else "UNEXPECTED BINDING INTACT"
    else:
        status = "UNKNOWN"

    return qc_pass, status


def process_spr_benchmark_file(csv_path: str | Path) -> pd.DataFrame:
    """Process an independent SPR benchmark CSV containing multi-construct sensorgrams."""
    df = pd.read_csv(csv_path)
    required_cols = {"construct", "replicate", "phase", "time_sec", "response_RU", "conc_nM"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in SPR dataset: {required_cols - set(df.columns)}")

    results = []
    groups = df.groupby(["construct", "replicate"])

    for (construct, rep), g in groups:
        assoc_data = g[g["phase"] == "association"]
        dissoc_data = g[g["phase"] == "dissociation"]

        t_a = assoc_data["time_sec"].to_numpy()
        r_a = assoc_data["response_RU"].to_numpy()
        t_d = dissoc_data["time_sec"].to_numpy()
        r_d = dissoc_data["response_RU"].to_numpy()
        conc_molar = float(assoc_data["conc_nM"].iloc[0]) * 1e-9

        fit = fit_spr_cycle(t_a, r_a, t_d, r_d, conc_molar)
        qc_pass, status = evaluate_spr_candidate(construct, fit["Kd_nM"], fit["Chi2"], fit["R_max"])

        results.append({
            "construct": construct,
            "replicate": rep,
            "k_on_M_s": fit["k_on"],
            "k_off_s": fit["k_off"],
            "Kd_nM": fit["Kd_nM"],
            "R_max_RU": fit["R_max"],
            "Chi2": fit["Chi2"],
            "QC_Pass": qc_pass,
            "Evaluation_Verdict": status,
        })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="SPR Kinetic Analysis Fitting Pipeline")
    parser.add_argument("--input", type=str, default="data/benchmarks/spr_independent_benchmark.csv", help="Input CSV path")
    parser.add_argument("--output", type=str, default="results/experimental_validation/spr_benchmark_results.csv", help="Output CSV path")
    args = parser.parse_args()

    in_path = Path(args.input)
    if in_path.exists():
        print(f"Processing SPR dataset: {in_path}")
        res_df = process_spr_benchmark_file(in_path)
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        res_df.to_csv(out_path, index=False)
        print(f"Results written to: {out_path}")
        print("\n--- Summary Table ---")
        print(res_df[["construct", "replicate", "Kd_nM", "Chi2", "QC_Pass", "Evaluation_Verdict"]].to_string(index=False))
    else:
        print(f"File not found: {in_path}. Running internal demonstration.")
        t_a = np.linspace(0, 180, 90)
        t_d = np.linspace(0, 600, 120)
        sim_wt_a = langmuir_1_to_1_association(t_a, 4.2e4, 0.0063, 65.0, 240e-9) + np.random.normal(0, 0.3, len(t_a))
        sim_wt_d = langmuir_1_to_1_dissociation(t_d, sim_wt_a[-1], 0.0063) + np.random.normal(0, 0.3, len(t_d))
        fit_wt = fit_spr_cycle(t_a, sim_wt_a, t_d, sim_wt_d, 240e-9)
        _qc_wt, stat_wt = evaluate_spr_candidate("TERT_WT", fit_wt["Kd_nM"], fit_wt["Chi2"], fit_wt["R_max"])
        print(f"TERT WT: Kd={fit_wt['Kd_nM']:.1f} nM -> {stat_wt}")


if __name__ == "__main__":
    main()
