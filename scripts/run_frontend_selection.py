"""Joint selection of TIA transimpedance and ADC input range."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lidar_sim.adc import ADC, quantize_voltage, sample_waveform
from lidar_sim.config import build_models, load_config
from lidar_sim.detector import APD, Background, background_power_w, multiplied_signal_current_a, shot_noise_rms_a
from lidar_sim.models import received_waveform, round_trip_delay_s
from lidar_sim.ranging import matched_filter_scores, sampled_peak_time_s, time_error_to_range_error_m
from lidar_sim.tia import TIA, clip_output_voltage, first_order_lowpass_transimpedance, total_output_noise_rms_v

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "results" / "data"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
TRANSIMPEDANCES_KOHM = (5.0, 10.0, 20.0)
ADC_POSITIVE_FULL_SCALE_V = (0.25, 0.50, 1.00)
WEAK_RANGES_M = np.arange(50.0, 251.0, 10.0)
STRONG_RANGES_M = np.arange(20.0, 151.0, 10.0)
SCREENING_TRIALS = 200
BASELINE_TRIALS = 500
RANDOM_SEED = 20261107


def build_fixed_components(config):
    laser, target, channel, receiver = build_models(config)
    a = config["apd"]
    apd = APD(a["responsivity_a_w"], a["gain"], a["ionization_ratio"], a["primary_dark_current_na"] * 1e-9, a["noise_bandwidth_mhz"] * 1e6)
    b = config["background"]
    background = Background(b["spectral_radiance_w_m2_sr_nm"], b["receiver_full_fov_mrad"] * 1e-3, b["filter_bandwidth_nm"])
    background_w = background_power_w(background, receiver.aperture_area_m2, receiver.optical_efficiency)
    return laser, target, channel, receiver, apd, background_w


def make_tia(config, transimpedance_kohm):
    t = config["tia"]
    return TIA(transimpedance_kohm * 1e3, t["bandwidth_mhz"] * 1e6, t["input_current_noise_pa_sqrt_hz"] * 1e-12, t["output_voltage_noise_nv_sqrt_hz"] * 1e-9, t["output_swing_v"], t["temperature_k"])


def make_adc(config, positive_full_scale_v):
    d = config["adc"]
    return ADC(d["sample_rate_ghz"] * 1e9, d["resolution_bits"], -0.2 * positive_full_scale_v, positive_full_scale_v)


def analog_return(laser, target, channel, receiver, apd, tia):
    tof_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(tof_s - 20e-9, tof_s + 30e-9, 10e-12)
    optical_w = received_waveform(time_s, laser, target, channel, receiver)
    current_a = multiplied_signal_current_a(optical_w, apd)
    raw_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], tia)
    clipped_v, tia_clip_fraction = clip_output_voltage(raw_v, tia.output_swing_v)
    return time_s, clipped_v, tof_s, float(np.max(optical_w)), tia_clip_fraction


def evaluate_condition(rng, dense_time, analog_v, true_tof_s, noise_rms_v, adc, template_v, calibration_offset_s, trials):
    errors_m = []
    adc_clipped_trials = 0
    threshold = 5.0 * noise_rms_v
    for _ in range(trials):
        phase = float(rng.random())
        sample_time, clean_v = sample_waveform(dense_time, analog_v, adc, phase)
        noisy_v = clean_v + rng.normal(0.0, noise_rms_v, clean_v.size)
        _, quantized_v, clipped_fraction = quantize_voltage(noisy_v, adc)
        if clipped_fraction > 0.0:
            adc_clipped_trials += 1
        scores = matched_filter_scores(quantized_v, template_v)
        estimated_s = sampled_peak_time_s(sample_time, scores, interpolate=True) - calibration_offset_s
        if float(np.max(scores)) >= threshold and abs(estimated_s - true_tof_s) <= 5e-9:
            errors_m.append(time_error_to_range_error_m(estimated_s - true_tof_s))
    errors = np.asarray(errors_m)
    return {
        "detection_rate": errors.size / trials,
        "rmse_m": float(np.sqrt(np.mean(errors**2))) if errors.size else np.nan,
        "adc_clip_rate": adc_clipped_trials / trials,
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, base_target, channel, receiver, apd, background_w = build_fixed_components(config)
    rng = np.random.default_rng(RANDOM_SEED)
    summary = []

    for transimpedance_kohm in TRANSIMPEDANCES_KOHM:
        tia = make_tia(config, transimpedance_kohm)
        for full_scale_v in ADC_POSITIVE_FULL_SCALE_V:
            adc = make_adc(config, full_scale_v)
            base_time, base_v, base_tof, base_peak_optical, _ = analog_return(laser, base_target, channel, receiver, apd, tia)
            base_sample_time, base_samples = sample_waveform(base_time, base_v, adc, 0.0)
            _, base_q, _ = quantize_voltage(base_samples, adc)
            peak_index = int(np.argmax(base_q))
            template_v = base_q[peak_index - 10 : peak_index + 11]
            base_scores = matched_filter_scores(base_q, template_v)
            calibration_offset_s = sampled_peak_time_s(base_sample_time, base_scores, interpolate=True) - base_tof

            weak_results = []
            baseline_result = None
            baseline_noise_rms_v = None
            for distance_m in WEAK_RANGES_M:
                target = replace(base_target, range_m=float(distance_m), reflectivity=0.10)
                dense_time, analog_v, tof_s, peak_optical, _ = analog_return(laser, target, channel, receiver, apd, tia)
                apd_density = shot_noise_rms_a(peak_optical, background_w, replace(apd, noise_bandwidth_hz=1.0))
                noise_rms_v = total_output_noise_rms_v(apd_density, tia)
                trials = BASELINE_TRIALS if distance_m == 50.0 else SCREENING_TRIALS
                result = evaluate_condition(rng, dense_time, analog_v, tof_s, noise_rms_v, adc, template_v, calibration_offset_s, trials)
                weak_results.append((distance_m, result))
                if distance_m == 50.0:
                    baseline_result = result
                    baseline_noise_rms_v = noise_rms_v

            reliable_ranges = [distance for distance, result in weak_results if result["detection_rate"] >= 0.95 and result["rmse_m"] <= 0.05]
            reliable_range = max(reliable_ranges) if reliable_ranges else np.nan
            reliable_boundary_observed = not (reliable_ranges and reliable_range == float(WEAK_RANGES_M[-1]))

            unsaturated_ranges = []
            for distance_m in STRONG_RANGES_M:
                target = replace(base_target, range_m=float(distance_m), reflectivity=0.80)
                dense_time, analog_v, tof_s, peak_optical, tia_clip = analog_return(laser, target, channel, receiver, apd, tia)
                apd_density = shot_noise_rms_a(peak_optical, background_w, replace(apd, noise_bandwidth_hz=1.0))
                noise_rms_v = total_output_noise_rms_v(apd_density, tia)
                result = evaluate_condition(rng, dense_time, analog_v, tof_s, noise_rms_v, adc, template_v, calibration_offset_s, SCREENING_TRIALS)
                if tia_clip == 0.0 and result["adc_clip_rate"] <= 0.01:
                    unsaturated_ranges.append(distance_m)
            minimum_unsaturated_range = min(unsaturated_ranges) if unsaturated_ranges else np.nan

            quantization_rms_v = adc.lsb_v / np.sqrt(12.0)
            summary.append((transimpedance_kohm, full_scale_v, adc.lsb_v, quantization_rms_v, baseline_noise_rms_v, quantization_rms_v / baseline_noise_rms_v, baseline_result["rmse_m"], baseline_result["detection_rate"], reliable_range, reliable_boundary_observed, minimum_unsaturated_range))

    path = DATA_DIR / "frontend_selection_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["transimpedance_kohm", "adc_positive_full_scale_v", "adc_lsb_v", "ideal_quantization_rms_v", "baseline_analog_noise_rms_v", "quantization_to_analog_noise_ratio", "baseline_50m_10pct_rmse_m", "baseline_50m_10pct_detection_rate", "10pct_reliable_range_m", "far_boundary_observed", "80pct_minimum_unsaturated_range_m"])
        writer.writerows(summary)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
    figure, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))
    colors = ("#1368aa", "#238636", "#c73e1d")
    for full_scale_v, color in zip(ADC_POSITIVE_FULL_SCALE_V, colors):
        subset = [row for row in summary if row[1] == full_scale_v]
        gains = [row[0] for row in subset]
        axes[0].plot(gains, [row[8] for row in subset], "o-", color=color, label=f"FS +{full_scale_v:g} V")
        axes[1].plot(gains, [row[10] for row in subset], "o-", color=color)
        axes[2].plot(gains, [100.0 * row[6] for row in subset], "o-", color=color)
    axes[0].set(title="10% target reliable far range", xlabel="TIA transimpedance (kOhm)", ylabel="Reliable range (m)", xticks=list(TRANSIMPEDANCES_KOHM))
    axes[0].legend()
    axes[1].set(title="80% target strong-return limit", xlabel="TIA transimpedance (kOhm)", ylabel="Minimum unsaturated range (m)", xticks=list(TRANSIMPEDANCES_KOHM))
    axes[2].set(title="50 m, 10% ranging precision", xlabel="TIA transimpedance (kOhm)", ylabel="Matched-filter RMSE (cm)", xticks=list(TRANSIMPEDANCES_KOHM))
    figure.suptitle("TIA gain and ADC range joint selection", fontsize=13)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"09_frontend_selection.{suffix}")
    plt.close(figure)

    print(f"Frontend selection: screening={SCREENING_TRIALS}, baseline={BASELINE_TRIALS}, seed={RANDOM_SEED}")
    print("Rf(kOhm)  ADC +FS(V)  LSB(mV)  Q/analog(%)  50m RMSE(cm)  10% far(m)  80% min unsat(m)")
    for row in summary:
        far_label = f">={row[8]:.0f}" if not row[9] else f"{row[8]:.0f}"
        unsat_label = "none" if np.isnan(row[10]) else f"{row[10]:.0f}"
        print(f"{row[0]:8.1f}  {row[1]:10.2f}  {row[2] * 1e3:7.3f}  {row[5] * 100:11.3f}  {row[6] * 100:12.3f}  {far_label:>10s}  {unsat_label:>16s}")
    print("Saved frontend_selection_summary.csv and figure 09_frontend_selection (PNG/PDF).")


if __name__ == "__main__":
    main()
