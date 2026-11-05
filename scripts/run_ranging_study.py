"""Monte Carlo comparison of bottom-level dToF ranging estimators."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lidar_sim.adc import ADC, quantize_voltage, sample_waveform
from lidar_sim.config import build_models, load_config
from lidar_sim.detector import APD, Background, background_power_w, multiplied_signal_current_a, shot_noise_rms_a
from lidar_sim.models import received_waveform, round_trip_delay_s
from lidar_sim.ranging import first_threshold_crossing_s, matched_filter_scores, sampled_peak_time_s, time_error_to_range_error_m
from lidar_sim.tia import TIA, first_order_lowpass_transimpedance, total_output_noise_rms_v

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "results" / "data"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
TRIALS = 2000
RANDOM_SEED = 20261105


def build_chain(config):
    laser, target, channel, receiver = build_models(config)
    a = config["apd"]
    apd = APD(a["responsivity_a_w"], a["gain"], a["ionization_ratio"], a["primary_dark_current_na"] * 1e-9, a["noise_bandwidth_mhz"] * 1e6)
    t = config["tia"]
    tia = TIA(t["transimpedance_kohm"] * 1e3, t["bandwidth_mhz"] * 1e6, t["input_current_noise_pa_sqrt_hz"] * 1e-12, t["output_voltage_noise_nv_sqrt_hz"] * 1e-9, t["output_swing_v"], t["temperature_k"])
    d = config["adc"]
    adc = ADC(d["sample_rate_ghz"] * 1e9, d["resolution_bits"], d["input_min_v"], d["input_max_v"])
    b = config["background"]
    background = Background(b["spectral_radiance_w_m2_sr_nm"], b["receiver_full_fov_mrad"] * 1e-3, b["filter_bandwidth_nm"])
    background_w = background_power_w(background, receiver.aperture_area_m2, receiver.optical_efficiency)
    return laser, target, channel, receiver, apd, tia, adc, background_w


def analog_waveform(laser, target, channel, receiver, apd, tia):
    tof_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(tof_s - 20e-9, tof_s + 30e-9, 10e-12)
    optical_w = received_waveform(time_s, laser, target, channel, receiver)
    current_a = multiplied_signal_current_a(optical_w, apd)
    voltage_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], tia)
    return time_s, voltage_v, tof_s, float(np.max(optical_w))


def quadratic_time(time_s, values):
    return sampled_peak_time_s(time_s, values, interpolate=True)


def estimate_all(sample_time_s, sample_v, template_v, threshold_v):
    threshold_time = first_threshold_crossing_s(sample_time_s, sample_v, threshold_v)
    integer_time = sampled_peak_time_s(sample_time_s, sample_v, interpolate=False)
    quadratic_peak_time = quadratic_time(sample_time_s, sample_v)
    correlation = matched_filter_scores(sample_v, template_v)
    matched_time = quadratic_time(sample_time_s, correlation)
    return threshold_time, integer_time, quadratic_peak_time, matched_time, float(np.max(correlation))


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, target, channel, receiver, apd, tia, adc, background_w = build_chain(config)
    dense_time_s, analog_v, true_tof_s, peak_optical_w = analog_waveform(laser, target, channel, receiver, apd, tia)

    apd_noise_density = shot_noise_rms_a(peak_optical_w, background_w, APD(apd.responsivity_a_w, apd.gain, apd.ionization_ratio, apd.primary_dark_current_a, 1.0))
    noise_rms_v = total_output_noise_rms_v(apd_noise_density, tia)
    threshold_v = 5.0 * noise_rms_v

    reference_time_s, reference_v = sample_waveform(dense_time_s, analog_v, adc, 0.0)
    _, reference_q, _ = quantize_voltage(reference_v, adc)
    peak_index = int(np.argmax(reference_q))
    half_width = 10
    template_v = reference_q[peak_index - half_width : peak_index + half_width + 1]
    if template_v.size != 2 * half_width + 1:
        raise RuntimeError("Reference waveform is too short for matched-filter template")

    reference_estimates = estimate_all(reference_time_s, reference_q, template_v, threshold_v)
    calibration_offsets = []
    for estimate in reference_estimates[:4]:
        if estimate is None:
            raise RuntimeError("Noiseless baseline calibration failed")
        calibration_offsets.append(estimate - true_tof_s)
    matched_detection_threshold = 5.0 * noise_rms_v

    rng = np.random.default_rng(RANDOM_SEED)
    algorithm_names = ("Fixed threshold", "Integer peak", "Quadratic peak", "Matched filter")
    errors_m = {name: [] for name in algorithm_names}
    detections = {name: 0 for name in algorithm_names}

    for _ in range(TRIALS):
        phase = float(rng.random())
        sample_time_s, clean_samples_v = sample_waveform(dense_time_s, analog_v, adc, phase)
        noisy_samples_v = clean_samples_v + rng.normal(0.0, noise_rms_v, clean_samples_v.size)
        _, quantized_v, _ = quantize_voltage(noisy_samples_v, adc)
        estimates = estimate_all(sample_time_s, quantized_v, template_v, threshold_v)
        matched_score = estimates[4]
        for index, name in enumerate(algorithm_names):
            estimate = estimates[index]
            detected = estimate is not None
            if name == "Matched filter":
                detected = detected and matched_score >= matched_detection_threshold
            if not detected:
                continue
            detections[name] += 1
            corrected_time_s = estimate - calibration_offsets[index]
            errors_m[name].append(time_error_to_range_error_m(corrected_time_s - true_tof_s))

    summary_rows = []
    for name in algorithm_names:
        values = np.asarray(errors_m[name])
        detection_rate = detections[name] / TRIALS
        bias_m = float(np.mean(values)) if values.size else np.nan
        rmse_m = float(np.sqrt(np.mean(values**2))) if values.size else np.nan
        std_m = float(np.std(values)) if values.size else np.nan
        p95_abs_m = float(np.percentile(np.abs(values), 95)) if values.size else np.nan
        summary_rows.append((name, detection_rate, bias_m, rmse_m, std_m, p95_abs_m))

    with (DATA_DIR / "ranging_monte_carlo_summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["algorithm", "detection_rate", "bias_m", "rmse_m", "std_m", "p95_absolute_error_m"])
        writer.writerows(summary_rows)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    colors = ("#1368aa", "#777777", "#238636", "#c73e1d")
    for name, color in zip(algorithm_names, colors):
        values_cm = np.asarray(errors_m[name]) * 100.0
        if values_cm.size:
            axes[0].hist(values_cm, bins=45, density=True, histtype="step", linewidth=1.7, color=color, label=name)
    axes[0].set(title=f"Range error at {target.range_m:.0f} m, {target.reflectivity:.0%} reflectivity", xlabel="Range error (cm)", ylabel="Probability density")
    axes[0].legend(fontsize=8)

    rmse_cm = np.array([row[3] for row in summary_rows]) * 100.0
    bias_cm = np.array([row[2] for row in summary_rows]) * 100.0
    x = np.arange(len(algorithm_names))
    width = 0.36
    axes[1].bar(x - width / 2, rmse_cm, width, color="#1368aa", label="RMSE")
    axes[1].bar(x + width / 2, np.abs(bias_cm), width, color="#f0a202", label="Absolute bias")
    axes[1].set(title=f"Monte Carlo performance ({TRIALS} trials)", ylabel="Range error (cm)", xticks=x, xticklabels=("Threshold", "Integer", "Quadratic", "Matched"))
    axes[1].legend()
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"06_ranging_monte_carlo.{suffix}")
    plt.close(figure)

    print(f"Ranging Monte Carlo: {TRIALS} trials, seed {RANDOM_SEED}")
    print(f"Baseline: {target.range_m:.1f} m, reflectivity {target.reflectivity:.0%}, noise RMS {noise_rms_v * 1e3:.3f} mV")
    print(f"Fixed threshold: {threshold_v * 1e3:.3f} mV (5 sigma)")
    print("Algorithm          Detection    Bias (cm)   RMSE (cm)    Std (cm)  P95 abs (cm)")
    for row in summary_rows:
        print(f"{row[0]:17s} {row[1]:9.3%} {row[2] * 100:12.3f} {row[3] * 100:11.3f} {row[4] * 100:11.3f} {row[5] * 100:13.3f}")
    print("Saved ranging_monte_carlo_summary.csv and figure 06_ranging_monte_carlo (PNG/PDF).")


if __name__ == "__main__":
    main()
