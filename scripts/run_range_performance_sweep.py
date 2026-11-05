"""Sweep matched-filter ranging performance over range and reflectivity."""

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
RANGES_M = np.arange(20.0, 251.0, 10.0)
REFLECTIVITIES = (0.10, 0.20, 0.50, 0.80)
TRIALS_PER_CONDITION = 500
RANDOM_SEED = 20261106
MINIMUM_DETECTION_RATE = 0.95
MAXIMUM_RMSE_M = 0.05
MAXIMUM_VALID_TIME_ERROR_S = 5e-9


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


def make_analog_waveform(laser, target, channel, receiver, apd, tia):
    tof_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(tof_s - 20e-9, tof_s + 30e-9, 10e-12)
    optical_w = received_waveform(time_s, laser, target, channel, receiver)
    current_a = multiplied_signal_current_a(optical_w, apd)
    raw_voltage_v = first_order_lowpass_transimpedance(
        current_a, time_s[1] - time_s[0], tia
    )
    voltage_v, clipped_sample_fraction = clip_output_voltage(
        raw_voltage_v, tia.output_swing_v
    )
    return (
        time_s,
        voltage_v,
        tof_s,
        float(np.max(optical_w)),
        clipped_sample_fraction,
    )


def reliable_range_m(rows, reflectivity):
    valid_ranges = [
        row[0]
        for row in rows
        if row[1] == reflectivity
        and row[4] >= MINIMUM_DETECTION_RATE
        and row[6] <= MAXIMUM_RMSE_M
    ]
    return max(valid_ranges) if valid_ranges else np.nan


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, baseline_target, channel, receiver, apd, tia, adc, background_w = build_chain(config)

    (
        baseline_time,
        baseline_analog,
        baseline_tof,
        baseline_peak_optical,
        _,
    ) = make_analog_waveform(laser, baseline_target, channel, receiver, apd, tia)
    reference_time, reference_samples = sample_waveform(baseline_time, baseline_analog, adc, 0.0)
    _, reference_q, _ = quantize_voltage(reference_samples, adc)
    peak_index = int(np.argmax(reference_q))
    template_half_width = 10
    template_v = reference_q[peak_index - template_half_width : peak_index + template_half_width + 1]
    reference_scores = matched_filter_scores(reference_q, template_v)
    reference_estimate = sampled_peak_time_s(reference_time, reference_scores, interpolate=True)
    calibration_offset_s = reference_estimate - baseline_tof

    baseline_apd_density = shot_noise_rms_a(
        baseline_peak_optical,
        background_w,
        replace(apd, noise_bandwidth_hz=1.0),
    )
    baseline_noise_rms_v = total_output_noise_rms_v(baseline_apd_density, tia)
    detection_threshold_v = 5.0 * baseline_noise_rms_v
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []

    for reflectivity in REFLECTIVITIES:
        for distance_m in RANGES_M:
            target = replace(
                baseline_target,
                range_m=float(distance_m),
                reflectivity=reflectivity,
            )
            (
                dense_time,
                analog_v,
                true_tof_s,
                peak_optical_w,
                tia_clipped_sample_fraction,
            ) = make_analog_waveform(laser, target, channel, receiver, apd, tia)
            apd_density = shot_noise_rms_a(
                peak_optical_w,
                background_w,
                replace(apd, noise_bandwidth_hz=1.0),
            )
            noise_rms_v = total_output_noise_rms_v(apd_density, tia)
            peak_snr = float(np.max(analog_v) / noise_rms_v)
            conditional_error_values_m = []
            all_error_values_m = []
            clipped_trials = 0

            

            for _ in range(TRIALS_PER_CONDITION):
                phase = float(rng.random())

                sample_time, clean_samples = sample_waveform(
                    dense_time,
                    analog_v,
                    adc,
                    phase,
                )

                noisy_samples = clean_samples + rng.normal(
                    0.0,
                    noise_rms_v,
                    clean_samples.size,
                )

                _, quantized_v, clipped_fraction = quantize_voltage(
                    noisy_samples,
                    adc,
                )

                if clipped_fraction > 0.0:
                    clipped_trials += 1

                scores = matched_filter_scores(
                    quantized_v,
                    template_v,
                )

                estimate_s = sampled_peak_time_s(
                    sample_time,
                    scores,
                    interpolate=True,
                )

                corrected_estimate_s = (
                    estimate_s - calibration_offset_s
                )

                range_error_m = time_error_to_range_error_m(
                    corrected_estimate_s - true_tof_s
                )

                # 无条件误差：所有试验均计入
                all_error_values_m.append(range_error_m)

                valid_detection = (
                    float(np.max(scores)) >= detection_threshold_v
                    and abs(corrected_estimate_s - true_tof_s)
                    <= MAXIMUM_VALID_TIME_ERROR_S
                )

                # 条件误差：只统计满足检出条件的试验
                if valid_detection:
                    conditional_error_values_m.append(
                        range_error_m
                    )

            conditional_errors_m = np.asarray(
                conditional_error_values_m
            )
            all_errors_m = np.asarray(
                all_error_values_m
            )

            detection_rate = (
                conditional_errors_m.size
                / TRIALS_PER_CONDITION
            )

            conditional_bias_m = (
                float(np.mean(conditional_errors_m))
                if conditional_errors_m.size
                else np.nan
            )

            conditional_rmse_m = (
                float(
                    np.sqrt(
                        np.mean(conditional_errors_m**2)
                    )
                )
                if conditional_errors_m.size
                else np.nan
            )

            conditional_p95_m = (
                float(
                    np.percentile(
                        np.abs(conditional_errors_m),
                        95,
                    )
                )
                if conditional_errors_m.size
                else np.nan
            )

            unconditional_rmse_m = float(
                np.sqrt(
                    np.mean(all_errors_m**2)
                )
            )

            rows.append(
                (
                    float(distance_m),
                    reflectivity,
                    float(np.max(analog_v)),
                    peak_snr,
                    detection_rate,
                    conditional_bias_m,
                    conditional_rmse_m,
                    conditional_p95_m,
                    clipped_trials / TRIALS_PER_CONDITION,
                    tia_clipped_sample_fraction,
                    unconditional_rmse_m,
                )
            )

    output_path = DATA_DIR / "range_reflectivity_performance.csv"

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "range_m",
                "reflectivity",
                "analog_peak_v",
                "analog_peak_snr",
                "detection_rate",
                "conditional_bias_m",
                "conditional_rmse_m",
                "conditional_p95_absolute_error_m",
                "adc_clipped_trial_fraction",
                "tia_clipped_sample_fraction",
                "unconditional_rmse_m",
            ]
        )

        writer.writerows(rows)

    plt.rcParams.update(

        {
            "font.family": "sans-serif",
            "font.size": 9,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "legend.frameon": False,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(10.4, 7.5), sharex=True)
    colors = ("#1368aa", "#238636", "#f0a202", "#c73e1d")
    for reflectivity, color in zip(REFLECTIVITIES, colors):
        subset = [row for row in rows if row[1] == reflectivity]
        ranges = np.array([row[0] for row in subset])
        axes[0, 0].plot(ranges, [20.0 * np.log10(row[3]) for row in subset], "o-", color=color, label=f"{reflectivity:.0%}")
        axes[0, 1].plot(ranges, [100.0 * row[4] for row in subset], "o-", color=color)
        axes[1, 0].plot(ranges, [100.0 * row[6] for row in subset], "o-", color=color,label=f"{reflectivity:.0%} conditional")
        axes[1, 0].plot(
    ranges,[100.0 * row[10] for row in subset],"--",color=color,alpha=0.75,label=f"{reflectivity:.0%} all-trial",)

        axes[1, 1].plot(ranges, [100.0 * row[7] for row in subset], "o-", color=color)
    axes[0, 0].set(title="Analog peak SNR", ylabel="SNR (dB)")
    axes[0, 0].legend(title="Reflectivity")
    axes[0, 1].axhline(95.0, color="#333333", linestyle="--", linewidth=1, label="95% criterion")
    axes[0, 1].set(title="Valid detection rate", ylabel="Detection rate (%)", ylim=(0.0, 105.0))
    axes[0, 1].legend()
    axes[1, 0].axhline(5.0, color="#333333", linestyle=":", linewidth=1.2, label="5 cm criterion")
    axes[1, 0].set(title="Conditional and all-trial ranging RMSE", xlabel="Target range (m)", ylabel="RMSE (cm)")
    axes[1, 0].legend(fontsize=7, ncol=2,)
    axes[1, 1].set(title="Conditional P95 absolute error", xlabel="Target range (m)", ylabel="P95 absolute error (cm)")
    figure.suptitle("Matched-filter performance versus range and reflectivity", fontsize=13)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"07_range_reflectivity_performance.{suffix}")
    plt.close(figure)

    saturation_figure, saturation_axes = plt.subplots(1, 2, figsize=(10.4, 4.0), sharex=True)
    for reflectivity, color in zip(REFLECTIVITIES, colors):
        subset = [row for row in rows if row[1] == reflectivity]
        ranges = np.array([row[0] for row in subset])
        saturation_axes[0].plot(
            ranges,
            [100.0 * row[9] for row in subset],
            "o-",
            color=color,
            label=f"{reflectivity:.0%}",
        )
        saturation_axes[1].plot(
            ranges,
            [100.0 * row[8] for row in subset],
            "o-",
            color=color,
            label=f"{reflectivity:.0%}",
        )
    saturation_axes[0].set(
        title="TIA saturation within the analog pulse",
        xlabel="Target range (m)",
        ylabel="Clipped dense samples (%)",
    )
    saturation_axes[1].set(
        title="ADC saturation across Monte Carlo trials",
        xlabel="Target range (m)",
        ylabel="Trials with clipping (%)",
    )
    saturation_axes[0].legend(title="Reflectivity")
    saturation_axes[1].legend(title="Reflectivity")
    saturation_figure.suptitle("Receiver dynamic-range limits", fontsize=13)
    saturation_figure.tight_layout()
    for suffix in ("png", "pdf"):
        saturation_figure.savefig(FIGURE_DIR / f"08_receiver_dynamic_range.{suffix}")
    plt.close(saturation_figure)

    print(
        f"Range/reflectivity sweep: {TRIALS_PER_CONDITION} trials per condition, "
        f"seed {RANDOM_SEED}"
    )
    print(
        f"Criteria: detection >= {MINIMUM_DETECTION_RATE:.0%}, "
        f"conditional RMSE <= {MAXIMUM_RMSE_M * 100:.1f} cm"
    )
    print("Reflectivity  Reliable range (m)  Boundary observed")
    for reflectivity in REFLECTIVITIES:
        reliable = reliable_range_m(rows, reflectivity)
        if np.isnan(reliable):
            label = "none"
            boundary = "yes"
        elif reliable == float(RANGES_M[-1]):
            label = f">={reliable:.0f}"
            boundary = "no"
        else:
            label = f"{reliable:.0f}"
            boundary = "yes"
        print(f"{reflectivity:11.0%}  {label:>18s}  {boundary:>17s}")

    saturated_conditions = [
        row for row in rows if row[8] > 0.0 or row[9] > 0.0
    ]
    print(f"Conditions with TIA or ADC clipping: {len(saturated_conditions)}")
    for row in saturated_conditions:
        print(
            f"  range={row[0]:.0f} m, rho={row[1]:.0%}, "
            f"ADC clipped trials={row[8]:.1%}, "
            f"TIA clipped samples={row[9]:.3%}"
        )
    print(f"Baseline receiver noise RMS: {baseline_noise_rms_v * 1e3:.3f} mV")
    print(f"Matched-filter threshold:    {detection_threshold_v * 1e3:.3f} mV")
    print("Saved range_reflectivity_performance.csv and figure 07 (PNG/PDF).")


if __name__ == "__main__":
    main()
