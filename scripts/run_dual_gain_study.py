"""Validate a 20 kOhm / 2 kOhm dual-gain LiDAR receiver."""

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
REFLECTIVITIES = (0.10, 0.80)
HIGH_GAIN_KOHM = 20.0
LOW_GAIN_KOHM = 2.0
TRIALS = 500
RANDOM_SEED = 20261108


def fixed_components(config):
    laser, target, channel, receiver = build_models(config)
    a = config["apd"]
    apd = APD(a["responsivity_a_w"], a["gain"], a["ionization_ratio"], a["primary_dark_current_na"] * 1e-9, a["noise_bandwidth_mhz"] * 1e6)
    b = config["background"]
    background = Background(b["spectral_radiance_w_m2_sr_nm"], b["receiver_full_fov_mrad"] * 1e-3, b["filter_bandwidth_nm"])
    background_w = background_power_w(background, receiver.aperture_area_m2, receiver.optical_efficiency)
    d = config["adc"]
    adc = ADC(d["sample_rate_ghz"] * 1e9, d["resolution_bits"], -0.2, 1.0)
    return laser, target, channel, receiver, apd, background_w, adc


def make_tia(config, gain_kohm):
    t = config["tia"]
    return TIA(gain_kohm * 1e3, t["bandwidth_mhz"] * 1e6, t["input_current_noise_pa_sqrt_hz"] * 1e-12, t["output_voltage_noise_nv_sqrt_hz"] * 1e-9, t["output_swing_v"], t["temperature_k"])


def analog_return(laser, target, channel, receiver, apd, tia):
    tof_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(tof_s - 20e-9, tof_s + 30e-9, 10e-12)
    optical_w = received_waveform(time_s, laser, target, channel, receiver)
    current_a = multiplied_signal_current_a(optical_w, apd)
    raw_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], tia)
    voltage_v, tia_clip_fraction = clip_output_voltage(raw_v, tia.output_swing_v)
    return time_s, voltage_v, tof_s, float(np.max(optical_w)), tia_clip_fraction


def calibrate_path(laser, target, channel, receiver, apd, tia, adc):
    time_s, voltage_v, tof_s, _, tia_clip = analog_return(laser, target, channel, receiver, apd, tia)
    sample_time, samples = sample_waveform(time_s, voltage_v, adc, 0.0)
    _, quantized, adc_clip = quantize_voltage(samples, adc)
    if tia_clip > 0.0 or adc_clip > 0.0:
        raise RuntimeError("Calibration target saturates the selected channel")
    peak_index = int(np.argmax(quantized))
    template = quantized[peak_index - 10 : peak_index + 11]
    scores = matched_filter_scores(quantized, template)
    offset_s = sampled_peak_time_s(sample_time, scores, interpolate=True) - tof_s
    return template, offset_s


def process_channel(
    sample_time,
    clean_v,
    shared_apd_noise_v,
    independent_tia_noise_v,
    noise_rms_v,
    tia_clipped,
    adc,
    template,
    offset_s,
):
    noisy_v = clean_v + shared_apd_noise_v + independent_tia_noise_v
    _, quantized_v, adc_clip_fraction = quantize_voltage(noisy_v, adc)
    clipped = tia_clipped or adc_clip_fraction > 0.0
    scores = matched_filter_scores(quantized_v, template)
    estimate_s = sampled_peak_time_s(sample_time, scores, interpolate=True) - offset_s
    detected = float(np.max(scores)) >= 5.0 * noise_rms_v
    return detected, clipped, estimate_s


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, base_target, channel, receiver, apd, background_w, adc = fixed_components(config)
    high_tia = make_tia(config, HIGH_GAIN_KOHM)
    low_tia = make_tia(config, LOW_GAIN_KOHM)

    high_template, high_offset = calibrate_path(laser, replace(base_target, range_m=50.0, reflectivity=0.10), channel, receiver, apd, high_tia, adc)
    low_template, low_offset = calibrate_path(laser, replace(base_target, range_m=20.0, reflectivity=0.80), channel, receiver, apd, low_tia, adc)
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []

    for reflectivity in REFLECTIVITIES:
        for distance_m in RANGES_M:
            target = replace(base_target, range_m=float(distance_m), reflectivity=reflectivity)
            high_time, high_v, tof_s, peak_optical, high_tia_clip = analog_return(laser, target, channel, receiver, apd, high_tia)
            low_time, low_v, _, _, low_tia_clip = analog_return(laser, target, channel, receiver, apd, low_tia)
            apd_density = shot_noise_rms_a(peak_optical, background_w, replace(apd, noise_bandwidth_hz=1.0))
            high_noise = total_output_noise_rms_v(apd_density, high_tia)
            low_noise = total_output_noise_rms_v(apd_density, low_tia)
            high_apd_noise = (
                apd_density
                * np.sqrt(high_tia.equivalent_noise_bandwidth_hz)
                * high_tia.transimpedance_ohm
            )
            low_apd_noise = (
                apd_density
                * np.sqrt(low_tia.equivalent_noise_bandwidth_hz)
                * low_tia.transimpedance_ohm
            )
            high_tia_noise = np.sqrt(max(0.0, high_noise**2 - high_apd_noise**2))
            low_tia_noise = np.sqrt(max(0.0, low_noise**2 - low_apd_noise**2))
            errors_m = []
            high_used = 0
            low_used = 0
            both_clipped = 0
            false_time_selections = 0

            for _ in range(TRIALS):
                phase = float(rng.random())
                high_sample_time, high_clean_v = sample_waveform(
                    high_time, high_v, adc, phase
                )
                low_sample_time, low_clean_v = sample_waveform(
                    low_time, low_v, adc, phase
                )
                if high_clean_v.size != low_clean_v.size:
                    raise RuntimeError("Dual-gain channels produced unequal sample counts")
                shared_standard = rng.normal(0.0, 1.0, high_clean_v.size)
                high_shared_noise = shared_standard * high_apd_noise
                low_shared_noise = shared_standard * low_apd_noise
                high_independent_noise = rng.normal(
                    0.0, high_tia_noise, high_clean_v.size
                )
                low_independent_noise = rng.normal(
                    0.0, low_tia_noise, low_clean_v.size
                )
                high_valid, high_clipped, high_estimate = process_channel(
                    high_sample_time,
                    high_clean_v,
                    high_shared_noise,
                    high_independent_noise,
                    high_noise,
                    high_tia_clip > 0.0,
                    adc,
                    high_template,
                    high_offset,
                )
                low_valid, low_clipped, low_estimate = process_channel(
                    low_sample_time,
                    low_clean_v,
                    low_shared_noise,
                    low_independent_noise,
                    low_noise,
                    low_tia_clip > 0.0,
                    adc,
                    low_template,
                    low_offset,
                )
                if high_valid and not high_clipped:
                    estimate_s = high_estimate
                    high_used += 1
                elif low_valid and not low_clipped:
                    estimate_s = low_estimate
                    low_used += 1
                else:
                    if high_clipped and low_clipped:
                        both_clipped += 1
                    continue
                time_error_s = estimate_s - tof_s
                if abs(time_error_s) <= 5e-9:
                    errors_m.append(time_error_to_range_error_m(time_error_s))
                else:
                    false_time_selections += 1

            errors = np.asarray(errors_m)
            selected_trials = high_used + low_used
            rows.append(
                (
                    float(distance_m),
                    reflectivity,
                    errors.size / TRIALS,
                    float(np.sqrt(np.mean(errors**2))) if errors.size else np.nan,
                    float(np.percentile(np.abs(errors), 95)) if errors.size else np.nan,
                    high_used / TRIALS,
                    low_used / TRIALS,
                    both_clipped / TRIALS,
                    false_time_selections / TRIALS,
                    selected_trials / TRIALS,
                )
            )

    with (DATA_DIR / "dual_gain_performance.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["range_m", "reflectivity", "correct_detection_rate", "conditional_rmse_m", "conditional_p95_absolute_error_m", "high_gain_usage_rate", "low_gain_usage_rate", "both_channels_clipped_rate", "false_time_selection_rate", "any_channel_selected_rate"])
        writer.writerows(rows)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), sharex=True)
    colors = {0.10: "#1368aa", 0.80: "#c73e1d"}
    for reflectivity in REFLECTIVITIES:
        subset = [row for row in rows if row[1] == reflectivity]
        ranges = [row[0] for row in subset]
        color = colors[reflectivity]
        axes[0, 0].plot(ranges, [100.0 * row[2] for row in subset], "o-", color=color, label=f"{reflectivity:.0%}")
        axes[0, 1].plot(ranges, [100.0 * row[3] for row in subset], "o-", color=color, label=f"{reflectivity:.0%}")
        axes[1, 0].plot(ranges, [100.0 * row[5] for row in subset], "o-", color=color, label=f"High, {reflectivity:.0%}")
        axes[1, 0].plot(ranges, [100.0 * row[6] for row in subset], "s--", color=color, label=f"Low, {reflectivity:.0%}")
        axes[1, 1].plot(
            ranges,
            [100.0 * row[8] for row in subset],
            "o-",
            color=color,
            label=f"False time, {reflectivity:.0%}",
        )
        axes[1, 1].plot(
            ranges,
            [100.0 * row[7] for row in subset],
            "s--",
            color=color,
            label=f"Both clipped, {reflectivity:.0%}",
        )
    axes[0, 0].axhline(95.0, color="#333333", linestyle="--", linewidth=1)
    axes[0, 0].set(title="Dual-gain detection rate", ylabel="Detection rate (%)", ylim=(0.0, 105.0))
    axes[0, 1].axhline(5.0, color="#333333", linestyle="--", linewidth=1)
    axes[0, 1].set(title="Conditional ranging RMSE", ylabel="RMSE (cm)")
    axes[1, 0].set(title="Selected channel usage", xlabel="Target range (m)", ylabel="Usage rate (%)")
    axes[1, 1].set(title="Fusion failure modes", xlabel="Target range (m)", ylabel="Trial rate (%)")
    for axis in axes.flat:
        axis.legend(fontsize=8)
    figure.suptitle("20 kOhm / 2 kOhm dual-gain receiver", fontsize=13)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"10_dual_gain_receiver.{suffix}")
    plt.close(figure)

    print(f"Dual-gain study: high={HIGH_GAIN_KOHM:.0f} kOhm, low={LOW_GAIN_KOHM:.0f} kOhm, trials={TRIALS}, seed={RANDOM_SEED}")
    print("Reflectivity  Reliable interval (m)  High-to-low transition (m)  Both-clipped near limit")
    for reflectivity in REFLECTIVITIES:
        subset = [row for row in rows if row[1] == reflectivity]
        reliable = [row[0] for row in subset if row[2] >= 0.95 and row[3] <= 0.05]
        interval = "none" if not reliable else f"{min(reliable):.0f}-{max(reliable):.0f}"
        low_used = [row[0] for row in subset if row[6] > row[5]]
        transition = "none" if not low_used else f"<= {max(low_used):.0f}"
        clipped = [row[0] for row in subset if row[7] > 0.01]
        clipped_limit = "none" if not clipped else f"<= {max(clipped):.0f} m"
        print(f"{reflectivity:11.0%}  {interval:>21s}  {transition:>26s}  {clipped_limit:>23s}")
    maximum_false_selection = max(row[8] for row in rows)
    print(f"Maximum false-time selection rate: {maximum_false_selection:.3%}")
    print("Saved dual_gain_performance.csv and figure 10_dual_gain_receiver (PNG/PDF).")


if __name__ == "__main__":
    main()
