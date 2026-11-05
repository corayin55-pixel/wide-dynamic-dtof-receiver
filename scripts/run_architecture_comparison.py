"""Compare single-gain and dual-gain LiDAR receiver architectures."""

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
NEAR_RANGES_M = np.array(
    [0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
)

FAR_RANGES_M = np.arange(
    20.0,
    251.0,
    10.0,
)

RANGES_M = np.concatenate(
    [NEAR_RANGES_M, FAR_RANGES_M]
)

REFLECTIVITIES = (
    0.10,
    0.20,
    0.50,
    0.80,
)

HIGH_GAIN_KOHM = 20.0
LOW_GAIN_KOHM = 2.0
TRIALS = 500
RANDOM_SEED = 20261109


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

def calculate_metrics(
    errors_m,
    false_time_count,
    trials,
):
    """Calculate detection and ranging metrics for one architecture."""
    errors = np.asarray(errors_m)

    correct_detection_rate = (
        errors.size / trials
    )

    false_time_rate = (
        false_time_count / trials
    )

    conditional_rmse_m = (
        float(
            np.sqrt(
                np.mean(errors**2)
            )
        )
        if errors.size
        else np.nan
    )

    conditional_p95_m = (
        float(
            np.percentile(
                np.abs(errors),
                95,
            )
        )
        if errors.size
        else np.nan
    )

    reliable = (
        correct_detection_rate >= 0.95
        and np.isfinite(conditional_rmse_m)
        and conditional_rmse_m <= 0.05
    )

    return (
        correct_detection_rate,
        false_time_rate,
        conditional_rmse_m,
        conditional_p95_m,
        reliable,
    )

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
    comparison_rows = []

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
            # 三种接收架构的正确测距误差
            high_only_errors_m = []
            low_only_errors_m = []
            dual_errors_m = []

            # 各方案的错误时间选择次数
            high_only_false_time = 0
            low_only_false_time = 0
            dual_false_time = 0

            # 双增益通道使用及饱和统计
            high_used = 0
            low_used = 0
            both_clipped = 0


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
                # 方案1：单20 kΩ高增益通道
                if high_valid and not high_clipped:
                    high_time_error_s = (
                        high_estimate - tof_s
                    )

                    if abs(high_time_error_s) <= 5e-9:
                        high_only_errors_m.append(
                            time_error_to_range_error_m(
                                high_time_error_s
                            )
                        )
                    else:
                        high_only_false_time += 1

                # 方案2：单2 kΩ低增益通道
                if low_valid and not low_clipped:
                    low_time_error_s = (
                        low_estimate - tof_s
                    )

                    if abs(low_time_error_s) <= 5e-9:
                        low_only_errors_m.append(
                            time_error_to_range_error_m(
                                low_time_error_s
                            )
                        )
                    else:
                        low_only_false_time += 1

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
                dual_time_error_s = estimate_s - tof_s

                if abs(dual_time_error_s) <= 5e-9:
                    dual_errors_m.append(
                        time_error_to_range_error_m(
                            dual_time_error_s
                        )
                    )
                else:
                    dual_false_time += 1
            (
                high_detection_rate,
                high_false_time_rate,
                high_rmse_m,
                high_p95_m,
                high_reliable,
            ) = calculate_metrics(
                high_only_errors_m,
                high_only_false_time,
                TRIALS,
            )
            (
                low_detection_rate,
                low_false_time_rate,
                low_rmse_m,
                low_p95_m,
                low_reliable,
            ) = calculate_metrics(
                low_only_errors_m,
                low_only_false_time,
                TRIALS,
            )
            (
                dual_detection_rate,
                dual_false_time_rate,
                dual_rmse_m,
                dual_p95_m,
                dual_reliable,
            ) = calculate_metrics(
                dual_errors_m,
                dual_false_time,
                TRIALS,
            )
            comparison_rows.append(
                (
                    float(distance_m),
                    reflectivity,
                    "single_20kohm",
                    high_detection_rate,
                    high_false_time_rate,
                    high_rmse_m,
                    high_p95_m,
                    high_reliable,
                )
            )

            comparison_rows.append(
                (
                    float(distance_m),
                    reflectivity,
                    "single_2kohm",
                    low_detection_rate,
                    low_false_time_rate,
                    low_rmse_m,
                    low_p95_m,
                    low_reliable,
                )
            )

            comparison_rows.append(
                (
                    float(distance_m),
                    reflectivity,
                    "dual_20kohm_2kohm",
                    dual_detection_rate,
                    dual_false_time_rate,
                    dual_rmse_m,
                    dual_p95_m,
                    dual_reliable,
                )
            )


    # 保存三种接收架构的逐条件结果
    output_path = (
        DATA_DIR / "architecture_comparison.csv"
    )

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
                "architecture",
                "correct_detection_rate",
                "false_time_rate",
                "conditional_rmse_m",
                "conditional_p95_absolute_error_m",
                "reliable",
            ]
        )

        writer.writerows(comparison_rows)

    architecture_names = (
        "single_20kohm",
        "single_2kohm",
        "dual_20kohm_2kohm",
    )

    architecture_labels = {
        "single_20kohm": "Single 20 kOhm",
        "single_2kohm": "Single 2 kOhm",
        "dual_20kohm_2kohm": "20/2 kOhm dual gain",
    }

    architecture_colors = {
        "single_20kohm": "#1368aa",
        "single_2kohm": "#f0a202",
        "dual_20kohm_2kohm": "#238636",
    }

    architecture_markers = {
        "single_20kohm": "o",
        "single_2kohm": "s",
        "dual_20kohm_2kohm": "^",
    }

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

    # 图12：10%与80%目标下的详细性能对照
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(10.5, 7.4),
        sharex=True,
    )

    comparison_reflectivities = (
        0.10,
        0.80,
    )

    for column, reflectivity in enumerate(
        comparison_reflectivities
    ):
        for architecture in architecture_names:
            subset = [
                row
                for row in comparison_rows
                if row[1] == reflectivity
                and row[2] == architecture
            ]

            ranges = np.array(
                [row[0] for row in subset]
            )

            detection_rates = np.array(
                [100.0 * row[3] for row in subset]
            )

            rmse_cm = np.array(
                [100.0 * row[5] for row in subset]
            )

            label = architecture_labels[
                architecture
            ]

            color = architecture_colors[
                architecture
            ]

            marker = architecture_markers[
                architecture
            ]

            axes[0, column].plot(
                ranges,
                detection_rates,
                marker=marker,
                color=color,
                linewidth=1.6,
                markersize=4,
                label=label,
            )

            axes[1, column].plot(
                ranges,
                rmse_cm,
                marker=marker,
                color=color,
                linewidth=1.6,
                markersize=4,
                label=label,
            )

        axes[0, column].axhline(
            95.0,
            color="#333333",
            linestyle="--",
            linewidth=1.0,
            label="95% criterion",
        )

        axes[1, column].axhline(
            5.0,
            color="#333333",
            linestyle="--",
            linewidth=1.0,
            label="5 cm criterion",
        )

        axes[0, column].set(
            title=(
                f"{reflectivity:.0%} target: "
                "correct detection rate"
            ),
            ylabel="Detection rate (%)",
            ylim=(-3.0, 103.0),
        )

        axes[1, column].set(
            title=(
                f"{reflectivity:.0%} target: "
                "conditional RMSE"
            ),
            xlabel="Target range (m)",
            ylabel="RMSE (cm)",
        )

        for row_index in range(2):
            axes[row_index, column].set_xscale(
                "log"
            )
            axes[row_index, column].set_xlim(
                0.18,
                300.0,
            )
            axes[row_index, column].set_xticks(
                [
                    0.2,
                    0.5,
                    1,
                    2,
                    5,
                    10,
                    20,
                    50,
                    100,
                    250,
                ]
            )
            axes[row_index, column].set_xticklabels(
                [
                    "0.2",
                    "0.5",
                    "1",
                    "2",
                    "5",
                    "10",
                    "20",
                    "50",
                    "100",
                    "250",
                ]
            )

    axes[0, 0].legend(fontsize=8)
    axes[1, 0].legend(fontsize=8)

    figure.suptitle(
        "Single-gain and dual-gain receiver comparison",
        fontsize=13,
    )

    figure.tight_layout()

    for suffix in ("png", "pdf"):
        figure.savefig(
            FIGURE_DIR
            / f"12_architecture_performance_comparison.{suffix}"
        )

    plt.close(figure)

    # 将离散可靠点分割为连续区间
    def get_reliable_intervals(
        reflectivity,
        architecture,
    ):
        reliable_indices = [
            index
            for index, distance in enumerate(RANGES_M)
            if any(
                row[0] == float(distance)
                and row[1] == reflectivity
                and row[2] == architecture
                and row[7]
                for row in comparison_rows
            )
        ]

        if not reliable_indices:
            return []

        index_groups = [
            [reliable_indices[0]]
        ]

        for index in reliable_indices[1:]:
            if index == index_groups[-1][-1] + 1:
                index_groups[-1].append(index)
            else:
                index_groups.append([index])

        return [
            (
                float(RANGES_M[group[0]]),
                float(RANGES_M[group[-1]]),
            )
            for group in index_groups
        ]

    # 图13：三种方案的可靠距离区间 vs 反射率
    interval_figure, interval_axis = plt.subplots(
        figsize=(8.2, 5.2)
    )

    reflectivity_percent = np.array(
        [
            100.0 * value
            for value in REFLECTIVITIES
        ]
    )

    architecture_offsets = {
        "single_20kohm": -1.5,
        "single_2kohm": 0.0,
        "dual_20kohm_2kohm": 1.5,
    }

    for architecture in architecture_names:
        first_segment = True

        for reflectivity, reflectivity_x in zip(
            REFLECTIVITIES,
            reflectivity_percent,
        ):
            intervals = get_reliable_intervals(
                reflectivity,
                architecture,
            )

            x_position = (
                reflectivity_x
                + architecture_offsets[architecture]
            )

            for lower_range, upper_range in intervals:
                interval_axis.vlines(
                    x_position,
                    lower_range,
                    upper_range,
                    color=architecture_colors[
                        architecture
                    ],
                    linewidth=6,
                    alpha=0.85,
                    label=(
                        architecture_labels[
                            architecture
                        ]
                        if first_segment
                        else None
                    ),
                )

                interval_axis.scatter(
                    [x_position, x_position],
                    [lower_range, upper_range],
                    color=architecture_colors[
                        architecture
                    ],
                    marker=architecture_markers[
                        architecture
                    ],
                    s=32,
                    zorder=3,
                )

                first_segment = False

    interval_axis.set_yscale("log")
    interval_axis.set_ylim(0.18, 300.0)

    interval_axis.set_xticks(
        reflectivity_percent
    )

    interval_axis.set_xticklabels(
        [
            f"{value:.0%}"
            for value in REFLECTIVITIES
        ]
    )

    interval_axis.set_yticks(
        [
            0.2,
            0.5,
            1,
            2,
            5,
            10,
            20,
            50,
            100,
            250,
        ]
    )

    interval_axis.set_yticklabels(
        [
            "0.2",
            "0.5",
            "1",
            "2",
            "5",
            "10",
            "20",
            "50",
            "100",
            "250",
        ]
    )

    interval_axis.set(
        title=(
        "Reliable range comparison: "
        "single-gain versus dual-gain"
        ),
        xlabel="Target reflectivity",
        ylabel="Reliable target range (m)",
    )

    interval_axis.legend()
    interval_figure.tight_layout()

    for suffix in ("png", "pdf"):
        interval_figure.savefig(
            FIGURE_DIR
            / f"13_architecture_reliable_intervals.{suffix}"
        )

    plt.close(interval_figure)

    # 终端输出可靠区间摘要
    print(
        "Architecture comparison: "
        f"{TRIALS} trials per condition, "
        f"seed={RANDOM_SEED}"
    )

    print(
        "Architecture              "
        "Reflectivity  Reliable intervals (m)"
    )

    for architecture in architecture_names:
        for reflectivity in REFLECTIVITIES:
            intervals = get_reliable_intervals(
                reflectivity,
                architecture,
            )

            if intervals:
                interval_text = ", ".join(
                    f"{lower:g}-{upper:g}"
                    for lower, upper in intervals
                )
            else:
                interval_text = "none"

            print(
                f"{architecture_labels[architecture]:26s}"
                f"{reflectivity:12.0%}  "
                f"{interval_text}"
            )

    print(
        "Saved architecture_comparison.csv, "
        "figure 12, and figure 13."
    )


if __name__ == "__main__":
    main()
