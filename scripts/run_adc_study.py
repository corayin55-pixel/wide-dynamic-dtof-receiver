"""Evaluate ADC sample-rate and bit-depth choices after the 250 MHz TIA."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lidar_sim.adc import ADC, quantization_error_rms_v, quantize_voltage, sample_waveform
from lidar_sim.config import build_models, load_config
from lidar_sim.detector import APD, multiplied_signal_current_a
from lidar_sim.models import received_waveform, round_trip_delay_s
from lidar_sim.tia import TIA, first_order_lowpass_transimpedance

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "results" / "data"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
SAMPLE_RATES_GHZ = (0.25, 0.50, 1.0, 2.0)
RESOLUTIONS_BITS = (8, 10, 12)
PHASES = np.linspace(0.0, 1.0, 101, endpoint=False)


def build_frontend_waveform(config):
    laser, target, channel, receiver = build_models(config)
    apd_values = config["apd"]
    apd = APD(
        apd_values["responsivity_a_w"],
        apd_values["gain"],
        apd_values["ionization_ratio"],
        apd_values["primary_dark_current_na"] * 1e-9,
        apd_values["noise_bandwidth_mhz"] * 1e6,
    )
    tia_values = config["tia"]
    tia = TIA(
        tia_values["transimpedance_kohm"] * 1e3,
        tia_values["bandwidth_mhz"] * 1e6,
        tia_values["input_current_noise_pa_sqrt_hz"] * 1e-12,
        tia_values["output_voltage_noise_nv_sqrt_hz"] * 1e-9,
        tia_values["output_swing_v"],
        tia_values["temperature_k"],
    )
    delay_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(delay_s - 20e-9, delay_s + 30e-9, 10e-12)
    optical_power_w = received_waveform(time_s, laser, target, channel, receiver)
    current_a = multiplied_signal_current_a(optical_power_w, apd)
    voltage_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], tia)
    return time_s, voltage_v, delay_s


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    time_s, analog_v, tof_s = build_frontend_waveform(config)
    adc_values = config["adc"]
    true_peak_v = float(np.max(analog_v))

    sample_rate_rows = []
    for sample_rate_ghz in SAMPLE_RATES_GHZ:
        adc = ADC(
            sample_rate_ghz * 1e9,
            adc_values["resolution_bits"],
            adc_values["input_min_v"],
            adc_values["input_max_v"],
        )
        peak_errors = []
        for phase in PHASES:
            _, sampled_v = sample_waveform(time_s, analog_v, adc, float(phase))
            _, quantized_v, clipped_fraction = quantize_voltage(sampled_v, adc)
            peak_errors.append((true_peak_v - float(np.max(quantized_v))) / true_peak_v)
            if clipped_fraction > 0.0:
                raise RuntimeError("Unexpected ADC clipping in sample-rate study")
        sample_rate_rows.append(
            (
                sample_rate_ghz,
                adc.sample_period_s * 1e9,
                adc.raw_range_bin_m,
                float(np.mean(peak_errors) * 100.0),
                float(np.max(peak_errors) * 100.0),
            )
        )

    bit_depth_rows = []
    for resolution_bits in RESOLUTIONS_BITS:
        adc = ADC(
            adc_values["sample_rate_ghz"] * 1e9,
            resolution_bits,
            adc_values["input_min_v"],
            adc_values["input_max_v"],
        )
        _, sampled_v = sample_waveform(time_s, analog_v, adc, 0.0)
        _, quantized_v, clipped_fraction = quantize_voltage(sampled_v, adc)
        bit_depth_rows.append(
            (
                resolution_bits,
                adc.lsb_v * 1e3,
                quantization_error_rms_v(sampled_v, quantized_v) * 1e3,
                clipped_fraction,
            )
        )

    with (DATA_DIR / "adc_sample_rate_summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["sample_rate_ghz", "sample_period_ns", "raw_range_bin_m", "mean_peak_error_percent", "worst_peak_error_percent"])
        writer.writerows(sample_rate_rows)
    with (DATA_DIR / "adc_bit_depth_summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["resolution_bits", "lsb_mv", "quantization_rms_mv", "clipped_fraction"])
        writer.writerows(bit_depth_rows)

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
    figure, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    axes[0].plot((time_s - tof_s) * 1e9, analog_v * 1e3, color="#333333", linewidth=2, label="Analog TIA output")
    sample_colors = ("#1368aa", "#238636", "#f0a202", "#c73e1d")
    for sample_rate_ghz, color in zip(SAMPLE_RATES_GHZ, sample_colors):
        adc = ADC(sample_rate_ghz * 1e9, adc_values["resolution_bits"], adc_values["input_min_v"], adc_values["input_max_v"])
        sampled_time_s, sampled_v = sample_waveform(time_s, analog_v, adc, 0.35)
        _, quantized_v, _ = quantize_voltage(sampled_v, adc)
        mask = ((sampled_time_s - tof_s) * 1e9 >= -7.0) & ((sampled_time_s - tof_s) * 1e9 <= 12.0)
        axes[0].plot((sampled_time_s[mask] - tof_s) * 1e9, quantized_v[mask] * 1e3, "o", markersize=3.5, color=color, label=f"{sample_rate_ghz:g} GHz")
    axes[0].set(title="ADC sampling of the 250 MHz TIA pulse", xlabel="Time relative to optical ToF (ns)", ylabel="Voltage (mV)", xlim=(-7.0, 12.0))
    axes[0].legend(fontsize=8)

    rates = np.array([row[0] for row in sample_rate_rows])
    worst_error = np.array([row[4] for row in sample_rate_rows])
    axes[1].plot(rates, worst_error, "o-", color="#1368aa", label="Worst phase peak error")
    axes[1].set(title="Sampling-rate sensitivity", xlabel="ADC sample rate (GHz)", ylabel="Worst peak error (%)", xticks=list(SAMPLE_RATES_GHZ))
    second_axis = axes[1].twinx()
    second_axis.plot(rates, [row[2] * 100 for row in sample_rate_rows], "s--", color="#c73e1d", label="Raw range interval")
    second_axis.set_ylabel("Raw range interval (cm)")
    lines = axes[1].get_lines() + second_axis.get_lines()
    axes[1].legend(lines, [line.get_label() for line in lines], loc="best", fontsize=8)
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"05_adc_sampling_tradeoff.{suffix}")
    plt.close(figure)

    print("ADC sample-rate study (10 bit, phase sweep)")
    print("Rate (GHz)  Period (ns)  Raw bin (cm)  Mean peak error (%)  Worst peak error (%)")
    for row in sample_rate_rows:
        print(f"{row[0]:10.2f}  {row[1]:11.3f}  {row[2] * 100:12.3f}  {row[3]:19.3f}  {row[4]:20.3f}")
    print("\nADC bit-depth study (1 GHz, phase = 0)")
    print("Bits  LSB (mV)  Quantization RMS (mV)  Clipped")
    for row in bit_depth_rows:
        print(f"{row[0]:4d}  {row[1]:8.4f}  {row[2]:21.4f}  {row[3]:7.3%}")
    print(f"\nAnalog TIA peak: {true_peak_v * 1e3:.3f} mV")
    print("Saved ADC CSV data and figure 05_adc_sampling_tradeoff (PNG/PDF).")


if __name__ == "__main__":
    main()
