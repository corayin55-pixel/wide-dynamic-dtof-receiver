"""Compare TIA bandwidth choices for the baseline APD return waveform."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lidar_sim.config import build_models, load_config
from lidar_sim.detector import (
    APD,
    Background,
    background_power_w,
    multiplied_signal_current_a,
    shot_noise_rms_a,
)
from lidar_sim.models import received_waveform, round_trip_delay_s
from lidar_sim.tia import (
    TIA,
    clip_output_voltage,
    first_order_lowpass_transimpedance,
    pulse_fwhm_s,
    total_output_noise_rms_v,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "results" / "data"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
BANDWIDTHS_MHZ = (100.0, 250.0, 500.0)


def make_apd(config) -> APD:
    values = config["apd"]
    return APD(
        values["responsivity_a_w"],
        values["gain"],
        values["ionization_ratio"],
        values["primary_dark_current_na"] * 1e-9,
        values["noise_bandwidth_mhz"] * 1e6,
    )


def make_tia(config, bandwidth_mhz: float) -> TIA:
    values = config["tia"]
    return TIA(
        values["transimpedance_kohm"] * 1e3,
        bandwidth_mhz * 1e6,
        values["input_current_noise_pa_sqrt_hz"] * 1e-12,
        values["output_voltage_noise_nv_sqrt_hz"] * 1e-9,
        values["output_swing_v"],
        values["temperature_k"],
    )


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, target, channel, receiver = build_models(config)
    apd = make_apd(config)

    background_cfg = config["background"]
    background = Background(
        background_cfg["spectral_radiance_w_m2_sr_nm"],
        background_cfg["receiver_full_fov_mrad"] * 1e-3,
        background_cfg["filter_bandwidth_nm"],
    )
    background_w = background_power_w(
        background, receiver.aperture_area_m2, receiver.optical_efficiency
    )

    delay_s = round_trip_delay_s(target.range_m)
    time_s = np.arange(delay_s - 20e-9, delay_s + 40e-9, 10e-12)
    optical_power_w = received_waveform(time_s, laser, target, channel, receiver)
    signal_current_a = multiplied_signal_current_a(optical_power_w, apd)
    peak_signal_power_w = float(np.max(optical_power_w))

    # One-hertz integration gives the white APD output-noise density in A/sqrt(Hz).
    apd_noise_density = shot_noise_rms_a(
        peak_signal_power_w, background_w, replace(apd, noise_bandwidth_hz=1.0)
    )

    colors = ("#1368aa", "#238636", "#c73e1d")
    summary_rows = []
    waveforms = {}
    for bandwidth_mhz in BANDWIDTHS_MHZ:
        tia = make_tia(config, bandwidth_mhz)
        raw_voltage_v = first_order_lowpass_transimpedance(
            signal_current_a, time_s[1] - time_s[0], tia
        )
        voltage_v, clipped_fraction = clip_output_voltage(
            raw_voltage_v, tia.output_swing_v
        )
        noise_rms_v = total_output_noise_rms_v(apd_noise_density, tia)
        peak_index = int(np.argmax(voltage_v))
        peak_voltage_v = float(voltage_v[peak_index])
        measured_fwhm_s = pulse_fwhm_s(time_s, voltage_v)
        peak_delay_s = time_s[peak_index] - delay_s
        peak_snr = peak_voltage_v / noise_rms_v
        waveforms[bandwidth_mhz] = voltage_v
        summary_rows.append(
            (
                bandwidth_mhz,
                peak_voltage_v,
                noise_rms_v,
                peak_snr,
                20.0 * np.log10(peak_snr),
                measured_fwhm_s * 1e9,
                peak_delay_s * 1e9,
                clipped_fraction,
            )
        )

    with (DATA_DIR / "tia_bandwidth_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "bandwidth_mhz",
                "peak_voltage_v",
                "noise_rms_v",
                "peak_snr_linear",
                "peak_snr_db",
                "output_fwhm_ns",
                "peak_delay_ns",
                "clipped_sample_fraction",
            ]
        )
        writer.writerows(summary_rows)

    with (DATA_DIR / "tia_waveforms.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            ["time_relative_to_tof_ns"]
            + [f"output_v_{int(value)}mhz" for value in BANDWIDTHS_MHZ]
        )
        for index, relative_time in enumerate((time_s - delay_s) * 1e9):
            writer.writerow(
                [relative_time]
                + [waveforms[value][index] for value in BANDWIDTHS_MHZ]
            )

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "legend.frameon": False,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    relative_time_ns = (time_s - delay_s) * 1e9
    for bandwidth_mhz, color in zip(BANDWIDTHS_MHZ, colors):
        axes[0].plot(
            relative_time_ns,
            waveforms[bandwidth_mhz] * 1e3,
            color=color,
            linewidth=1.8,
            label=f"{bandwidth_mhz:.0f} MHz",
        )
    axes[0].set(
        title="TIA pulse shaping",
        xlabel="Time relative to optical ToF (ns)",
        ylabel="Output voltage (mV)",
        xlim=(-10.0, 20.0),
    )
    axes[0].legend()

    bandwidth_array = np.array([row[0] for row in summary_rows])
    snr_db = np.array([row[4] for row in summary_rows])
    fwhm_ns = np.array([row[5] for row in summary_rows])
    axes[1].plot(
        bandwidth_array, snr_db, "o-", color="#1368aa", label="Peak SNR (dB)"
    )
    axes[1].set(
        title="Bandwidth trade-off",
        xlabel="TIA bandwidth (MHz)",
        ylabel="Peak SNR (dB)",
        xticks=list(BANDWIDTHS_MHZ),
    )
    second_axis = axes[1].twinx()
    second_axis.plot(
        bandwidth_array,
        fwhm_ns,
        "s--",
        color="#c73e1d",
        label="Output FWHM (ns)",
    )
    second_axis.set_ylabel("Output FWHM (ns)")
    lines = axes[1].get_lines() + second_axis.get_lines()
    axes[1].legend(lines, [line.get_label() for line in lines], loc="best")
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"04_tia_bandwidth_tradeoff.{suffix}")
    plt.close(figure)

    print("TIA bandwidth study")
    print("BW (MHz)  Peak (mV)  Noise RMS (mV)  SNR (dB)  FWHM (ns)  Peak delay (ns)  Clipped")
    for row in summary_rows:
        print(
            f"{row[0]:8.0f}  {row[1] * 1e3:9.3f}  {row[2] * 1e3:14.3f}"
            f"  {row[4]:8.3f}  {row[5]:9.3f}  {row[6]:15.3f}  {row[7]:7.3%}"
        )
    print(f"APD output-noise density: {apd_noise_density * 1e12:.3f} pA/sqrt(Hz)")
    print("Saved CSV data and figure 04_tia_bandwidth_tradeoff (PNG/PDF).")


if __name__ == "__main__":
    main()
