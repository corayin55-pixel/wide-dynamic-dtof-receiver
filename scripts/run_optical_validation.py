"""Generate validation data and report-ready figures for the optical link model."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from lidar_sim.config import build_models, load_config
from lidar_sim.models import (
    gaussian_power,
    received_pulse_energy_j,
    received_waveform,
    round_trip_delay_s,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "results" / "data"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures"
REFLECTIVITIES = (0.10, 0.20, 0.50, 0.80)


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.7,
            "legend.frameon": False,
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def save_figure(figure: plt.Figure, filename: str) -> None:
    for suffix in ("png", "pdf"):
        figure.savefig(FIGURE_DIR / f"{filename}.{suffix}")
    plt.close(figure)


def write_csv(path: Path, header: list[str], rows: list[tuple[float, ...]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(rows)


def plot_pulse_waveforms(laser, target, channel, receiver) -> None:
    relative_time_s = np.linspace(-15e-9, 15e-9, 2401)
    transmit_w = gaussian_power(
        relative_time_s, laser.pulse_energy_j, laser.pulse_fwhm_s
    )
    delay_s = round_trip_delay_s(target.range_m)
    receive_time_s = delay_s + relative_time_s
    receive_w = received_waveform(receive_time_s, laser, target, channel, receiver)

    write_csv(
        DATA_DIR / "pulse_waveforms.csv",
        [
            "transmit_relative_time_ns",
            "transmit_power_w",
            "receive_absolute_time_ns",
            "receive_power_w",
        ],
        list(
            zip(
                relative_time_s * 1e9,
                transmit_w,
                receive_time_s * 1e9,
                receive_w,
            )
        ),
    )

    figure, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))
    axes[0].plot(relative_time_s * 1e9, transmit_w, color="#1368aa", linewidth=2)
    axes[0].axvspan(
        -laser.pulse_fwhm_s * 0.5e9,
        laser.pulse_fwhm_s * 0.5e9,
        color="#f0a202",
        alpha=0.18,
        label=f"FWHM = {laser.pulse_fwhm_s * 1e9:.1f} ns",
    )
    axes[0].set(
        title="Transmitted Gaussian pulse",
        xlabel="Time relative to pulse center (ns)",
        ylabel="Optical power (W)",
    )
    axes[0].legend()

    axes[1].plot(
        receive_time_s * 1e9, receive_w * 1e6, color="#238636", linewidth=2
    )
    axes[1].axvline(
        delay_s * 1e9,
        color="#c73e1d",
        linestyle="--",
        linewidth=1.3,
        label=f"ToF = {delay_s * 1e9:.2f} ns",
    )
    axes[1].set(
        title=f"Received pulse at {target.range_m:.0f} m",
        xlabel="Absolute time (ns)",
        ylabel="Received power (uW)",
    )
    axes[1].legend()
    figure.tight_layout()
    save_figure(figure, "01_pulse_waveforms")


def plot_energy_sweep(laser, target, channel, receiver) -> float:
    ranges_m = np.geomspace(1.0, 150.0, 220)
    energy_by_reflectivity = {}
    for reflectivity in REFLECTIVITIES:
        energy_by_reflectivity[reflectivity] = np.array(
            [
                received_pulse_energy_j(
                    laser,
                    replace(target, range_m=float(distance), reflectivity=reflectivity),
                    channel,
                    receiver,
                )
                for distance in ranges_m
            ]
        )

    zero_atmosphere = replace(channel, atmospheric_extinction_per_m=0.0)
    ideal_energy_j = np.array(
        [
            received_pulse_energy_j(
                laser,
                replace(target, range_m=float(distance), reflectivity=0.10),
                zero_atmosphere,
                receiver,
            )
            for distance in ranges_m
        ]
    )
    slope, intercept = np.polyfit(np.log10(ranges_m), np.log10(ideal_energy_j), 1)
    fitted_energy_j = 10 ** (intercept + slope * np.log10(ranges_m))

    header = ["range_m"] + [
        f"received_energy_j_rho_{int(rho * 100):02d}" for rho in REFLECTIVITIES
    ]
    rows = [
        tuple(
            [float(distance)]
            + [float(energy_by_reflectivity[rho][index]) for rho in REFLECTIVITIES]
        )
        for index, distance in enumerate(ranges_m)
    ]
    write_csv(DATA_DIR / "received_energy_vs_range.csv", header, rows)

    colors = ("#1368aa", "#238636", "#f0a202", "#c73e1d")
    figure, axis = plt.subplots(figsize=(6.6, 4.6))
    for reflectivity, color in zip(REFLECTIVITIES, colors):
        axis.loglog(
            ranges_m,
            energy_by_reflectivity[reflectivity] * 1e15,
            color=color,
            linewidth=2,
            label=f"Reflectivity = {reflectivity:.0%}",
        )
    axis.loglog(
        ranges_m,
        fitted_energy_j * 1e15,
        color="#333333",
        linestyle="--",
        linewidth=1.3,
        label=f"No-atmosphere fit: slope = {slope:.3f}",
    )
    axis.set(
        title="Received pulse energy versus target range",
        xlabel="Target range (m)",
        ylabel="Received pulse energy (fJ)",
    )
    axis.legend()
    figure.tight_layout()
    save_figure(figure, "02_received_energy_vs_range")
    return float(slope)


def plot_time_of_flight() -> None:
    ranges_m = np.linspace(1.0, 150.0, 300)
    delays_ns = np.array(
        [round_trip_delay_s(float(distance)) for distance in ranges_m]
    ) * 1e9
    write_csv(
        DATA_DIR / "time_of_flight_vs_range.csv",
        ["range_m", "round_trip_delay_ns"],
        list(zip(ranges_m, delays_ns)),
    )

    marked_ranges_m = np.array([50.0, 100.0, 150.0])
    marked_delays_ns = np.array(
        [round_trip_delay_s(float(distance)) for distance in marked_ranges_m]
    ) * 1e9

    figure, axis = plt.subplots(figsize=(6.6, 4.4))
    axis.plot(ranges_m, delays_ns, color="#1368aa", linewidth=2)
    axis.scatter(marked_ranges_m, marked_delays_ns, color="#c73e1d", zorder=3)
    for distance, delay in zip(marked_ranges_m, marked_delays_ns):
        axis.annotate(
            f"{distance:.0f} m: {delay:.2f} ns",
            xy=(distance, delay),
            xytext=(-8, 10),
            textcoords="offset points",
            ha="right",
            fontsize=9,
        )
    axis.set(
        title="Round-trip time of flight versus target range",
        xlabel="Target range (m)",
        ylabel="Round-trip delay (ns)",
        xlim=(0.0, 157.0),
        ylim=(0.0, 1060.0),
    )
    figure.tight_layout()
    save_figure(figure, "03_time_of_flight_vs_range")


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_plot_style()

    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, target, channel, receiver = build_models(config)
    plot_pulse_waveforms(laser, target, channel, receiver)
    fitted_slope = plot_energy_sweep(laser, target, channel, receiver)
    plot_time_of_flight()

    print("Optical validation completed")
    print(f"No-atmosphere log-log slope: {fitted_slope:.6f} (expected: -2.0)")
    print(f"CSV data: {DATA_DIR}")
    print(f"Figures:  {FIGURE_DIR}")
    print("Generated 3 PNG figures, 3 PDF figures, and 3 CSV files")


if __name__ == "__main__":
    main()
