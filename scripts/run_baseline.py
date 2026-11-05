from pathlib import Path

import numpy as np

from lidar_sim.config import build_models, load_config
from lidar_sim.models import (
    received_pulse_energy_j,
    received_waveform,
    round_trip_delay_s,
    spot_diameter_m,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, target, channel, receiver = build_models(config)
    simulation = config["simulation"]

    delay_s = round_trip_delay_s(target.range_m)
    start_s = delay_s - simulation["window_before_ns"] * 1e-9
    stop_s = delay_s + simulation["window_after_ns"] * 1e-9
    step_s = simulation["time_step_ps"] * 1e-12
    time_s = np.arange(start_s, stop_s, step_s)

    waveform_w = received_waveform(time_s, laser, target, channel, receiver)
    received_energy_j = received_pulse_energy_j(laser, target, channel, receiver)

    print("Baseline 905 nm dToF optical link")
    print(f"Target range:          {target.range_m:10.3f} m")
    print(f"Round-trip delay:      {delay_s * 1e9:10.3f} ns")
    print(f"Spot diameter:         {spot_diameter_m(laser, target.range_m) * 1e3:10.3f} mm")
    print(f"Receiver aperture:     {receiver.aperture_diameter_m * 1e3:10.3f} mm")
    print(f"Received pulse energy: {received_energy_j * 1e15:10.3f} fJ")
    print(f"Peak received power:   {np.max(waveform_w) * 1e6:10.3f} uW")
    print(f"Numerical energy:      {np.trapezoid(waveform_w, time_s) * 1e15:10.3f} fJ")


if __name__ == "__main__":
    main()
