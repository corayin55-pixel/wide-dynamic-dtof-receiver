"""Load engineering-unit YAML parameters and convert them to SI units."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .models import Channel, Laser, Receiver, Target


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as file:
        return yaml.safe_load(file)


def build_models(config: dict[str, Any]) -> tuple[Laser, Target, Channel, Receiver]:
    laser_cfg = config["laser"]
    target_cfg = config["target"]
    channel_cfg = config["channel"]
    receiver_cfg = config["receiver"]

    laser = Laser(
        wavelength_m=laser_cfg["wavelength_nm"] * 1e-9,
        pulse_energy_j=laser_cfg["pulse_energy_nj"] * 1e-9,
        pulse_fwhm_s=laser_cfg["pulse_fwhm_ns"] * 1e-9,
        initial_beam_diameter_m=laser_cfg["initial_beam_diameter_mm"] * 1e-3,
        full_angle_divergence_rad=laser_cfg["full_angle_divergence_mrad"] * 1e-3,
        transmitter_efficiency=laser_cfg["transmitter_efficiency"],
    )
    target = Target(
        range_m=target_cfg["range_m"],
        reflectivity=target_cfg["reflectivity"],
        incidence_angle_rad=np.deg2rad(target_cfg["incidence_angle_deg"]),
        projected_area_m2=target_cfg["projected_area_m2"],
    )
    channel = Channel(
        atmospheric_extinction_per_m=channel_cfg["atmospheric_extinction_per_m"]
    )
    receiver = Receiver(
        aperture_diameter_m=receiver_cfg["aperture_diameter_mm"] * 1e-3,
        optical_efficiency=receiver_cfg["optical_efficiency"],
    )
    return laser, target, channel, receiver
