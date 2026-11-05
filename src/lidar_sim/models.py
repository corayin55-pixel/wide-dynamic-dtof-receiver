"""Deterministic optical models for one pulsed dToF LiDAR channel."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

SPEED_OF_LIGHT_M_S = 299_792_458.0


@dataclass(frozen=True)
class Laser:
    wavelength_m: float
    pulse_energy_j: float
    pulse_fwhm_s: float
    initial_beam_diameter_m: float
    full_angle_divergence_rad: float
    transmitter_efficiency: float


@dataclass(frozen=True)
class Target:
    range_m: float
    reflectivity: float
    incidence_angle_rad: float
    projected_area_m2: float


@dataclass(frozen=True)
class Channel:
    atmospheric_extinction_per_m: float


@dataclass(frozen=True)
class Receiver:
    aperture_diameter_m: float
    optical_efficiency: float

    @property
    def aperture_area_m2(self) -> float:
        return np.pi * (self.aperture_diameter_m / 2.0) ** 2


def gaussian_sigma_from_fwhm(fwhm_s: float) -> float:
    if fwhm_s <= 0.0:
        raise ValueError("Pulse FWHM must be positive")
    return fwhm_s / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def gaussian_power(
    time_s: NDArray[np.float64],
    energy_j: float,
    fwhm_s: float,
    center_s: float = 0.0,
) -> NDArray[np.float64]:
    """Return a Gaussian pulse whose time integral equals energy_j."""
    if energy_j < 0.0:
        raise ValueError("Pulse energy cannot be negative")
    sigma_s = gaussian_sigma_from_fwhm(fwhm_s)
    peak_power_w = energy_j / (sigma_s * np.sqrt(2.0 * np.pi))
    return peak_power_w * np.exp(-0.5 * ((time_s - center_s) / sigma_s) ** 2)


def round_trip_delay_s(range_m: float) -> float:
    if range_m <= 0.0:
        raise ValueError("Range must be positive")
    return 2.0 * range_m / SPEED_OF_LIGHT_M_S


def spot_diameter_m(laser: Laser, range_m: float) -> float:
    """Approximate spot diameter using the full-angle beam divergence."""
    return laser.initial_beam_diameter_m + range_m * laser.full_angle_divergence_rad


def target_fill_factor(laser: Laser, target: Target) -> float:
    """Fraction of the beam footprint intercepted by the target."""
    spot_area_m2 = np.pi * (spot_diameter_m(laser, target.range_m) / 2.0) ** 2
    return min(1.0, target.projected_area_m2 / spot_area_m2)


def received_pulse_energy_j(
    laser: Laser,
    target: Target,
    channel: Channel,
    receiver: Receiver,
) -> float:
    """Received pulse energy for a coaxial far-field Lambertian model."""
    if not 0.0 <= target.reflectivity <= 1.0:
        raise ValueError("Target reflectivity must be within [0, 1]")
    if not 0.0 <= target.incidence_angle_rad < np.pi / 2.0:
        raise ValueError("Incidence angle must be within [0, pi/2)")

    atmospheric_two_way = np.exp(
        -2.0 * channel.atmospheric_extinction_per_m * target.range_m
    )
    geometric_coupling = receiver.aperture_area_m2 / (np.pi * target.range_m**2)

    return (
        laser.pulse_energy_j
        * laser.transmitter_efficiency
        * target_fill_factor(laser, target)
        * target.reflectivity
        * np.cos(target.incidence_angle_rad)
        * geometric_coupling
        * atmospheric_two_way
        * receiver.optical_efficiency
    )


def received_waveform(
    time_s: NDArray[np.float64],
    laser: Laser,
    target: Target,
    channel: Channel,
    receiver: Receiver,
) -> NDArray[np.float64]:
    energy_j = received_pulse_energy_j(laser, target, channel, receiver)
    delay_s = round_trip_delay_s(target.range_m)
    return gaussian_power(time_s, energy_j, laser.pulse_fwhm_s, delay_s)
