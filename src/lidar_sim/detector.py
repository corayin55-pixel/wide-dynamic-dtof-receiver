"""Photon conversion and Si APD receiver models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PLANCK_CONSTANT_J_S = 6.626_070_15e-34
SPEED_OF_LIGHT_M_S = 299_792_458.0
ELEMENTARY_CHARGE_C = 1.602_176_634e-19


@dataclass(frozen=True)
class APD:
    responsivity_a_w: float
    gain: float
    ionization_ratio: float
    primary_dark_current_a: float
    noise_bandwidth_hz: float

    @property
    def excess_noise_factor(self) -> float:
        """McIntyre model for electron-initiated multiplication."""
        k = self.ionization_ratio
        multiplication = self.gain
        return k * multiplication + (1.0 - k) * (2.0 - 1.0 / multiplication)


@dataclass(frozen=True)
class Background:
    spectral_radiance_w_m2_sr_nm: float
    receiver_full_fov_rad: float
    filter_bandwidth_nm: float

    @property
    def receiver_solid_angle_sr(self) -> float:
        half_angle = self.receiver_full_fov_rad / 2.0
        return 2.0 * np.pi * (1.0 - np.cos(half_angle))


def photon_energy_j(wavelength_m: float) -> float:
    if wavelength_m <= 0.0:
        raise ValueError("Wavelength must be positive")
    return PLANCK_CONSTANT_J_S * SPEED_OF_LIGHT_M_S / wavelength_m


def photons_per_pulse(pulse_energy_j: float, wavelength_m: float) -> float:
    if pulse_energy_j < 0.0:
        raise ValueError("Pulse energy cannot be negative")
    return pulse_energy_j / photon_energy_j(wavelength_m)


def quantum_efficiency(responsivity_a_w: float, wavelength_m: float) -> float:
    """Convert responsivity to external quantum efficiency."""
    efficiency = responsivity_a_w * photon_energy_j(wavelength_m) / ELEMENTARY_CHARGE_C
    if not 0.0 <= efficiency <= 1.0:
        raise ValueError("Responsivity implies a quantum efficiency outside [0, 1]")
    return efficiency


def mean_primary_electrons(
    pulse_energy_j: float, wavelength_m: float, responsivity_a_w: float
) -> float:
    return photons_per_pulse(pulse_energy_j, wavelength_m) * quantum_efficiency(
        responsivity_a_w, wavelength_m
    )


def background_power_w(
    background: Background,
    receiver_aperture_area_m2: float,
    receiver_optical_efficiency: float,
) -> float:
    """Background power from spectral radiance and receiver etendue."""
    return (
        background.spectral_radiance_w_m2_sr_nm
        * receiver_aperture_area_m2
        * background.receiver_solid_angle_sr
        * background.filter_bandwidth_nm
        * receiver_optical_efficiency
    )


def primary_photocurrent_a(optical_power_w: float, responsivity_a_w: float) -> float:
    if optical_power_w < 0.0:
        raise ValueError("Optical power cannot be negative")
    return optical_power_w * responsivity_a_w


def multiplied_signal_current_a(optical_power_w: float, apd: APD) -> float:
    return apd.gain * primary_photocurrent_a(optical_power_w, apd.responsivity_a_w)


def shot_noise_rms_a(
    signal_power_w: float,
    background_power_w_value: float,
    apd: APD,
) -> float:
    """APD output-referred shot-noise RMS over the configured bandwidth."""
    signal_current = primary_photocurrent_a(signal_power_w, apd.responsivity_a_w)
    background_current = primary_photocurrent_a(
        background_power_w_value, apd.responsivity_a_w
    )
    total_primary_current = (
        signal_current + background_current + apd.primary_dark_current_a
    )
    variance_a2 = (
        2.0
        * ELEMENTARY_CHARGE_C
        * apd.noise_bandwidth_hz
        * apd.gain**2
        * apd.excess_noise_factor
        * total_primary_current
    )
    return float(np.sqrt(variance_a2))


def shot_noise_limited_peak_snr(
    signal_power_w: float,
    background_power_w_value: float,
    apd: APD,
) -> float:
    noise_rms_a = shot_noise_rms_a(signal_power_w, background_power_w_value, apd)
    if noise_rms_a == 0.0:
        return np.inf
    return multiplied_signal_current_a(signal_power_w, apd) / noise_rms_a
