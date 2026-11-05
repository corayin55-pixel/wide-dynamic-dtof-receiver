import numpy as np
import pytest

from lidar_sim.detector import (
    APD,
    Background,
    background_power_w,
    mean_primary_electrons,
    photon_energy_j,
    photons_per_pulse,
    quantum_efficiency,
    shot_noise_rms_a,
)


def test_905_nm_photon_energy():
    assert photon_energy_j(905e-9) == pytest.approx(2.195e-19, rel=5e-4)


def test_responsivity_quantum_efficiency_round_trip():
    efficiency = quantum_efficiency(0.50, 905e-9)
    assert efficiency == pytest.approx(0.685, rel=2e-3)
    photons = photons_per_pulse(1e-15, 905e-9)
    electrons = mean_primary_electrons(1e-15, 905e-9, 0.50)
    assert electrons / photons == pytest.approx(efficiency)


def test_background_power_scales_with_filter_bandwidth():
    narrow = Background(0.30, 2e-3, 5.0)
    wide = Background(0.30, 2e-3, 10.0)
    aperture_area_m2 = np.pi * (25e-3 / 2.0) ** 2
    narrow_power = background_power_w(narrow, aperture_area_m2, 0.75)
    wide_power = background_power_w(wide, aperture_area_m2, 0.75)
    assert wide_power / narrow_power == pytest.approx(2.0)


def test_apd_noise_is_positive_and_increases_with_background():
    apd = APD(0.50, 50.0, 0.02, 5e-9, 250e6)
    low_noise = shot_noise_rms_a(1e-7, 0.0, apd)
    high_noise = shot_noise_rms_a(1e-7, 1e-8, apd)
    assert low_noise > 0.0
    assert high_noise > low_noise
