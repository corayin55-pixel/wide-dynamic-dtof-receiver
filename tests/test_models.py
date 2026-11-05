from dataclasses import replace

import numpy as np
import pytest

from lidar_sim.models import (
    Channel,
    Laser,
    Receiver,
    Target,
    gaussian_power,
    received_pulse_energy_j,
    round_trip_delay_s,
)


@pytest.fixture
def baseline():
    laser = Laser(905e-9, 500e-9, 5e-9, 3e-3, 1e-3, 0.85)
    target = Target(50.0, 0.1, 0.0, 10.0)
    channel = Channel(0.0)
    receiver = Receiver(25e-3, 0.75)
    return laser, target, channel, receiver


def test_gaussian_pulse_integrates_to_requested_energy():
    time_s = np.linspace(-20e-9, 20e-9, 100_001)
    requested_energy_j = 500e-9
    power_w = gaussian_power(time_s, requested_energy_j, 5e-9)
    assert np.trapezoid(power_w, time_s) == pytest.approx(requested_energy_j, rel=1e-8)


def test_round_trip_delay_at_approximately_150_m_is_one_microsecond():
    assert round_trip_delay_s(149.896229) == pytest.approx(1e-6)


def test_large_target_obeys_inverse_square_law(baseline):
    laser, target, channel, receiver = baseline
    near_energy = received_pulse_energy_j(laser, target, channel, receiver)
    far_energy = received_pulse_energy_j(
        laser, replace(target, range_m=100.0), channel, receiver
    )
    assert near_energy / far_energy == pytest.approx(4.0)


def test_received_energy_scales_with_reflectivity(baseline):
    laser, target, channel, receiver = baseline
    low = received_pulse_energy_j(laser, target, channel, receiver)
    high = received_pulse_energy_j(
        laser, replace(target, reflectivity=0.2), channel, receiver
    )
    assert high / low == pytest.approx(2.0)
