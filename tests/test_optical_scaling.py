from dataclasses import replace

import numpy as np
import pytest

from lidar_sim.models import Channel, Laser, Receiver, Target, received_pulse_energy_j


def test_no_atmosphere_range_exponent_is_minus_two():
    laser = Laser(905e-9, 500e-9, 5e-9, 3e-3, 1e-3, 0.85)
    target = Target(50.0, 0.1, 0.0, 10.0)
    channel = Channel(0.0)
    receiver = Receiver(25e-3, 0.75)
    ranges_m = np.geomspace(1.0, 150.0, 100)
    energies_j = np.array(
        [
            received_pulse_energy_j(
                laser, replace(target, range_m=float(distance)), channel, receiver
            )
            for distance in ranges_m
        ]
    )
    slope, _ = np.polyfit(np.log10(ranges_m), np.log10(energies_j), 1)
    assert slope == pytest.approx(-2.0, abs=1e-12)
