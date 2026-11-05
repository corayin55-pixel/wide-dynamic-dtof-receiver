import numpy as np
import pytest

from lidar_sim.tia import (
    TIA,
    clip_output_voltage,
    first_order_lowpass_transimpedance,
    pulse_fwhm_s,
    total_output_noise_rms_v,
)


def make_tia(bandwidth_hz=250e6):
    return TIA(20e3, bandwidth_hz, 3e-12, 5e-9, 1.0, 300.0)


def test_constant_current_converges_to_dc_transimpedance():
    dt_s = 10e-12
    current_a = np.full(5000, 2e-6)
    output_v = first_order_lowpass_transimpedance(current_a, dt_s, make_tia())
    assert output_v[-1] == pytest.approx(0.04, rel=1e-8)


def test_lower_bandwidth_broadens_gaussian_pulse():
    time_s = np.linspace(-20e-9, 40e-9, 6001)
    sigma_s = 5e-9 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    current_a = np.exp(-0.5 * (time_s / sigma_s) ** 2)
    low_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], make_tia(100e6))
    high_v = first_order_lowpass_transimpedance(current_a, time_s[1] - time_s[0], make_tia(500e6))
    assert pulse_fwhm_s(time_s, low_v) > pulse_fwhm_s(time_s, high_v)


def test_noise_increases_with_square_root_of_bandwidth():
    low = total_output_noise_rms_v(1e-12, make_tia(100e6))
    high = total_output_noise_rms_v(1e-12, make_tia(400e6))
    assert high / low == pytest.approx(2.0)


def test_output_clipping_reports_saturated_samples():
    waveform_v = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    clipped_v, fraction = clip_output_voltage(waveform_v, 1.0)
    assert np.max(np.abs(clipped_v)) == pytest.approx(1.0)
    assert fraction == pytest.approx(0.4)
