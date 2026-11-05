import numpy as np
import pytest

from lidar_sim.ranging import (
    first_threshold_crossing_s,
    matched_filter_time_s,
    quadratic_peak_offset,
    sampled_peak_time_s,
    time_error_to_range_error_m,
)


def test_threshold_crossing_uses_linear_interpolation():
    time_s = np.array([0.0, 1.0, 2.0]) * 1e-9
    voltage_v = np.array([0.0, 0.4, 1.0])
    crossing_s = first_threshold_crossing_s(time_s, voltage_v, 0.7)
    assert crossing_s == pytest.approx(1.5e-9)


def test_threshold_returns_none_when_not_detected():
    time_s = np.arange(5, dtype=float) * 1e-9
    voltage_v = np.zeros(5)
    assert first_threshold_crossing_s(time_s, voltage_v, 0.1) is None


def test_quadratic_interpolation_recovers_subsample_peak():
    sample_positions = np.arange(-3.0, 4.0)
    true_peak = 0.25
    parabola = -(sample_positions - true_peak) ** 2 + 10.0
    peak_index = int(np.argmax(parabola))
    assert quadratic_peak_offset(parabola, peak_index) == pytest.approx(0.25)
    time_s = sample_positions * 1e-9
    assert sampled_peak_time_s(time_s, parabola) == pytest.approx(0.25e-9)


def test_matched_filter_finds_shifted_template_center():
    time_s = np.arange(101, dtype=float) * 1e-9
    template_x = np.arange(-5.0, 6.0)
    template = np.exp(-0.5 * (template_x / 2.0) ** 2)
    waveform = np.zeros(101)
    waveform[55:66] = template
    estimate_s = matched_filter_time_s(time_s, waveform, template, interpolate=False)
    assert estimate_s == pytest.approx(60e-9)


def test_one_nanosecond_time_error_is_about_15_cm():
    assert time_error_to_range_error_m(1e-9) == pytest.approx(0.149896229)
