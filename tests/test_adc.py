import numpy as np
import pytest

from lidar_sim.adc import ADC, quantization_error_rms_v, quantize_voltage, sample_waveform


def make_adc(sample_rate_hz=1e9, resolution_bits=8):
    return ADC(sample_rate_hz, resolution_bits, -1.0, 1.0)


def test_adc_lsb_and_raw_range_bin():
    adc = make_adc(1e9, 8)
    assert adc.lsb_v == pytest.approx(2.0 / 255.0)
    assert adc.raw_range_bin_m == pytest.approx(0.149896229)


def test_uniform_sampling_respects_rate_and_phase():
    time_s = np.linspace(0.0, 20e-9, 2001)
    voltage_v = time_s * 1e8
    adc = make_adc(1e9, 10)
    sampled_time_s, sampled_voltage_v = sample_waveform(
        time_s, voltage_v, adc, phase_fraction=0.5
    )
    assert sampled_time_s[0] == pytest.approx(0.5e-9)
    assert np.diff(sampled_time_s) == pytest.approx(np.full(sampled_time_s.size - 1, 1e-9))
    assert sampled_voltage_v == pytest.approx(sampled_time_s * 1e8)


def test_quantizer_clips_and_uses_valid_codes():
    adc = make_adc(1e9, 8)
    voltage_v = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    codes, reconstructed_v, clipped_fraction = quantize_voltage(voltage_v, adc)
    assert codes[0] == 0
    assert codes[-1] == 255
    assert np.all(reconstructed_v >= -1.0)
    assert np.all(reconstructed_v <= 1.0)
    assert clipped_fraction == pytest.approx(0.4)


def test_quantization_error_is_bounded_by_half_lsb_without_clipping():
    adc = make_adc(1e9, 10)
    voltage_v = np.linspace(-0.99, 0.99, 1000)
    _, reconstructed_v, clipped_fraction = quantize_voltage(voltage_v, adc)
    assert clipped_fraction == 0.0
    assert np.max(np.abs(reconstructed_v - voltage_v)) <= adc.lsb_v / 2.0 + 1e-15
    assert quantization_error_rms_v(voltage_v, reconstructed_v) < adc.lsb_v / 2.0
