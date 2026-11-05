import numpy as np
import pytest


def quantization_noise_rms_v(lsb_v):
    return lsb_v / np.sqrt(12.0)


def minimum_unsaturated_range(ranges_m, tia_clip, adc_clip, maximum_adc_clip=0.01):
    valid = (tia_clip == 0.0) & (adc_clip <= maximum_adc_clip)
    return float(np.min(ranges_m[valid])) if np.any(valid) else np.nan


def test_ideal_quantization_noise_is_lsb_over_root_twelve():
    assert quantization_noise_rms_v(1e-3) == pytest.approx(1e-3 / np.sqrt(12.0))


def test_minimum_unsaturated_range_uses_both_clipping_limits():
    ranges = np.array([20.0, 30.0, 40.0, 50.0])
    tia_clip = np.array([0.1, 0.0, 0.0, 0.0])
    adc_clip = np.array([1.0, 0.2, 0.005, 0.0])
    assert minimum_unsaturated_range(ranges, tia_clip, adc_clip) == 40.0
