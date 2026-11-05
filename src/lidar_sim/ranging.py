"""Time-of-flight estimators for sampled LiDAR return waveforms."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def first_threshold_crossing_s(
    time_s: NDArray[np.float64], voltage_v: NDArray[np.float64], threshold_v: float
) -> float | None:
    """Return the rising-edge threshold crossing with linear interpolation."""
    if time_s.size != voltage_v.size or time_s.size < 2:
        raise ValueError("Time and voltage arrays must have equal length >= 2")
    crossings = np.flatnonzero(
        (voltage_v[:-1] < threshold_v) & (voltage_v[1:] >= threshold_v)
    )
    if crossings.size == 0:
        return None
    index = int(crossings[0])
    delta_v = voltage_v[index + 1] - voltage_v[index]
    if delta_v == 0.0:
        return float(time_s[index])
    fraction = (threshold_v - voltage_v[index]) / delta_v
    return float(time_s[index] + fraction * (time_s[index + 1] - time_s[index]))


def quadratic_peak_offset(samples: NDArray[np.float64], peak_index: int) -> float:
    """Estimate a peak offset in samples from three adjacent points."""
    if peak_index <= 0 or peak_index >= samples.size - 1:
        return 0.0
    left = samples[peak_index - 1]
    center = samples[peak_index]
    right = samples[peak_index + 1]
    denominator = left - 2.0 * center + right
    if denominator == 0.0:
        return 0.0
    offset = 0.5 * (left - right) / denominator
    return float(np.clip(offset, -1.0, 1.0))


def sampled_peak_time_s(
    time_s: NDArray[np.float64], samples: NDArray[np.float64], interpolate: bool = True
) -> float:
    if time_s.size != samples.size or time_s.size < 3:
        raise ValueError("Time and sample arrays must have equal length >= 3")
    peak_index = int(np.argmax(samples))
    if not interpolate or peak_index == 0 or peak_index == samples.size - 1:
        return float(time_s[peak_index])
    offset = quadratic_peak_offset(samples, peak_index)
    sample_period_s = time_s[1] - time_s[0]
    return float(time_s[peak_index] + offset * sample_period_s)


def matched_filter_scores(
    voltage_v: NDArray[np.float64], template_v: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Cross-correlate a waveform with an odd-length, zero-mean template."""
    if template_v.size < 3 or template_v.size % 2 == 0:
        raise ValueError("Matched-filter template must have odd length >= 3")
    if voltage_v.size < template_v.size:
        raise ValueError("Waveform must not be shorter than the template")
    centered_template = template_v - np.mean(template_v)
    norm = np.linalg.norm(centered_template)
    if norm == 0.0:
        raise ValueError("Matched-filter template cannot be constant")
    return np.correlate(voltage_v, centered_template / norm, mode="same")


def matched_filter_time_s(
    time_s: NDArray[np.float64],
    voltage_v: NDArray[np.float64],
    template_v: NDArray[np.float64],
    interpolate: bool = True,
) -> float:
    scores = matched_filter_scores(voltage_v, template_v)
    return sampled_peak_time_s(time_s, scores, interpolate=interpolate)


def time_error_to_range_error_m(time_error_s: float) -> float:
    return 0.5 * 299_792_458.0 * time_error_s
