"""Uniform sampling and ideal quantization models for the LiDAR ADC."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class ADC:
    sample_rate_hz: float
    resolution_bits: int
    input_min_v: float
    input_max_v: float

    @property
    def levels(self) -> int:
        return 2**self.resolution_bits

    @property
    def lsb_v(self) -> float:
        return (self.input_max_v - self.input_min_v) / (self.levels - 1)

    @property
    def sample_period_s(self) -> float:
        return 1.0 / self.sample_rate_hz

    @property
    def raw_range_bin_m(self) -> float:
        return 299_792_458.0 / (2.0 * self.sample_rate_hz)

    def validate(self) -> None:
        if self.sample_rate_hz <= 0.0:
            raise ValueError("ADC sample rate must be positive")
        if self.resolution_bits < 1:
            raise ValueError("ADC resolution must be at least one bit")
        if self.input_max_v <= self.input_min_v:
            raise ValueError("ADC input maximum must exceed input minimum")


def sample_waveform(
    source_time_s: NDArray[np.float64],
    source_voltage_v: NDArray[np.float64],
    adc: ADC,
    phase_fraction: float = 0.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Uniformly sample a dense waveform using linear interpolation."""
    adc.validate()
    if source_time_s.size != source_voltage_v.size or source_time_s.size < 2:
        raise ValueError("Source time and voltage arrays must have equal length >= 2")
    if not 0.0 <= phase_fraction < 1.0:
        raise ValueError("Sampling phase fraction must be within [0, 1)")
    if np.any(np.diff(source_time_s) <= 0.0):
        raise ValueError("Source time must be strictly increasing")

    start_s = source_time_s[0] + phase_fraction * adc.sample_period_s
    sample_time_s = np.arange(start_s, source_time_s[-1], adc.sample_period_s)
    sample_voltage_v = np.interp(sample_time_s, source_time_s, source_voltage_v)
    return sample_time_s, sample_voltage_v


def quantize_voltage(
    voltage_v: NDArray[np.float64], adc: ADC
) -> tuple[NDArray[np.int64], NDArray[np.float64], float]:
    """Clip and quantize voltage; return codes, reconstructed voltage, clipping fraction."""
    adc.validate()
    voltage_array = np.asarray(voltage_v, dtype=float)
    clipped = np.clip(voltage_array, adc.input_min_v, adc.input_max_v)
    codes = np.rint((clipped - adc.input_min_v) / adc.lsb_v).astype(np.int64)
    codes = np.clip(codes, 0, adc.levels - 1)
    reconstructed_v = adc.input_min_v + codes * adc.lsb_v
    clipped_fraction = float(
        np.mean((voltage_array < adc.input_min_v) | (voltage_array > adc.input_max_v))
    )
    return codes, reconstructed_v, clipped_fraction


def quantization_error_rms_v(
    original_v: NDArray[np.float64], reconstructed_v: NDArray[np.float64]
) -> float:
    if original_v.shape != reconstructed_v.shape:
        raise ValueError("Original and reconstructed arrays must have the same shape")
    return float(np.sqrt(np.mean((reconstructed_v - original_v) ** 2)))
