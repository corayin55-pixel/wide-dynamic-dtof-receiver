"""First-order transimpedance amplifier model for the LiDAR receiver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

BOLTZMANN_CONSTANT_J_K = 1.380_649e-23


@dataclass(frozen=True)
class TIA:
    transimpedance_ohm: float
    bandwidth_hz: float
    input_current_noise_a_sqrt_hz: float
    output_voltage_noise_v_sqrt_hz: float
    output_swing_v: float
    temperature_k: float = 300.0

    @property
    def equivalent_noise_bandwidth_hz(self) -> float:
        """Equivalent noise bandwidth of a one-pole low-pass filter."""
        return 0.5 * np.pi * self.bandwidth_hz

    @property
    def feedback_resistor_noise_a_sqrt_hz(self) -> float:
        return np.sqrt(
            4.0 * BOLTZMANN_CONSTANT_J_K * self.temperature_k
            / self.transimpedance_ohm
        )


def first_order_lowpass_transimpedance(
    current_a: NDArray[np.float64], sample_interval_s: float, tia: TIA
) -> NDArray[np.float64]:
    """Apply an exact discrete-time one-pole transimpedance response."""
    if sample_interval_s <= 0.0:
        raise ValueError("Sample interval must be positive")
    if tia.bandwidth_hz <= 0.0 or tia.transimpedance_ohm <= 0.0:
        raise ValueError("TIA bandwidth and transimpedance must be positive")

    decay = np.exp(-2.0 * np.pi * tia.bandwidth_hz * sample_interval_s)
    output_v = np.empty_like(current_a, dtype=float)
    state_v = 0.0
    for index, sample_a in enumerate(current_a):
        state_v = decay * state_v + (1.0 - decay) * tia.transimpedance_ohm * sample_a
        output_v[index] = state_v
    return output_v


def clip_output_voltage(
    voltage_v: NDArray[np.float64], output_swing_v: float
) -> tuple[NDArray[np.float64], float]:
    """Clip to a symmetric output swing and return the clipped sample fraction."""
    if output_swing_v <= 0.0:
        raise ValueError("Output swing must be positive")
    clipped_v = np.clip(voltage_v, -output_swing_v, output_swing_v)
    clipped_fraction = float(np.mean(np.abs(voltage_v) > output_swing_v))
    return clipped_v, clipped_fraction


def total_output_noise_rms_v(apd_noise_a_sqrt_hz: float, tia: TIA) -> float:
    """Integrate independent APD and TIA white-noise sources over one-pole ENBW."""
    total_input_density_squared = (
        apd_noise_a_sqrt_hz**2
        + tia.input_current_noise_a_sqrt_hz**2
        + tia.feedback_resistor_noise_a_sqrt_hz**2
        + (tia.output_voltage_noise_v_sqrt_hz / tia.transimpedance_ohm) ** 2
    )
    input_noise_rms_a = np.sqrt(
        total_input_density_squared * tia.equivalent_noise_bandwidth_hz
    )
    return float(input_noise_rms_a * tia.transimpedance_ohm)


def pulse_fwhm_s(time_s: NDArray[np.float64], waveform: NDArray[np.float64]) -> float:
    """Estimate pulse FWHM by linearly interpolating both half-height crossings."""
    if time_s.size != waveform.size or time_s.size < 3:
        raise ValueError("Time and waveform arrays must have equal length >= 3")
    peak_index = int(np.argmax(waveform))
    half_height = 0.5 * waveform[peak_index]

    left_candidates = np.flatnonzero(waveform[:peak_index] < half_height)
    right_candidates = np.flatnonzero(waveform[peak_index + 1 :] < half_height)
    if left_candidates.size == 0 or right_candidates.size == 0:
        raise ValueError("Waveform window does not contain both half-height crossings")

    left_low = int(left_candidates[-1])
    left_high = left_low + 1
    right_high = peak_index + 1 + int(right_candidates[0])
    right_low = right_high - 1

    left_time = np.interp(
        half_height,
        [waveform[left_low], waveform[left_high]],
        [time_s[left_low], time_s[left_high]],
    )
    right_time = np.interp(
        half_height,
        [waveform[right_high], waveform[right_low]],
        [time_s[right_high], time_s[right_low]],
    )
    return float(right_time - left_time)
