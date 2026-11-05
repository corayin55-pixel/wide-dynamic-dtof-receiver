"""Evaluate the baseline optical return at the Si APD detector."""

from pathlib import Path

import numpy as np

from lidar_sim.config import build_models, load_config
from lidar_sim.detector import (
    APD,
    Background,
    background_power_w,
    mean_primary_electrons,
    multiplied_signal_current_a,
    photons_per_pulse,
    quantum_efficiency,
    shot_noise_limited_peak_snr,
    shot_noise_rms_a,
)
from lidar_sim.models import gaussian_sigma_from_fwhm, received_pulse_energy_j

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    config = load_config(PROJECT_ROOT / "configs" / "baseline.yaml")
    laser, target, channel, receiver = build_models(config)
    apd_cfg = config["apd"]
    background_cfg = config["background"]

    apd = APD(
        responsivity_a_w=apd_cfg["responsivity_a_w"],
        gain=apd_cfg["gain"],
        ionization_ratio=apd_cfg["ionization_ratio"],
        primary_dark_current_a=apd_cfg["primary_dark_current_na"] * 1e-9,
        noise_bandwidth_hz=apd_cfg["noise_bandwidth_mhz"] * 1e6,
    )
    background = Background(
        spectral_radiance_w_m2_sr_nm=background_cfg[
            "spectral_radiance_w_m2_sr_nm"
        ],
        receiver_full_fov_rad=background_cfg["receiver_full_fov_mrad"] * 1e-3,
        filter_bandwidth_nm=background_cfg["filter_bandwidth_nm"],
    )

    return_energy_j = received_pulse_energy_j(laser, target, channel, receiver)
    pulse_sigma_s = gaussian_sigma_from_fwhm(laser.pulse_fwhm_s)
    peak_return_power_w = return_energy_j / (pulse_sigma_s * np.sqrt(2.0 * np.pi))
    background_w = background_power_w(
        background, receiver.aperture_area_m2, receiver.optical_efficiency
    )
    received_photons = photons_per_pulse(return_energy_j, laser.wavelength_m)
    primary_electrons = mean_primary_electrons(
        return_energy_j, laser.wavelength_m, apd.responsivity_a_w
    )
    signal_peak_a = multiplied_signal_current_a(peak_return_power_w, apd)
    noise_rms_a = shot_noise_rms_a(peak_return_power_w, background_w, apd)
    peak_snr = shot_noise_limited_peak_snr(peak_return_power_w, background_w, apd)

    print("Baseline Si APD detector model")
    print(f"Wavelength:                  {laser.wavelength_m * 1e9:12.3f} nm")
    print(f"Received pulse energy:       {return_energy_j * 1e15:12.3f} fJ")
    print(f"Mean received photons:       {received_photons:12.1f}")
    print(f"APD quantum efficiency:      {quantum_efficiency(apd.responsivity_a_w, laser.wavelength_m):12.3%}")
    print(f"Mean primary electrons:      {primary_electrons:12.1f}")
    print(f"APD excess noise factor:     {apd.excess_noise_factor:12.3f}")
    print(f"Background optical power:    {background_w * 1e9:12.3f} nW")
    print(f"Multiplied peak current:     {signal_peak_a * 1e6:12.3f} uA")
    print(f"Shot-noise RMS current:      {noise_rms_a * 1e6:12.3f} uA")
    print(f"Shot-noise-limited peak SNR: {peak_snr:12.3f}")
    print()
    print("Model boundary: APD shot noise only; TIA thermal noise and ADC are excluded.")


if __name__ == "__main__":
    main()
