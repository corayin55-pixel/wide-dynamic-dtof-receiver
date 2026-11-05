# 905 nm Pulsed dToF LiDAR System Simulator

A reproducible system-level simulation of one outdoor 905 nm pulsed direct
time-of-flight (dToF) LiDAR channel. The project connects optical-link physics,
Si APD detection, TIA/ADC design, bottom-level ranging algorithms, Monte Carlo
performance evaluation, and receiver dynamic-range architecture.

> Status: simulation validated by automated consistency tests. Results are model
> predictions under documented assumptions, not hardware measurements or product
> specifications.

## Engineering question

How should the optical receiver, detector gain, analog bandwidth, ADC, and
ranging algorithm be configured to detect low-reflectivity distant targets while
avoiding saturation on high-reflectivity near targets?

## Signal chain

```text
905 nm Gaussian pulse
  -> atmospheric propagation and Lambertian target
  -> receiver aperture and optical efficiency
  -> Si APD photon conversion and multiplication noise
  -> first-order TIA pulse shaping and analog noise
  -> ADC sampling, clipping, and quantization
  -> threshold / peak / matched-filter ranging
  -> single-gain and dual-gain receiver evaluation
```

## Baseline design

| Parameter | Baseline |
|---|---:|
| Wavelength | 905 nm |
| Pulse energy / FWHM | 500 nJ / 5 ns |
| Beam divergence | 1 mrad full angle |
| Receiver aperture | 25 mm |
| APD responsivity / gain | 0.50 A/W / 50 |
| TIA bandwidth | 250 MHz |
| ADC | 1 GHz, 10 bit |
| Dual-gain transimpedance | 20 kOhm / 2 kOhm |
| Dual-gain ADC range | -0.2 to 1.0 V |

These are simulation design parameters. They are not all tied to selected
commercial part numbers yet.

## Key results

- The deterministic optical model reproduces the inverse-square range law with a
  fitted log-log exponent of `-2.000` when atmospheric extinction is disabled.
- At 50 m and 10% reflectivity, the baseline return contains about `8,986`
  photons and produces about `6,155` primary photoelectrons.
- A 250 MHz TIA provides a practical trade-off: `5.189 ns` output FWHM,
  `26.409 dB` peak SNR, and no clipping for the 50 m / 10% baseline.
- A 1 GHz, 10-bit ADC has a `14.99 cm` raw range interval, but matched filtering
  and sub-sample interpolation reduce baseline Monte Carlo RMSE to `1.173 cm` in
  the initial 50 m experiment.
- Single-gain reception cannot cover both near high-reflectivity returns and far
  low-reflectivity returns without saturation or sensitivity loss.
- The corrected 20 kOhm / 2 kOhm dual-gain model gives a reliable simulated
  interval of `20-160 m` for 10% reflectivity and `20-250 m` for 80%
  reflectivity on the tested grid, with no simultaneous channel clipping and a
  maximum false-time selection rate of `1.6%`.

A condition is called reliable only when detection rate is at least 95% and
conditional ranging RMSE is at most 5 cm. A range ending at 250 m means only that
the condition passed at the simulation boundary; it is not a measured maximum.

## Selected figures

### Optical-link validation

![Pulse waveforms](results/figures/01_pulse_waveforms.png)

![Received energy versus range](results/figures/02_received_energy_vs_range.png)

### Receiver and ADC trade-offs

![TIA bandwidth trade-off](results/figures/04_tia_bandwidth_tradeoff.png)

![ADC sampling trade-off](results/figures/05_adc_sampling_tradeoff.png)

### Ranging and architecture

![Ranging Monte Carlo](results/figures/06_ranging_monte_carlo.png)

![Receiver dynamic range](results/figures/08_receiver_dynamic_range.png)

![Frontend selection](results/figures/09_frontend_selection.png)

![Dual-gain receiver](results/figures/10_dual_gain_receiver.png)

## Reproduce

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -v

python scripts/run_baseline.py
python scripts/run_optical_validation.py
python scripts/run_detector_baseline.py
python scripts/run_tia_bandwidth_study.py
python scripts/run_adc_study.py
python scripts/run_ranging_study.py
python scripts/run_range_performance_sweep.py
python scripts/run_frontend_selection.py
python scripts/run_dual_gain_study.py
```

Monte Carlo scripts use fixed random seeds. CSV data are written to
`results/data/`; PNG and PDF figures are written to `results/figures/`.

## Repository structure

```text
configs/                System design parameters
src/lidar_sim/          Optical, detector, TIA, ADC, and ranging models
scripts/                Reproducible experiments
results/data/           Numerical outputs
results/figures/        Report figures
 tests/                  Physical and numerical consistency tests
 docs/                   Model derivations, assumptions, and report
```

## Model boundaries

- Lambertian flat target, coaxial far-field receiver, and fully captured return
- Parameterized spectral background rather than measured 100 klx spectrum
- White Gaussian sample-domain noise approximation
- First-order TIA rather than transistor/op-amp or SPICE model
- Ideal ADC apart from clipping, sampling, and quantization
- No clock jitter, temperature drift, PCB parasitics, crosstalk, rain/fog,
  multipath, eye-safety certification, or hardware validation

See `docs/technical_report_cn.md` and `docs/assumptions_and_sources.md` before
interpreting the numerical results.
