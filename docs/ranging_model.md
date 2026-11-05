# Bottom-level ranging algorithms

## Scope

This stage estimates time of flight from the 1 GHz, 10-bit ADC waveform after
the 250 MHz TIA. It compares fixed-threshold crossing, sampled peak detection,
quadratic sub-sample peak interpolation, and matched filtering.

A static timing offset is calibrated from a noiseless baseline waveform. This
represents production calibration of deterministic analog delay. Random timing
error and amplitude-dependent threshold walk remain in the reported results.

## Algorithms

Fixed threshold uses the first rising-edge crossing and linear interpolation
between adjacent ADC samples. It is computationally cheap but its trigger time
changes with pulse amplitude.

Sampled peak selects the largest ADC sample. Its timing is constrained to the
sample grid.

Quadratic interpolation fits a parabola through the peak and its two neighbors:

```text
delta = 0.5 (y[n-1] - y[n+1]) / (y[n-1] - 2 y[n] + y[n+1])
```

Matched filtering cross-correlates the waveform with a known pulse template.
The implementation removes the template mean and normalizes template energy.
The correlation peak is refined with the same quadratic interpolation.

## Monte Carlo assumptions

The first experiment adds independent Gaussian voltage noise at the ADC input
using the TIA model's RMS noise. It randomizes ADC sampling phase and quantizes
to 10 bits. Noise is treated as sample-independent; analog noise correlation,
clock jitter, non-white interference, and false-alarm optimization are deferred.

The raw range interval `c/(2 fs)` is not equivalent to final precision. A
sub-sample estimator may provide lower RMSE when waveform SNR and calibration
support it.
