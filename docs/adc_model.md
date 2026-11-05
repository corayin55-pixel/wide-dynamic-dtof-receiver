# ADC sampling and quantization model

## Scope

This stage applies ideal uniform sampling and quantization to the deterministic
250 MHz TIA output. It evaluates sample-rate and bit-depth effects separately.
Clock jitter, aperture uncertainty, ADC thermal noise, ENOB degradation, and
nonlinearity are intentionally excluded from this first ADC model.

## Quantizer

For input range `[Vmin, Vmax]` and `N` bits:

```text
levels = 2^N
LSB = (Vmax - Vmin) / (levels - 1)
code = round((clip(V) - Vmin) / LSB)
Vq = Vmin + code * LSB
```

The endpoint convention gives exact codes for both rails. It differs slightly
from some ADC datasheet conventions that define one LSB as full scale divided by
`2^N`; the convention must be stated when comparing numerical results.

## Sampling

Uniform sample spacing is `Ts = 1/fs`. The corresponding raw ToF range interval
is

```text
Delta R = c / (2 fs)
```

This is a sample-grid interval, not the final ranging accuracy. Sub-sample timing
estimators can achieve finer precision when waveform SNR and bandwidth permit.

## Baseline design

The baseline ADC range is 0 to 0.25 V because the modeled TIA output is a
single-polarity pulse with approximately 0.178 V peak. This uses the converter
range more effectively than a symmetric +/-1 V range while retaining headroom.
