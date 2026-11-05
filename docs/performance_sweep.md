# Range and reflectivity performance sweep

## Purpose

This experiment evaluates the selected 250 MHz TIA, 1 GHz ADC, and matched-filter
ranging chain over target range and Lambertian reflectivity. It reports analog
peak SNR, detection rate, conditional ranging RMSE, and conditional 95th-percentile
absolute error.

## Reliable-range criterion

A condition is classified as reliable only when both requirements are met:

```text
detection rate >= 95%
conditional ranging RMSE <= 5 cm
```

The reported reliable range is the farthest simulated grid point satisfying both
requirements. It is not an interpolated hardware specification.

## Important interpretation

RMSE and P95 are conditional on successful detection. Detection rate must always
be reported beside them; otherwise a method that rejects difficult returns can
appear artificially precise.

The experiment uses independent Gaussian ADC-input noise samples, a fixed pulse
template, ideal ADC timing, and a simplified solar-background radiance. Results
are model predictions under stated assumptions, not measured product performance.

The ADC input range is changed from 0--0.25 V to -0.05--0.25 V so zero-mean
receiver noise is not clipped at the lower rail. This slightly changes the LSB
while keeping quantization noise far below the analog receiver noise.
