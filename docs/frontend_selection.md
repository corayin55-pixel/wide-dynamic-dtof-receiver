# TIA gain and ADC range selection

## Design question

A large transimpedance improves weak-return voltage but reduces strong-return
headroom. A wider ADC range reduces clipping but increases the LSB for a fixed
bit depth. This experiment evaluates both parameters jointly rather than
optimizing maximum range alone.

## Fixed parameters

- TIA bandwidth: 250 MHz
- ADC sample rate: 1 GHz
- ADC resolution: 10 bit
- APD gain: 50
- Matched-filter ranging with quadratic peak interpolation

## Candidate grid

- TIA transimpedance: 5, 10, and 20 kOhm
- ADC positive full scale: 0.25, 0.5, and 1.0 V
- ADC negative rail: -20% of positive full scale

The asymmetric ADC range retains negative noise headroom while allocating most
codes to the unipolar return pulse.

## Reported metrics

- 10% target farthest reliable range using detection >= 95% and conditional
  RMSE <= 5 cm;
- 80% target minimum range with no TIA clipping and <= 1% ADC-clipped trials;
- 50 m, 10% target matched-filter RMSE;
- ADC LSB and ideal quantization RMS;
- quantization RMS divided by analog receiver noise RMS.

A single-gain design should not be called acceptable unless both weak-return and
strong-return constraints are stated. If no candidate gives sufficient near and
far coverage, dual-gain or variable-gain reception is the engineering outcome.
