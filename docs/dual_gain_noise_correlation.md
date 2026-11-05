# Dual-gain noise-correlation correction

Both gain channels observe the same APD current, so APD shot noise is a shared
physical source. The initial screening generated independent total-noise samples
and independent ADC sampling phases for each path, which could overstate the
benefit of selecting between channels.

The corrected simulation uses:

- one ADC sampling phase shared by both channels in each trial;
- one standard Gaussian APD-noise sequence shared by both channels and scaled by
  each channel's transimpedance;
- independent TIA-internal noise sequences for the two analog paths;
- each path's total RMS noise for its own detection threshold.

The APD/TIA split is computed from RMS variances. This remains a white-noise,
ADC-sample-domain approximation; analog temporal correlation and crosstalk are
not modeled. It is nevertheless physically more consistent for dual-channel
architecture comparison than independent total noise.
