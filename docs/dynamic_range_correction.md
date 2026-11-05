# Dynamic-range correction

The initial 50--150 m sweep did not apply the configured TIA output swing to the
analog waveform. This is immaterial for the 50 m, 10% baseline but invalid for
strong near-range returns. For example, scaling the baseline return to 80%
reflectivity predicts an unclipped voltage above the TIA and ADC ranges.

The corrected sweep applies TIA clipping before ADC sampling and reports:

- fraction of dense analog samples clipped by the TIA;
- fraction of Monte Carlo trials containing ADC-clipped samples;
- reliable range only when detection and conditional RMSE criteria both pass;
- `>= limit` when performance still passes at the farthest simulated range.

The distance grid is expanded to 20--250 m. The resulting upper range is still a
model boundary, not a hardware specification. Saturated near-range conditions
must be discussed separately from weak-return far-range limits because both can
fail despite high nominal optical SNR.
