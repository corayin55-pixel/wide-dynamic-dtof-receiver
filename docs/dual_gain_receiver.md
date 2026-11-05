# Dual-gain receiver validation

## Motivation

The single-channel sweep shows that 20 kOhm transimpedance preserves weak-return
range but saturates on strong returns, while 5 kOhm improves headroom but still
cannot cover a 20 m, 80% target. A 2 kOhm low-gain branch is therefore evaluated
with the 20 kOhm high-gain branch.

## Architecture

Both channels receive the same APD current:

```text
APD current -> 20 kOhm high-gain TIA -> 10-bit ADC (-0.2 to 1.0 V)
            ->  2 kOhm low-gain TIA  -> 10-bit ADC (-0.2 to 1.0 V)
```

Both TIAs retain 250 MHz bandwidth in this system-level study. A practical
implementation would require separate noise, loading, power, matching, and timing
characterization for each analog path.

## Channel selection

The fusion logic selects the high-gain channel when neither TIA nor ADC is
clipped. If the high-gain path clips, it selects the low-gain channel when that
path is valid. If both paths clip, the return is marked saturated. If neither
path passes detection and timing gates, it is a missed detection.

A noiseless per-channel timing calibration removes deterministic path delay.
Templates are generated independently for both channels.

## Interpretation

This simulation tests dynamic-range coverage, not implementation cost. A dual
channel requires additional analog circuitry, ADC resources or multiplexing,
calibration, and channel-switching logic. The recommendation is justified only
if the coverage improvement outweighs those costs.
