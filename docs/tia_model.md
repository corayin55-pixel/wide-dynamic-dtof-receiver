# TIA receiver model

## Scope

This stage converts the multiplied APD signal current into voltage through a
one-pole transimpedance amplifier. It models bandwidth-dependent pulse shaping,
APD shot noise, TIA input-current noise, feedback-resistor thermal noise,
output-voltage noise, and symmetric output clipping.

The model does not yet include APD/TIA capacitance explicitly, op-amp stability,
non-white noise, recovery from overload, ADC sampling, or PCB parasitics.

## Signal response

The analog transfer function is

```text
H(s) = Rf / (1 + s / (2 pi fc))
```

The implementation uses the exact sampled recursion for a zero-order-held input:

```text
a = exp(-2 pi fc dt)
y[n] = a y[n-1] + (1-a) Rf x[n]
```

## Noise model

For a one-pole low-pass response, equivalent noise bandwidth is

```text
BENBW = pi fc / 2
```

Feedback-resistor input-current noise density is

```text
i_R = sqrt(4 k T / Rf)
```

Independent white-noise densities are summed in quadrature and integrated over
`BENBW`. APD shot-noise density includes signal, background, dark current,
multiplication gain, and excess noise factor.

## Interpretation

A lower bandwidth reduces integrated white noise approximately with `sqrt(B)`,
but broadens the pulse and shifts its peak. A higher bandwidth better preserves
timing information but admits more noise. The best bandwidth must therefore be
selected using ranging performance, not peak SNR alone.
