# Si APD detector model

## Scope

This stage converts received optical energy and peak power into photon count,
primary photoelectrons, multiplied APD current, and shot-noise-limited peak SNR.
The TIA, ADC, temperature dependence, saturation, and bandwidth-shaped waveform
are intentionally excluded.

## Equations

Photon energy and pulse photon count:

```text
E_ph = h c / lambda
N_ph = E_rx / E_ph
```

Responsivity and external quantum efficiency:

```text
eta_q = R_lambda h c / (q lambda)
N_e = eta_q N_ph
```

For APD gain `M`, the multiplied signal current is:

```text
I_sig = M R_lambda P_sig
```

The McIntyre excess noise factor is approximated by:

```text
F(M) = k M + (1 - k) (2 - 1/M)
```

The current shot-noise variance over equivalent noise bandwidth `B` is:

```text
sigma_i^2 = 2 q B M^2 F(M) (I_sig,primary + I_bg,primary + I_dark)
```

Background power uses the receiver etendue approximation:

```text
P_bg = L_lambda A_rx Omega_FOV Delta_lambda eta_rx
```

## Baseline assumptions

The baseline spectral radiance is a design assumption for model bring-up, not a
conversion from illuminance in klx and not a measured solar spectrum. APD dark
current is treated as entirely multiplied primary current. Both assumptions must
be replaced or refined during component selection.
