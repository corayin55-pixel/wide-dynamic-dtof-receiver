# Assumptions and parameter provenance

## Classification

Every model parameter should be tagged as one of:

- **Design**: chosen system parameter to be optimized;
- **Datasheet typical**: typical value from a named part and document revision;
- **Datasheet limit**: guaranteed minimum or maximum;
- **Literature**: value or equation from a cited publication;
- **Calibration**: estimated from measured data;
- **Simplification**: modeling assumption without a hardware guarantee.

## Current parameter status

| Parameter | Current value | Current class | Required final source |
|---|---:|---|---|
| Laser wavelength | 905 nm | Design | Candidate laser datasheet |
| Pulse energy | 500 nJ | Design | Laser + eye-safety budget |
| Pulse FWHM | 5 ns | Design | Candidate laser datasheet |
| Beam divergence | 1 mrad | Design | Transmit-optics budget |
| Receiver aperture | 25 mm | Design | Mechanical/optical design |
| Optical efficiencies | 0.85 / 0.75 | Simplification | Coating/filter/lens budget |
| Atmospheric extinction | 0.0001 1/m | Simplification | Visibility model or literature |
| Spectral background radiance | 0.30 W/m2/sr/nm | Simplification | Solar spectrum and scene model |
| APD responsivity | 0.50 A/W | Design baseline | Named APD datasheet |
| APD gain | 50 | Design baseline | Named APD datasheet + bias |
| APD ionization ratio | 0.02 | Literature baseline | APD technology reference |
| Primary dark current | 5 nA | Design baseline | Named APD datasheet |
| TIA noise densities | 3 pA/sqrtHz, 5 nV/sqrtHz | Design baseline | Named TIA/AFE datasheet |
| ADC resolution/rate | 10 bit / 1 GHz | Design | Named ADC datasheet |

## Rules for the final report

1. Do not call a baseline parameter a selected component specification.
2. Record manufacturer, full part number, document revision, page, test
   conditions, and whether a value is typical or guaranteed.
3. Preserve parameter ranges when temperature or bias changes the value.
4. Do not convert illuminance in lux directly to 905 nm radiance without an
   explicit source spectrum and photopic weighting calculation.
5. Do not claim IEC 60825-1 compliance from pulse energy alone.

## Recommended component-screening fields

### Pulsed laser

Wavelength, pulse energy, pulse width, peak current, repetition rate, divergence,
temperature drift, package, monitor diode, eye-safety implications, and cost.

### Si APD

Responsivity at 905 nm, multiplication gain versus bias and temperature, excess
noise, dark current, capacitance, active area, bandwidth, breakdown voltage,
package, and automotive/industrial qualification.

### TIA or AFE

Transimpedance, bandwidth versus detector capacitance, input-current noise,
input-voltage noise, overload recovery, output swing, power, stability, and gain
switching support.

### ADC

Sample rate, nominal bits, ENOB at the relevant input frequency, input range,
input bandwidth, clock jitter, latency, power, interface data rate, and package.
