# Dual-gain fusion correction

The first dual-gain screening implementation included the simulated true ToF in
the per-channel validity decision. That makes channel selection unavailable to a
real receiver and can cause misleading fallback behavior.

The corrected online decision uses only observable quantities:

- TIA/ADC clipping state;
- matched-filter correlation peak relative to its threshold.

The high-gain channel remains preferred when valid and unsaturated. The low-gain
channel is used only when the high-gain channel is clipped or below its detection
threshold.

Ground-truth ToF is now used only after channel selection for offline scoring. A
selected estimate farther than 5 ns from truth is counted as a false-time
selection rather than a correct detection. The CSV and figure report that failure
rate separately. This distinction prevents an algorithm from using evaluation
ground truth in its runtime logic.
