import numpy as np


def reliable_range(ranges_m, detection_rates, rmse_m, minimum_detection=0.95, maximum_rmse_m=0.05):
    valid = (detection_rates >= minimum_detection) & (rmse_m <= maximum_rmse_m)
    return float(np.max(ranges_m[valid])) if np.any(valid) else np.nan


def test_reliable_range_uses_both_detection_and_rmse_requirements():
    ranges = np.array([10.0, 20.0, 30.0, 40.0])
    detection = np.array([1.0, 0.99, 0.94, 0.99])
    rmse = np.array([0.01, 0.04, 0.03, 0.06])
    assert reliable_range(ranges, detection, rmse) == 20.0


def test_reliable_range_is_nan_when_no_condition_passes():
    ranges = np.array([10.0, 20.0])
    detection = np.array([0.5, 0.4])
    rmse = np.array([0.1, 0.2])
    assert np.isnan(reliable_range(ranges, detection, rmse))
