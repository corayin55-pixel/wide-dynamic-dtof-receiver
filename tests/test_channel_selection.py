

def select_channel(high_valid, high_clipped, low_valid, low_clipped):
    if high_valid and not high_clipped:
        return "high"
    if low_valid and not low_clipped:
        return "low"
    return None


def test_fusion_prefers_unsaturated_high_gain():
    assert select_channel(True, False, True, False) == "high"


def test_fusion_falls_back_to_low_gain_when_high_clips():
    assert select_channel(True, True, True, False) == "low"


def test_fusion_rejects_when_both_channels_clip():
    assert select_channel(True, True, True, True) is None
