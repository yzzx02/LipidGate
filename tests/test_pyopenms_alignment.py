from __future__ import annotations

import pytest

pyopenms = pytest.importorskip("pyopenms")

from lipidbench.runners.run_pyopenms import _feature_rt_bounds


def test_feature_rt_bounds_fall_back_to_detected_feature_width() -> None:
    feature = pyopenms.Feature()
    feature.setRT(300.0)
    feature.setWidth(12.0)

    assert _feature_rt_bounds(feature) == (294.0, 306.0)
