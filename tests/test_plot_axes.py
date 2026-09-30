import pytest

from lipidgate.gui.plot_axes import (
    MIN_VIEW_SPAN, clamp_window, decimal_tick_step, scaled_window,
    tick_label, tick_values,
)


def test_navigation_mz_ticks_use_clean_decimal_intervals():
    step = decimal_tick_step(244, 1226, 6)
    assert step == 100
    assert list(tick_values(244, 1226, step)) == list(range(300, 1201, 100))
    assert tick_label(300, step) == "300"
    assert decimal_tick_step(10.00012, 10.00082, 6) == pytest.approx(0.0001)
    assert tick_label(10.0003, 0.0001) == "10.0003"


def test_axis_zoom_stops_at_data_bounds_and_reaches_high_precision():
    bounds = (200.0, 800.0)
    assert clamp_window(-4000, 5000, *bounds) == bounds
    assert scaled_window(200, 800, 700, 1.25, *bounds) == bounds
    low, high = scaled_window(200, 800, 700, 1e-8, *bounds)
    assert high - low == pytest.approx(MIN_VIEW_SPAN)
    assert bounds[0] <= low < high <= bounds[1]
