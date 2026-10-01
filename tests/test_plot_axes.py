import pytest

from lipidgate.gui.plot_axes import (
    MIN_VIEW_SPAN, axis_drag_factor, clamp_window, decimal_tick_step,
    intensity_tick_label, scaled_from_minimum, scaled_window, tick_label, tick_values,
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


@pytest.mark.parametrize("axis,start,end", [
    ("x", 100, 140), ("x", 60, 80), ("x", 80, 120),
    ("y", 100, 60), ("y", 150, 120), ("y", 120, 80),
])
def test_axis_drag_direction_does_not_change_across_center(axis, start, end):
    factor = axis_drag_factor(start, end, axis)
    assert 0 < factor < 1
    assert axis_drag_factor(end, start, axis) > 1
    assert factor * axis_drag_factor(end, start, axis) == pytest.approx(1)


@pytest.mark.parametrize("bounds", [(0, 105), (244, 1226)])
def test_vertical_scale_keeps_original_minimum_even_after_pan(bounds):
    low, high = scaled_from_minimum(bounds[0] + 20, bounds[1] - 20, 0.5, *bounds)
    assert low == bounds[0]
    assert high - low == pytest.approx((bounds[1] - bounds[0] - 40) * 0.5)
    assert scaled_from_minimum(low, high, 1000, *bounds) == bounds


@pytest.mark.parametrize("value,step,label", [
    (0, 1000, "0"), (1000, 1000, "1000"), (10000, 2000, "10000"),
    (100000, 20000, "1e5"), (2500000, 500000, "2.5e6"), (0.0005, 0.0001, "5e-4"),
])
def test_intensity_ticks_show_numbers_without_abbreviated_units(value, step, label):
    assert intensity_tick_label(value, step) == label
