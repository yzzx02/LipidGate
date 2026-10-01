"""Shared numeric bounds for interactive result plots."""

from __future__ import annotations

import math


MIN_VIEW_SPAN = 1e-4


def axis_drag_factor(previous: float, current: float, axis: str) -> float:
    """Right/up zooms in, left/down zooms out, including across the center."""
    if axis not in {"x", "y"}:
        raise ValueError("Axis must be x or y")
    delta = current - previous
    exponent = delta * (-0.012 if axis == "x" else 0.012)
    return math.exp(max(-3.0, min(3.0, exponent)))


def decimal_tick_step(low: float, high: float, target_count: float = 6) -> float:
    """Choose a power-of-ten interval: 100, 10, 1, 0.1, ..."""
    span = max(float(high) - float(low), MIN_VIEW_SPAN)
    raw = span / max(float(target_count), 1)
    return 10.0 ** math.floor(math.log10(raw) + 0.5)


def tick_label(value: float, step: float) -> str:
    decimals = max(0, min(8, -int(math.floor(math.log10(step) + 1e-12))))
    rendered = f"{value:.{decimals}f}"
    return "0" if rendered.startswith("-0") and abs(value) < step * 1e-6 else rendered


def intensity_tick_label(value: float, step: float) -> str:
    """Use actual intensity values, with compact scientific notation if needed."""
    if value and (abs(value) >= 100_000 or abs(value) < 0.001):
        exponent = math.floor(math.log10(abs(value)))
        return f"{value / 10.0 ** exponent:.3g}e{exponent}"
    return tick_label(value, step)


def tick_values(low: float, high: float, step: float):
    first = math.ceil(low / step - 1e-10)
    last = math.floor(high / step + 1e-10)
    for index in range(first, min(last, first + 200) + 1):
        yield index * step


def clamp_window(low: float, high: float, bound_low: float, bound_high: float,
                 min_span: float = MIN_VIEW_SPAN) -> tuple[float, float]:
    """Keep a view within its data range while preserving its requested span."""
    values = (low, high, bound_low, bound_high)
    if not all(math.isfinite(value) for value in values) or bound_high <= bound_low:
        raise ValueError("Plot bounds must be finite and increasing")
    capacity = bound_high - bound_low
    span = min(capacity, max(min_span, high - low))
    if span >= capacity:
        return bound_low, bound_high
    start = max(bound_low, min(low, bound_high - span))
    return start, start + span


def scaled_window(low: float, high: float, anchor: float, factor: float,
                  bound_low: float, bound_high: float,
                  min_span: float = MIN_VIEW_SPAN) -> tuple[float, float]:
    if factor <= 0 or not math.isfinite(factor):
        raise ValueError("Zoom factor must be positive")
    span = high - low
    fraction = (anchor - low) / span if span else 0.5
    new_span = max(min_span, span * factor)
    start = anchor - fraction * new_span
    return clamp_window(start, start + new_span, bound_low, bound_high, min_span)


def scaled_from_minimum(low: float, high: float, factor: float,
                        bound_low: float, bound_high: float) -> tuple[float, float]:
    """Change vertical scale with the original lowest value pinned in place."""
    return scaled_window(bound_low, bound_low + high - low, bound_low, factor,
                         bound_low, bound_high)
