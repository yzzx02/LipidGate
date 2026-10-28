"""Device-pixel-aware measured MS/MS spectrum for the result browser."""

from __future__ import annotations

import html
import math

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .plot_axes import MIN_VIEW_SPAN, axis_drag_factor, clamp_window, decimal_tick_step, intensity_tick_label, scaled_from_minimum, scaled_window, tick_label, tick_values
from .result_plots import FRAGMENT_COLORS, FRAGMENT_LABELS, display_role
from .ui_icons import plot_cursor, reset_view_button


def _nice_step(span: float) -> float:
    raw = max(span / 5.0, 1e-9)
    power = 10 ** math.floor(math.log10(raw))
    for multiplier in (1, 2, 2.5, 5, 10):
        if raw <= multiplier * power:
            return multiplier * power
    return 10 * power


class _SpectrumCanvas(QtWidgets.QWidget):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.setMinimumSize(260, 185)
        self.setMouseTracking(True)
        self._press = None
        self._hover_peak = None
        self._right_axis = None
        self._left_axis = None
        self._drag_button = None
        self._last = None
        self._axis_anchor = None

    def plot_rect(self):
        return QtCore.QRectF(84, 17, max(1, self.width() - 99), max(1, self.height() - 76))

    def _pixel(self, mz: float, intensity: float):
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        return (
            rect.left() + (mz - x0) / (x1 - x0) * rect.width(),
            rect.bottom() - (intensity - y0) / (y1 - y0) * rect.height(),
        )

    def _data(self, point):
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        return (
            x0 + (point.x() - rect.left()) / rect.width() * (x1 - x0),
            y0 + (rect.bottom() - point.y()) / rect.height() * (y1 - y0),
        )

    def _font(self, size: float, bold: bool = False):
        font = self.font()
        font.setPointSizeF(size)
        font.setBold(bold)
        return font

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing)
        painter.fillRect(self.rect(), QtGui.QColor("#ffffff"))
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        painter.setFont(self._font(10.2))

        # Draw ticks at native screen resolution; leave the data field ungridded.
        tick_pen = QtGui.QPen(QtGui.QColor("#657487"), 1.2)
        tick_text = QtGui.QColor("#4b5864")
        y_step = _nice_step(y1 - y0)
        intensity = math.ceil(y0 / y_step) * y_step
        while intensity <= y1 + 1e-9:
            y = self._pixel(x0, intensity)[1]
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(rect.left() - 5, y), QtCore.QPointF(rect.left(), y))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(30, y - 10, rect.left() - 39, 20), QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, intensity_tick_label(intensity, y_step))
            intensity += y_step

        step = decimal_tick_step(x0, x1, max(3.0, rect.width() / 90.0))
        for tick in tick_values(x0, x1, step):
            x = self._pixel(tick, y0)[0]
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(x, rect.bottom()), QtCore.QPointF(x, rect.bottom() + 5))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(x - 33, rect.bottom() + 5, 66, 21), QtCore.Qt.AlignmentFlag.AlignCenter, tick_label(tick, step))

        painter.setPen(QtGui.QPen(QtGui.QColor("#596b7d"), 1.55))
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())
        painter.drawLine(rect.bottomLeft(), rect.topLeft())
        painter.setPen(QtGui.QColor("#243746"))
        painter.setFont(self._font(10.5, bold=True))
        painter.drawText(QtCore.QRectF(rect.left(), rect.bottom() + 27, rect.width(), 22), QtCore.Qt.AlignmentFlag.AlignCenter, "m/z")
        painter.save()
        painter.translate(23, rect.center().y())
        painter.rotate(-90)
        painter.drawText(QtCore.QRectF(-rect.height() / 2, -14, rect.height(), 28), QtCore.Qt.AlignmentFlag.AlignCenter, self.owner.y_axis_label)
        painter.restore()

        if not self.owner.has_peaks:
            painter.setPen(QtGui.QColor("#75838e"))
            painter.setFont(self._font(10.5))
            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter | QtCore.Qt.TextFlag.TextWordWrap,
                             self.owner.status)
            return

        painter.save()
        painter.setClipRect(rect)
        baseline = self._pixel(x0, 0)[1]
        # Draw unassigned peaks first so a weak matched peak remains visible.
        visible = self.owner.visible_peaks(max(1, int(rect.width())))
        for emphasized in (False, True):
            for peak in visible:
                matched = bool(peak["fragment"])
                precursor = self.owner.is_precursor_peak(peak)
                if (matched or precursor) != emphasized:
                    continue
                x, y = self._pixel(peak["mz"], peak["relative"])
                if x < rect.left() or x > rect.right():
                    continue
                role = display_role(peak["fragment"]) if matched else "unmatched"
                color = QtGui.QColor(("#c55243" if precursor else "#525b62")
                                     if self.owner.ms1_spectrum is not None else FRAGMENT_COLORS[role])
                pen = QtGui.QPen(color, 2.1 if matched else 1.6 if precursor else 1.2)
                pen.setCapStyle(QtCore.Qt.PenCapStyle.FlatCap)
                painter.setPen(pen)
                painter.drawLine(QtCore.QPointF(x, baseline), QtCore.QPointF(x, y))
                if matched:
                    painter.setPen(QtCore.Qt.PenStyle.NoPen)
                    painter.setBrush(color)
                    # The stick retains its true height; the dot marks even a
                    # sub-1% match that would otherwise disappear at the axis.
                    painter.drawEllipse(QtCore.QPointF(x, min(y, baseline - 4)), 3.2, 3.2)
                    painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.restore()

    def _axis_at(self, position):
        rect = self.plot_rect()
        if rect.left() <= position.x() <= rect.right() and rect.bottom() < position.y() <= self.height():
            return "x"
        if 0 <= position.x() < rect.left() and rect.top() <= position.y() <= rect.bottom():
            return "y"
        return None

    def _cursor_at(self, position):
        axis = (self._right_axis or self._left_axis) if self._press is not None else self._axis_at(position)
        self.setCursor(plot_cursor(axis) if axis else QtGui.QCursor(QtCore.Qt.CursorShape.ArrowCursor))

    def mousePressEvent(self, event):
        position = event.position()
        axis = self._axis_at(position)
        if axis is None or event.button() not in {
            QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.RightButton,
        }:
            return
        self._right_axis = axis if event.button() == QtCore.Qt.MouseButton.RightButton else None
        self._left_axis = axis if event.button() == QtCore.Qt.MouseButton.LeftButton else None
        self._drag_button = event.button()
        self._press = position
        self._last = position
        if self._right_axis is not None:
            x, y = self._data(position)
            self._axis_anchor = x if self._right_axis == "x" else y
        self._cursor_at(position)
        event.accept()

    def mouseMoveEvent(self, event):
        self._cursor_at(event.position())
        if self._press is not None:
            position = event.position()
            if self._right_axis is not None:
                coordinate = "x" if self._right_axis == "x" else "y"
                current = position.x() if coordinate == "x" else position.y()
                previous = self._last.x() if coordinate == "x" else self._last.y()
                delta = current - previous
                if delta:
                    factor = axis_drag_factor(previous, current, coordinate)
                    self.owner.zoom_axis(self._right_axis, self._axis_anchor, factor)
            elif self._left_axis is not None:
                x0, x1, y0, y1 = self.owner.limits
                rect = self.plot_rect()
                if self._left_axis == "x":
                    shift = (position.x() - self._last.x()) / rect.width() * (x1 - x0)
                    self.owner.set_limits((x0 - shift, x1 - shift, y0, y1))
                else:
                    shift = (position.y() - self._last.y()) / rect.height() * (y1 - y0)
                    self.owner.set_limits((x0, x1, y0 + shift, y1 + shift))
            self._last = position
            return
        self._cursor_at(event.position())
        self._show_peak(event)

    def mouseReleaseEvent(self, event):
        if self._press is not None and event.button() == self._drag_button:
            self._press = None
            self._last = None
            self._right_axis = None
            self._left_axis = None
            self._drag_button = None
            self._axis_anchor = None
            self._cursor_at(event.position())
            event.accept()

    def mouseDoubleClickEvent(self, event):
        if self.plot_rect().contains(event.position()) or self._axis_at(event.position()):
            self.owner.reset_view()

    def wheelEvent(self, event):
        if not self.plot_rect().contains(event.position()):
            return
        mz, intensity = self._data(event.position())
        factor = 0.8 if event.angleDelta().y() > 0 else 1.25
        self.owner.zoom_axis("x", mz, factor)
        self.owner.zoom_axis("y", intensity, factor)
        event.accept()

    def _show_peak(self, event):
        if not self.plot_rect().contains(event.position()):
            self._hover_peak = None
            QtWidgets.QToolTip.hideText()
            return
        nearby = []
        mz, _ = self._data(event.position())
        tolerance = 6 / self.plot_rect().width() * (self.owner.limits[1] - self.owner.limits[0])
        for peak in self.owner.nearby_peaks(mz, tolerance):
            x, y = self._pixel(peak["mz"], peak["relative"])
            dx = abs(x - event.position().x())
            if dx <= 6 and y - 6 <= event.position().y() <= self.plot_rect().bottom() + 5:
                nearby.append((dx, 0 if peak["fragment"] else 1, abs(y - event.position().y()), peak))
        if not nearby:
            self._hover_peak = None
            QtWidgets.QToolTip.hideText()
            return
        peak = min(nearby, key=lambda item: item[:3])[3]
        if peak == self._hover_peak:
            return
        self._hover_peak = peak
        parts = [
            f"实测 m/z：{peak['mz']:.4f}",
            f"相对强度：{peak['relative']:.2f}%",
        ]
        if self.owner.ms1_spectrum is not None:
            parts.insert(1, f"强度：{peak['intensity']:,.0f}")
        fragment = peak["fragment"]
        if fragment:
            parts.append(f"匹配碎片：{html.escape(str(fragment.get('name') or '—'))}")
            parts.append(f"理论 m/z：{float(fragment['theoretical_mz']):.4f}")
            if fragment.get("error_ppm") is not None:
                parts.append(f"误差：{float(fragment['error_ppm']):+.2f} ppm")
            parts.append(f"类型：{FRAGMENT_LABELS[display_role(fragment)]}")
        elif self.owner.ms1_spectrum is None:
            parts.append("未匹配")
        QtWidgets.QToolTip.showText(event.globalPosition().toPoint(), "<br>".join(parts), self)


class SpectrumPlot(QtWidgets.QWidget):
    """One measured spectrum; no equal-height theoretical sticks."""

    y_axis_label = "Intensity (%)"

    def __init__(self):
        super().__init__()
        self.peaks = []
        self.ms1_spectrum = None
        self.target_mz = None
        self.target_peak_mz = None
        self.status = "没有保存可显示的实验谱图"
        self._survey_base = 1.0
        self._visible_cache = None
        self.limits = (0.0, 1.0, 0.0, 110.0)
        self.home_limits = self.limits
        self.canvas = _SpectrumCanvas(self)
        self.reset_button = reset_view_button(self, self.reset_view)
        controls = QtWidgets.QHBoxLayout()
        controls.setContentsMargins(0, 0, 0, 0)
        controls.addStretch(1)
        controls.addWidget(self.reset_button)
        self.controls = controls
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        layout.addLayout(controls)
        layout.addWidget(self.canvas, 1)

    @property
    def has_peaks(self):
        return bool(self.peaks) or (self.ms1_spectrum is not None and bool(len(self.ms1_spectrum.peaks)))

    def _survey_slice(self, lower, upper):
        peaks = self.ms1_spectrum.peaks
        start = int(np.searchsorted(peaks[:, 0], lower, side="left"))
        stop = int(np.searchsorted(peaks[:, 0], upper, side="right"))
        return peaks[start:stop]

    def _survey_point(self, peak):
        return dict(mz=float(peak[0]), intensity=float(peak[1]),
                    relative=float(peak[1]) / self._survey_base * 100.0, fragment=None)

    def visible_peaks(self, columns):
        if self.ms1_spectrum is None:
            return self.peaks
        key = self.limits[:2], columns
        if self._visible_cache is None or self._visible_cache[0] != key:
            x0, x1 = self.limits[:2]
            peaks = self._survey_slice(x0, x1)
            if len(peaks) > columns * 2:
                bins = np.clip(((peaks[:, 0] - x0) / (x1 - x0) * columns).astype(int), 0, columns - 1)
                edges = np.r_[0, np.flatnonzero(np.diff(bins)) + 1, len(peaks)]
                indices = [left + int(np.argmax(peaks[left:right, 1]))
                           for left, right in zip(edges[:-1], edges[1:])]
                # Retain the actual precursor even when a taller peak in its
                # pixel column wins display thinning in the full-scan view.
                if self.target_peak_mz is not None:
                    target = int(np.searchsorted(peaks[:, 0], self.target_peak_mz))
                    if target < len(peaks) and peaks[target, 0] == self.target_peak_mz:
                        indices = np.unique(np.r_[indices, target])
                peaks = peaks[indices]
            self._visible_cache = key, [self._survey_point(peak) for peak in peaks]
        return self._visible_cache[1]

    def nearby_peaks(self, mz, tolerance):
        if self.ms1_spectrum is None:
            return self.peaks
        return (self._survey_point(peak) for peak in self._survey_slice(mz - tolerance, mz + tolerance))

    def set_status(self, message):
        self.set_evidence(None)
        self.status = message
        self.canvas.update()

    def is_precursor_peak(self, peak):
        return (self.ms1_spectrum is not None and self.target_peak_mz is not None
                and peak["mz"] == self.target_peak_mz)

    def set_ms1_spectrum(self, spectrum, target_mz=None, *, tolerance_ppm=10.):
        self.set_evidence(None)
        self.ms1_spectrum = spectrum
        self.target_mz = float(target_mz) if target_mz is not None and math.isfinite(target_mz) and target_mz > 0 else None
        if spectrum is not None and len(spectrum.peaks):
            peaks = spectrum.peaks
            self._survey_base = float(peaks[:, 1].max()) or 1.0
            left, right = float(peaks[0, 0]), float(peaks[-1, 0])
            step = decimal_tick_step(left, right, 6) if right > left else 1.0
            self.home_limits = (max(0.0, math.floor(left / step) * step),
                                max(right + 1.0, math.ceil(right / step) * step), 0.0, 105.0)
            if self.target_mz is not None:
                tolerance = self.target_mz * tolerance_ppm * 1e-6
                nearby = self._survey_slice(self.target_mz - tolerance, self.target_mz + tolerance)
                nearby = nearby[nearby[:, 1] > 0]
                if len(nearby):
                    self.target_peak_mz = float(nearby[np.argmin(abs(nearby[:, 0] - self.target_mz)), 0])
        else:
            self.status = "所选时间没有可显示的 MS1 峰"
        self.reset_view()
        self.focus_precursor()

    def focus_precursor(self):
        if self.ms1_spectrum is None or self.target_mz is None:
            return
        lower, upper = self.target_mz - 3.0, self.target_mz + 3.0
        peaks = self._survey_slice(lower, upper)
        height = float(peaks[:, 1].max()) / self._survey_base * 100 if len(peaks) else 1.0
        self.set_limits((lower, upper, 0.0, max(height * 1.08, MIN_VIEW_SPAN)))

    def set_limits(self, limits):
        x0, x1, y0, y1 = (float(value) for value in limits)
        if not all(map(math.isfinite, (x0, x1, y0, y1))) or x1 <= x0 or y1 <= y0:
            return
        bx0, bx1, by0, by1 = self.home_limits
        x0, x1 = clamp_window(x0, x1, bx0, bx1)
        y0, y1 = clamp_window(y0, y1, by0, by1)
        self.limits = (x0, x1, y0, y1)
        self.canvas.update()

    def zoom_axis(self, axis, anchor, factor):
        x0, x1, y0, y1 = self.limits
        bx0, bx1, by0, by1 = self.home_limits
        if axis == "x":
            x0, x1 = scaled_window(x0, x1, anchor, factor, bx0, bx1)
        elif axis == "y":
            y0, y1 = scaled_from_minimum(y0, y1, factor, by0, by1)
        self.set_limits((x0, x1, y0, y1))

    def reset_view(self):
        self.limits = self.home_limits
        self.canvas.update()

    def set_evidence(self, evidence, experimental=None):
        self.ms1_spectrum = None
        self.target_mz = None
        self.target_peak_mz = None
        self._visible_cache = None
        self.status = "没有保存可显示的实验谱图"
        fragments = evidence.get("fragments", []) if evidence else []
        spectrum = experimental if experimental is not None else (evidence or {}).get("spectrum", [])
        observed = {}
        priority = {"headgroup": 0, "chain": 1, "neutral_loss": 2, "precursor": 3, "other": 4}
        for fragment in fragments:
            mz = fragment.get("observed_mz")
            if mz is not None:
                key = float(mz)
                if key not in observed or priority.get(display_role(fragment), 5) < priority.get(display_role(observed[key]), 5):
                    observed[key] = fragment
        valid = []
        for mz, intensity in spectrum or []:
            try:
                mz, intensity = float(mz), float(intensity)
            except (TypeError, ValueError):
                continue
            if math.isfinite(mz) and math.isfinite(intensity) and intensity >= 0:
                valid.append((mz, intensity))
        base = max((intensity for _, intensity in valid), default=0) or 1
        self.peaks = [
            dict(mz=mz, relative=intensity / base * 100.0, fragment=observed.get(mz))
            for mz, intensity in valid
        ]
        if valid:
            left = min(mz for mz, _ in valid)
            right = max(mz for mz, _ in valid)
            step = decimal_tick_step(left, right, 6) if right > left else 1.0
            bound_left = max(0.0, math.floor(left / step) * step)
            bound_right = math.ceil(right / step) * step
            if bound_right - bound_left < MIN_VIEW_SPAN:
                bound_left = max(0.0, left - 1.0)
                bound_right = right + 1.0
            self.home_limits = (bound_left, bound_right, 0.0, 105.0)
        else:
            self.home_limits = (0.0, 1.0, 0.0, 110.0)
        self.canvas._hover_peak = None
        self.reset_view()
