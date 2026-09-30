"""Device-pixel-aware measured MS/MS spectrum for the result browser."""

from __future__ import annotations

import html
import math

from PySide6 import QtCore, QtGui, QtWidgets

from .plot_axes import MIN_VIEW_SPAN, clamp_window, decimal_tick_step, scaled_window, tick_label, tick_values
from .result_plots import FRAGMENT_COLORS, FRAGMENT_LABELS, display_role
from .ui_icons import reset_view_button


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
            painter.drawText(QtCore.QRectF(30, y - 10, rect.left() - 39, 20), QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, tick_label(intensity, y_step))
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
        painter.drawText(QtCore.QRectF(-rect.height() / 2, -14, rect.height(), 28), QtCore.Qt.AlignmentFlag.AlignCenter, "相对强度 (%)")
        painter.restore()

        if not self.owner.peaks:
            painter.setPen(QtGui.QColor("#75838e"))
            painter.setFont(self._font(10.5))
            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter, "没有保存可显示的实验谱图")
            return

        painter.save()
        painter.setClipRect(rect)
        baseline = self._pixel(x0, 0)[1]
        # Draw unassigned peaks first so a weak matched peak remains visible.
        for matched in (False, True):
            for peak in self.owner.peaks:
                if bool(peak["fragment"]) != matched:
                    continue
                x, y = self._pixel(peak["mz"], peak["relative"])
                if x < rect.left() or x > rect.right():
                    continue
                role = display_role(peak["fragment"]) if matched else "unmatched"
                color = QtGui.QColor(FRAGMENT_COLORS[role])
                pen = QtGui.QPen(color, 2.1 if matched else 1.2)
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
        axis = self._axis_at(position)
        cursor = (
            QtCore.Qt.CursorShape.SizeHorCursor if axis == "x" else
            QtCore.Qt.CursorShape.SizeVerCursor if axis == "y" else
            QtCore.Qt.CursorShape.ArrowCursor
        )
        self.setCursor(cursor)

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
        event.accept()

    def mouseMoveEvent(self, event):
        if self._press is not None:
            position = event.position()
            if self._right_axis is not None:
                coordinate = "x" if self._right_axis == "x" else "y"
                current = position.x() if coordinate == "x" else position.y()
                previous = self._last.x() if coordinate == "x" else self._last.y()
                delta = current - previous
                if delta:
                    factor = math.exp(max(-3.0, min(3.0, delta * 0.012)))
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
        for peak in self.owner.peaks:
            x, y = self._pixel(peak["mz"], peak["relative"])
            dx = abs(x - event.position().x())
            if dx <= 6 and y - 6 <= event.position().y() <= self.plot_rect().bottom() + 5:
                nearby.append((dx, 0 if peak["fragment"] else 1, abs(y - event.position().y()), peak))
        if not nearby:
            self._hover_peak = None
            QtWidgets.QToolTip.hideText()
            return
        peak = min(nearby, key=lambda item: item[:3])[3]
        if peak is self._hover_peak:
            return
        self._hover_peak = peak
        parts = [
            f"实测 m/z：{peak['mz']:.4f}",
            f"相对强度：{peak['relative']:.2f}%",
        ]
        fragment = peak["fragment"]
        if fragment:
            parts.append(f"匹配碎片：{html.escape(str(fragment.get('name') or '—'))}")
            parts.append(f"理论 m/z：{float(fragment['theoretical_mz']):.4f}")
            if fragment.get("error_ppm") is not None:
                parts.append(f"误差：{float(fragment['error_ppm']):+.2f} ppm")
            parts.append(f"类型：{FRAGMENT_LABELS[display_role(fragment)]}")
        else:
            parts.append("未匹配")
        QtWidgets.QToolTip.showText(event.globalPosition().toPoint(), "<br>".join(parts), self)


class SpectrumPlot(QtWidgets.QWidget):
    """One measured spectrum; no equal-height theoretical sticks."""

    def __init__(self):
        super().__init__()
        self.peaks = []
        self.limits = (0.0, 1.0, 0.0, 110.0)
        self.home_limits = self.limits
        self.canvas = _SpectrumCanvas(self)
        self.canvas.setToolTip("坐标轴上右键拖动缩放、左键拖动平移；图内滚轮缩放")
        self.reset_button = reset_view_button(self, self.reset_view)
        controls = QtWidgets.QHBoxLayout()
        controls.setContentsMargins(0, 0, 0, 0)
        controls.addStretch(1)
        controls.addWidget(self.reset_button)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        layout.addLayout(controls)
        layout.addWidget(self.canvas, 1)

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
            y0, y1 = scaled_window(y0, y1, anchor, factor, by0, by1)
        self.set_limits((x0, x1, y0, y1))

    def reset_view(self):
        self.limits = self.home_limits
        self.canvas.update()

    def set_evidence(self, evidence, experimental=None):
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
