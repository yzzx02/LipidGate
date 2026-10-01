"""Compact native MS1 EIC view for linked and unlinked annotations."""

from __future__ import annotations

import math

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .plot_axes import axis_drag_factor, clamp_window, intensity_tick_label, scaled_from_minimum, scaled_window, tick_label, tick_values
from .eic_trace import EIC_HALF_WINDOW_MIN
from .ui_icons import plot_cursor


def _step(raw: float) -> float:
    power = 10.0 ** math.floor(math.log10(max(raw, 1e-9)))
    return next(multiplier * power for multiplier in (1, 2, 5, 10) if raw <= multiplier * power)


def display_indices(times, intensities, x0, x1, columns):
    """Keep minima/maxima per pixel column so narrow peaks survive display thinning."""
    if len(times) <= 2 * columns:
        return np.arange(len(times))
    bins = np.clip(((times - x0) / (x1 - x0) * columns).astype(int), 0, columns - 1)
    edges = np.r_[0, np.flatnonzero(np.diff(bins)) + 1, len(times)]
    indices = [0, len(times) - 1]
    for left, right in zip(edges[:-1], edges[1:]):
        values = intensities[left:right]
        indices.extend((left + int(np.argmin(values)), left + int(np.argmax(values))))
    return np.unique(indices)


class EICPlot(QtWidgets.QWidget):
    y_axis_label = "Intensity"

    def __init__(self):
        super().__init__()
        self.setMinimumSize(260, 185)
        self.setMouseTracking(True)
        self.times = np.empty(0)
        self.intensities = np.empty(0)
        self.feature_rt = None
        self.mz = None
        self.ppm = None
        self.half_window_min = EIC_HALF_WINDOW_MIN
        self._trace_bounds = None
        self.limits = None
        self._path_cache = None
        self.peak_bounds = None
        self.status = "选择特征后查看 EIC"
        self.signal_status = ""
        self._clear_drag()

    def _clear_drag(self):
        self._press = self._last = self._drag_axis = self._drag_button = self._axis_anchor = None

    def set_status(self, message):
        self.status = message
        self.setToolTip("")
        self.times = np.empty(0)
        self.intensities = np.empty(0)
        self._trace_bounds = None
        self.limits = None
        self._path_cache = None
        self.peak_bounds = None
        self.signal_status = ""
        self._clear_drag()
        self.unsetCursor()
        self.update()

    def set_eic(self, trace, *, peak_bounds=None):
        self.set_trace(trace.times, trace.intensities, feature_rt=trace.center_rt,
                       mz=trace.target_mz, ppm=trace.ppm, half_window_min=trace.half_window_min,
                       peak_bounds=peak_bounds)

    def set_trace(self, times, intensities, *, feature_rt, mz, ppm,
                  half_window_min=EIC_HALF_WINDOW_MIN, peak_bounds=None):
        self.times = np.asarray(times, dtype=np.float64)
        self.intensities = np.asarray(intensities, dtype=np.float64)
        self.feature_rt = float(feature_rt)
        self.mz = float(mz)
        self.ppm = float(ppm)
        self.half_window_min = half_window_min
        self.peak_bounds = peak_bounds
        self.status = "" if len(self.times) else "所选窗口没有 MS1 扫描"
        x0, x1 = max(0.0, self.feature_rt - half_window_min), self.feature_rt + half_window_min
        ymax = max(float(self.intensities.max()), 1.0) if len(self.intensities) else 1.0
        self._trace_bounds = x0, x1, 0.0, ymax * 1.08
        self.limits = self._trace_bounds
        self._path_cache = None
        self._clear_drag()
        self.signal_status = ("有 MS1 信号" if np.any(self.intensities > 0) else
                              "目标 m/z 无信号" if len(self.times) else "窗口内无 MS1 扫描")
        self.setToolTip(
            f"MS1 EIC · RT ±{half_window_min:g} min · {self.ppm:g} ppm · 原始强度\n"
            "轴上右键：右拖/上拖放大，左拖/下拖缩小；强度从 0 开始；双击复位。"
            if len(self.times) else "所选时间范围没有 MS1 扫描；可在 MS2 标签查看 m/z 谱图。"
        )
        self.update()

    def _plot_rect(self):
        return QtCore.QRectF(84, 48, max(1, self.width() - 100), max(1, self.height() - 94))

    def _bounds(self):
        return self.limits

    def set_limits(self, limits):
        if self._trace_bounds is None:
            return
        x0, x1, y0, y1 = limits
        bx0, bx1, by0, by1 = self._trace_bounds
        self.limits = (*clamp_window(x0, x1, bx0, bx1), *clamp_window(y0, y1, by0, by1))
        self.update()

    def zoom_axis(self, axis, anchor, factor):
        if self.limits is None:
            return
        x0, x1, y0, y1 = self.limits
        bx0, bx1, by0, by1 = self._trace_bounds
        if axis == "x":
            x0, x1 = scaled_window(x0, x1, anchor, factor, bx0, bx1)
        elif axis == "y":
            y0, y1 = scaled_from_minimum(y0, y1, factor, by0, by1)
        self.set_limits((x0, x1, y0, y1))

    def reset_view(self):
        self.limits = self._trace_bounds
        self.update()

    def _trace_path(self, rect, bounds):
        size = self.width(), self.height(), bounds
        if self._path_cache is None or self._path_cache[0] != size:
            start = max(0, int(np.searchsorted(self.times, bounds[0])) - 1)
            stop = min(len(self.times), int(np.searchsorted(self.times, bounds[1], side="right")) + 1)
            indices = start + display_indices(self.times[start:stop], self.intensities[start:stop],
                                              *bounds[:2], max(1, int(rect.width())))
            path = QtGui.QPainterPath()
            for position, index in enumerate(indices):
                x, y = self._pixel(float(self.times[index]), float(self.intensities[index]), bounds)
                if position == 0:
                    path.moveTo(x, y)
                else:
                    path.lineTo(x, y)
            fill = QtGui.QPainterPath(path)
            if len(indices):
                fill.lineTo(*self._pixel(float(self.times[indices[-1]]), 0.0, bounds))
                fill.lineTo(*self._pixel(float(self.times[indices[0]]), 0.0, bounds))
                fill.closeSubpath()
            self._path_cache = size, path, fill
        return self._path_cache[1:]

    def _pixel(self, rt, intensity, bounds):
        rect = self._plot_rect()
        x0, x1, y0, y1 = bounds
        return (
            rect.left() + (rt - x0) / (x1 - x0) * rect.width(),
            rect.bottom() - (intensity - y0) / (y1 - y0) * rect.height(),
        )

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing)
        painter.fillRect(self.rect(), QtGui.QColor("#ffffff"))
        rect = self._plot_rect()
        font = self.font()
        font.setPointSizeF(10.2)
        painter.setFont(font)
        if self._trace_bounds is None:
            painter.setPen(QtGui.QColor("#6e7379"))
            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter, self.status)
            return
        title = f"MS1 EIC · m/z {self.mz:.4f}"
        subtitle = f"RT {self.feature_rt:.3f} min · ±{self.half_window_min:g} min · {self.signal_status}"
        painter.setPen(QtGui.QColor("#42474d"))
        metrics = painter.fontMetrics()
        painter.drawText(QtCore.QRectF(8, 1, self.width() - 16, 22), QtCore.Qt.AlignmentFlag.AlignLeft,
                         metrics.elidedText(title, QtCore.Qt.TextElideMode.ElideRight, self.width() - 16))
        font.setPointSizeF(9.2)
        painter.setFont(font)
        painter.setPen(QtGui.QColor("#777c82"))
        painter.drawText(QtCore.QRectF(8, 23, self.width() - 16, 20), QtCore.Qt.AlignmentFlag.AlignLeft,
                         painter.fontMetrics().elidedText(subtitle, QtCore.Qt.TextElideMode.ElideRight, self.width() - 16))
        font.setPointSizeF(10.2)
        painter.setFont(font)
        bounds = self._bounds()
        x0, x1, y0, y1 = bounds
        x_step = _step((x1 - x0) / max(3.0, rect.width() / 90.0))
        y_step = _step((y1 - y0) / max(3.0, rect.height() / 60.0))
        tick_pen = QtGui.QPen(QtGui.QColor("#7a8086"), 1.0)
        tick_text = QtGui.QColor("#565c62")
        for rt in tick_values(x0, x1, x_step):
            x, _ = self._pixel(rt, 0, bounds)
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(x, rect.bottom()), QtCore.QPointF(x, rect.bottom() + 5))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(x - 33, rect.bottom() + 5, 66, 20), QtCore.Qt.AlignmentFlag.AlignCenter, tick_label(rt, x_step))
        for intensity in tick_values(y0, y1, y_step):
            _, y = self._pixel(x0, intensity, bounds)
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(rect.left() - 5, y), QtCore.QPointF(rect.left(), y))
            label = intensity_tick_label(intensity, y_step)
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(26, y - 10, rect.left() - 35, 20), QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, label)
        painter.setPen(QtGui.QPen(QtGui.QColor("#777d83"), 1.1))
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())
        painter.drawLine(rect.bottomLeft(), rect.topLeft())
        painter.setPen(QtGui.QColor("#42474d"))
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.drawText(QtCore.QRectF(rect.left(), self.height() - 25, rect.width(), 20), QtCore.Qt.AlignmentFlag.AlignCenter, "RT (min)")
        painter.save()
        painter.translate(16, rect.center().y())
        painter.rotate(-90)
        painter.drawText(QtCore.QRectF(-rect.height() / 2, -14, rect.height(), 28), QtCore.Qt.AlignmentFlag.AlignCenter, self.y_axis_label)
        painter.restore()

        painter.save()
        painter.setClipRect(rect)
        if x0 <= self.feature_rt <= x1:
            marker_x, _ = self._pixel(self.feature_rt, 0, bounds)
            painter.setPen(QtGui.QPen(QtGui.QColor("#a3a7ab"), 1.0, QtCore.Qt.PenStyle.DashLine))
            painter.drawLine(QtCore.QPointF(marker_x, rect.top()), QtCore.QPointF(marker_x, rect.bottom()))
        if not len(self.times):
            painter.setPen(QtGui.QColor("#6e7379"))
            painter.drawText(rect.adjusted(8, 0, -8, 0), QtCore.Qt.AlignmentFlag.AlignCenter | QtCore.Qt.TextFlag.TextWordWrap, self.status)
            painter.restore()
            return
        path, fill = self._trace_path(rect, bounds)
        if self.peak_bounds is not None:
            left, right = self.peak_bounds
            left = max(left, x0)
            right = min(right, x1)
            if right > left:
                lx, _ = self._pixel(left, 0, bounds)
                rx, _ = self._pixel(right, 0, bounds)
                painter.save()
                painter.setClipRect(QtCore.QRectF(lx, rect.top(), rx - lx, rect.height()), QtCore.Qt.ClipOperation.IntersectClip)
                painter.fillPath(fill, QtGui.QColor(189, 200, 207, 110))
                painter.restore()
                painter.setPen(QtGui.QPen(QtGui.QColor("#b2bcc3"), 0.9, QtCore.Qt.PenStyle.DotLine))
                for edge in self.peak_bounds:
                    if x0 <= edge <= x1:
                        ex, _ = self._pixel(edge, 0, bounds)
                        painter.drawLine(QtCore.QPointF(ex, rect.top()), QtCore.QPointF(ex, rect.bottom()))
        painter.setPen(QtGui.QPen(QtGui.QColor("#3e444a"), 1.6))
        painter.drawPath(path)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor("#3e444a"))
        for index in [int(np.argmax(self.intensities))]:
            px, py = self._pixel(float(self.times[index]), float(self.intensities[index]), bounds)
            painter.drawEllipse(QtCore.QPointF(px, py), 2.6, 2.6)
        painter.restore()

    def _data(self, position):
        rect = self._plot_rect()
        x0, x1, y0, y1 = self.limits
        return (x0 + (position.x() - rect.left()) / rect.width() * (x1 - x0),
                y0 + (rect.bottom() - position.y()) / rect.height() * (y1 - y0))

    def _axis_at(self, position):
        if self.limits is None:
            return None
        rect = self._plot_rect()
        if rect.left() <= position.x() <= rect.right() and rect.bottom() < position.y() <= self.height():
            return "x"
        if 0 <= position.x() < rect.left() and rect.top() <= position.y() <= rect.bottom():
            return "y"
        return None

    def _cursor_at(self, position):
        axis = self._drag_axis if self._press is not None else self._axis_at(position)
        self.setCursor(plot_cursor(axis) if axis else QtGui.QCursor(QtCore.Qt.CursorShape.ArrowCursor))

    def mousePressEvent(self, event):
        axis = self._axis_at(event.position())
        if axis is None or event.button() not in {QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.RightButton}:
            return
        self._drag_axis = axis
        self._drag_button = event.button()
        self._press = self._last = event.position()
        x, y = self._data(event.position())
        self._axis_anchor = x if axis == "x" else y
        self._cursor_at(event.position())
        event.accept()

    def mouseMoveEvent(self, event):
        position = event.position()
        self._cursor_at(position)
        if self._press is not None:
            QtWidgets.QToolTip.hideText()
            current = position.x() if self._drag_axis == "x" else position.y()
            previous = self._last.x() if self._drag_axis == "x" else self._last.y()
            if self._drag_button == QtCore.Qt.MouseButton.RightButton:
                self.zoom_axis(self._drag_axis, self._axis_anchor,
                               axis_drag_factor(previous, current, self._drag_axis))
            else:
                x0, x1, y0, y1 = self.limits
                rect = self._plot_rect()
                if self._drag_axis == "x":
                    shift = (current - previous) / rect.width() * (x1 - x0)
                    self.set_limits((x0 - shift, x1 - shift, y0, y1))
                else:
                    shift = (current - previous) / rect.height() * (y1 - y0)
                    self.set_limits((x0, x1, y0 + shift, y1 + shift))
            self._last = position
            return
        if not len(self.times) or not self._plot_rect().contains(event.position()):
            QtWidgets.QToolTip.hideText()
            return
        bounds = self._bounds()
        x0, x1 = bounds[:2]
        rect = self._plot_rect()
        rt = x0 + (event.position().x() - rect.left()) / rect.width() * (x1 - x0)
        position = int(np.searchsorted(self.times, rt))
        choices = [max(0, min(len(self.times) - 1, i)) for i in (position - 1, position)]
        index = min(choices, key=lambda i: abs(self.times[i] - rt))
        QtWidgets.QToolTip.showText(
            event.globalPosition().toPoint(),
            f"MS1 EIC<br>RT {self.times[index]:.4f} min<br>强度 {self.intensities[index]:,.0f}<br>m/z {self.mz:.4f} ± {self.ppm:g} ppm",
            self,
        )

    def mouseReleaseEvent(self, event):
        if self._press is not None and event.button() == self._drag_button:
            self._clear_drag()
            self._cursor_at(event.position())
            event.accept()

    def mouseDoubleClickEvent(self, event):
        if self.limits is not None and (self._plot_rect().contains(event.position()) or self._axis_at(event.position())):
            self.reset_view()

    def wheelEvent(self, event):
        if self.limits is not None and self._plot_rect().contains(event.position()):
            x, y = self._data(event.position())
            factor = 0.8 if event.angleDelta().y() > 0 else 1.25
            self.zoom_axis("x", x, factor)
            self.zoom_axis("y", y, factor)
            event.accept()
