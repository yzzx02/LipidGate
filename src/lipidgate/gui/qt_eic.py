"""Compact native MS1 EIC view for a selected result feature."""

from __future__ import annotations

import math

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .plot_axes import tick_label, tick_values


def _step(raw: float) -> float:
    power = 10.0 ** math.floor(math.log10(max(raw, 1e-9)))
    return next(multiplier * power for multiplier in (1, 2, 5, 10) if raw <= multiplier * power)


class EICPlot(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setMinimumSize(260, 185)
        self.setMouseTracking(True)
        self.times = np.empty(0)
        self.intensities = np.empty(0)
        self.feature_rt = None
        self.mz = None
        self.ppm = None
        self.status = "选择特征后查看 MS1 EIC"

    def set_status(self, message):
        self.status = message
        self.setToolTip("")
        self.times = np.empty(0)
        self.intensities = np.empty(0)
        self.update()

    def set_trace(self, times, intensities, *, feature_rt, mz, ppm):
        self.times = np.asarray(times, dtype=np.float64)
        self.intensities = np.asarray(intensities, dtype=np.float64)
        self.feature_rt = float(feature_rt)
        self.mz = float(mz)
        self.ppm = float(ppm)
        self.status = "" if len(self.times) else "所选 RT 范围没有 MS1 扫描"
        self.setToolTip("" if len(self.times) else "该 mzML 的所选时间窗口没有一级扫描；MS2 证据仍可在 MS2 标签查看。")
        self.update()

    def _plot_rect(self):
        return QtCore.QRectF(68, 26, max(1, self.width() - 84), max(1, self.height() - 72))

    def _bounds(self):
        x0 = float(self.times.min())
        x1 = float(self.times.max())
        if x1 <= x0:
            x0 -= 0.01
            x1 += 0.01
        ymax = max(float(self.intensities.max()), 1.0)
        return x0, x1, 0.0, ymax * 1.08

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
        if not len(self.times):
            painter.setPen(QtGui.QColor("#66788a"))
            painter.drawText(rect, QtCore.Qt.AlignmentFlag.AlignCenter, self.status)
            return
        bounds = self._bounds()
        x0, x1, _, y1 = bounds
        x_step = _step((x1 - x0) / max(3.0, rect.width() / 90.0))
        y_step = _step(y1 / max(3.0, rect.height() / 60.0))
        tick_pen = QtGui.QPen(QtGui.QColor("#657487"), 1.2)
        tick_text = QtGui.QColor("#4b5864")
        for rt in tick_values(x0, x1, x_step):
            x, _ = self._pixel(rt, 0, bounds)
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(x, rect.bottom()), QtCore.QPointF(x, rect.bottom() + 5))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(x - 33, rect.bottom() + 5, 66, 20), QtCore.Qt.AlignmentFlag.AlignCenter, tick_label(rt, x_step))
        for intensity in tick_values(0, y1, y_step):
            _, y = self._pixel(x0, intensity, bounds)
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(rect.left() - 5, y), QtCore.QPointF(rect.left(), y))
            label = f"{intensity / 1000:g}k" if intensity >= 1000 else tick_label(intensity, y_step)
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(2, y - 10, 57, 20), QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, label)
        painter.setPen(QtGui.QPen(QtGui.QColor("#596b7d"), 1.55))
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())
        painter.drawLine(rect.bottomLeft(), rect.topLeft())
        painter.setPen(QtGui.QColor("#243746"))
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.drawText(QtCore.QRectF(rect.left(), self.height() - 25, rect.width(), 20), QtCore.Qt.AlignmentFlag.AlignCenter, "RT (min)")
        painter.save()
        painter.translate(16, rect.center().y())
        painter.rotate(-90)
        painter.drawText(QtCore.QRectF(-rect.height() / 2, -14, rect.height(), 28), QtCore.Qt.AlignmentFlag.AlignCenter, "Intensity")
        painter.restore()

        painter.save()
        painter.setClipRect(rect)
        if x0 <= self.feature_rt <= x1:
            marker_x, _ = self._pixel(self.feature_rt, 0, bounds)
            painter.setPen(QtGui.QPen(QtGui.QColor("#b76555"), 1.2, QtCore.Qt.PenStyle.DashLine))
            painter.drawLine(QtCore.QPointF(marker_x, rect.top()), QtCore.QPointF(marker_x, rect.bottom()))
        path = QtGui.QPainterPath()
        for index, (rt, intensity) in enumerate(zip(self.times, self.intensities)):
            x, y = self._pixel(float(rt), float(intensity), bounds)
            if index == 0:
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
        painter.setPen(QtGui.QPen(QtGui.QColor("#2f6581"), 2.0))
        painter.drawPath(path)
        peak = int(np.argmax(self.intensities))
        px, py = self._pixel(float(self.times[peak]), float(self.intensities[peak]), bounds)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QColor("#2f6581"))
        painter.drawEllipse(QtCore.QPointF(px, py), 3.2, 3.2)
        painter.restore()

    def mouseMoveEvent(self, event):
        if not len(self.times) or not self._plot_rect().contains(event.position()):
            QtWidgets.QToolTip.hideText()
            return
        bounds = self._bounds()
        x0, x1 = bounds[:2]
        rect = self._plot_rect()
        rt = x0 + (event.position().x() - rect.left()) / rect.width() * (x1 - x0)
        index = int(np.argmin(np.abs(self.times - rt)))
        QtWidgets.QToolTip.showText(
            event.globalPosition().toPoint(),
            f"RT {self.times[index]:.4f} min<br>强度 {self.intensities[index]:,.0f}<br>m/z {self.mz:.4f} ± {self.ppm:g} ppm",
            self,
        )
