"""High DPI, responsive RT–m/z navigator drawn directly with Qt."""

from __future__ import annotations

import html
import math

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .plot_axes import MIN_VIEW_SPAN, clamp_window, decimal_tick_step, scaled_window, tick_label, tick_values
from .result_data import FAMILY_COLORS, family_for, number, text
from .ui_icons import plot_tool_icon, reset_view_button


def _readable_step(raw: float) -> float:
    power = 10.0 ** math.floor(math.log10(max(raw, MIN_VIEW_SPAN)))
    return next(multiplier * power for multiplier in (1, 2, 5, 10) if raw <= multiplier * power)


class _Canvas(QtWidgets.QWidget):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.setMouseTracking(True)
        self.setMinimumSize(280, 145)
        self._press = None
        self._last = None
        self._rubber = None
        self._drag_button = None
        self._right_axis = None
        self._left_axis = None
        self._dragged = False
        self._pixel_points = np.empty((0, 2))
        self._pixel_rows = []
        self._hover_key = None

    def plot_rect(self):
        width = max(1, self.width() - 74)
        columns = max(1, int(width // 56))
        legend_rows = math.ceil(len(self.owner.legend) / columns)
        top = 11 + 18 * legend_rows
        return QtCore.QRectF(62, top, width, max(1, self.height() - top - 46))

    def _to_pixel(self, x, y):
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        return rect.left() + (x - x0) / (x1 - x0) * rect.width(), rect.bottom() - (y - y0) / (y1 - y0) * rect.height()

    def _to_data(self, point):
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        return x0 + (point.x() - rect.left()) / rect.width() * (x1 - x0), y0 + (rect.bottom() - point.y()) / rect.height() * (y1 - y0)

    def _nearest(self, point):
        if not len(self._pixel_points):
            return None
        delta = self._pixel_points - np.array([point.x(), point.y()])
        distances = np.einsum("ij,ij->i", delta, delta)
        nearest = int(np.argmin(distances))
        if distances[nearest] > 9 * 9:
            return None
        close = np.flatnonzero(distances <= distances[nearest] + 2.0)
        for position in close:
            if self._pixel_rows[int(position)]["_key"] == self.owner.key:
                return self._pixel_rows[int(position)]
        return self._pixel_rows[nearest]

    def paintEvent(self, event):
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        painter.setRenderHint(QtGui.QPainter.RenderHint.TextAntialiasing)
        painter.fillRect(self.rect(), QtGui.QColor("#ffffff"))
        rect = self.plot_rect()
        x0, x1, y0, y1 = self.owner.limits
        font = self.font()
        font.setPointSizeF(9.5)
        painter.setFont(font)
        tick_pen = QtGui.QPen(QtGui.QColor("#657487"), 1.2)
        tick_text = QtGui.QColor("#526477")
        x_step = _readable_step((x1 - x0) / max(2.0, rect.width() / 95.0))
        for value in tick_values(x0, x1, x_step):
            x = self._to_pixel(value, y0)[0]
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(x, rect.bottom()), QtCore.QPointF(x, rect.bottom() + 5))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(x - 34, rect.bottom() + 4, 68, 19), QtCore.Qt.AlignmentFlag.AlignCenter, tick_label(value, x_step))
        y_step = decimal_tick_step(y0, y1, max(4.0, rect.height() / 48.0))
        for value in tick_values(y0, y1, y_step):
            y = self._to_pixel(x0, value)[1]
            painter.setPen(tick_pen)
            painter.drawLine(QtCore.QPointF(rect.left() - 5, y), QtCore.QPointF(rect.left(), y))
            painter.setPen(tick_text)
            painter.drawText(QtCore.QRectF(0, y - 10, rect.left() - 10, 20), QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, tick_label(value, y_step))
        painter.setPen(QtGui.QPen(QtGui.QColor("#596b7d"), 1.55))
        painter.drawLine(rect.bottomLeft(), rect.bottomRight())
        painter.drawLine(rect.bottomLeft(), rect.topLeft())
        painter.setPen(QtGui.QColor("#374b61"))
        font.setPointSizeF(10.5)
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        painter.setFont(font)
        painter.drawText(QtCore.QRectF(rect.left(), rect.bottom() + 23, rect.width(), 21), QtCore.Qt.AlignmentFlag.AlignCenter, "RT (min)")
        painter.save()
        painter.translate(17, rect.center().y())
        painter.rotate(-90)
        painter.drawText(QtCore.QRectF(-rect.height()/2, -13, rect.height(), 26), QtCore.Qt.AlignmentFlag.AlignCenter, "m/z")
        painter.restore()

        font.setPointSizeF(9.2)
        font.setWeight(QtGui.QFont.Weight.Medium)
        painter.setFont(font)
        columns = max(1, int(rect.width() // 56))
        for index, (family, color) in enumerate(self.owner.legend):
            row, column = divmod(index, columns)
            items_in_row = min(columns, len(self.owner.legend) - row * columns)
            legend_x = rect.right() - items_in_row * 56 + column * 56
            legend_y = 7 + row * 19
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor(color))
            painter.drawEllipse(QtCore.QPointF(legend_x + 5, legend_y + 7), 3.5, 3.5)
            painter.setPen(QtGui.QColor("#536477"))
            painter.drawText(QtCore.QRectF(legend_x + 13, legend_y - 2, 43, 18), QtCore.Qt.AlignmentFlag.AlignVCenter, family)

        painter.save()
        painter.setClipRect(rect)
        if self.owner.dim:
            background = []
            for row in self.owner.background_rows:
                x, y = self._to_pixel(row["_rt"], row["_mz"])
                if rect.contains(QtCore.QPointF(x, y)):
                    background.append(QtCore.QPointF(x, y))
            pen = QtGui.QPen(QtGui.QColor(148, 163, 184, 50), 4.5)
            pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawPoints(QtGui.QPolygonF(background))
        grouped = {}
        pixel_points = []
        pixel_rows = []
        for row in self.owner.rows:
            x, y = self._to_pixel(row["_rt"], row["_mz"])
            if not rect.contains(QtCore.QPointF(x, y)):
                continue
            grouped.setdefault(family_for(row["compound_class"]), []).append(QtCore.QPointF(x, y))
            pixel_points.append((x, y))
            pixel_rows.append(row)
        # Paint large families first so rare classes remain visible in dense areas.
        for family, points in sorted(grouped.items(), key=lambda item: -len(item[1])):
            color = QtGui.QColor(FAMILY_COLORS[family])
            color.setAlphaF(0.70)
            pen = QtGui.QPen(color, 6.2)
            pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawPoints(QtGui.QPolygonF(points))
        self._pixel_points = np.asarray(pixel_points, dtype=float).reshape(-1, 2)
        self._pixel_rows = pixel_rows
        selected = next((row for row in self.owner.rows if row["_key"] == self.owner.key), None)
        if selected is not None:
            x, y = self._to_pixel(selected["_rt"], selected["_mz"])
            if rect.contains(QtCore.QPointF(x, y)):
                painter.setPen(QtCore.Qt.PenStyle.NoPen)
                painter.setBrush(QtGui.QColor("#ffffff"))
                painter.drawEllipse(QtCore.QPointF(x, y), 11, 11)
                painter.setPen(QtGui.QPen(QtGui.QColor("#182b3d"), 2.2))
                painter.setBrush(QtGui.QColor(FAMILY_COLORS[family_for(selected["compound_class"])]))
                painter.drawEllipse(QtCore.QPointF(x, y), 7.2, 7.2)
        if self._rubber is not None:
            painter.setPen(QtGui.QPen(QtGui.QColor("#2563eb"), 1.4, QtCore.Qt.PenStyle.DashLine))
            painter.setBrush(QtGui.QColor(37, 99, 235, 28))
            painter.drawRect(self._rubber.normalized())
        painter.restore()

    def wheelEvent(self, event):
        if not self.plot_rect().contains(event.position()):
            return
        x, y = self._to_data(event.position())
        factor = 0.8 if event.angleDelta().y() > 0 else 1.25
        self.owner.zoom_axis("x", x, factor)
        self.owner.zoom_axis("y", y, factor)
        event.accept()

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
            QtCore.Qt.CursorShape.CrossCursor if self.owner.mode == "box" and self.plot_rect().contains(position) else
            QtCore.Qt.CursorShape.OpenHandCursor if self.plot_rect().contains(position) else
            QtCore.Qt.CursorShape.ArrowCursor
        )
        self.setCursor(cursor)

    def mousePressEvent(self, event):
        position = event.position()
        if event.button() == QtCore.Qt.MouseButton.RightButton:
            self._right_axis = self._axis_at(position)
            if self._right_axis is None:
                return
            x, y = self._to_data(position)
            self._axis_anchor = x if self._right_axis == "x" else y
            self._press = position
            self._last = position
            self._drag_button = event.button()
            self._dragged = False
        elif event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._left_axis = self._axis_at(position)
            if self._left_axis is None and not self.plot_rect().contains(position):
                return
        else:
            return
        self._press = position
        self._last = position
        self._drag_button = event.button()
        self._dragged = False
        if self._left_axis is None:
            if self.owner.mode == "box":
                self._rubber = QtCore.QRectF(position, position)
            else:
                self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)
        event.accept()

    def mouseMoveEvent(self, event):
        position = event.position()
        if self._press is not None and self._drag_button == QtCore.Qt.MouseButton.RightButton:
            coordinate = "x" if self._right_axis == "x" else "y"
            current = position.x() if coordinate == "x" else position.y()
            previous = self._last.x() if coordinate == "x" else self._last.y()
            delta = current - previous
            if delta:
                factor = math.exp(max(-3.0, min(3.0, delta * 0.012)))
                self.owner.zoom_axis(self._right_axis, self._axis_anchor, factor, record=False)
                self._dragged = True
            self._last = position
            return
        if self._press is not None and self._left_axis is not None:
            if not self._dragged and (position - self._press).manhattanLength() < 4:
                return
            x0, x1, y0, y1 = self.owner.limits
            rect = self.plot_rect()
            if self._left_axis == "x":
                shift = (position.x() - self._last.x()) / rect.width() * (x1 - x0)
                self.owner._set_limits((x0 - shift, x1 - shift, y0, y1), record=False)
            else:
                shift = (position.y() - self._last.y()) / rect.height() * (y1 - y0)
                self.owner._set_limits((x0, x1, y0 + shift, y1 + shift), record=False)
            self._last = position
            self._dragged = True
            return
        if self._press is not None and self.owner.mode == "box":
            self._rubber = QtCore.QRectF(self._press, position).intersected(self.plot_rect())
            self.update()
            return
        if self._press is not None and self._drag_button == QtCore.Qt.MouseButton.LeftButton:
            if not self._dragged and (position - self._press).manhattanLength() < 4:
                return
            self._dragged = True
            x0, x1, y0, y1 = self.owner.limits
            rect = self.plot_rect()
            dx = (position.x() - self._last.x()) / rect.width() * (x1-x0)
            dy = (position.y() - self._last.y()) / rect.height() * (y1-y0)
            self.owner._set_limits((x0-dx, x1-dx, y0+dy, y1+dy), record=False)
            self._last = position
            self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)
            return
        self._cursor_at(position)
        if not self.plot_rect().contains(position):
            QtWidgets.QToolTip.hideText()
            self._hover_key = None
            return
        row = self._nearest(position)
        key = row["_key"] if row else None
        if key == self._hover_key:
            return
        self._hover_key = key
        if row:
            tooltip = f"<b>{html.escape(text(row['_feature']))}</b><br>{html.escape(text(row['matched_name']))}<br>m/z {number(row['_mz']):.4f} · RT {number(row['_rt']):.3f} min<br>Score {number(row['final_score']):.2f}"
            QtWidgets.QToolTip.showText(event.globalPosition().toPoint(), tooltip, self)
        else:
            QtWidgets.QToolTip.hideText()

    def mouseReleaseEvent(self, event):
        if self._press is None:
            return
        if self._drag_button == QtCore.Qt.MouseButton.RightButton:
            if self._dragged:
                self.owner._record_limits()
        elif self._left_axis is not None:
            if self._dragged:
                self.owner._record_limits()
        elif self.owner.mode == "box" and self._rubber is not None:
            area = self._rubber.normalized()
            if area.width() >= 5 and area.height() >= 5:
                x0, y1 = self._to_data(area.topLeft())
                x1, y0 = self._to_data(area.bottomRight())
                count = sum(x0 <= row["_rt"] <= x1 and y0 <= row["_mz"] <= y1 for row in self.owner.rows)
                self.owner.selection_count.setText(f"框选 {count:,} 个特征")
                self.owner.region_changed.emit((x0, x1, y0, y1))
        elif self._dragged:
            self.owner._record_limits()
        else:
            row = self._nearest(event.position())
            if row:
                self.owner.selected.emit(row["_key"])
        self._press = None
        self._last = None
        self._rubber = None
        self._drag_button = None
        self._right_axis = None
        self._left_axis = None
        self._dragged = False
        self._cursor_at(event.position())
        self.update()


class NavigationPlot(QtWidgets.QWidget):
    selected = QtCore.Signal(str)
    region_changed = QtCore.Signal(object)

    def __init__(self):
        super().__init__()
        self.rows = []
        self.background_rows = []
        self.legend = []
        self.dim = True
        self.key = None
        self.mode = "select"
        self.limits = (0.0, 1.0, 0.0, 1.0)
        self.home_limits = self.limits
        self._history = [self.limits]
        self._history_position = 0
        self.canvas = _Canvas(self)
        toolbar = QtWidgets.QHBoxLayout()
        toolbar.setContentsMargins(0, 0, 0, 2)
        toolbar.setSpacing(3)
        self.home_button = reset_view_button(self, self.reset_view)
        self.home_button.setIconSize(QtCore.QSize(24, 24))
        self.home_button.setFixedSize(36, 36)
        self.home_button.setStyleSheet(
            "QToolButton { border:0; border-radius:6px; background:transparent; padding:0; }"
            "QToolButton:hover { background:#e9f1fc; }"
        )
        toolbar.addWidget(self.home_button)
        self.back_button = self._button(toolbar, "back", "上一个视图", self._back)
        self.forward_button = self._button(toolbar, "forward", "下一个视图", self._forward)
        self.box_button = self._button(toolbar, "select", "框选结果区域", self._mode, checked=True)
        self.save_button = self._button(toolbar, "save", "保存导航图", self._save)
        toolbar.addStretch(1)
        self.selection_count = QtWidgets.QLabel("")
        self.selection_count.setStyleSheet("color: #64748b; font-size: 10pt;")
        toolbar.addWidget(self.selection_count)
        self._update_history_buttons()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(toolbar)
        layout.addWidget(self.canvas, 1)

    @staticmethod
    def _button(layout, icon_name, tooltip, callback, checked=False):
        button = QtWidgets.QToolButton()
        button.setIcon(plot_tool_icon(icon_name))
        button.setIconSize(QtCore.QSize(23, 23))
        button.setFixedSize(36, 36)
        button.setAutoRaise(True)
        button.setToolTip(tooltip)
        button.setCheckable(checked)
        button.setStyleSheet("QToolButton { border:0; border-radius:6px; background:transparent; padding:0; }"
                             "QToolButton:hover, QToolButton:checked { background:#e9f1fc; }"
                             "QToolButton:disabled { border:0; background:transparent; }")
        (button.toggled if checked else button.clicked).connect(callback)
        layout.addWidget(button)
        return button

    def _mode(self):
        self.mode = "box" if self.box_button.isChecked() else "select"
        self.canvas._cursor_at(self.canvas.mapFromGlobal(QtGui.QCursor.pos()))

    def _record_limits(self):
        if self.limits == self._history[self._history_position]:
            return
        self._history = self._history[:self._history_position + 1] + [self.limits]
        self._history_position += 1
        self._update_history_buttons()

    def _update_history_buttons(self):
        self.back_button.setEnabled(self._history_position > 0)
        self.forward_button.setEnabled(self._history_position + 1 < len(self._history))

    def _set_limits(self, limits, record=True):
        if limits[1] <= limits[0] or limits[3] <= limits[2]:
            return
        bx0, bx1, by0, by1 = self.home_limits
        x0, x1 = clamp_window(float(limits[0]), float(limits[1]), bx0, bx1)
        y0, y1 = clamp_window(float(limits[2]), float(limits[3]), by0, by1)
        self.limits = (x0, x1, y0, y1)
        if record:
            self._record_limits()
        self.canvas.update()

    def zoom_axis(self, axis, anchor, factor, record=True):
        x0, x1, y0, y1 = self.limits
        bx0, bx1, by0, by1 = self.home_limits
        if axis == "x":
            x0, x1 = scaled_window(x0, x1, anchor, factor, bx0, bx1)
        elif axis == "y":
            y0, y1 = scaled_window(y0, y1, anchor, factor, by0, by1)
        self._set_limits((x0, x1, y0, y1), record=record)

    def _back(self):
        if self._history_position:
            self._history_position -= 1
            self.limits = self._history[self._history_position]
            self.canvas.update()
            self._update_history_buttons()

    def _forward(self):
        if self._history_position + 1 < len(self._history):
            self._history_position += 1
            self.limits = self._history[self._history_position]
            self.canvas.update()
            self._update_history_buttons()

    def reset_view(self):
        self._set_limits(self.home_limits)
        self.selection_count.clear()
        self.region_changed.emit(None)

    def set_data(self, all_rows, visible, dim=True):
        valid = all_rows.dropna(subset=["_rt", "_mz"])
        selected = visible.dropna(subset=["_rt", "_mz"])
        was_empty = not self.rows
        self.dim = dim
        self.rows = selected.to_dict("records")
        self.background_rows = valid.loc[~valid._key.isin(selected._key)].to_dict("records") if dim else []
        families = {}
        for row in self.rows:
            family = family_for(row["compound_class"])
            families.setdefault(family, FAMILY_COLORS[family])
        self.legend = [(name, families[name]) for name in ("GP", "SP", "GL", "FA", "ST", "Other") if name in families]
        if not valid.empty:
            x0, x1 = float(valid._rt.min()), float(valid._rt.max())
            y0, y1 = float(valid._mz.min()), float(valid._mz.max())
            padx = max((x1-x0)*0.03, 0.1)
            pady = max((y1-y0)*0.05, 1.0)
            self.home_limits = (max(0.0, x0-padx), x1+padx, max(0.0, y0-pady), y1+pady)
            if was_empty:
                self.limits = self.home_limits
                self._history = [self.limits]
                self._history_position = 0
                self._update_history_buttons()
            else:
                self._set_limits(self.limits, record=False)
        self.canvas.update()

    def select_key(self, key):
        self.key = key
        self.canvas.update()

    def _save(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "保存导航图", "LipidGate_RT_mz.png", "PNG (*.png)")
        if path:
            self.canvas.grab().save(path)
