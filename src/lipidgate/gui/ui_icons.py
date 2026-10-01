"""Small native toolbar icons shared by result plots."""

from pathlib import Path
import sys

from PySide6 import QtCore, QtGui, QtWidgets


_PLOT_CURSORS = {}


def plot_cursor(kind: str) -> QtGui.QCursor:
    """Thin 18 px cursors, independent of oversized native resize/hand cursors."""
    if kind in _PLOT_CURSORS:
        return _PLOT_CURSORS[kind]
    pixmap = QtGui.QPixmap(36, 36)
    pixmap.setDevicePixelRatio(2)
    pixmap.fill(QtCore.Qt.GlobalColor.transparent)
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
    pen = QtGui.QPen(QtGui.QColor("#434a51"), 1.0)
    pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    if kind in {"x", "y"}:
        if kind == "y":
            painter.translate(9, 9)
            painter.rotate(90)
            painter.translate(-9, -9)
        painter.drawLine(QtCore.QPointF(2.5, 9), QtCore.QPointF(15.5, 9))
        for tip, direction in ((2.5, 1), (15.5, -1)):
            painter.drawLine(QtCore.QPointF(tip, 9), QtCore.QPointF(tip + direction * 3, 6))
            painter.drawLine(QtCore.QPointF(tip, 9), QtCore.QPointF(tip + direction * 3, 12))
    elif kind in {"grab", "grabbing"}:
        path = QtGui.QPainterPath()
        if kind == "grab":
            path.moveTo(5.5, 9)
            path.lineTo(5.5, 4)
            path.cubicTo(5.5, 2.4, 7.2, 2.4, 7.2, 4)
            path.lineTo(7.2, 7.5)
            path.lineTo(7.2, 2.8)
            path.cubicTo(7.2, 1.2, 9, 1.2, 9, 2.8)
            path.lineTo(9, 7.2)
            path.lineTo(9, 3.4)
            path.cubicTo(9, 1.9, 10.8, 1.9, 10.8, 3.4)
            path.lineTo(10.8, 7.8)
            path.lineTo(10.8, 5)
            path.cubicTo(10.8, 3.5, 12.6, 3.5, 12.6, 5)
            path.lineTo(12.6, 10)
        else:
            path.moveTo(5.5, 9)
            path.lineTo(5.5, 6.5)
            path.cubicTo(5.5, 5, 7.3, 5, 7.3, 6.5)
            path.cubicTo(7.3, 4.8, 9, 4.8, 9, 6.5)
            path.cubicTo(9, 5.1, 10.8, 5.1, 10.8, 6.8)
            path.cubicTo(10.8, 5.8, 12.6, 5.8, 12.6, 7.3)
            path.lineTo(12.6, 10)
        path.cubicTo(12.6, 12.6, 11, 13.4, 11, 15)
        path.lineTo(6.6, 15)
        path.cubicTo(6.6, 13.6, 3.9, 12, 3, 10)
        path.cubicTo(2.1, 8.1, 3.2, 7.4, 4.1, 8.3)
        path.lineTo(5.5, 9.8)
        path.closeSubpath()
        painter.setBrush(QtGui.QColor("#ffffff"))
        painter.drawPath(path)
    elif kind == "cross":
        painter.drawLine(3, 9, 15, 9)
        painter.drawLine(9, 3, 9, 15)
    else:
        painter.end()
        raise ValueError(f"Unknown plot cursor: {kind}")
    painter.end()
    cursor = QtGui.QCursor(pixmap, 9, 9)
    _PLOT_CURSORS[kind] = cursor
    return cursor


def icon_path(name: str) -> Path:
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[3]))
    return root / "assets" / "icons" / name


def plot_tool_icon(kind: str) -> QtGui.QIcon:
    """Small vector-like icons painted at 2x resolution for the plot controls."""
    pixmap = QtGui.QPixmap(48, 48)
    pixmap.fill(QtCore.Qt.GlobalColor.transparent)
    pixmap.setDevicePixelRatio(2)
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
    pen = QtGui.QPen(QtGui.QColor("#526477"), 1.8)
    pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
    painter.setPen(pen)
    if kind in {"back", "forward"}:
        painter.drawLine(5, 12, 19, 12)
        direction = -1 if kind == "back" else 1
        center = 5 if direction < 0 else 19
        painter.drawLine(center, 12, center - 6 * direction, 6)
        painter.drawLine(center, 12, center - 6 * direction, 18)
    elif kind == "select":
        pen.setStyle(QtCore.Qt.PenStyle.DashLine)
        painter.setPen(pen)
        painter.drawRoundedRect(QtCore.QRectF(4.5, 5.5, 15, 13), 1.4, 1.4)
    elif kind == "save":
        painter.drawRoundedRect(QtCore.QRectF(5, 4, 14, 16), 1.5, 1.5)
        painter.drawLine(8, 4.5, 8, 10)
        painter.drawLine(8, 10, 16, 10)
        painter.drawLine(16, 10, 16, 4.5)
        painter.drawLine(8, 16, 16, 16)
    else:
        painter.end()
        raise ValueError(f"Unknown plot icon: {kind}")
    painter.end()
    return QtGui.QIcon(pixmap)


def reset_view_button(parent, callback):
    icon = QtGui.QIcon(str(icon_path("view_reset.svg")))
    icon.addFile(str(icon_path("view_reset_active.svg")), mode=QtGui.QIcon.Mode.Active)
    button = QtWidgets.QToolButton(parent)
    button.setIcon(icon)
    button.setIconSize(QtCore.QSize(22, 22))
    button.setFixedSize(32, 32)
    button.setAutoRaise(True)
    button.setToolTip("恢复完整视图")
    button.clicked.connect(callback)
    return button
