"""Small native toolbar icons shared by result plots."""

from pathlib import Path
import sys

from PySide6 import QtCore, QtGui, QtWidgets


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
