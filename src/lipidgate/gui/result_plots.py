"""Interactive canvases for recorded feature and fragment evidence."""

import html
import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector

from .result_data import class_color, family_for, number, text


FRAGMENT_COLORS = {
    "headgroup": "#2F6581",
    "chain": "#C05D4A",
    "neutral_loss": "#81709A",
    "precursor": "#3F8074",
    "other": "#65717D",
    "unmatched": "#B7C1C9",
}
FRAGMENT_LABELS = {
    "headgroup": "头基",
    "chain": "脂肪链",
    "neutral_loss": "中性丢失",
    "precursor": "母离子",
    "other": "其他碎片",
    "unmatched": "未匹配",
}


def display_role(fragment):
    role = fragment.get("role", "other")
    if role in {"headgroup", "chain"}:
        return role
    fragment_type = str(fragment.get("fragment_type", "")).strip().lower().replace("_", " ")
    if "loss" in fragment_type:
        return "neutral_loss"
    if fragment_type == "precursor ion":
        return "precursor"
    return role if role in FRAGMENT_COLORS else "other"


def style_axis(ax):
    ax.set_facecolor("white")
    ax.tick_params(labelsize=9.5, colors="#526171")
    for spine in ax.spines.values():
        spine.set_color("#dce2e9")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25, color="#dce2e9", linewidth=0.5)


class NavigationPlot(QtWidgets.QWidget):
    selected = QtCore.Signal(str)
    region_changed = QtCore.Signal(object)

    def __init__(self):
        super().__init__()
        self.figure = Figure(figsize=(6, 3), facecolor="white", layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)
        self.canvas.setMinimumSize(280, 145)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar.setIconSize(QtCore.QSize(16, 16))
        self.toolbar.setMinimumHeight(38)
        self.toolbar.setMaximumHeight(40)
        self.select_button = QtWidgets.QToolButton()
        self.select_button.setText("框选")
        self.select_button.setCheckable(True)
        self.reset_button = QtWidgets.QToolButton()
        self.reset_button.setText("恢复视图")
        controls = QtWidgets.QHBoxLayout()
        controls.setContentsMargins(0, 0, 0, 0)
        controls.addWidget(self.toolbar, 1)
        controls.addWidget(self.select_button)
        controls.addWidget(self.reset_button)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        layout.addLayout(controls)
        layout.addWidget(self.canvas, 1)
        self.scatter = None
        self.highlight = None
        self.rows = []
        self.key = None
        self.selector = RectangleSelector(
            self.ax,
            self._region,
            useblit=False,
            button=[1],
            minspanx=5,
            minspany=5,
            spancoords="pixels",
            interactive=False,
        )
        self.selector.set_active(False)
        self.select_button.toggled.connect(self._selection_mode)
        for name in ("pan", "zoom"):
            self.toolbar._actions[name].triggered.connect(
                lambda: self.select_button.setChecked(False)
            )
        self.reset_button.clicked.connect(self.reset_view)
        self.canvas.mpl_connect("pick_event", self._pick)
        self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.canvas.mpl_connect("scroll_event", self._zoom)

    def _selection_mode(self, on):
        if on and self.toolbar.mode:
            if "pan" in str(self.toolbar.mode).lower():
                self.toolbar.pan()
            else:
                self.toolbar.zoom()
        self.selector.set_active(on)

    def reset_view(self):
        self.region_changed.emit(None)
        if hasattr(self, "limits"):
            self.ax.set_xlim(self.limits[0])
            self.ax.set_ylim(self.limits[1])
            self.canvas.draw_idle()

    def _region(self, a, b):
        if None in (a.xdata, a.ydata, b.xdata, b.ydata):
            return
        self.region_changed.emit(
            (*sorted([a.xdata, b.xdata]), *sorted([a.ydata, b.ydata]))
        )

    def set_data(self, all_rows, visible, dim=True):
        self.selector.disconnect_events()
        self.selector.set_active(False)
        old_limits = (self.ax.get_xlim(), self.ax.get_ylim()) if self.rows else None
        self.ax.clear()
        style_axis(self.ax)
        self.ax.set_xlabel("RT (min)", fontsize=9)
        self.ax.set_ylabel("m/z", fontsize=9)
        self.highlight = None
        valid = all_rows.dropna(subset=["_rt", "_mz"])
        selected = visible.dropna(subset=["_rt", "_mz"])
        if dim:
            background = valid.loc[~valid._key.isin(selected._key)]
            self.ax.scatter(
                background._rt,
                background._mz,
                c="#939ba5",
                s=16,
                alpha=0.13,
                edgecolors="none",
                rasterized=True,
            )
        self.rows = selected.to_dict("records")
        self.scatter = self.ax.scatter(
            selected._rt,
            selected._mz,
            c=[class_color(v) for v in selected.compound_class],
            s=23,
            alpha=1,
            linewidths=0.2,
            edgecolors="white",
            picker=5,
            rasterized=True,
        )
        if valid.empty:
            self.ax.text(
                0.5,
                0.5,
                "No RT / m/z data",
                transform=self.ax.transAxes,
                ha="center",
                color="#7b8798",
                fontsize=10,
            )
        for family, label in [
            ("GP", "GP"),
            ("SP", "SP"),
            ("GL", "GL"),
            ("FA", "FA"),
            ("ST", "ST"),
            ("Other", "Other"),
        ]:
            examples = [
                r["compound_class"]
                for r in self.rows
                if family_for(r["compound_class"]) == family
            ]
            if examples:
                self.ax.scatter([], [], c=class_color(examples[0]), s=20, label=label)
        if self.rows:
            self.ax.legend(
                loc="upper left",
                ncol=6,
                frameon=False,
                fontsize=7,
                handletextpad=0.3,
                columnspacing=0.7,
            )
        self.limits = (self.ax.get_xlim(), self.ax.get_ylim())
        if old_limits:
            self.ax.set_xlim(old_limits[0])
            self.ax.set_ylim(old_limits[1])
        self.selector = RectangleSelector(
            self.ax,
            self._region,
            useblit=False,
            button=[1],
            minspanx=5,
            minspany=5,
            spancoords="pixels",
            interactive=False,
        )
        self.selector.set_active(self.select_button.isChecked())
        self.select_key(self.key)
        self.canvas.draw_idle()

    def select_key(self, key):
        self.key = key
        if self.highlight is not None:
            self.highlight.remove()
            self.highlight = None
        row = next((r for r in self.rows if r["_key"] == key), None)
        if row:
            self.highlight = self.ax.scatter(
                [row["_rt"]],
                [row["_mz"]],
                s=82,
                facecolors="none",
                edgecolors="#253447",
                linewidths=1.4,
                zorder=10,
            )
        self.canvas.draw_idle()

    def _pick(self, event):
        if (
            event.artist is self.scatter
            and len(event.ind)
            and not self.select_button.isChecked()
            and not self.toolbar.mode
        ):
            self.selected.emit(self.rows[int(event.ind[0])]["_key"])

    def _hover(self, event):
        if event.inaxes is not self.ax or self.scatter is None:
            return
        hit, info = self.scatter.contains(event)
        if not hit:
            QtWidgets.QToolTip.hideText()
            return
        row = self.rows[int(info["ind"][0])]
        content = f"<b>{html.escape(text(row['_feature']))}</b><br>{html.escape(text(row['matched_name']))}<br>m/z {row['_mz']:.4f}<br>RT {row['_rt']:.3f} min<br>Score {number(row['final_score']):.2f}<br>ECN {html.escape(row['_ecn'])}"
        QtWidgets.QToolTip.showText(QtGui.QCursor.pos(), content, self.canvas)

    def _zoom(self, event):
        if event.inaxes is not self.ax or event.xdata is None:
            return
        factor = 0.8 if event.button == "up" else 1.25
        for get, setter, center in [
            (self.ax.get_xlim, self.ax.set_xlim, event.xdata),
            (self.ax.get_ylim, self.ax.set_ylim, event.ydata),
        ]:
            low, high = get()
            setter(center + (low - center) * factor, center + (high - center) * factor)
        self.canvas.draw_idle()


class SpectrumPlot(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.figure = Figure(figsize=(4, 3.5), dpi=145, facecolor="white", layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.subplots(2, 1, sharex=True)
        self.canvas.setMinimumSize(260, 185)
        self.note = QtWidgets.QLabel("选择特征查看已保存的谱图证据")
        self.note.setWordWrap(True)
        self.note.setObjectName("resultMuted")
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(3)
        layout.addWidget(self.note)
        layout.addWidget(self.canvas, 1)
        self.hovers = {}
        self.canvas.mpl_connect("motion_notify_event", self._hover)
        self.set_evidence(None)

    def set_evidence(self, evidence, experimental=None):
        self.hovers = {}
        top, bottom = self.axes
        for ax in self.axes:
            ax.clear()
            style_axis(ax)
        top.set_title("Theoretical fragment positions", fontsize=11, loc="left", pad=5)
        top.set_ylim(0, 1.25)
        top.set_yticks([])
        bottom.set_title("Experimental MS/MS", fontsize=11, loc="left", pad=5)
        bottom.set_ylabel("Relative intensity (%)", fontsize=10)
        bottom.set_xlabel("m/z", fontsize=10)
        bottom.set_ylim(0, 110)
        if not evidence:
            top.text(
                0.5,
                0.45,
                "No theoretical snapshot",
                transform=top.transAxes,
                ha="center",
                fontsize=10,
                color="#7b8798",
            )
            bottom.text(
                0.5,
                0.45,
                "No measured spectrum saved",
                transform=bottom.transAxes,
                ha="center",
                fontsize=10,
                color="#7b8798",
            )
            self.note.setText(
                "此结果未保存完整谱图快照；仅在上方显示已有碎片文字，不推测理论位置或峰强度。"
            )
            self.canvas.draw_idle()
            return
        self.note.setText("上方等高棒线只表示理论位置；下方为实际参与匹配的实验峰。")
        fragments = evidence.get("fragments", [])
        spectrum = (
            experimental if experimental is not None else evidence.get("spectrum", [])
        )
        observed = {}
        importance = {"headgroup": 0, "chain": 1, "precursor": 2, "neutral_loss": 3, "other": 4}
        for fragment in sorted(fragments, key=lambda item: importance.get(display_role(item), 5), reverse=True):
            if fragment.get("observed_mz") is not None:
                observed[float(fragment["observed_mz"])] = fragment
        for role in FRAGMENT_COLORS:
            if role == "unmatched":
                continue
            for matched, opacity, width in ((False, 0.38, 1.2), (True, 1.0, 2.2)):
                group = [f for f in fragments if display_role(f) == role and (f.get("observed_mz") is not None) == matched]
                if group:
                    top.vlines([f["theoretical_mz"] for f in group], 0, 1, colors=FRAGMENT_COLORS[role], linewidths=width, alpha=opacity)
        self.hovers[top] = [(f["theoretical_mz"], 1, f) for f in fragments]
        valid = [
            (float(m), float(i))
            for m, i in spectrum
            if np.isfinite(m) and np.isfinite(i) and i >= 0
        ]
        scale = max((i for _, i in valid), default=0) or 1
        items = []
        for mz, intensity in valid:
            f = observed.get(mz, dict(role="unmatched", name="", observed_mz=mz))
            items.append((mz, intensity / scale * 100, f))
        for role, color in FRAGMENT_COLORS.items():
            group = [p for p in items if display_role(p[2]) == role]
            if group:
                bottom.vlines(
                    [p[0] for p in group],
                    0,
                    [p[1] for p in group],
                    colors=color,
                    linewidths=0.8 if role == "unmatched" else 1.8,
                )
        self.hovers[bottom] = items
        positions = [p[0] for p in valid] + [f["theoretical_mz"] for f in fragments]
        if positions:
            bottom.set_xlim(max(0, min(positions) - 10), max(positions) + 15)
        self.canvas.draw_idle()

    def _hover(self, event):
        entries = self.hovers.get(event.inaxes, [])
        if not entries or event.x is None:
            return
        # Hit-test both m/z and the actual stick height in display coordinates.
        points = event.inaxes.transData.transform(
            np.array([[m, h] for m, h, _ in entries])
        )
        baseline = event.inaxes.transData.transform((0, 0))[1]
        distances = np.abs(points[:, 0] - event.x)
        distances[(event.y < baseline - 5) | (event.y > points[:, 1] + 6)] = np.inf
        index = int(np.argmin(distances))
        if distances[index] > 7:
            QtWidgets.QToolTip.hideText()
            return
        _, _, f = entries[index]
        parts = []
        for label, key in [
            ("Experimental m/z", "observed_mz"),
            ("Theoretical m/z", "theoretical_mz"),
            ("Error (ppm)", "error_ppm"),
        ]:
            if f.get(key) is not None:
                parts.append(f"{label}: {float(f[key]):.4f}")
        parts += [
            f"Fragment: {html.escape(text(f.get('name'))) or '—'}",
            f"Type: {FRAGMENT_LABELS.get(display_role(f), '其他碎片')}",
        ]
        QtWidgets.QToolTip.showText(
            QtGui.QCursor.pos(), "<br>".join(parts), self.canvas
        )
