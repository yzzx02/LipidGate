"""Visible run status for the desktop workflow."""

from __future__ import annotations

from collections import deque
import re
import time

from PySide6 import QtCore, QtWidgets


_FILE_PROGRESS = re.compile(r"^(检查|MS1|MS2) 已完成 (\d+)/(\d+)：(.+)$")
_STAGES = ("检查输入", "MS1 特征", "MS2 匹配", "整理导出")


def _duration(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


class RunProgressPanel(QtWidgets.QFrame):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("runProgressPanel")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self._started: float | None = None
        self._file_times: deque[float] = deque(maxlen=6)
        self._library_ready_at: float | None = None
        self._stage = 0
        self._total = 0
        self._done = 0
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(1000)
        self._timer.timeout.connect(self._tick)

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 18)
        layout.setSpacing(12)
        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QLabel("分析进度")
        title.setObjectName("runProgressTitle")
        self.badge = QtWidgets.QLabel("准备中")
        self.badge.setObjectName("runProgressBadge")
        header.addWidget(title)
        header.addStretch()
        header.addWidget(self.badge)
        layout.addLayout(header)

        self.stage_labels = []
        stages = QtWidgets.QHBoxLayout()
        stages.setSpacing(8)
        for name in _STAGES:
            stage = QtWidgets.QLabel(name)
            stage.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            stage.setObjectName("runStage")
            stages.addWidget(stage, 1)
            self.stage_labels.append(stage)
        layout.addLayout(stages)

        self.bar = QtWidgets.QProgressBar()
        self.bar.setObjectName("runProgressBar")
        self.bar.setTextVisible(False)
        self.bar.setFixedHeight(12)
        layout.addWidget(self.bar)

        self.current = QtWidgets.QLabel("等待开始")
        self.current.setObjectName("runProgressCurrent")
        self.current.setWordWrap(True)
        layout.addWidget(self.current)

        timing = QtWidgets.QHBoxLayout()
        self.count = QtWidgets.QLabel("文件 0/0")
        self.elapsed = QtWidgets.QLabel("已用 00:00:00")
        self.remaining = QtWidgets.QLabel("预计剩余：计算中")
        for item in (self.count, self.elapsed, self.remaining):
            item.setObjectName("runProgressMeta")
            timing.addWidget(item)
        timing.addStretch()
        layout.addLayout(timing)
        self.setStyleSheet("""
            QFrame#runProgressPanel { background: #eff6ff; border: 1px solid #bfdbfe; border-radius: 12px; }
            QLabel#runProgressTitle { color: #172554; font-size: 17px; font-weight: 700; }
            QLabel#runProgressBadge { color: #1d4ed8; background: #dbeafe; border-radius: 10px; padding: 5px 10px; font-weight: 600; }
            QLabel#runStage { color: #64748b; background: #e2e8f0; border-radius: 7px; padding: 7px 3px; font-size: 12px; }
            QLabel#runStage[active="true"] { color: #ffffff; background: #2563eb; font-weight: 600; }
            QLabel#runStage[complete="true"] { color: #166534; background: #dcfce7; }
            QProgressBar#runProgressBar { background: #dbeafe; border: none; border-radius: 6px; }
            QProgressBar#runProgressBar::chunk { background: #2563eb; border-radius: 6px; }
            QLabel#runProgressCurrent { color: #1e3a8a; font-weight: 600; }
            QLabel#runProgressMeta { color: #475569; font-size: 12px; margin-right: 16px; }
        """)
        self._paint_stages()

    def start(self, total_files: int, ms1_enabled: bool = True) -> None:
        self._started = time.monotonic()
        self._file_times.clear()
        self._library_ready_at = None
        self._total = total_files
        self._done = 0
        self._stage = 0
        self._ms1_enabled = ms1_enabled
        self.count.setText(f"文件 0/{total_files}")
        self.elapsed.setText("已用 00:00:00")
        self.remaining.setText("预计剩余：计算中")
        self.current.setText("正在启动分析进程…")
        self.badge.setText("运行中")
        self.bar.setRange(0, 0)
        self._paint_stages()
        self.show()
        self._timer.start()

    def update_message(self, message: str) -> None:
        if self._started is None:
            return
        self.current.setText(message)
        if message.startswith("检查输入"):
            self._set_stage(0)
        elif message.startswith("提取并对齐 MS1"):
            self._set_stage(1)
        elif message.startswith("匹配 MS2"):
            self._set_stage(2)
        elif message.startswith("MS2 正在准备谱库缓存"):
            self.bar.setRange(0, 0)
            self.remaining.setText("预计剩余：谱库准备中")
        elif message.startswith("MS2 谱库已就绪"):
            self._library_ready_at = time.monotonic()
            self.remaining.setText("预计剩余：计算中")
        elif message.startswith("应用分数与 ECN"):
            self._set_stage(3)
            self.bar.setRange(0, 0)
            self.remaining.setText("预计剩余：正在整理结果")
        elif message.startswith(("MS1 正在", "MS2 正在")):
            self.bar.setRange(0, 0)
            self.remaining.setText("预计剩余：正在整理结果")
        match = _FILE_PROGRESS.match(message)
        if match:
            kind, done, total, _ = match.groups()
            stage = {"检查": 0, "MS1": 1, "MS2": 2}[kind]
            self._set_stage(stage)
            done, total = int(done), int(total)
            if stage != self._stage or total <= 0:
                return
            if done > self._done:
                self._file_times.append(time.monotonic())
            self._done = done
            self._total = total
            self.count.setText(f"文件 {done}/{total}")
            self.bar.setRange(0, total)
            self.bar.setValue(done)
            self._update_eta()
        self._tick()

    def finish(self) -> None:
        self._timer.stop()
        self._tick()
        self._started = None
        self._stage = len(_STAGES)
        self._paint_stages()
        self.bar.setRange(0, 1)
        self.bar.setValue(1)
        self.current.setText("分析完成，结果已保存")
        self.badge.setText("已完成")
        self.remaining.setText("预计剩余：00:00:00")

    def fail(self) -> None:
        self._timer.stop()
        self._tick()
        self._started = None
        self.badge.setText("运行失败")
        self.remaining.setText("预计剩余：—")

    def _set_stage(self, stage: int) -> None:
        if stage == self._stage:
            return
        self._stage = stage
        self._done = 0
        self._file_times.clear()
        self._library_ready_at = None
        self.count.setText(f"文件 0/{self._total}")
        self.remaining.setText("预计剩余：计算中")
        self.bar.setRange(0, 0)
        self._paint_stages()

    def _paint_stages(self) -> None:
        for index, label in enumerate(self.stage_labels):
            label.setProperty("active", index == self._stage)
            label.setProperty("complete", index < self._stage)
            label.style().unpolish(label)
            label.style().polish(label)

    def _update_eta(self) -> None:
        if self._done >= self._total:
            self.remaining.setText("预计剩余：等待下一阶段")
        elif len(self._file_times) >= 2:
            intervals = [b - a for a, b in zip(self._file_times, list(self._file_times)[1:])]
            per_file = sum(intervals) / len(intervals)
            self.remaining.setText(f"本阶段预计剩余：约 {_duration(per_file * (self._total - self._done))}")
        elif self._stage == 2 and self._file_times and self._library_ready_at is not None:
            per_file = self._file_times[0] - self._library_ready_at
            self.remaining.setText(f"本阶段预计剩余：约 {_duration(per_file * (self._total - self._done))}")
        else:
            self.remaining.setText("预计剩余：计算中")

    def _tick(self) -> None:
        if self._started is not None:
            self.elapsed.setText(f"已用 {_duration(time.monotonic() - self._started)}")
