"""Run the unchanged analysis pipeline outside the Qt application's DLL space."""

import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from types import SimpleNamespace
from collections import deque

PREFIX = "LIPIDGATE_EVENT "


def run_isolated(project_path, settings, progress):
    env = os.environ.copy()
    # The GUI adds Qt DLL directories to PATH; those are not analysis dependencies.
    env["PATH"] = os.pathsep.join(
        p
        for p in env.get("PATH", "").split(os.pathsep)
        if not any(v in p.lower() for v in ("pyside6", "shiboken6"))
    )
    env["PYTHONIOENCODING"] = "utf-8"
    options = dict(
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )
    if os.name == "nt":
        options["creationflags"] = subprocess.CREATE_NO_WINDOW
    command = (
        [sys.executable, "--worker", str(project_path)]
        if getattr(sys, "frozen", False)
        else [sys.executable, "-m", "lipidgate.gui.project_worker", str(project_path)]
    )
    process = subprocess.Popen(command, **options)
    tail = deque(maxlen=30)
    result = None
    error = None
    try:
        process.stdin.write(json.dumps(settings, ensure_ascii=True))
        process.stdin.close()
        for line in process.stdout:
            tail.append(line.rstrip())
            if not line.startswith(PREFIX):
                continue
            event = json.loads(line[len(PREFIX) :])
            if event["type"] == "progress":
                progress(event["message"])
            elif event["type"] == "result":
                result = event
            elif event["type"] == "error":
                error = event["message"]
        code = process.wait()
        if error:
            raise RuntimeError(error)
        if code or result is None:
            raise RuntimeError("分析进程未正常完成：\n" + "\n".join(tail))
        return SimpleNamespace(
            row_count=result["row_count"],
            xlsx_path=Path(result["xlsx_path"]) if result["xlsx_path"] else None,
            csv_path=Path(result["csv_path"]),
        )
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait()
        process.stdout.close()


def main(project_path=None):
    from lipidgate.pipeline import run_project

    settings = json.load(sys.stdin)

    def emit(event):
        print(PREFIX + json.dumps(event, ensure_ascii=True), flush=True)

    try:
        _, _, result = run_project(
            project_path or sys.argv[1],
            settings,
            progress=lambda message: emit(dict(type="progress", message=message)),
        )
    except Exception as exc:
        emit(dict(type="error", message=str(exc)))
        traceback.print_exc()
        return 1
    emit(
        dict(
            type="result",
            row_count=len(result.passed_data),
            xlsx_path=str(result.xlsx_path) if result.xlsx_path else None,
            csv_path=str(result.csv_path),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
