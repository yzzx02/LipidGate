from __future__ import annotations

from pathlib import Path

import pandas as pd

from lipidgate.ecn_filter import ECNFilterResult
from lipidgate.ms2 import MS2SearchResult


def test_cli_ms2_search_prints_only_ms2_outputs(tmp_path: Path, monkeypatch, capsys) -> None:
    import lipidgate.cli as cli

    def fake_run_ms2_search_result(**kwargs) -> MS2SearchResult:
        csv_path = tmp_path / "ms2_results.csv"
        csv_path.write_text("scan_id\n", encoding="utf-8")
        return MS2SearchResult(
            data=pd.DataFrame(),
            csv_path=csv_path,
            xlsx_path=None,
            mode="negative",
            library_path=tmp_path / "tiny.msp",
            output_dir=tmp_path,
            row_count=0,
        )

    monkeypatch.setattr(cli, "run_ms2_search_result", fake_run_ms2_search_result)

    code = cli.main(["ms2-search", "--mzml", str(tmp_path / "sample.mzML"), "--output", str(tmp_path)])

    output = capsys.readouterr().out
    assert code == 0
    assert "csv:" in output
    assert "passed_csv:" not in output
    assert "model_summary_csv:" not in output


def test_cli_ecn_filter_prints_passed_and_model_outputs(tmp_path: Path, monkeypatch, capsys) -> None:
    import lipidgate.cli as cli

    def fake_run_ecn_filter_result(**kwargs) -> ECNFilterResult:
        csv_path = tmp_path / "ecn_filter_results.csv"
        passed_path = tmp_path / "ecn_passed_results.csv"
        summary_path = tmp_path / "ecn_model_summary.csv"
        for path in [csv_path, passed_path, summary_path]:
            path.write_text("x\n", encoding="utf-8")
        return ECNFilterResult(
            data=pd.DataFrame({"x": [1]}),
            csv_path=csv_path,
            xlsx_path=None,
            output_dir=tmp_path,
            row_count=1,
            passed_data=pd.DataFrame({"x": [1]}),
            model_summary=pd.DataFrame({"subclass": ["PC"]}),
            passed_csv_path=passed_path,
            model_summary_csv_path=summary_path,
        )

    monkeypatch.setattr(cli, "run_ecn_filter_result", fake_run_ecn_filter_result)

    code = cli.main(["ecn-filter", "--input", str(tmp_path / "in.csv"), "--output", str(tmp_path)])

    output = capsys.readouterr().out
    assert code == 0
    assert "passed_csv:" in output
    assert "model_summary_csv:" in output
    assert "passed_rows: 1" in output
