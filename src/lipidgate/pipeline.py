"""One-mode project execution with immutable settings and isolated run outputs."""

from datetime import datetime
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from lipidgate.ms1 import run_feature_detection_result
from lipidgate.ms2 import run_ms2_feature_annotation_result
from lipidgate.final_results import export_final_results
from lipidgate.project import Project


def check_polarity(path, mode):
    if mode not in {"positive", "negative"}:
        raise ValueError("无效的离子模式")
    expected = "MS:1000130" if mode == "positive" else "MS:1000129"
    has_ms1 = False
    in_spectrum = False
    # Parse the whole stream, releasing binary arrays as soon as consumed.
    for event, element in ET.iterparse(path, events=("start", "end")):
        name = element.tag.rsplit("}", 1)[-1]
        if event == "start":
            if name == "spectrum":
                in_spectrum = True
            continue
        if name == "cvParam" and element.get("accession") in {
            "MS:1000130",
            "MS:1000129",
        }:
            if element.get("accession") != expected:
                raise ValueError(
                    f"{Path(path).name} 与项目离子模式不一致；一个项目仅处理一个模式"
                )
        if in_spectrum and name == "cvParam" and element.get("accession") == "MS:1000511":
            has_ms1 |= element.get("value") == "1"
        if name == "spectrum":
            in_spectrum = False
        element.clear()
    return has_ms1


def run_project(project_path, settings, progress=lambda message: None):
    import pandas as pd
    from lipidgate.ms2.provenance import sha256, code_fingerprint
    from lipidgate.ms1.detection import default_config
    from lipidgate.paths import default_positive_msp, default_negative_msp

    project = Project.open(project_path)
    files = project.mzml_files()
    settings = json.loads(json.dumps(settings))
    mode = settings["ms2"]["mode"]
    if mode not in {"positive", "negative"}:
        raise ValueError("请选择正或负模式")
    options = dict(settings["ms2"])
    library = Path(
        options.get("library_path")
        or (default_positive_msp() if mode == "positive" else default_negative_msp())
    ).resolve()
    if not library.is_file():
        raise FileNotFoundError(f"谱库不存在：{library}")
    if library.stat().st_size < 1024 and library.read_bytes().startswith(
        b"version https://git-lfs.github.com/spec/v1"
    ):
        raise ValueError("谱库尚未下载，请在项目目录运行 git lfs pull")
    from lipidgate.ms2.file_search import validate_parallel_capacity

    validate_parallel_capacity(options.get("workers", 1), len(files), library)
    ms1 = settings["ms1"]
    if ms1.get("enabled", True):
        params = ms1.get("params", {})
        widths = params.get(
            "peakwidth", [params.get("min_fwhm", 5), params.get("max_fwhm", 60)]
        )
        if len(widths) != 2 or not 0 < widths[0] <= widths[1]:
            raise ValueError("峰宽范围无效")
        if (
            ms1["algo"] in {"msdial", "ms-dial"}
            and not Path(ms1.get("msdial_table") or "").is_file()
        ):
            raise ValueError("请选择有效的 MS-DIAL 特征表")
    progress("检查输入文件和离子模式…")
    ms1_files = []
    for index, file in enumerate(files, 1):
        if check_polarity(file, mode):
            ms1_files.append(file)
        progress(f"检查 已完成 {index}/{len(files)}：{file.name}")
    out = project.root / "runs" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    out.mkdir(parents=True)
    (out / "run_settings.json").write_text(
        json.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (out / "provenance.json").write_text(
        json.dumps(
            {
                "product_code_sha256": code_fingerprint(Path(__file__).parent),
                "ms1_runtime_sha256": code_fingerprint(
                    Path(__file__).parent.parent / "lipidbench"
                ),
                "inputs": [
                    {
                        "path": str(p),
                        "bytes": p.stat().st_size,
                        "mtime_ns": p.stat().st_mtime_ns,
                    }
                    for p in files
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (out / "library.json").write_text(
        json.dumps({"path": str(library), "sha256": sha256(library)}, indent=2),
        encoding="utf-8",
    )
    project.settings = settings
    project.save()
    feature = None
    external_ms1 = ms1.get("algo") in {"msdial", "ms-dial"}
    if ms1.get("enabled", True) and (ms1_files or external_ms1):
        progress("提取并对齐 MS1 特征…")
        if not external_ms1 and len(ms1_files) < len(files):
            progress(f"MS1 跳过 {len(files) - len(ms1_files)} 个不含 MS1 扫描的文件；这些文件的结果仍可进行 MS2 匹配。")
        config = default_config()
        config["parameters"]["asari"]["mode"] = "pos" if mode == "positive" else "neg"
        config["parameters"]["xcms"]["polarity"] = mode
        feature = run_feature_detection_result(
            algo=ms1["algo"],
            input_path=files if external_ms1 else ms1_files,
            output_dir=out / "ms1",
            params=ms1["params"],
            msdial_table=ms1.get("msdial_table"),
            config=config,
            progress=progress,
        )
    elif ms1.get("enabled", True):
        progress("所有文件均无 MS1 扫描，跳过特征提取并继续 MS2 匹配。")
    progress("匹配 MS2 谱库…")
    options["library_path"] = library
    feature_table = (
        (getattr(feature, "native_table_path", None) or feature.table_path)
        if feature
        else None
    )
    result = run_ms2_feature_annotation_result(
        mzml_input=files,
        feature_table=feature_table,
        output_dir=out / "ms2",
        export_csv=True,
        export_xlsx=False,
        map_to_features=feature is not None,
        min_total_score=0.0,
        progress=progress,
        **options,
    )
    try:
        audit = pd.read_csv(out / "ms2/audit/ms2_candidates.csv")
    except pd.errors.EmptyDataError:
        audit = pd.DataFrame()
    audit["mode"] = mode
    if feature is None:
        audit["ms1_support_status"] = "MS2-only"
        audit["ms1_support_reason"] = (
            "no_ms1_scans" if ms1.get("enabled", True) else "feature_detection_disabled"
        )
    progress("应用分数与 ECN 筛选，导出结果…")
    filter_settings = {key: value for key, value in settings["filter"].items()
                       if key in {"use_ecn", "use_score", "min_score", "rt_tolerance", "retain_unmodeled"}}
    final = export_final_results(audit, out / "results", export_xlsx=False,
                                 export_plots=False, **filter_settings)
    progress("完成")
    return feature, result, final
