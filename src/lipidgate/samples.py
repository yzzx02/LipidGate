"""Sample roles and experimental groups are independent project metadata."""
from __future__ import annotations

SAMPLE_TYPES = ("Sample", "QC", "Blank")


def normalize_sample_metadata(value=None, *, order=1):
    if value is not None and not isinstance(value, dict):
        raise ValueError("样本信息必须为对象。")
    value = value or {}
    sample_type = value.get("sample_type", "Sample")
    aliases = {"sample": "Sample", "qc": "QC", "blank": "Blank", "bk": "Blank"}
    sample_type = aliases.get(str(sample_type).strip().lower())
    if sample_type is None:
        raise ValueError("样本类型须为 Sample、QC 或 Blank。")
    if value.get("group", "1") is None:
        raise ValueError("实验分组不能为空；默认分组为 1。")
    group = str(value.get("group", "1")).strip()
    if not group:
        raise ValueError("实验分组不能为空；默认分组为 1。")
    result = {"sample_type": sample_type, "group": group}
    for key, default in (("batch", 1), ("injection_order", order)):
        number = value.get(key, default)
        if isinstance(number, bool):
            raise ValueError("批次和进样顺序须为正整数。")
        try:
            parsed = int(number)
            if float(number) != parsed or parsed < 1:
                raise ValueError
        except (TypeError, ValueError, OverflowError):
            raise ValueError("批次和进样顺序须为正整数。") from None
        result[key] = parsed
    included = value.get("included", True)
    if not isinstance(included, bool):
        raise ValueError("是否参与分析须为布尔值。")
    result["included"] = included
    return result
