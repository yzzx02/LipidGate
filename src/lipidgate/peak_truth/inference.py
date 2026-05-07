from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from lipidgate.paths import default_peak_truth_model_dir

from .convnext_fusion import AttrScaler, FusionModel, build_eval_transform


DEFAULT_ATTR_COLUMNS = ["SNR", "CV", "GS", "TPAS", "H2B", "ZZ", "DZZ", "PCC", "SKEW", "DENT", "DM", "ENT", "JAG"]


def _load_attr_scaler(obj: dict | None, scaler_path: Path | None = None) -> AttrScaler | None:
    if obj is None and scaler_path is not None and scaler_path.exists():
        obj = json.loads(scaler_path.read_text(encoding="utf-8"))
    if not obj:
        return None
    return AttrScaler(
        fill={str(k): float(v) for k, v in obj["fill"].items()},
        mean={str(k): float(v) for k, v in obj["mean"].items()},
        std={str(k): float(v) for k, v in obj["std"].items()},
    )


def _resolve_image(row: pd.Series, image_root: Path, mzml_stem: str | None) -> Path:
    if "image" in row and pd.notna(row["image"]):
        p = Path(str(row["image"]))
        return p if p.is_absolute() else image_root / p
    feature_id = str(row.get("Feature_ID", "")).strip()
    if not feature_id:
        raise ValueError("Missing Feature_ID and image column")
    candidates = []
    if mzml_stem:
        candidates.append(image_root / mzml_stem / f"{feature_id}.png")
    candidates.append(image_root / f"{feature_id}.png")
    candidates.extend(sorted(image_root.rglob(f"{feature_id}.png")))
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"EIC image not found for Feature_ID={feature_id} under {image_root}")


def _attr_tensor(row: pd.Series, attr_columns: list[str], scaler: AttrScaler | None):
    import torch

    vals = []
    for c in attr_columns:
        v = pd.to_numeric(row.get(c), errors="coerce")
        if pd.isna(v):
            v = float(scaler.fill[c]) if scaler is not None and c in scaler.fill else 0.0
        value = float(v)
        if scaler is not None and c in scaler.mean and c in scaler.std:
            value = (value - float(scaler.mean[c])) / float(scaler.std[c])
        vals.append(value)
    return torch.tensor(vals, dtype=torch.float32)


def predict_peak_truth(
    *,
    attributes_csv: str | Path,
    image_root: str | Path,
    model_dir: str | Path | None = None,
    output_csv: str | Path | None = None,
    mzml_stem: str | None = None,
    threshold: float = 0.5,
    batch_size: int = 16,
    cpu: bool = False,
) -> pd.DataFrame:
    import torch
    from PIL import Image

    model_path = Path(model_dir).resolve() if model_dir else default_peak_truth_model_dir()
    checkpoint_path = model_path / "best_model.pth"
    if not checkpoint_path.exists():
        raise FileNotFoundError(checkpoint_path)

    attrs_df = pd.read_csv(attributes_csv).copy()
    device = torch.device("cuda" if torch.cuda.is_available() and not cpu else "cpu")
    try:
        ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    except TypeError:
        ckpt = torch.load(checkpoint_path, map_location=device)
    ckpt_args = ckpt.get("args", {})
    attr_columns = list(ckpt.get("attr_columns") or DEFAULT_ATTR_COLUMNS)
    scaler = _load_attr_scaler(ckpt.get("attr_scaler"), model_path / "attr_scaler.json")

    model_mode = str(ckpt_args.get("model_mode", "gated_fusion"))
    input_size = int(ckpt_args.get("input_size", 224))
    transform = build_eval_transform(input_size=input_size)
    model = FusionModel(
        attr_dim=len(attr_columns),
        out_dim=1,
        dropout=float(ckpt_args.get("dropout", 0.2)),
        pretrained=False,
        vision_backbone=str(ckpt_args.get("vision_backbone", "convnext_tiny")),
        lwga_depth=int(ckpt_args.get("lwga_depth", 2)),
        lwga_groups=int(ckpt_args.get("lwga_groups", 8)),
        lwga_mlp_ratio=float(ckpt_args.get("lwga_mlp_ratio", 2.0)),
        lwga_dropout=float(ckpt_args.get("lwga_dropout", 0.0)),
        model_mode=model_mode,
    ).to(device)
    model.load_state_dict(ckpt["model_state"])
    model.eval()

    image_root = Path(image_root).resolve()
    probabilities: list[float] = []
    image_paths: list[str] = []

    with torch.no_grad():
        for start in range(0, len(attrs_df), int(batch_size)):
            batch = attrs_df.iloc[start : start + int(batch_size)]
            images = []
            attr_tensors = []
            for _, row in batch.iterrows():
                image_path = _resolve_image(row, image_root, mzml_stem)
                image_paths.append(str(image_path.relative_to(image_root) if image_path.is_relative_to(image_root) else image_path))
                image = Image.open(image_path).convert("RGB")
                images.append(transform(image))
                attr_tensors.append(_attr_tensor(row, attr_columns, scaler))
            image_tensor = torch.stack(images).to(device)
            attr_tensor = torch.stack(attr_tensors).to(device)
            logits = model(image_tensor, attr_tensor).squeeze(1)
            probabilities.extend(torch.sigmoid(logits).cpu().tolist())

    out = attrs_df.iloc[: len(probabilities)].copy()
    out["image"] = image_paths
    out["prob_true_peak"] = probabilities
    out["pred_true_peak"] = (out["prob_true_peak"] >= float(threshold)).astype(int)
    if output_csv is not None:
        out_path = Path(output_csv).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, index=False)
    return out
