from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Iterator, Optional

import pandas as pd


def _read_acmop_json_fragment(path: Path) -> list[tuple[str, object]]:
    """Read ACMOP's archive JSON while preserving duplicate top-level keys."""
    text = path.read_text(encoding="utf-8-sig").strip()
    if not text:
        return []
    if text.startswith("{"):
        payload = text
    elif text.startswith(","):
        payload = "{" + text[1:] + "}"
    else:
        raise ValueError(f"Unsupported JSON archive format: {path}")

    data = json.loads(payload, object_pairs_hook=list)
    return data if isinstance(data, list) else list(data.items())


def _pairs_to_plain(value):
    if isinstance(value, list):
        if all(isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str) for item in value):
            return {key: _pairs_to_plain(item_value) for key, item_value in value}
        return [_pairs_to_plain(item) for item in value]
    return value


def _read_settings(settings_path: Path) -> tuple[Optional[str], Optional[str]]:
    if not settings_path.exists():
        return None, None
    parts = settings_path.read_text(encoding="utf-8", errors="ignore").split("|", 1)
    select_spec = parts[0].strip() if parts else None
    select_fea_config = parts[1].strip() if len(parts) > 1 else None
    return select_spec or None, select_fea_config or None


def find_acmop_archives(path: str | Path) -> list[Path]:
    """Find top-level ACMOP result archives, excluding jsonpickle object dumps."""
    root = Path(path)
    if root.is_file():
        return [root]

    archives: list[Path] = []
    if (root / "acmop-settings.txt").exists():
        select_spec, _ = _read_settings(root / "acmop-settings.txt")
        if select_spec:
            expected = root / f"{select_spec}.json"
            if expected.exists():
                return [expected]
        return sorted(p for p in root.glob("*.json") if p.parent.name != "jsonpickle")

    for settings_path in root.rglob("acmop-settings.txt"):
        folder = settings_path.parent
        select_spec, _ = _read_settings(settings_path)
        if select_spec and (folder / f"{select_spec}.json").exists():
            archives.append(folder / f"{select_spec}.json")
        else:
            archives.extend(p for p in folder.glob("*.json") if p.parent.name != "jsonpickle")

    archives.extend(p for p in root.glob("*.json") if p.parent.name != "jsonpickle")
    return sorted(set(archives))


def _flatten_geometric_parameters(items: list[dict]) -> dict:
    row = {}
    for item in items or []:
        for key, value in item.items():
            if isinstance(value, dict):
                row[f"gp.{key}.type"] = value.get("type")
                row[f"gp.{key}.name"] = value.get("name")
                row[f"gp.{key}.value"] = value.get("value")
                bounds = value.get("bounds") or [None, None]
                row[f"gp.{key}.min"] = bounds[0] if len(bounds) > 0 else None
                row[f"gp.{key}.max"] = bounds[1] if len(bounds) > 1 else None
            else:
                row[f"gp.{key}.value"] = value
    return row


def iter_archive_records(archive_path: str | Path) -> Iterator[dict]:
    archive_path = Path(archive_path)
    data = _read_acmop_json_fragment(archive_path)
    select_spec, select_fea_config = _read_settings(archive_path.parent / "acmop-settings.txt")

    archive_record_no = 0
    for run_key, run_payload in data:
        if not isinstance(run_payload, list):
            continue
        for design_key, design_payload_raw in run_payload:
            archive_record_no += 1
            design_payload = _pairs_to_plain(design_payload_raw)
            if not isinstance(design_payload, dict) or "Performance" not in design_payload:
                continue
            row_select_spec = select_spec or design_key.split("-gen", 1)[0]

            row = {
                "archive_record_no": archive_record_no,
                "source_file": str(archive_path),
                "source_folder": str(archive_path.parent),
                "archive_run_key": run_key,
                "design_key": design_key,
                "select_spec": row_select_spec,
                "select_fea_config": select_fea_config,
            }

            for key, value in (design_payload.get("Spec inputs") or {}).items():
                row[f"spec.{key}"] = value
            for key, value in (design_payload.get("x_denorm_dict") or {}).items():
                row[f"x.{key}"] = value
            row.update(_flatten_geometric_parameters(design_payload.get("Geometric parameters") or []))
            for key, value in (design_payload.get("Excitations") or {}).items():
                if key != "wily":
                    row[f"ex.{key}"] = value
            for key, value in (design_payload.get("Performance") or {}).items():
                if key in row:
                    row[f"perf.{key}"] = value
                else:
                    row[key] = value
            if row.get("select_fea_config") is None and row.get("select_fea_config_dict") is not None:
                row["select_fea_config"] = row["select_fea_config_dict"]

            yield row


def load_acmop_metrics(paths: str | Path | Iterable[str | Path]) -> pd.DataFrame:
    if isinstance(paths, (str, Path)):
        archive_paths = find_acmop_archives(paths)
    else:
        archive_paths = []
        for path in paths:
            archive_paths.extend(find_acmop_archives(path))

    records = []
    for archive_path in archive_paths:
        records.extend(iter_archive_records(archive_path))

    df = pd.DataFrame(records)
    if df.empty:
        return df

    _add_derived_metric_columns(df)
    return df


def _to_numeric(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series([pd.NA] * len(df), index=df.index, dtype="Float64")
    return pd.to_numeric(df[column], errors="coerce")


def _add_derived_metric_columns(df: pd.DataFrame) -> None:
    if "rated_efficiency" not in df.columns and "f2" in df.columns:
        df["rated_efficiency"] = -_to_numeric(df, "f2")
    if "rated_efficiency_pct" not in df.columns:
        df["rated_efficiency_pct"] = _to_numeric(df, "rated_efficiency") * 100
    if "normalized_torque_ripple_pct" not in df.columns:
        df["normalized_torque_ripple_pct"] = _to_numeric(df, "normalized_torque_ripple") * 100
    if "normalized_force_error_magnitude_pct" not in df.columns:
        df["normalized_force_error_magnitude_pct"] = _to_numeric(df, "normalized_force_error_magnitude") * 100
    if "TRV_kNm_per_m3" not in df.columns:
        df["TRV_kNm_per_m3"] = _to_numeric(df, "TRV") / 1000
    if "valid_metrics" not in df.columns:
        df["valid_metrics"] = (
            (_to_numeric(df, "Cost") > 0)
            & (_to_numeric(df, "rated_stack_length_mm") > 0)
            & (_to_numeric(df, "TRV") > 0)
            & (_to_numeric(df, "rated_efficiency") > 0)
        )


def summarize_metric_archive(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    rows = []
    for select_spec, group in df.groupby("select_spec", dropna=False):
        valid = group[group["valid_metrics"]] if "valid_metrics" in group else group
        if valid.empty:
            valid = group

        def best(metric: str, ascending: bool = True):
            series = _to_numeric(valid, metric)
            if series.dropna().empty:
                return None, None
            idx = series.idxmin() if ascending else series.idxmax()
            return valid.loc[idx, "design_key"], valid.loc[idx, metric]

        low_cost_key, low_cost = best("Cost", True)
        high_eta_key, high_eta = best("rated_efficiency", False)
        low_ripple_key, low_ripple = best("f3", True)

        rows.append(
            {
                "select_spec": select_spec,
                "records": len(group),
                "valid_records": len(valid),
                "low_cost_design": low_cost_key,
                "low_cost": low_cost,
                "high_efficiency_design": high_eta_key,
                "high_efficiency": high_eta,
                "low_ripple_design": low_ripple_key,
                "low_ripple_f3": low_ripple,
            }
        )
    return pd.DataFrame(rows)


def export_metrics_csv(df: pd.DataFrame, output_path: str | Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")
    return output_path


def _main() -> None:
    parser = argparse.ArgumentParser(description="Extract ACMOP post-process metrics from result JSON archives.")
    parser.add_argument("path", help="A result archive, result folder, or project root such as _default.")
    parser.add_argument("--csv", dest="csv_path", help="Optional CSV export path.")
    args = parser.parse_args()

    df = load_acmop_metrics(args.path)
    summary = summarize_metric_archive(df)
    print(summary.to_string(index=False) if not summary.empty else "No ACMOP metric records found.")
    if args.csv_path:
        export_metrics_csv(df, args.csv_path)
        print(f"Exported {len(df)} records to {args.csv_path}")


if __name__ == "__main__":
    _main()
