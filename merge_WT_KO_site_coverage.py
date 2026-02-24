import argparse
import ast
from pathlib import Path

import numpy as np
import pandas as pd


KEY_COLS = ["Chromosome", "Start", "End", "Flag"]
NUM_COLS = [
    "count_interval_WT",
    "count_reads_WT",
    "count_interval_KO",
    "count_reads_KO",
    "Coverage_WT",
    "Coverage_KO",
]
LIST_COLS = ["Count_Read_Interval_WT", "Count_Read_Interval_KO"]
OUT_COL_ORDER = [
    "Chromosome",
    "Start",
    "End",
    "count_interval_WT",
    "count_reads_WT",
    "Flag",
    "Count_Read_Interval_WT",
    "count_interval_KO",
    "count_reads_KO",
    "Count_Read_Interval_KO",
    "Coverage_WT",
    "Coverage_KO",
    "m6A_level_WT",
    "m6A_level_KO",
]


def _parse_list_cell(x) -> list:
    if x is None:
        return []
    if isinstance(x, float) and np.isnan(x):
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, (int, np.integer)):
        return [] if x == 0 else [int(x)]
    if isinstance(x, (float, np.floating)):
        return [] if x == 0.0 else [float(x)]

    s = str(x).strip()
    if s in {"", "0", "0.0", "nan", "None"}:
        return []
    if s.startswith("[") and s.endswith("]"):
        try:
            v = ast.literal_eval(s)
        except Exception:
            return []
        if isinstance(v, list):
            return v
        return [] if v in {0, 0.0, None} else [v]
    try:
        v = ast.literal_eval(s)
        if isinstance(v, list):
            return v
        return [] if v in {0, 0.0, None} else [v]
    except Exception:
        return []


def _list_to_cell(v: list) -> str:
    if not v:
        return "0"
    return str(v)


def read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": "string", "Flag": "string"}, low_memory=False)

    missing = [c for c in (KEY_COLS + NUM_COLS + LIST_COLS) if c not in df.columns]
    if missing:
        raise ValueError(f"文件 {path} 缺少列: {missing}")

    df["Start"] = pd.to_numeric(df["Start"], errors="raise").astype("int64")
    df["End"] = pd.to_numeric(df["End"], errors="raise").astype("int64")

    for c in NUM_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    for c in LIST_COLS:
        df[c] = df[c].map(_parse_list_cell)

    return df


def merge_and_aggregate(df_a: pd.DataFrame, df_b: pd.DataFrame, eps: float = 1e-10) -> pd.DataFrame:
    df = pd.concat([df_a, df_b], ignore_index=True, copy=False)

    agg_spec = {c: "sum" for c in NUM_COLS}
    for c in LIST_COLS:
        agg_spec[c] = lambda s: [item for sub in s for item in sub]

    out = df.groupby(KEY_COLS, dropna=False, sort=False, as_index=False).agg(agg_spec)

    out["m6A_level_WT"] = out["count_reads_WT"] / (out["Coverage_WT"] + eps)
    out["m6A_level_KO"] = out["count_reads_KO"] / (out["Coverage_KO"] + eps)

    for c in LIST_COLS:
        out[c] = out[c].map(_list_to_cell)

    # 保持输出列顺序（若输入里还有额外列，这里不输出；按需求只输出这些列）
    out = out[OUT_COL_ORDER]
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="合并两个 site_coverage TSV：按 Chromosome/Start/End/Flag 聚合求和并拼接 list 列。"
    )
    parser.add_argument(
        "--a",
        default="sampled_merged_xgb_WT3_KO1_site_coverage.tsv",
        help="第一个输入 TSV（默认：sampled_merged_xgb_WT3_KO1_site_coverage.tsv）",
    )
    parser.add_argument(
        "--b",
        default="sampled_merged_xgb_WT1_KO1_site_coverage.tsv",
        help="第二个输入 TSV（默认：sampled_merged_xgb_WT1_KO1_site_coverage.tsv）",
    )
    parser.add_argument(
        "--out",
        default="merged_xgb_WT_KO_site_coverage.tsv",
        help="输出 TSV（默认：merged_xgb_WT_KO_site_coverage.tsv）",
    )
    args = parser.parse_args()

    path_a = Path(args.a)
    path_b = Path(args.b)
    out_path = Path(args.out)

    df_a = read_tsv(path_a)
    df_b = read_tsv(path_b)
    merged = merge_and_aggregate(df_a, df_b)
    merged.to_csv(out_path, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

