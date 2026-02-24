"""
将 GLORI CSV（1-based 单碱基位点）转为 BED（0-based）。
- BED name 列 = 输入 Gene 列
- BED score = 三列 m6A_level 的平均值
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="CSV(1-based) 转 BED(0-based)")
    parser.add_argument(
        "input",
        nargs="?",
        default="GLORI_v2.0_50ng.csv",
        help="输入 CSV 路径（默认: GLORI_v2.0_50ng.csv）",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="输出 BED 路径（默认: 输入文件名.bed）",
    )
    args = parser.parse_args()

    inp = Path(args.input)
    out = Path(args.output) if args.output else inp.with_suffix(".bed")

    df = pd.read_csv(inp, dtype={"Chr": str, "Gene": str})

    # 三列 m6A_level 列名（可能带空格和 %)
    m6a_cols = [c for c in df.columns if "m6A_level" in c and "(%" in c]
    if len(m6a_cols) < 3:
        m6a_cols = [c for c in df.columns if "m6A_level" in c]
    if not m6a_cols:
        raise SystemExit("未找到 m6A_level 列")

    # 1-based 单碱基 -> 0-based BED: chromStart = Site-1, chromEnd = Site
    bed = pd.DataFrame({
        "chrom": df["Chr"],
        "chromStart": df["Site"].astype(int) - 1,
        "chromEnd": df["Site"].astype(int),
        "name": df["Gene"].fillna("").astype(str),
        "score": df[m6a_cols].astype(float).mean(axis=1).round(2),
        "strand": df["Strand"].fillna(".").astype(str),
    })

    bed.to_csv(out, sep="\t", index=False, header=False)
    print(f"已写入 {out}，共 {len(bed)} 行")


if __name__ == "__main__":
    main()
