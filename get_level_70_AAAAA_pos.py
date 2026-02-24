"""
在 get_AAAAA_pos.py 基础上：仅输出三列 m6A_level_rep 平均值 > 70 的位点，
并输出 5mer（AAAAA/TTTTT）起止位置的标准 BED 文件。
"""
import argparse
import csv
import sys
from typing import Optional, Tuple
import os

import pysam

# 仅输出 level 平均值大于此阈值的位点
MIN_MEAN_LEVEL = 70
LEVEL_KEYS = ["m6A_level_rep1 (%)", "m6A_level_rep2 (%)", "m6A_level_rep3 (%)"]


def normalize_chrom_name(chrom: str, fasta: pysam.FastaFile) -> Optional[str]:
    """
    直接使用 CSV 里的 Chr 名字与 fasta 里的 contig 名匹配。
    假设两个文件中的染色体命名完全一致（例如都为 'chr1', 'chr2', ...）。
    """
    chrom = chrom.strip()
    if chrom in fasta.references:
        return chrom
    return None


def infer_coordinate_system(
    fasta: pysam.FastaFile,
    chrom: str,
    site: int,
    strand: str,
) -> Optional[str]:
    """
    对单个记录尝试判断 Site 是 1-based 还是 0-based。
    返回: '1-based' / '0-based' / None (无法判断)
    """
    expected_base = "A" if strand == "+" else "T"

    base_1based = ""
    if site - 1 >= 0:
        base_1based = fasta.fetch(chrom, site - 1, site).upper()

    base_0based = fasta.fetch(chrom, site, site + 1).upper()

    matches_1 = (base_1based == expected_base)
    matches_0 = (base_0based == expected_base)

    if matches_1 and not matches_0:
        return "1-based"
    if matches_0 and not matches_1:
        return "0-based"
    return None


def get_window_and_center_base(
    fasta: pysam.FastaFile,
    chrom: str,
    site: int,
    coord_system: str,
    flank: int = 4,
) -> Tuple[str, str, int]:
    """
    给定 site 和坐标体系，返回：
    - 9bp 窗口序列（左右各 flank, 共 2*flank + 1 bp）
    - 中心碱基（site 对应的碱基）
    - 窗口起始位置（0-based）
    """
    if coord_system == "1-based":
        center0 = site - 1
    elif coord_system == "0-based":
        center0 = site
    else:
        raise ValueError(f"未知坐标体系: {coord_system}")

    if center0 < 0:
        raise ValueError(f"中心位点小于 0，site={site}, coord_system={coord_system}")

    start = max(0, center0 - flank)
    end = center0 + flank + 1

    window_seq = fasta.fetch(chrom, start, end).upper()
    center_base = fasta.fetch(chrom, center0, center0 + 1).upper()
    return window_seq, center_base, start


def mean_m6a_level(row: dict) -> Optional[float]:
    """计算三列 m6A_level_rep 的平均值；若缺列或非数字则返回 None。"""
    vals = []
    for k in LEVEL_KEYS:
        v = row.get(k, "")
        if v == "":
            return None
        try:
            vals.append(float(v))
        except (ValueError, TypeError):
            return None
    if len(vals) != 3:
        return None
    return sum(vals) / 3.0


def main():
    parser = argparse.ArgumentParser(
        description="过滤 AAAAA/TTTTT 模体位点，仅输出 m6A level 平均值 > 70 的位点，并输出 5mer 起止位置的 BED"
    )
    parser.add_argument("--input", "-i", required=True, help="输入 CSV 文件，例如 GLORI_v2.0_50ng.csv")
    parser.add_argument("--ref", "-r", required=True, help="参考基因组 fasta，例如 hg38.fa")
    parser.add_argument("--output", "-o", required=True, help="输出 CSV 文件名")
    parser.add_argument(
        "--min-level",
        type=float,
        default=MIN_MEAN_LEVEL,
        help=f"三列 m6A_level_rep 平均值的下限，仅输出大于此值的位点（默认 {MIN_MEAN_LEVEL}）",
    )
    parser.add_argument(
        "--max_infer_rows",
        type=int,
        default=100,
        help="用于推断坐标体系的最大行数（默认 100）",
    )
    args = parser.parse_args()

    try:
        fasta = pysam.FastaFile(args.ref)
    except Exception as e:
        print(f"无法打开参考基因组 {args.ref}: {e}", file=sys.stderr)
        sys.exit(1)

    # ---------- 第一次遍历：推断坐标体系 ----------
    coord_system: Optional[str] = None
    count_1based = 0
    count_0based = 0
    infer_tried = 0

    with open(args.input, newline="", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)
        row_idx = 0
        for row in reader:
            if infer_tried >= args.max_infer_rows:
                break
            row_idx += 1
            try:
                chrom_raw = row["Chr"]
                site = int(row["Site"])
                strand = row["Strand"].strip()
            except Exception:
                continue

            chrom = normalize_chrom_name(chrom_raw, fasta)
            if chrom is None or strand not in {"+", "-"}:
                continue

            tmp = infer_coordinate_system(fasta, chrom, site, strand)
            if tmp == "1-based":
                count_1based += 1
                infer_tried += 1
            elif tmp == "0-based":
                count_0based += 1
                infer_tried += 1

    if count_1based > 0 and count_0based == 0:
        coord_system = "1-based"
        print(f"在前 {infer_tried} 条可用记录中，统一推断为 1-based 坐标体系。", file=sys.stderr)
    elif count_0based > 0 and count_1based == 0:
        coord_system = "0-based"
        print(f"在前 {infer_tried} 条可用记录中，统一推断为 0-based 坐标体系。", file=sys.stderr)
    else:
        coord_system = "1-based"
        print(
            f"警告: 在前 {infer_tried} 条可用记录中，无法得到一致的坐标体系推断结果，默认按 1-based 处理，请人工确认。",
            file=sys.stderr,
        )

    # ---------- 第二次遍历：筛选 AAAAA/TTTTT 且 level 平均值 > 阈值，写 CSV 与 BED ----------
    bed_path = os.path.splitext(args.output)[0] + ".bed"
    with open(args.input, newline="", encoding="utf-8") as f_in, \
            open(args.output, "w", newline="", encoding="utf-8") as f_out, \
            open(bed_path, "w", encoding="utf-8") as bed_out:
        reader = csv.DictReader(f_in)
        fieldnames = ["Chr", "Site", "Strand",
                      "m6A_level_rep1 (%)", "m6A_level_rep2 (%)", "m6A_level_rep3 (%)"]
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        row_idx = 0
        for row in reader:
            row_idx += 1

            try:
                chrom_raw = row["Chr"]
                site = int(row["Site"])
                strand = row["Strand"].strip()
            except Exception as e:
                print(f"第 {row_idx} 行解析 Chr/Site/Strand 出错: {e}", file=sys.stderr)
                continue

            mean_level = mean_m6a_level(row)
            if mean_level is None or mean_level <= args.min_level:
                continue

            chrom = normalize_chrom_name(chrom_raw, fasta)
            if chrom is None:
                continue

            try:
                window_seq, center_base, window_start = get_window_and_center_base(
                    fasta, chrom, site, coord_system, flank=4
                )
            except Exception as e:
                print(f"第 {row_idx} 行获取序列窗口失败: {e}", file=sys.stderr)
                continue

            if strand == "+":
                if center_base != "A":
                    continue
                motif = "AAAAA"
                if motif not in window_seq:
                    continue
                idx = window_seq.find(motif)
                motif_start = window_start + idx
                motif_end = motif_start + len(motif)

                out_row = {
                    "Chr": chrom_raw,
                    "Site": site,
                    "Strand": strand,
                    "m6A_level_rep1 (%)": row.get("m6A_level_rep1 (%)", ""),
                    "m6A_level_rep2 (%)": row.get("m6A_level_rep2 (%)", ""),
                    "m6A_level_rep3 (%)": row.get("m6A_level_rep3 (%)", ""),
                }
                writer.writerow(out_row)
                # BED6: chrom, start, end, name, score, strand（5mer 起止位置，0-based）
                bed_out.write(f"{chrom}\t{motif_start}\t{motif_end}\t{site}\t0\t{strand}\n")

            elif strand == "-":
                if center_base != "T":
                    continue
                motif = "TTTTT"
                if motif not in window_seq:
                    continue
                idx = window_seq.find(motif)
                motif_start = window_start + idx
                motif_end = motif_start + len(motif)

                out_row = {
                    "Chr": chrom_raw,
                    "Site": site,
                    "Strand": strand,
                    "m6A_level_rep1 (%)": row.get("m6A_level_rep1 (%)", ""),
                    "m6A_level_rep2 (%)": row.get("m6A_level_rep2 (%)", ""),
                    "m6A_level_rep3 (%)": row.get("m6A_level_rep3 (%)", ""),
                }
                writer.writerow(out_row)
                bed_out.write(f"{chrom}\t{motif_start}\t{motif_end}\t{site}\t0\t{strand}\n")

    fasta.close()
    print(f"处理完成。仅输出 level 平均值 > {args.min_level} 的位点；CSV 与 BED 已写入。", file=sys.stderr)


if __name__ == "__main__":
    main()
