import argparse
import csv
import sys
from typing import Optional, Tuple
import os

import pysam


def normalize_chrom_name(chrom: str, fasta: pysam.FastaFile) -> Optional[str]:
    """
    直接使用 CSV 里的 Chr 名字与 fasta 里的 contig 名匹配。
    假设两个文件中的染色体命名完全一致（例如都为 'chr1', 'chr2', ...）。
    """
    chrom = chrom.strip()
    if chrom in fasta.references:
        return chrom
    # 若完全一致仍找不到，返回 None，由调用处报错信息
    return None


def infer_coordinate_system(
    fasta: pysam.FastaFile,
    chrom: str,
    site: int,
    strand: str,
) -> Optional[str]:
    """
    对单个记录尝试判断 Site 是 1-based 还是 0-based。
    返回:
        '1-based' / '0-based' / None (无法判断)
    """
    expected_base = "A" if strand == "+" else "T"

    # 1-based: 中心碱基在 0-based 下的位置是 site-1
    base_1based = ""
    if site - 1 >= 0:
        base_1based = fasta.fetch(chrom, site - 1, site).upper()

    # 0-based: 中心碱基在 0-based 下的位置是 site
    base_0based = ""
    base_0based = fasta.fetch(chrom, site, site + 1).upper()

    matches_1 = (base_1based == expected_base)
    matches_0 = (base_0based == expected_base)

    if matches_1 and not matches_0:
        return "1-based"
    if matches_0 and not matches_1:
        return "0-based"

    # 两者都符合或者都不符合，无法仅靠这一条判断
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
    """
    if coord_system == "1-based":
        center0 = site - 1   # 转成 0-based
    elif coord_system == "0-based":
        center0 = site       # 已是 0-based
    else:
        raise ValueError(f"未知坐标体系: {coord_system}")

    if center0 < 0:
        raise ValueError(f"中心位点小于 0，site={site}, coord_system={coord_system}")

    start = max(0, center0 - flank)
    end = center0 + flank + 1  # pysam.fetch 右开区间

    window_seq = fasta.fetch(chrom, start, end).upper()
    center_base = fasta.fetch(chrom, center0, center0 + 1).upper()
    return window_seq, center_base, start


def main():
    parser = argparse.ArgumentParser(description="根据 hg38.fa 过滤 GLORI_v2.0_50ng.csv 中满足 AAAAA/TTTTT 模体的位点")
    parser.add_argument("--input", "-i", required=True, help="输入 CSV 文件，例如 GLORI_v2.0_50ng.csv")
    parser.add_argument("--ref", "-r", required=True, help="参考基因组 fasta，例如 hg38.fa")
    parser.add_argument("--output", "-o", required=True, help="输出 CSV 文件名")
    parser.add_argument(
        "--max_infer_rows",
        type=int,
        default=100,
        help="用于推断坐标体系的最大行数（默认 100）"
    )
    args = parser.parse_args()

    # 打开参考基因组
    try:
        fasta = pysam.FastaFile(args.ref)
    except Exception as e:
        print(f"无法打开参考基因组 {args.ref}: {e}", file=sys.stderr)
        sys.exit(1)

    # ---------- 第一次遍历：使用前 max_infer_rows 条可用记录推断坐标体系 ----------
    coord_system: Optional[str] = None  # '1-based' or '0-based'
    count_1based = 0
    count_0based = 0
    infer_tried = 0

    with open(args.input, newline="") as f_in:
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
            if chrom is None:
                continue

            if strand not in {"+", "-"}:
                continue

            tmp = infer_coordinate_system(fasta, chrom, site, strand)
            if tmp == "1-based":
                count_1based += 1
                infer_tried += 1
            elif tmp == "0-based":
                count_0based += 1
                infer_tried += 1
            # tmp 为 None 的行不计入 infer_tried

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

    # ---------- 第二次遍历：根据已经确定的 coord_system 进行实际筛选 ----------
    bed_path = os.path.splitext(args.output)[0] + ".bed"
    with open(args.input, newline="") as f_in, \
            open(args.output, "w", newline="") as f_out, \
            open(bed_path, "w") as bed_out:
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

            chrom = normalize_chrom_name(chrom_raw, fasta)
            if chrom is None:
                print(f"第 {row_idx} 行: 在参考基因组中找不到染色体 {chrom_raw}", file=sys.stderr)
                continue

            # 获取 9bp 窗口和中心碱基
            try:
                window_seq, center_base, window_start = get_window_and_center_base(
                    fasta, chrom, site, coord_system, flank=4
                )
            except Exception as e:
                print(f"第 {row_idx} 行获取序列窗口失败: {e}", file=sys.stderr)
                continue

            # Strand 为 '+' 时，中心碱基应该为 A，并要求窗口中含有 'AAAAA'
            if strand == "+":
                if center_base != "A":
                    print(
                        f"第 {row_idx} 行 Strand=+，期望中心碱基为 A，但实际为 {center_base} "
                        f"(Chr={chrom_raw}, Site={site}, coord={coord_system})",
                        file=sys.stderr,
                    )
                    # 一般这种情况建议跳过，不写入结果
                    continue

                motif = "AAAAA"
                if motif in window_seq:
                    # 计算 AAAAA 片段在基因组中的起止坐标（0-based, BED 风格，end 为开区间）
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

                    # BED6: chrom, start, end, name(这里用Site), score(0), strand
                    bed_out.write(f"{chrom}\t{motif_start}\t{motif_end}\t{site}\t0\t{strand}\n")

            # Strand 为 '-' 时，中心碱基应该为 T，并要求窗口中含有 'TTTTT'
            elif strand == "-":
                if center_base != "T":
                    print(
                        f"第 {row_idx} 行 Strand=-，期望中心碱基为 T，但实际为 {center_base} "
                        f"(Chr={chrom_raw}, Site={site}, coord={coord_system})",
                        file=sys.stderr,
                    )
                    # 同样建议跳过
                    continue

                motif = "TTTTT"
                if motif in window_seq:
                    # 计算 TTTTT 片段在基因组中的起止坐标（0-based, BED 风格，end 为开区间）
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

                    # BED6: chrom, start, end, name(这里用Site), score(0), strand
                    bed_out.write(f"{chrom}\t{motif_start}\t{motif_end}\t{site}\t0\t{strand}\n")

            else:
                print(
                    f"第 {row_idx} 行 Strand 既不是 '+' 也不是 '-'，为 '{strand}'，已跳过。",
                    file=sys.stderr,
                )

    fasta.close()
    print("处理完成。请检查 stderr 输出中关于坐标体系和碱基匹配的提示。", file=sys.stderr)


if __name__ == "__main__":
    main()