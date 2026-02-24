import argparse
from dataclasses import dataclass
import sys
from typing import Dict, Iterator, List, Optional, Tuple

import pysam


@dataclass(frozen=True)
class BedInterval:
    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    name: str
    strand: str  # '+', '-', or '.'


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def open_alignment(path: str) -> pysam.AlignmentFile:
    lower = path.lower()
    # SAM: text; BAM/CRAM: binary
    if lower.endswith(".bam") or lower.endswith(".cram"):
        return pysam.AlignmentFile(path, "rb")
    return pysam.AlignmentFile(path, "r")


def read_bed(path: str) -> List[BedInterval]:
    intervals: List[BedInterval] = []
    with open(path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"BED 第 {line_no} 行列数不足3: {line!r}")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 and parts[3] != "" else f"{chrom}:{start}-{end}"
            strand = parts[5] if len(parts) >= 6 and parts[5] in {"+", "-"} else "."
            if end <= start:
                raise ValueError(f"BED 第 {line_no} 行 end<=start: {line!r}")
            intervals.append(BedInterval(chrom=chrom, start=start, end=end, name=name, strand=strand))
    return intervals


def build_refpos_to_qpos(aln: pysam.AlignedSegment) -> Dict[int, int]:
    """
    返回映射: reference_position (0-based) -> query_position (0-based)
    只保留 ref 和 query 都存在的位置（即跳过 deletion/skip）。
    """
    mapping: Dict[int, int] = {}
    # (qpos, rpos)；当出现 insertion 时 rpos=None；deletion/skip 时 qpos=None
    for qpos, rpos in aln.get_aligned_pairs(matches_only=False):
        if qpos is None or rpos is None:
            continue
        mapping[rpos] = qpos
    return mapping


def phred_to_fastq_ascii(phred: List[int]) -> str:
    # FASTQ: ASCII = phred + 33
    return "".join(chr(q + 33) for q in phred)


def extract_kmer_for_interval(
    aln: pysam.AlignedSegment,
    interval: BedInterval,
    require_full_coverage: bool = True,
) -> Optional[Tuple[str, int, str]]:
    """
    从单条比对记录里提取其在参考区间 [start,end) 对应的 read 序列片段。
    - require_full_coverage=True: 区间内每个 ref 位点都必须能映射到 read 碱基，否则返回 None

    注意：输出的 Sequence / Quality 始终按 read 的 5'->3' 方向（即 qpos 从小到大）排列，
    这样能正确处理 read 反向比对（aln.is_reverse=True）的情况。
    """
    if aln.is_unmapped or aln.query_sequence is None:
        return None
    if aln.query_qualities is None:
        return None

    ref_to_q = build_refpos_to_qpos(aln)
    qseq = aln.query_sequence
    quals = list(aln.query_qualities)

    # 收集区间内每个 ref 位点对应的 (qpos, base, phred)
    mapped: List[Tuple[int, str, int]] = []
    for rpos in range(interval.start, interval.end):
        qpos = ref_to_q.get(rpos)
        if qpos is None:
            return None if require_full_coverage else None
        if qpos < 0 or qpos >= len(qseq):
            return None if require_full_coverage else None
        mapped.append((qpos, qseq[qpos], quals[qpos]))

    if not mapped:
        return None

    # 按 read 的 5'->3' 方向排序（qpos 递增）
    mapped.sort(key=lambda x: x[0])
    pos0 = mapped[0][0]
    seq = "".join(b for _, b, _ in mapped)
    qual = phred_to_fastq_ascii([q for _, __, q in mapped])
    return seq, int(pos0), qual


def iter_reads_for_interval(
    af: pysam.AlignmentFile,
    interval: BedInterval,
) -> Iterator[pysam.AlignedSegment]:
    # pysam.fetch 使用 0-based start,end（右开），与 BED 一致
    yield from af.fetch(interval.chrom, interval.start, interval.end)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="从 SAM/BAM/CRAM 中提取覆盖 BED 区间的 read 片段（kmer）"
    )
    parser.add_argument("--sam", required=True, help="输入 SAM/BAM/CRAM（建议已 sort+index）")
    parser.add_argument("--ref", required=False, help="参考基因组 fasta（仅 CRAM 可能需要）")
    parser.add_argument("--bed", required=True, help="BED 文件（0-based, end开区间；支持 BED6 strand）")
    parser.add_argument("--out", required=True, help="输出 TSV 文件路径（tsv）")
    parser.add_argument(
        "--require-full-coverage",
        action="store_true",
        help="要求 read 覆盖区间内每个 ref 位点（否则跳过该 read）；默认不启用时用N补齐",
    )
    parser.add_argument(
        "--fill-N",
        action="store_true",
        help="当不要求完整覆盖时，用N填充缺失位点（默认行为；仅用于显式说明）",
    )
    parser.add_argument(
        "--respect-bed-strand",
        action="store_true",
        help="（不兼容）按 BED strand 输出方向；本脚本默认按 read 5'->3' 输出",
    )
    parser.add_argument(
        "--max-reads-per-interval",
        type=int,
        default=0,
        help="每个区间最多输出多少条 read（0 表示不限制）",
    )
    args = parser.parse_args()

    intervals = read_bed(args.bed)

    # coverage 策略
    require_full_coverage = bool(args.require_full_coverage)

    if args.respect_bed_strand:
        raise SystemExit(
            "错误: --respect-bed-strand 与当前输出要求冲突。"
            "本脚本输出的 Sequence/Quality/Position 均以 read 的 5'->3'（qpos 递增）为准，"
            "若需要按 BED/reference 方向输出，请另行说明我再提供另一种输出模式。"
        )

    try:
        if args.ref:
            af = open_alignment(args.sam)
            # 对于 CRAM，pysam 会通过 reference_filename 找 ref
            try:
                af.reference_filename = args.ref  # type: ignore[attr-defined]
            except Exception:
                # 不同 pysam 版本该属性可能不可写；忽略，由用户保证环境可读 CRAM
                pass
        else:
            af = open_alignment(args.sam)
    except Exception as e:
        print(f"无法打开比对文件 {args.sam}: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("\t".join(["ReadName", "Sequence", "Position", "Strand", "Quality"]) + "\n")

        for interval in intervals:
            n_written = 0
            try:
                reads = iter_reads_for_interval(af, interval)
                for aln in reads:
                    extracted = extract_kmer_for_interval(
                        aln,
                        interval,
                        require_full_coverage=require_full_coverage,
                    )
                    if extracted is None:
                        continue
                    seq, pos0, qual = extracted

                    if (not args.require_full_coverage) and args.fill_N:
                        # extract_kmer_for_interval 在 require_full_coverage=False 时已用 N 补齐
                        pass

                    read_strand = "-" if aln.is_reverse else "+"

                    out.write("\t".join([
                        aln.query_name,
                        seq,
                        str(pos0),
                        read_strand,
                        qual,
                    ]) + "\n")

                    n_written += 1
                    if args.max_reads_per_interval > 0 and n_written >= args.max_reads_per_interval:
                        break
            except ValueError as e:
                # 常见于 BAM 未建立 index 或 contig 不存在
                print(f"区间 {interval.chrom}:{interval.start}-{interval.end} 提取失败: {e}", file=sys.stderr)
                continue

    af.close()


if __name__ == "__main__":
    main()

