# -*- coding: utf-8 -*-
"""
从两个 BED 文件计算重叠位点，输出 CSV：Chr, pos, Strand, score_bed1, score_bed2
"""
import sys
import pandas as pd


def read_bed_with_score(file_path):
    """读取 BED（无表头）：chrom, start, end, name, score, strand -> (chrom, end, strand) 为位点，带 score"""
    df = pd.read_csv(file_path, sep='\t', header=None)
    # 标准 6 列 BED: 0=chrom, 1=start, 2=end, 3=name, 4=score, 5=strand；用 end 作为位点
    df = df.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 3: 'name', 4: 'score', 5: 'strand'})
    df['pos'] = df['end'].astype(int)
    df['chrom'] = df['chrom'].astype(str)
    df['strand'] = df['strand'].astype(str)
    return df[['chrom', 'pos', 'strand', 'score']]


def calculate_overlap_to_csv(bed1_path, bed2_path, output_csv=None):
    # 读取两个 BED，保留 score
    df1 = read_bed_with_score(bed1_path)
    df2 = read_bed_with_score(bed2_path)

    # 位点键 (chrom, pos, strand)
    df1['_key'] = list(zip(df1['chrom'], df1['pos'], df1['strand']))
    df2['_key'] = list(zip(df2['chrom'], df2['pos'], df2['strand']))

    set1 = set(df1['_key'])
    set2 = set(df2['_key'])
    overlap_keys = set1 & set2

    # 取重叠位点，并带上两个文件的 score
    score1 = df1.set_index('_key')['score']
    score2 = df2.set_index('_key')['score']

    rows = []
    for key in overlap_keys:
        chrom, pos, strand = key
        rows.append({
            'Chr': chrom,
            'pos': pos,
            'Strand': strand,
            'score_bed1': score1[key],
            'score_bed2': score2[key],
        })

    overlap_df = pd.DataFrame(rows)
    if output_csv is None:
        output_csv = bed1_path.replace('.bed', '_overlap.csv')
    overlap_df.to_csv(output_csv, index=False)

    print(f"BED1 位点数: {len(set1)}")
    print(f"BED2 位点数: {len(set2)}")
    print(f"重叠位点数: {len(overlap_keys)}")
    print(f"已保存: {output_csv}")
    return len(set1), len(set2), len(overlap_keys)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python calculate_overlap_from_bed_to_csv.py <bed1> <bed2> [output.csv]")
        sys.exit(1)
    bed1 = sys.argv[1]
    bed2 = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None
    calculate_overlap_to_csv(bed1, bed2, out)
