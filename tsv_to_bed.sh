#!/bin/bash

# tsv_to_bed.sh
# 将TSV文件转换为BED格式，使用m6A_level_WT作为第五列score值

# 设置输入输出文件
INPUT_FILE="wilcoxon_xgb_WT_KO_results_significant_3.tsv"
OUTPUT_FILE="predicted_m6A_sites_WT_KO.bed"

# 检查输入文件
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在"
    exit 1
fi

# 转换文件格式
awk 'BEGIN{OFS="\t"}
NR>1 {  # 跳过标题行
    # 提取染色体、起始、结束位置、m6A_level_WT值和正负链信息
    chr=$1
    start=$2
    end=$3
    strand=$4  # Flag列，+ 或 -
    wt_mean=$10  # m6A_level_WT
    
    # 生成name列（使用位置信息）
    name = chr ":" start "-" end
    
    # 输出BED格式：chr start end name score strand
    print chr, start, end, name, wt_mean, strand
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "转换完成。输出文件: $OUTPUT_FILE"
echo "前5行预览："
head -n 5 "$OUTPUT_FILE"
