#!/bin/bash

# bed_to_gene.sh
# 用途：从GTF提取基因位置并将BED文件中的位置对应到基因

# 设置错误时退出
set -e

# 定义输入文件
GTF_FILE="hg38.gtf"  # 请替换为你的GTF文件路径
BED_FILE="GLORI_v2.0_50ng.bed"       # 请替换为你的BED文件路径
OUTPUT_DIR="gene_results_GLORI_v2_50ng"

# 检查输入文件
if [ ! -f "$GTF_FILE" ]; then
    echo "错误: GTF文件 $GTF_FILE 不存在"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "错误: BED文件 $BED_FILE 不存在"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

echo "开始处理..."

# 从GTF文件提取基因信息，并添加chr前缀
echo "提取基因信息..."
awk -F"\t" '$3=="gene" {
    split($9,a,";");
    gene_name="";
    for(i in a) {
        if(a[i] ~ /gene_name/) {
            split(a[i],b,"\"");
            gene_name=b[2];
            break;
        }
    }
    if(gene_name!="") {
        # 添加chr前缀，除非已经有chr前缀或是特殊染色体名
        chr=$1;
        if(chr !~ /^chr/ && chr !~ /^MT/ && chr !~ /^M/ && chr !~ /^NC/) {
            chr="chr"chr;
        }
        print chr"\t"$4-1"\t"$5"\t"gene_name"\t.\t"$7
    }
}' "$GTF_FILE" > "$OUTPUT_DIR/genes.bed"

# 对基因文件排序
sort -k1,1 -k2,2n "$OUTPUT_DIR/genes.bed" > "$OUTPUT_DIR/genes.sorted.bed"

# 对输入BED文件排序
sort -k1,1 -k2,2n "$BED_FILE" > "$OUTPUT_DIR/input.sorted.bed"

# 使用bedtools intersect找到重叠
echo "查找基因重叠..."
bedtools intersect -a "$OUTPUT_DIR/input.sorted.bed" -b "$OUTPUT_DIR/genes.sorted.bed" -wa -wb > "$OUTPUT_DIR/overlaps.bed"

# 处理重叠结果，合并同一位点的多个基因名称
echo "更新基因名称..."
awk 'BEGIN{OFS="\t"} 
{
    key=$1":"$2"-"$3
    if(!(key in genes)) {
        chr[key]=$1
        start[key]=$2
        end[key]=$3
        score[key]=$5
        strand[key]=$6
        genes[key]=$10
    } else {
        genes[key]=genes[key]";"$10
    }
}
END{
    for(key in genes) {
        print chr[key], start[key], end[key], genes[key], score[key], strand[key]
    }
}' "$OUTPUT_DIR/overlaps.bed" > "$OUTPUT_DIR/final.bed"

# 统计结果
echo "统计结果..."
echo "输入位点数: $(wc -l < "$BED_FILE")"
echo "找到基因注释的位点数: $(wc -l < "$OUTPUT_DIR/final.bed")"
echo "基因总数: $(wc -l < "$OUTPUT_DIR/genes.bed")"

echo "处理完成。结果保存在 $OUTPUT_DIR/final.bed"

# 显示前几行结果
echo "结果示例（前5行）："
head -n 5 "$OUTPUT_DIR/final.bed"
