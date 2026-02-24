# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
#from matplotlib_venn import venn2
import matplotlib.pyplot as plt

def read_m6A_file(file_path):
	# 读取m6A位点文件并返回位点集合
	df = pd.read_csv(file_path, sep='\t')
	# 使用染色体、位置和链作为唯一标识
	sites = set(zip(df['Chromosome'], df['End'], df['Flag']))
	return sites

def read_bed_file(file_path):
	#读取BED格式文件并返回位点集合
	df = pd.read_csv(file_path, sep='\t', header=None)
	# 使用染色体、位置和链作为唯一标识
	sites = set(zip(df[0], df[2], df[5]))  # 假设列0是染色体，列2是end位置，列5是链
	return sites

def calculate_overlap(input_m6Apos_file, input_refpos_file):
	#计算m6A位点文件和参考位点bed文件之间的重叠
	
	# 读取位点文件
#    m6A_sites = read_m6A_file(input_m6Apos_file)
    m6A_sites = read_bed_file(input_m6Apos_file)
    ref_sites = read_bed_file(input_refpos_file)
	
	# 计算重叠
    overlap = m6A_sites.intersection(ref_sites)
	
	# 计算独特的位点
    only_m6A = m6A_sites - ref_sites
    only_ref = ref_sites - m6A_sites
	
	# 打印结果
    print(f"Total m6A sites: {len(m6A_sites)}")
    print(f"Total reference sites: {len(ref_sites)}")
    print(f"Overlapping sites: {len(overlap)}")
    print(f"m6A-specific sites: {len(only_m6A)}")
    print(f"Reference-specific sites: {len(only_ref)}")
	
	# 绘制Venn图
#	plt.figure(figsize=(10, 10))
#	venn2(subsets=(len(only_m6A), len(only_ref), len(overlap)),
#		  set_labels=('m6A sites', 'Reference sites'))
#	plt.title('Overlap between m6A sites and Reference sites')
	
	# 保存结果
    output_overlap = input_m6Apos_file.replace('.bed', '_overlap.tsv')
#	output_m6A_specific = input_m6Apos_file.replace('.tsv', '_specific.tsv')
	
	# 将结果转换回DataFrame并保存
    overlap_df = pd.DataFrame(list(overlap), columns=['chrom', 'pos', 'strand'])
#	m6A_specific_df = pd.DataFrame(list(only_m6A), columns=['chrom', 'pos', 'strand'])
	
	# 添加必要的BED格式列
#	for df in [overlap_df, m6A_specific_df]:
    overlap_df['start'] = overlap_df['pos'] - 1  # BED格式是0-based
    overlap_df['name'] = 'm6A_site'
    overlap_df['score'] = '.'
    overlap_df = overlap_df[['chrom', 'start', 'pos', 'name', 'score', 'strand']]
	
	# 保存文件
    overlap_df.to_csv(output_overlap, sep='\t', index=False, header=True)
#	m6A_specific_df.to_csv(output_m6A_specific, sep='\t', index=False, header=True)
	
#	plt.savefig(input_m6Apos_file.replace('.tsv', '_overlap.png'))
#	plt.close()
	
    return len(m6A_sites), len(ref_sites), len(overlap)

if __name__ == "__main__":
    import sys
	
    if len(sys.argv) != 3:
        print("Usage: python calculate_overlap.py <input_m6A_tsv> <input_ref_bed>")
        sys.exit(1)
	
    input_m6A_tsv = sys.argv[1]
    input_ref_bed = sys.argv[2]
    calculate_overlap(input_m6A_tsv, input_ref_bed)
