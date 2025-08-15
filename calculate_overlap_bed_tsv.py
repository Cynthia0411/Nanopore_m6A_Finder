import pandas as pd
#from scipy.stats import fisher_exact

def read_in_bed(filename):
	# 读取bed格式文件，列名为：chromosome, start, end, name, score, strand
	df = pd.read_csv(filename, sep="\t", 
					names=["chrom", "start", "end", "name", "score", "strand"],
					header=None)
	return df

def read_in_df(filename):
	df = pd.read_csv(filename, sep="\t", header=0)
	return df

def transform_to_list(df_pos):
	index_pos = df_pos.index.values
	pos_list = []
	for i in index_pos:
		chrom_cur = df_pos["Chromosome"][i]
		flag_cur = str(df_pos["Flag"][i])
		pos_cur = df_pos["End"][i]
		pos_list.append((chrom_cur, pos_cur, flag_cur))
	return pos_list

def main():
	# Input and output files
	input_refpos_file = "predicted_m6A_sites_HEK293.bed"  # 改为bed文件
	input_m6Apos_file = "merged_xgb_0.7_WT_ALKBH5_site_coverage.tsv"
	output_overlappos_file = "overlap_sites_of_xgb_0.7_WT_ALKBH5_HEK293.tsv"    
	
	# 读取输入文件
	df_ref = read_in_bed(input_refpos_file)  # 使用新的bed文件读取函数
	df_pos = read_in_df(input_m6Apos_file)
	
	# 创建参考集
	ref_list = []
	for i in df_ref.index.values:
		chromosome = df_ref["chrom"][i]
		flag = df_ref["strand"][i]  # bed文件中的strand列
		site = df_ref["end"][i]     # 使用bed文件中的end位置
		ref_list.append((chromosome, site, flag))
	ref_set = set(ref_list)
	
	# 创建位置集
	pos_list = transform_to_list(df_pos)
	pos_set = set(pos_list)
	
	# 找到交集
	intersection_result = ref_set & pos_set
	
	# 过滤df_pos以保留重叠的行
	overlapped_rows = []
	for result in intersection_result:
		mask = (df_pos["Chromosome"] == result[0]) & \
			   (df_pos["End"] == result[1]) & \
			   (df_pos["Flag"].astype(str) == result[2])
		overlapped_rows.append(df_pos[mask])
	
	# 合并所有重叠的行并保存到文件
	if overlapped_rows:
		final_df = pd.concat(overlapped_rows)
		final_df.to_csv(output_overlappos_file, sep='\t', index=False)

if __name__ == "__main__":
	main()
