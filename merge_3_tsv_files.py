import pandas as pd

# 读取三个TSV文件
df1 = pd.read_csv('merged_xgb_0.7_WT1_KO3_site_coverage.tsv', sep='\t')
df2 = pd.read_csv('merged_xgb_0.7_WT2_KO2_site_coverage.tsv', sep='\t')
df3 = pd.read_csv('merged_xgb_0.7_WT3_KO1_site_coverage.tsv', sep='\t')

# 合并三个数据框，保留在"Chromosome"、"Start"、"End"、"Flag"列相同的行
merged_df = pd.merge(df1, df2, on=['Chromosome', 'Start', 'End', 'Flag'], how='inner')
merged_df = pd.merge(merged_df, df3, on=['Chromosome', 'Start', 'End', 'Flag'], how='inner')

# 将合并后的数据保存为新的TSV文件
merged_df.to_csv('merged_xgb_0.7_WT1_WT2_WT3_coverage.tsv', sep='\t', index=False)


