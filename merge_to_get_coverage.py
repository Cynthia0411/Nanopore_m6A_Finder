import pandas as pd

col_names_WT = ["Chromosome","Start","End","N1","N2","Flag","Coverage_WT", "N3", "N4", "N5"]
col_names_KO = ["Chromosome","Start","End","N1","N2","Flag","Coverage_KO", "N3", "N4", "N5"]

WT_df = pd.read_csv('WT_coverage_xgb_0.7_WT.txt', delimiter = "\t", header = None, names = col_names_WT)
KO_df = pd.read_csv('ALKBH5_coverage_xgb_0.7_WT.txt', delimiter = "\t", header = None, names = col_names_KO)

df_WT = WT_df[["Chromosome","End", "Flag", "Coverage_WT"]]
df_KO = KO_df[["Chromosome","End", "Flag", "Coverage_KO"]]
df = pd.read_csv('merged_xgb_0.7_WT_ALKBH5.tsv', sep = "\t")

# Convert multiple columns to string
columns_to_convert = ['Chromosome']
df_WT[columns_to_convert] = df_WT[columns_to_convert].astype(str)
df_KO[columns_to_convert] = df_KO[columns_to_convert].astype(str)
df['Flag'] = df['Flag'].replace({0: '+', 16: '-'})
print("rename_done!")

# Merge dataframes
merged_df = pd.merge(df, df_WT, on=['Chromosome', 'End', 'Flag'], how='left')
merged_df = pd.merge(merged_df, df_KO, on=['Chromosome', 'End', 'Flag'], how='left')

# Fill NA values with 0
merged_df.fillna(0, inplace=True)

# Calculate m6A levels (Count_Reads/coverage)
EPSILON = 1e-10  # Small number to avoid division by zero
merged_df['m6A_level_WT'] = merged_df['count_reads_WT'] / (merged_df['Coverage_WT'] + EPSILON)
merged_df['m6A_level_KO'] = merged_df['count_reads_KO'] / (merged_df['Coverage_KO'] + EPSILON)

# Save the results
merged_df.to_csv('merged_xgb_0.7_WT_ALKBH5_site_coverage.tsv', sep='\t', index=False)

print("Processing complete! Added m6A level columns to the output file.")
