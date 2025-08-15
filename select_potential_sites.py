import pandas as pd

# Step 1: Read in the DataFrame from a CSV file
#filename = 'merged_WT_KO_site_coverage.tsv'
#filename = 'merged_xgb_0.7_WT3_KO1_site_coverage.tsv'
#filename = "merged_processed_xgb_0.7_WT1_WT2_WT3.tsv"
#filename = 'bootstrap_results_coverage10_WT1_KO1.tsv'
filename ='wilcoxon_xgb_0.7_WT_KO_results_all_3.tsv'

df = pd.read_csv(filename,sep = "\t")

# Step 2: Define multiple conditions for selecting rows
# Example condition: Select rows where the value in 'column1' is greater than 10,
# the value in 'column2' is equal to 'A', and the value in 'column3' is not null

#condition = (df['Coverage_WT1'] > 49) & (df['Coverage_WT2'] > 49) & (df['Coverage_WT3']> 49) & (df['Coverage_KO1']>49) & (df['Coverage_KO2']>49) & (df['Coverage_KO3']>49)

#condition = ((df["count_interval_WT1"]/df["count_read_WT1"])>1.2) & ((df["count_interval_WT2"]/df["count_read_WT2"])>1.2) & ((df["count_interval_WT3"]/df["count_read_WT3"])>1.2) & ((df["count_interval_KO1"]/df["count_read_KO1"])>1.2) & ((df["count_interval_KO2"]/df["count_read_KO2"])>1.2) & ((df["count_interval_KO3"]/df["count_read_KO3"])>1.2)

#condition =(df["Chromosome"] != "scplasm1") & (df["Coverage_WT"]>2) & (df["Coverage_KO"]>2)
condition = (df["Significant"]== False)

# Step 3: Apply the conditions to the DataFrame to filter rows
filtered_df = df[condition]

# Optional: Save the filtered DataFrame to a new CSV file
#output_filename = 'intersections_xgb_0.7_WT1_WT2_WT3_coverage_3.tsv'
output_filename= "wilcoxon_xgb_0.7_WT_KO_results_not_significant_3.tsv"
filtered_df.to_csv(output_filename, sep = "\t", index=False)

