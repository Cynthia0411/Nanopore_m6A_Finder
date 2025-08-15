import pandas as pd

def extract_and_transform_columns(input_file, output_file):
    """
    Extract specific columns and transform to BED format
    BED format: chromosome, start, end, name, score, strand
    """
    try:
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
        
        
        # Create BED format dataframe
        bed_df = pd.DataFrame()
        bed_df['chromosome'] = df['Chromosome']
        bed_df['start'] = df['Start']
        bed_df['end'] = df['End']
        bed_df['name'] = '.'  # Default name field
        bed_df['score'] = 0   # Default score field
        
        # Transform Flag values to strand
        bed_df['strand'] = df['Flag'].map({0: '+', 16: '-'})
        
        # Save as BED file (tab-separated, no header)
        bed_df.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"BED format file saved to: {output_file}")
        
    except Exception as e:
        print(f"Error during extraction and transformation: {str(e)}")

def merge_files(wt_file, ko_file, output_path):
    """
    Merge WT and KO files based on Chromosome, Start, End, and Flag columns
    """
    try:
        # Read the TSV files, skip the index column
        df_wt = pd.read_csv(wt_file, sep='\t', index_col=0)
        df_ko = pd.read_csv(ko_file, sep='\t', index_col=0)
        
        # Merge files based on Chromosome, Start, End, and Flag
        merged_df = pd.merge(
            df_wt,
            df_ko,
            on=['Chromosome', 'Start', 'End', 'Flag'],
            how='left',
            suffixes=('_WT', '_KO')
        )
        
        # Fill NaN values with 0
        merged_df = merged_df.fillna(0)
        
        # Save the merged result without index column
        merged_df.to_csv(output_path, sep='\t', index=False)
        print(f"Files successfully merged and saved to: {output_path}")
        
    except Exception as e:
        print(f"Error during merging: {str(e)}")

def merge_files_outer(wt_file, ko_file, output_path, bed_output=None):
    """
    Merge WT and KO files based on Chromosome, End, and Flag columns using outer join
    and optionally create a BED format file
    
    Args:
        wt_file: Path to WT file
        ko_file: Path to KO file
        output_path: Path for output merged file
        bed_output: Optional path for BED format output file
    """
    try:
        # Read the TSV files
        print(f"Reading WT file: {wt_file}")
        df_wt = pd.read_csv(wt_file, sep='\t')
        print(f"WT file contains {len(df_wt)} rows")
        
        print(f"Reading KO file: {ko_file}")
        df_ko = pd.read_csv(ko_file, sep='\t')
        print(f"KO file contains {len(df_ko)} rows")
        
        # Merge files based on Chromosome, End, and Flag
        print("\nMerging files...")
        merged_df = pd.merge(
            df_wt,
            df_ko,
            on=['Chromosome','Start', 'End', 'Flag'],
            how='outer',  # Use outer join to get all rows
            suffixes=('_WT', '_KO')
        )
        
        # Sort the merged dataframe
        merged_df = merged_df.sort_values(['Chromosome', 'Start','End', 'Flag'])
        
        # Fill NaN values with 0
        merged_df = merged_df.fillna(0)
        
        # Save the merged result
        print(f"\nSaving merged results to: {output_path}")
        merged_df.to_csv(output_path, sep='\t', index=False)
        print(f"Merged file contains {len(merged_df)} rows")
        
        # Create and save BED format file if requested
        if bed_output:
            print(f"\nCreating BED format file: {bed_output}")
            bed_df = pd.DataFrame()
            bed_df['chromosome'] = merged_df['Chromosome']
            
            bed_df['start'] = merged_df['Start']
                
            bed_df['end'] = merged_df['End']
            bed_df['name'] = '.'  # Default name field
            bed_df['score'] = 0   # Default score field
            bed_df['strand'] = merged_df['Flag'].map({0: '+', 16: '-'})
            
            # Save BED file (tab-separated, no header)
            bed_df.to_csv(bed_output, sep='\t', index=False, header=False)
            print(f"BED format file saved with {len(bed_df)} rows")
        
        # Print merge statistics
        print("\nMerge statistics:")
        print(f"Rows only in WT: {len(df_wt) - len(df_wt.merge(df_ko, on=['Chromosome', 'End', 'Flag']))}")
        print(f"Rows only in KO: {len(df_ko) - len(df_ko.merge(df_wt, on=['Chromosome', 'End', 'Flag']))}")
        print(f"Rows in both files: {len(df_wt.merge(df_ko, on=['Chromosome', 'End', 'Flag']))}")
        
    except Exception as e:
        print(f"Error during merging: {str(e)}")

# Usage
wt_file = "intersections_m6A_sites_xgb_0.7_WT.tsv"
ko_file = "intersections_m6A_sites_xgb_0.7_ALKBH5.tsv"
#output_path_outer = "merged_xgb_0.7_WT3_union_xgb_0.7_KO1.tsv"
output_path = "merged_xgb_0.7_WT_ALKBH5.tsv"
extracted_output = "extracted_xgb_0.7_WT.bed"
#bed_output = "xgb_0.7_WT3_union_xgb_0.7_KO1.bed"
#wt_file = "overlap_WT1_WT2_combined_1308.tsv"
#ko_file = "overlap_sites_of_xgb_WT3_union_xgb_KO1_1308_combined.tsv"
#output_path_outer = "overlap_WT1_WT2_WT3_combined_1308.tsv"
#bed_output = None


# Extract and transform columns from WT1 file
extract_and_transform_columns(wt_file, extracted_output)

# Original merge function remains available if needed
merge_files(wt_file, ko_file, output_path)

# Use the new outer merge function
#merge_files_outer(wt_file, ko_file, output_path_outer, bed_output=bed_output)
