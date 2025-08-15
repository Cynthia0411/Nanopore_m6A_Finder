import pandas as pd
import ast
import numpy as np
from scipy import stats
from tqdm import tqdm
import time
from config_bootstrap import (
    N_BOOTSTRAP,
    ALPHA,
    EPSILON,
    INPUT_FILE,
    OUTPUT_ALL_FILE,
    OUTPUT_SIG_FILE,
    COLUMNS,
    DISPLAY_COLUMNS,
    USE_WILCOXON
)

def estimate_runtime(df, n_bootstrap):
    """
    Estimate the runtime based on dataset size and number of bootstrap iterations
    
    Args:
        df: input DataFrame
        n_bootstrap: number of bootstrap iterations
    """
    # Estimate ~0.001 seconds per row per 1000 iterations (approximate baseline)
    estimated_seconds = len(df) * (n_bootstrap/1000) * 0.001
    
    # Convert to hours, minutes, seconds
    hours = int(estimated_seconds // 3600)
    minutes = int((estimated_seconds % 3600) // 60)
    seconds = int(estimated_seconds % 60)
    
    print(f"Estimated runtime: {hours}h {minutes}m {seconds}s")

def process_interval(interval_str, coverage):
    """
    Process interval string and pad with zeros up to coverage value.
    Convert any positive numbers to 1 before padding.
    """
    if interval_str == '0' or pd.isna(interval_str):
        return [0] * int(coverage)
    interval_list = ast.literal_eval(interval_str)
    interval_list = [int(x) for x in interval_list]
    # Convert positive numbers to 1
#    interval_list = [1 if x > 0 else 0 for x in interval_list]
    padding = int(coverage) - len(interval_list)
    if padding > 0:
        interval_list.extend([0] * padding)
    return interval_list

def bootstrap_test(wt1_list, ko1_list, n_iterations=10000):  # Changed from 1000 to 10000
    """
    Perform bootstrap test to compare means of two lists by resampling both groups together
    
    Args:
        wt1_list: list of values from WT1
        ko1_list: list of values from KO1
        n_iterations: number of bootstrap iterations (default: 10000)
    """
    # Calculate original difference in means
    orig_diff = np.mean(wt1_list) - np.mean(ko1_list)
    
    # Combine data for permutation
    combined = np.concatenate([wt1_list, ko1_list])
    n1 = len(wt1_list)
    n2 = len(ko1_list)
    
    # Perform bootstrap
    count_greater = 0
    for _ in range(n_iterations):
        np.random.shuffle(combined)
        perm_wt1 = combined[:n1]
        perm_ko1 = combined[n1:n1+n2]
        perm_diff = np.mean(perm_wt1) - np.mean(perm_ko1)
        if abs(perm_diff) >= abs(orig_diff):
            count_greater += 1
    
    return count_greater / n_iterations, orig_diff

def wilcoxon_test(wt1_list, ko1_list):
    """
    Perform Wilcoxon rank-sum test to compare two lists
    
    Args:
        wt1_list: list of values from WT1
        ko1_list: list of values from KO1
    Returns:
        p_value: two-sided p-value
        diff: difference in means (kept for consistency with bootstrap output)
    """
    # Calculate difference in means (kept for consistency with bootstrap output)
    diff = np.mean(wt1_list) - np.mean(ko1_list)
    
    # Perform Wilcoxon rank-sum test
    _, p_value = stats.ranksums(wt1_list, ko1_list)
    
    return p_value, diff

def analyze_count_intervals(input_file=INPUT_FILE, 
                          output_all_file=OUTPUT_ALL_FILE, 
                          output_sig_file=OUTPUT_SIG_FILE, 
                          n_bootstrap=N_BOOTSTRAP, 
                          alpha=ALPHA,
                          use_wilcoxon=USE_WILCOXON):  # Updated to use config parameter
    """
    Analyze count intervals using either bootstrap or Wilcoxon test to compare WT1 and KO1 with FDR correction
    """
    try:
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
        print(f"\nInput file contains {len(df)} rows")
        print(f"Using {'Wilcoxon' if use_wilcoxon else 'Bootstrap'} test method")
        
        # Estimate runtime only for bootstrap method
        if not use_wilcoxon:
            estimate_runtime(df, n_bootstrap)
        
        # Process and analyze each row
        results = []
        start_time = time.time()
        
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Analyzing intervals"):
            wt1_intervals = process_interval(str(row[COLUMNS['wt1_interval']]), row[COLUMNS['wt1_coverage']])
            ko1_intervals = process_interval(str(row[COLUMNS['ko1_interval']]), row[COLUMNS['ko1_coverage']])
            
            # Skip if both intervals are all zeros
            if sum(wt1_intervals) == 0 and sum(ko1_intervals) == 0:
                continue
            
            # Choose test method
            if use_wilcoxon:
                p_value, diff = wilcoxon_test(wt1_intervals, ko1_intervals)
            else:
                p_value, diff = bootstrap_test(wt1_intervals, ko1_intervals, n_bootstrap)
                
            wt1_mean = np.mean(wt1_intervals)
            ko1_mean = np.mean(ko1_intervals)
            
            # Calculate fold change
            fold_change = np.log2((wt1_mean + EPSILON) / (ko1_mean + EPSILON))
            
            results.append({
                'Chromosome': row[COLUMNS['chromosome']],
                'Start': row[COLUMNS['start']],
                'End': row[COLUMNS['end']],
                'Flag': row[COLUMNS['flag']],
                'WT_mean': wt1_mean,
                'KO_mean': ko1_mean,
                'Difference': diff,
                'Fold_Change': fold_change,
                'p_value': p_value,
                'm6A_level_WT': row[COLUMNS['wt1_m6a_level']],
                'm6A_level_KO': row[COLUMNS['ko1_m6a_level']]
            })
        
        # Calculate actual runtime
        actual_time = time.time() - start_time
        hours = int(actual_time // 3600)
        minutes = int((actual_time % 3600) // 60)
        seconds = int(actual_time % 60)
        print(f"\nActual runtime: {hours}h {minutes}m {seconds}s")
        
        # Process results
        results_df = pd.DataFrame(results)
        
        if len(results_df) > 0:
            # Sort by p-value
            results_df = results_df.sort_values('p_value')
            
            # Calculate FDR using Benjamini-Hochberg method
            p_values = results_df['p_value'].values
            n_tests = len(p_values)
            
            # Calculate FDR
            rank = np.arange(1, n_tests + 1)
            fdr_values = p_values * n_tests / rank
            
            # Ensure FDR values are monotonically increasing
            for i in range(n_tests-2, -1, -1):
                fdr_values[i] = min(fdr_values[i], fdr_values[i+1])
            
            # Add FDR to results
            results_df['FDR'] = fdr_values
            
            # Update significance based on FDR
            results_df['Significant'] = (results_df['p_value'] < alpha) & (results_df['Fold_Change'] > 0)
        
            # Save all results to TSV file
            results_df.to_csv(output_all_file, sep='\t', index=False)
            print(f"All results saved to: {output_all_file}")
            
            # Filter and save significant rows
            significant_rows = results_df[results_df['Significant']]
            significant_rows.to_csv(output_sig_file, sep='\t', index=False)
            print(f"Significant results saved to: {output_sig_file}")
            
            print(f"\nFound {len(significant_rows)} significant regions where WT1 > KO1 (p_value < {alpha})")
            
            # Print summary statistics
            print("\nSummary statistics:")
            print(f"Total tests: {n_tests}")
            print(f"Significant at p_value < {alpha}: {sum(results_df['Significant'])}")
            print(f"Minimum FDR: {min(fdr_values):.6f}")
            print(f"Median FDR: {np.median(fdr_values):.6f}")
            
            return results_df, significant_rows
        else:
            print("No results to process")
            return None, None
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        return None, None

# Usage
if __name__ == "__main__":
    # Use the configuration from config_bootstrap.py
    results_df, significant_df = analyze_count_intervals()

    # Display results
    if results_df is not None:
        print("\nFirst few rows of complete results:")
        print(results_df[DISPLAY_COLUMNS].head())
        
        if len(significant_df) > 0:
            print("\nFirst few significant regions where WT1 > KO1:")
            print(significant_df[DISPLAY_COLUMNS].head())
            
            # Print FDR distribution for significant results
            print("\nFDR distribution for significant results:")
            print(significant_df['FDR'].describe())
