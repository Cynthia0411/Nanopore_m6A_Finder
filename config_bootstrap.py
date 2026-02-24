# Analysis parameters
N_BOOTSTRAP = 10000
ALPHA = 0.05
EPSILON = 1e-10  # Small number to avoid division by zero

# Test method
USE_WILCOXON = True  # Set to True to use Wilcoxon test, False for bootstrap

# File paths
INPUT_FILE = "merged_xgb_WT_KO_site_coverage_3.tsv"
OUTPUT_ALL_FILE = "statistical_results_all.tsv"  # Generic name for both methods
OUTPUT_SIG_FILE = "statistical_results_significant.tsv"  # Generic name for both methods

# Column names
COLUMNS = {
    'wt1_interval': 'Count_Read_Interval_WT',
    'ko1_interval': 'Count_Read_Interval_KO',
    'wt1_coverage': 'Coverage_WT',
    'ko1_coverage': 'Coverage_KO',
    'chromosome': 'Chromosome',
    'start': 'Start',
    'end': 'End',
    'flag': 'Flag',
    'wt1_m6a_level': 'm6A_level_WT',
    'ko1_m6a_level': 'm6A_level_KO'
}

# Output columns for display
DISPLAY_COLUMNS = [
    'Chromosome', 'Start', 'End', 'Flag',
    'WT_mean', 'KO_mean', 'p_value', 'FDR', 'Significant',
    'm6A_level_WT', 'm6A_level_KO'
]

# Get output file names based on test method
def get_output_files():
    prefix = "wilcoxon" if USE_WILCOXON else "bootstrap"
    return (
        f"{prefix}_xgb_WT_KO_results_all_3.tsv",
        f"{prefix}_xgb_WT_KO_results_significant_3.tsv"
    )

OUTPUT_ALL_FILE, OUTPUT_SIG_FILE = get_output_files()
