import pysam
from pybedtools import BedTool
import pandas as pd

def get_sequence_and_sites(fasta_file, bed_file, window=20):
    """
    Get A and T sites within window bp of positions in bed file
    """
    # Read fasta file
    fasta = pysam.FastaFile(fasta_file)
    
    # Read original bed file
    original_bed = pd.read_csv(bed_file, sep='\t', header=None,
                             names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    new_sites = []
    
    # Process each position in bed file
    for _, row in original_bed.iterrows():
        chrom = row['chrom']
        center = row['start']
        
        # Get sequence window around site
        window_start = max(0, center - window)
        window_end = center + window
        seq = fasta.fetch(chrom, window_start, window_end)
        
        # Find A positions (+ strand)
        for i, base in enumerate(seq):
            if base.upper() == 'A':
                abs_pos = window_start + i
                new_sites.append({
                    'chrom': chrom,
                    'start': abs_pos,
                    'end': abs_pos + 1,
                    'name': 'A_site',
                    'score': '.',
                    'strand': '+'
                })
                
        # Find T positions (- strand)
        for i, base in enumerate(seq):
            if base.upper() == 'T':
                abs_pos = window_start + i
                new_sites.append({
                    'chrom': chrom,
                    'start': abs_pos,
                    'end': abs_pos + 1,
                    'name': 'T_site', 
                    'score': '.',
                    'strand': '-'
                })
    
    # Convert new sites to DataFrame
    new_sites_df = pd.DataFrame(new_sites)
    
    # Combine original and new sites
    combined_df = pd.concat([original_bed, new_sites_df], ignore_index=True)
    
    # Sort by chromosome and position
    combined_df = combined_df.sort_values(['chrom', 'start'])
    
    # Remove duplicates
    combined_df = combined_df.drop_duplicates(['chrom', 'start', 'end', 'strand'])
    
    return combined_df

# Run the function
fasta_file = "hg38.fa"  # Replace with your genome fasta file path
bed_file = "union.bed"
output_file = "combined_m6A_sites_union.bed"

combined_sites = get_sequence_and_sites(fasta_file, bed_file)
combined_sites.to_csv(output_file, sep='\t', header=False, index=False)
