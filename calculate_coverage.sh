#!/bin/bash

# Input files
SAM_FILE="ALKBH5_reads.sam"
#FIXED_BED_FILE="extracted_weighted_WT1.bed"
FIXED_BED_FILE="extracted_xgb_0.7_WT.bed"
BAM_FILE="ALKBH5_reads.bam"
SORTED_BAM_FILE="ALKBH5.sorted.bam"
COVERAGE_FILE="ALKBH5_coverage_xgb_0.7_WT.txt"

# Convert SAM to BAM
samtools view -S -b $SAM_FILE > $BAM_FILE

# Sort BAM file
samtools sort $BAM_FILE -o $SORTED_BAM_FILE

# Index the sorted BAM file
samtools index $SORTED_BAM_FILE

# Fix the BED file by adding placeholder columns
#awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", "0", $4}' $BED_FILE > $FIXED_BED_FILE

# Calculate strand-specific coverage
bedtools coverage -s -a $FIXED_BED_FILE -b $SORTED_BAM_FILE > $COVERAGE_FILE

echo "Strand-specific coverage written to $COVERAGE_FILE"

