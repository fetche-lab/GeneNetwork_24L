#!/usr/bin/env python3

## import necessary packages 
import pandas as pd
import numpy as np
import os 

# Load genotype data
genotype_file = "Celegans_genotyped_testing.csv"
genotype_data = pd.read_csv(genotype_file)

# Load SNP metadata
metadata_file = "SupTab2_Celegans.csv"
metadata_data = pd.read_csv(metadata_file)

# Merge genotype data with metadata on 'Locus' or equivalent column
merged_data = pd.merge(genotype_data, metadata_data, on="Locus", how="inner")

# Define smoothing function
def smooth_genotypes(row, window_size=3):
    """
    Smooth genotypes by considering neighboring SNPs within a window.
    """
    smoothed_row = row.copy()
    for i in range(len(row)):
        # Get window around current SNP
        start = max(0, i - window_size)
        end = min(len(row), i + window_size + 1)
        window = row[start:end]
        
        # Apply majority rule for smoothing
        if row[i] == 0.5:  # Heterozygous case
            counts = np.bincount(window.astype(int))
            smoothed_row[i] = np.argmax(counts) if counts.any() else row[i]
    return smoothed_row

# Apply smoothing to genotype columns
genotype_columns = [col for col in genotype_data.columns if col.startswith("WN")]
for col in genotype_columns:
    merged_data[col] = smooth_genotypes(merged_data[col].values)

# Save smoothed data to a new file
output_file = "Smoothed_Celegans_Genotypes.csv"
merged_data.to_csv(output_file, index=False)

print(f"Genotype smoothing completed. Smoothed data saved to {output_file}.")

