#!/usr/bin/env python3

import pandas as pd
import numpy as np

# Define output path
output_path = "/home/felixl/felixl/Project2024/HS_rats/to_experiment/genotypes/processed/HS_rats-bed-genotypes.csv"

# Load the SNP metadata (.bim file)
try:
    bim_data = pd.read_csv("/home/felixl/felixl/Project2024/HS_rats/to_experiment/genotypes/raw/r10.3.2.bim", 
                           delim_whitespace=True, header=None, 
                           names=["Chromosome", "Marker_ID", "Genetic_Distance", "Position", "Allele1", "Allele2"])
    print("Loaded r10.3.2.bim successfully.")
except FileNotFoundError:
    print("Error: r10.3.2.bim not found. Check file path.")
    exit(1)

# Initialize output DataFrame with SNP metadata
output_data = pd.DataFrame()
output_data["Chromosome"] = bim_data["Chromosome"]
output_data["Genetic_Distance"] = bim_data["Genetic_Distance"]
output_data["Position"] = bim_data["Position"]
output_data["Marker_ID"] = bim_data["Marker_ID"]

# Load .raw file to get sample data structure (just header and first few rows)
raw_header = pd.read_csv("/home/felixl/felixl/Project2024/HS_rats/to_experiment/genotypes/processed/r10.3.2.raw", 
                         delim_whitespace=True, nrows=1)
sample_ids = raw_header["IID"].tolist()
total_snps = len(raw_header.columns) - 6  # Exclude FID, IID, PAT, MAT, SEX, PHENOTYPE
print("Found {} samples and {} SNPs.".format(len(sample_ids), total_snps))

# Process .raw file in SNP chunks
chunk_size = 10000  # Process 10,000 SNPs at a time (adjust based on memory)
for start_idx in range(0, total_snps, chunk_size):
    end_idx = min(start_idx + chunk_size, total_snps)
    print("Processing SNPs from index {} to {}".format(start_idx, end_idx - 1))
    
    # Load specific SNP columns for all samples
    cols_to_load = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + \
                   list(raw_header.columns[6:][start_idx:end_idx])
    chunk_data = pd.read_csv("/home/felixl/felixl/Project2024/HS_rats/to_experiment/genotypes/processed/r10.3.2.raw", 
                             delim_whitespace=True, usecols=cols_to_load)
    
    # Extract genotypes for this chunk
    genotypes = chunk_data.iloc[:, 6:]  # Skip metadata columns
    
    # Add genotypes to output_data
    for i, sample_id in enumerate(sample_ids):
        output_data["Sample_{}".format(sample_id)] = genotypes.iloc[:, i]
    
    # Save incrementally
    if start_idx == 0:
        output_data.to_csv(output_path, index=False)
    else:
        output_data.to_csv(output_path, mode='a', header=False, index=False)

print("Genotype data with metadata saved to {}".format(output_path))