#!/usr/bin/env python3 

#import packages 
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from scipy.stats import chi2_contingency

# Load data
data = pd.read_csv('../data/Celegans_genotyped_testing.csv', low_memory=False)

# Extract genotype columns
genotype_cols = [col for col in data.columns if col.startswith('WN')]
genotypes = data[genotype_cols]

# 1. Quality Control
# Remove markers with high missing rate (>10%)
missing_rate = genotypes.isna().mean(axis=1)
data = data[missing_rate < 0.1]

# Calculate MAF and remove low MAF markers (<0.01)
def calculate_maf(row):
    alleles = row.value_counts(normalize=True)
    freqs = [alleles.get(x, 0) for x in [0.0, 0.5, 1.0]]
    p = (2 * freqs[0] + freqs[1]) / (2 * sum(freqs))
    q = 1 - p
    return min(p, q)

mafs = genotypes.apply(calculate_maf, axis=1)
data = data[mafs >= 0.005]

# Remove monomorphic markers
std = genotypes.std(axis=1)
data = data[std > 0]

# 2. Imputation
imputer = SimpleImputer(strategy='most_frequent')
genotypes_imputed = pd.DataFrame(imputer.fit_transform(genotypes), columns=genotype_cols)
data[genotype_cols] = genotypes_imputed

# 3. LD Pruning (simplified correlation-based)
def ld_pruning(genotypes, threshold=0.8, window=100):
    keep = []
    for i in range(0, len(genotypes), window):
        window_data = genotypes.iloc[i:i+window]
        corr = window_data.corr().abs()
        np.fill_diagonal(corr.values, 0)
        # Greedily select markers with low correlation
        selected = []
        for j in range(len(window_data)):
            if all(corr.iloc[j, k] < threshold for k in selected):
                selected.append(j)
        keep.extend(window_data.index[selected])
    return keep

keep_indices = ld_pruning(genotypes_imputed)
data = data.loc[keep_indices]

# 4. Save cleaned data
data.to_csv('../data/processed/plink_like_results/cleaned_celegans_genotype-r0.8_w100.csv', index=False)

