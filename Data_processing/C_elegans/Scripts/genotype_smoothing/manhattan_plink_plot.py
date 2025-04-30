#!/usr/bin/env python3 

#import packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the GWAS results
data = pd.read_csv("celegans_newPheno-lmm.assoc.txt", sep="\s+")

# Remove any rows with missing p-values
data = data.dropna(subset=['p_wald'])

# Calculate -log10(p_wald) for plotting
data['neg_log_p'] = -np.log10(data['p_wald'])

# Define chromosomes and sort them (1, 2, 3, 4, 5, X)
chromosomes = sorted(data['chr'].unique(), key=lambda x: int(x) if x != 'X' else 6)
chrom_sizes = data.groupby('chr')['ps'].max().reindex(chromosomes).fillna(0)

# Assign cumulative positions for plotting
cum_pos = []
ticks = []
current_pos = 0
for chrom in chromosomes:
    chrom_data = data[data['chr'] == chrom]
    cum_pos.extend(chrom_data['ps'] + current_pos)
    ticks.append(current_pos + chrom_sizes[chrom] / 2)
    current_pos += chrom_sizes[chrom]

data['cum_pos'] = cum_pos

# Set up colors for alternating chromosomes
colors = ['#1f77b4', '#ff7f0e']  # Blue and orange
chrom_colors = {chrom: colors[i % 2] for i, chrom in enumerate(chromosomes)}

# Create the Manhattan plot
plt.figure(figsize=(12, 6))
for chrom in chromosomes:
    chrom_data = data[data['chr'] == chrom]
    plt.scatter(chrom_data['cum_pos'], chrom_data['neg_log_p'], 
                c=chrom_colors[chrom], s=20, alpha=0.6, label=chrom if chrom == chromosomes[0] else "")

# Add significance threshold line (p_wald < 5e-8)
#threshold = -np.log10(5e-8)
threshold = 2 
plt.axhline(y=threshold, color='red', linestyle='--', linewidth=1, label='p = 0.01')

# Customize plot
plt.xlabel('Chromosome')
plt.ylabel('-log10(p-value)')
plt.title('Manhattan Plot of GWAS Results (C. elegans)')
plt.xticks(ticks, chromosomes)
plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig('./plots/manhattan_plotnewPhenoRaw.png', dpi=300)
plt.close()
