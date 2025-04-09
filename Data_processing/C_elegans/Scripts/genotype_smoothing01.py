#!/usr/bin/env python3

## import packages 
import pandas as pd

## define function to smoothen the genotype 
def smoothen_genotype(input_file, output_file):
    # Read input file
    df = pd.read_csv(input_file)
    df = df.reset_index()
    # Group by columns starting from the 5th (0-based index 4)
    group_cols = df.columns[5:].tolist()  # Columns 5+ contain genotype patterns
    df['group_id'] = df.groupby(group_cols).ngroup()

    # Select rows based on grouping criteria
    keep_indices = []
    for _, group in df.groupby('group_id'):
        size = len(group)
        if size > 2:
            keep_indices.extend([group.iloc[0]['index'], group.iloc[-1]['index']])
        elif size == 2:
            keep_indices.append(group.iloc[0]['index'])
        else:
            keep_indices.append(group.iloc[0]['index'])

    # Create smoothed dataset and maintain original order
    smoothed_df = (
        df[df['index'].isin(keep_indices)]
        .sort_values('index')
        .drop(columns=['index', 'group_id'])
    )
    
    # save the file 
    smoothed_df.to_csv(output_file, index=False, header=True)
    print(f"Smoothed file saved to {output_file}")

# Usage
smoothen_genotype('/home/fetche-lab/Desktop/gn_remote/C_elegans/Data/Raw/Celegans_dataset/Dataset2/Celegans_Tables-GN2/genotypes/Celegans_genotyped_testing.csv', '/home/fetche-lab/Desktop/gn_remote/C_elegans/Data/Raw/Celegans_dataset/Dataset2/Celegans_Tables-GN2/genotypes/smoothed_test1.csv')

