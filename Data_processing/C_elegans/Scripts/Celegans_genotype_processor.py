#!/usr/bin/env python3
# 
# import required modules 
import argparse
import pandas as pd 

def convert_geno_files(geno_file, marker_file, output_file):
    """
    Converts phenotype and marker files into GEMMA-compatible format. 

    Parameters:
    - geno_file (str): Path to the genotype file (e.g., file_map.txt) 
    - marker_file (str): Path to the marker file (e.g., file_marker.txt) 
    - output_file (str): Path to save the processed output file. 
    
    """
    # load in geno_file and marker_file 
    geno_df = pd.read_table(geno_file)
    marker_df = pd.read_table(marker_file)

    # Encode values in geno_file into gn2/gemmma-compatible format 
    for values in geno_df.columns: 
        geno_df[values].replace({-1:'A', 1:'B', 0.5:'H'}, inplace=True)

    # Fill empty cells in geno_file with `NA` string 
    geno_df.fillna('NA', inplace=True)

    # Create and modify columns (chr, locus, cM, Mb) from the marker dataframe 
    # Add cM and Mb columns to marker_df 
    marker_df['Mb'] = marker_df.iloc[:, 2] / 10**6
    marker_df['cM'] = marker_df.iloc[:, 2] / 10**6 

    # Rename and rearrange the dataframe columns 
    # Renaming columns by index
    marker_df.columns.values[0] = 'Locus'  # Rename the first column (name -> Locus)
    marker_df.columns.values[1] = 'Chr'   # Rename the second column (Chromosome -> Chr)

    # Dropping the column by index
    marker_df.drop(marker_df.columns[2], axis=1, inplace=True)  # Drop the 3rd column (Pos_WS258)

    # Reordering columns using indexing
    marker_df = marker_df[[marker_df.columns[1], marker_df.columns[0], 'cM', 'Mb']]  

    # Merge marker dataframe with genotype dataframe
    geno_df1 = geno_df.reset_index(drop=True)
    merged_df = pd.concat([marker_df, geno_df1], axis=1)

    # Save the output 
    merged_df.to_csv(output_file, index=False, header=True, sep='\t')

    # Exit message 
    print(f"Processing complete: File saved at: {output_file}")
def main(): 
    parser = argparse.ArgumentParser(
        prog="Celegans-geno-processor",
        description="Convert genotype data to GN2/GEMMA-compatible format"
        )
    parser.add_argument("geno_file", help="Path to the genotype file.")
    parser.add_argument("marker_file", help="Path to the marker file.")
    parser.add_argument("output_file", help="Path to the processed output file.")

    args = parser.parse_args() 
    convert_geno_files(args.geno_file, args.marker_file, args.output_file)

if __name__ == "__main__":
    main() 