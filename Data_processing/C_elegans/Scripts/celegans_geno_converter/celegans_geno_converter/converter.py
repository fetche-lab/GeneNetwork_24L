import pandas as pd 
import os 

def convert_genotype_files(geno_file, marker_file, output_file): 
    """
    
    Converts genotype and marker files into GEMMA-compatible format. 

    Parameters:
    - geno_file (str): Path to the genotype file (e.g., file_map.txt) 
    - marker_file (str): Path to the marker file (e.g., file_marker.txt) 
    - output_file (str): Path to save the processed output file. 

""" 

    # check if input files exist 
    if not os.path.exists(geno_file):
        print(f"Error: File not found: {geno_file}") 
        return 
    if not os.path.exists(marker_file):
        print(f"Error: File not found: {marker_file}") 
        return 

    try: 
        # Load genotype and marker files 
        geno_df = pd.read_table(geno_file) 
        marker_df = pd.read_table(geno_file) 

        # Encode genotype values to GEMMA format (-1 -> 'A', 1 -> 'B', 0.5 -> 'H')
        ## Validate expected values in genotype data 
        expected_values = {-1, 1, 0.5} 
        actual_values = set(geno_df.values.flatten()) # Get all the unique values in the dataframe 

        ## Check if all expected values are present 
        if not expected_values.issubset(actual_values): 
            missing_values = expected_values - actual_values 
            print(f"Warning: The following expected values are missing from the genotype data: {missing_values}") 
            print("This may indicate an issue with the input file format or data preprocessing") 
            ## Apply conditional rounding: 
            ## Values between 1 and 0 should be rounded to 0.5 
            geno_df = geno_df.applymap(lambda x: 0.5 if 0 < x < 1 else x) 

            # Proceed with encoding 
            geno_df.replace({-1: 'A', 1: 'B', 0.5: 'H'}, inplace=True)

            # Fill missing values with "NA" 
            geno_df.fillna("NA", inplace=True) 

        # Add cM and Mb columns to marker_df 
        marker_df['Mb'] = marker_df.iloc[:, 2] / 10**6
        marker_df['cM'] = marker_df.iloc[:, 2] / 10**6  # Placeholder calculation

        # Rename and rearrange the dataframe columns 
        # Renaming columns by index
        marker_df.columns.values[0] = 'Locus'  # Rename the first column (name -> Locus)
        marker_df.columns.values[1] = 'Chr'    # Rename the second column (Chromosome -> Chr)

        # Dropping the column by index
        marker_df.drop(marker_df.columns[2], axis=1, inplace=True)  # Drop the 3rd column (Pos_WS258)

        # Reordering columns using indexing
        marker_df = marker_df[[marker_df.columns[1], marker_df.columns[0], 'cM', 'Mb']]  

        # Merge marker dataframe with genotype dataframe
        geno_df1 = geno_df.reset_index(drop=True)
        merged_df = pd.concat([marker_df, geno_df1], axis=1)

        # Save output
        os.makedirs(os.path.dirname(output_file), exist_ok=True)  # Ensure output directory exists
        merged_df.to_csv(output_file, index=False, header=True, sep='\t')

        print(f"Processing complete: File saved at: {output_file}") 

    except Exception as e: 
        print(f"Error: {e}") 


     
