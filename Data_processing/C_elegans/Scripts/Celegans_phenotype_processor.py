#!/usr/bin/env python3 

## import the required libraries 
import argparse
import pandas as pd 

def convert_pheno_files(pheno_file, sheet_name, output_file):
    """
    Process phenotype data into GeneNetwork (GN) required formatting. 

    Parameters:
    - pheno_file (str): Path to the genotype file (e.g., file_map.txt) 
    - sheet_name (str): Path to the marker file (e.g., file_marker.txt) 
    - output_file (str): Path to save the processed output file. 
    
    """

    if pheno_file.split(".")[1] == "xlsx":
        ## load in the excel sheet with the experimental phenotypes 
        phenotype_df = pd.read_excel(pheno_file, sheet_name=sheet_name)

        ## drop the columns not needed 
        cols_to_drop = ['Seq_order', 'Ril_ID', 'Seqence_ID', 'rils']
        phenotype_df.drop(columns=cols_to_drop, axis=1, inplace=True)

        ## add the SE and N cols after each of the available columns, then add string x as the values for both added cols 
        #Create a new DataFrame with the desired structure
        new_cols = [] 

        for cols in phenotype_df.columns: 
            new_cols.append(cols)
            new_cols.append(f"SE")
            new_cols.append(f"N") 

        #Create a new DataFrame with the same number of rows and fill SE and N with 'x'
        phenotype01_df = pd.DataFrame(columns=new_cols)

        for col in phenotype_df.columns:

            # Assign values from original DataFrame
            phenotype01_df[col] = phenotype_df[col]

            # Assign 'x' to SE and N columns
            phenotype01_df[f'SE'] = 'x'
            phenotype01_df[f'N'] = 'x'

        ## Transpose the df and maintain the indexing 
        pheno_df_T = phenotype01_df.set_index('WN_ID').T 
        pheno_df_T = pheno_df_T.reset_index()
        pheno_df_T.rename(columns={'index':'Strain'}, inplace=True)

        ## drop the first two rows 'SE' and 'N' as they are not needed for the col headers
        pheno_df_T = pheno_df_T.drop(index=[0,1])
        pheno_df_T.reset_index(drop=True, inplace=True)

        ## round the values to at least 6 decimal places 
        pheno_df_T = pheno_df_T.map(lambda x: round(x, 6) if isinstance(x,(int, float))else x)

        ## Save your file 
        pheno_df_T.to_csv(output_file, index=None, header=True, sep='\t', float_format="%.6f")
    else: 
        ## load in the excel sheet with the experimental phenotypes 
        phenotype_df = pd.read_table(pheno_file)

        ## rearrange the column headers 
        phenotype_df.reset_index(inplace=True)
        phenotype_df.rename(columns={"index":"WN_ID"}, inplace=True)

        ## restructure rownames to be short and informative 
        phenotype_df['WN_IDt'] = phenotype_df['WN_ID'].str.split(';').apply(lambda x: f"{x[0]}_{x[1]}")
        phenotype_df.drop('WN_ID', axis=1, inplace=True)
        phenotype_df.insert(0, 'WN_ID', phenotype_df['WN_IDt'])
        phenotype_df.drop('WN_IDt', axis=1, inplace=True)

        ## Insert N and SE rows under each of the traits for this file 
        ###Create new dataframe with rows to be inserted 
        new_rows = pd.DataFrame({
            phenotype_df.columns[0]: ['SE', 'N'],
            **{col: ['x', 'x'] for col in phenotype_df.columns[1:]}
        })

        ###combine the original dataframe with the new rows 
        pheno_df2 = pd.DataFrame()
        for i in range(len(phenotype_df)):
            pheno_df2 = pd.concat([pheno_df2, phenotype_df.iloc[[i]], new_rows])

        ## round the values to at least 6 decimal places 
        pheno_df2 = pheno_df2.map(lambda x: round(x, 6) if isinstance(x,(int, float))else x)
        ## Save your file 
        pheno_df2.to_csv(output_file, index=None, header=True, sep='\t', float_format="%.6f")

    # Exit message 
    print(f"Processing complete: File saved at: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        prog="Celegans-Phenotype-Processor", 
        description="Processes Celegans experimental phenotype files into GN format."
        )
    parser.add_argument("pheno_file", help="Path to the phenotype file (excel or txt)")
    parser.add_argument("--sheet_name", help="Name of the sheet with your data in excel (use this if you provide an excel input file)")
    parser.add_argument("output_file", help="Path to save your output file")

    args = parser.parse_args()
    convert_pheno_files(args.pheno_file, args.sheet_name, args.output_file)

if __name__ == "__main__":
    main()