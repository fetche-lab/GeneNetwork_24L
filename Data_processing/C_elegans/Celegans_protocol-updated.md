## Caenorhabditis elegans (*C.elegans*) Data Processing and Entry Manual 

### This document contains many of the basic steps and hacks on how to process C.elegans datasets, and prepare them in the format that is fit to be uploaded into GeneNetwork.  

### In general, the following are file types GeneNetwork expects;  

01. Genotype file 
  - Raw file: Usually in `Variant Call Format (VCF)` file format, also could be customary to the data provider  
  - Processed file: in `*.geno and *BIMBAM.txt` format 

02. Classical Phenotype file 
  - Contain experimental phenotypes
  - Usually containing measurements on factors like; blood pressure, weight, height, growth rates, etc.  

03. Log2/Expression files (Non classical phenotypes) 
  - Contain log2 normalized gene expression values 
  - (Note yet applicable in the C.elegans case)

04. Annotation/Description files  
  - Contain descriptive information about the rest of the values in files mentioned above

### The following are general steps that can be used to process the files mentioned above to give out desired outputs for GN2 production server 

### 01 Processing genotype files 
- The script to process genotypes is written in Python  
- What is expected on the genotype input file(s); 
  - 2 files from which, genotype information is extracted 
    - `*_map.txt`
     - A simple matrix with individual names as column headers and first column containing marker ids
     - An example file: [Snoek_map.txt](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/example_files/Snoek_map.txt){:target="_blank"}
    - `*_marker.txt`
     - 3 columns of interest; 
       - column 01 to have marker ids 
       - column 02 to have chromosomes 
       - column 03 to have position values
     - An example file: [Snoek_marker.txt](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/example_files/Snoek_marker.txt){:target="_blank"}

- Before running the script below, make sure your input files adhere to the formatting explained above, with corresponding examples for each as per the links provided 

- NB; the `*_marker.txt` file usually has disarranged columns. Check and compare with the example provided above. If not in the expected order, you can use excel to create a new sheet with the correct format, then use the following script to convert it to the txt format as follows; 

```python 
#!/usr/bin/env python3 

## import packages 
import argparse
import pandas as pd

## the program 
def marker_cols_rearrange(marker_file, sheet_name, output_file):
    """ 
    Checks and rearrange the columns in marker file into a good order 
    
    Parameters: 
    - Marker_file(str): Path to the marker file (e.g marker_file.xlsx)
    - sheet_name(str): The name of your sheet in excel (e.g sheet01) 
    - Output_file(str): Path to save the processed output file  
    """
    ## load in the file 
    marker_df = pd.read_excel(marker_file, sheet_name=sheet_name)

    ## drop the first column 
    del marker_df["name"]

    ## rearrange the columns 
    marker_df1 = marker_df.set_axis(['name', 'Chromosome', 'Pos_WS258'], axis=1)

    ## save the file in txt format 
    marker_df1.to_csv(output_file, header = True, index = None, sep="\t") 

    # Exit message 
    print(f"Processing complete: File saved at: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Rearranges columns in marker file in the right order")
    parser.add_argument("marker_file", help="Path to the marker file {it should be in Excel format}")
    parser.add_argument("sheet_name", help="Name of the sheet with your data in excel")
    parser.add_argument("output_file", help="Path to save your output file")

    args = parser.parse_args()
    marker_cols_rearrange(args.marker_file, args.sheet_name, args.output_file)

if __name__ == "__main__":
    main()

```
- N.B; make sure `openpyxl` and `pandas` are installed in your system. 
    - you can use `pip install xyz` to install them, where `xyz` represents the package(s) to be installed 

- your output file should have a `.txt` extension 

- Then you can run the following script to process the genotype file for GN2 
 
```python 
#!/usr/bin/env python3
#
# import required modules 
import argparse
import pandas as pd 

def convert_geno_files(geno_file, marker_file, output_file):
    """
    Converts genotype and marker files into GEMMA-compatible format. 

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
    parser = argparse.ArgumentParser(description="Convert genotype data to GN2/GEMMA-compatible format")
    parser.add_argument("geno_file", help="Path to the genotype file.")
    parser.add_argument("marker_file", help="Path to the marker file.")
    parser.add_argument("output_file", help="Path to the processed output file.")

    args = parser.parse_args() 
    convert_geno_files(args.geno_file, args.marker_file, args.output_file)

if __name__ == "__main__":
    main() 

```
- NB; 
  - The script above assumes the genotype values in the raw file to be (-1, 1, 0.5)

- HOW TO USE IT: 
  01. Download the script from this link: [celegans_geno_processor](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/Scripts/Celegans_genotype_processor.py)
  ```sh 
  git clone "paste here the link to the script above and run" 
  ```
  02. To run the script, either: 
     - use python directly;
     ```sh 
     python3 Celegans_geno_processor.py geno_file marker_file output_file 
     ```
     - make it executable, then run it 
     ```sh
     #make it executable 
     chmod +x Celegans_geno_processor.py

     #run the script  
     ./Celegans_geno_processor.py geno_file marker_file output_file
     ``` 

### 02 Processing classical phenotypes 

- In this case, we will use python/pandas to process the classical phenotype file 

```python 
## import the necessary libraries 
import pandas as pd 
import numpy as np 

## load in the phenotype file into a dataframe 
celegans_pheno_df = pd.read_table('Snoek_2019_data.txt')

## inspect the file {use .head(), .columns, .shape}

## convert the index section into a column to be used as the rownames 
celegans_pheno_df.reset_index(inplace=True) #reset the index 
celegans_pheno_df.rename(columns={'index':'IDs'}, inplace=True) #rename the index header name 
celegans_pheno_df['IDs'].head() #inspect the new column 

## Extract the desirable part of the rownames 
celegans_pheno_df['Ids'] = celegans_pheno_df['IDs'].str.split(';').str[1] #take the desirable part of the col values, create a new col 
celegans_pheno_df.drop('IDs', axis=1, inplace=True) #drop the original rownames
celegans_pheno_df.insert(0, 'IDs', celegans_pheno_df['Ids']) #put the new col as the first col/rownames
celegans_pheno_df.drop('Ids', axis = 1, inplace=True) #drop the replica of the rownames at the end of the df

##Inspect your changes {use .head(), .columns, .shape}

## Perform log2 transformation on the df 
celegans_pheno_df1 = celegans_pheno_df.apply(pd.to_numeric, errors='coerce')#this helps prevent errors due to non numeric values vs log2 
log2_celegans_df = celegans_pheno_df1.apply(np.log2)#the actual log2 transformation 
log2_celegans_df.insert(0, 'Ids', celegans_pheno_df['IDs'])#insert the original IDs column to replace the transformed one (which has NaN)
log2_celegans_df.drop('IDs', axis=1, inplace=True)#drop the transformed IDs column 

## round the values to 6 decimal places 
log2_celegans_df1 = log2_celegans_df.round(6)

## fill the NaN with NA 
log2_celegans_df1.fillna('NA', inplace=True)

## Save your dataframe into a file 
log2_celegans_df1.to_csv('C_elegans_Snoek_log2_data01.tsv', index=None, header=True, sep='\t', float_format='%.6f')

```
### 03 Processing Experimental Phenotypes 
- On this part, we expect to have two files. One for the phenotypes and one for the corresponding descriptions

- processing the experimental phenotypes 

```python 
## import the required libraries 
import pandas as pd 

## load in the excel sheet with the experimental phenotypes 
phenotype_df = pd.read_excel('./Supplementary_materials/7843112/C-elegans_sup.xlsx', sheet_name='SupTab5')

## inspect the file {use .head(), .columns, .shape}

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

## test the new df and inspect it 
phenotype01_df.columns
phenotype01_df.head() 

## Transpose the df and maintain the indexing 
pheno_df_T = phenotype01_df.set_index('WN_ID').T 
pheno_df_T = pheno_df_T.reset_index()
pheno_df_T.rename(columns={'index':'Strain'}, inplace=True)

## drop the first two rows 'SE' and 'N' as they are not needed for the col headers
pheno_df_T = pheno_df_T.drop(index=[0,1])
pheno_df_T.reset_index(drop=True, inplace=True)

## Save your file 
pheno_df_T.to_csv('Snoek_Celegans_Av-phenotypes_updated01.tsv', index=None, header=True, sep='\t')

```

### 04 Processing Experimental Phenotypes Descriptions 
- This category of files contain description on the experimental phenotypes used in the study 
- Before generating the final descriptive file, it's important to consider the following; 
 - Necessary columns needed for the final descriptive file. In this case; 

   ```python 
    cols = ['Pubmed ID','Pre Publication Description','Post Publication Description','Original Description','Pre Publication Abbreviation','Post Publication Abbreviation','Lab Code','Submitter','Owner','Authorized Users','Authors','Title','Abstract','Journal','Volume','Pages','Month','Year','Units']

   ```
  - The availability of files/columns that fit the needed column structure above, and how to get them especially if they are found from different sources. 
  - This step is highly subjective, therefore it is important to consider the column names needed, as well as example files from gn2, or gn2 experts as guidance to make sure one is successful in generating the files. 
  - For now, this part is still under development, until we manage to clearly describe every bit of the columns necessary for this part.  

### 05 Strains and StrainXref 
- `TODO` 










