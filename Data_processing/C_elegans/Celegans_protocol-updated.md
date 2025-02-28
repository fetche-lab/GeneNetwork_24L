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

- NB; the `*_marker.txt` file usually has disarranged columns ([see example](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/example_files/Snoek_2019_marker01.xlsx)). Check and compare with the example provided above. If not in the expected order, you can use excel to create a new sheet with the correct format, then use the following script to convert it to the txt format as follows; 

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
    parser = argparse.ArgumentParser(
      prog="Marker-cols",
      description="Rearranges columns in marker file in the right order"
      )
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

- HOW TO USE IT: 
  01. Download the script from this link: [marker_cols](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/Scripts/marker_cols.py)
  ```sh 
  git clone "paste here the link to the script above and run" 
  ```
  02. To run the script, either: 
     - use python directly;
     ```sh 
     python3 marker_cols.py -h ## more information, how to use
     python3 marker_cols.py marker_file sheet_name output_file
     ```
     - make it executable, then run it 
     ```sh
     #make it executable 
     chmod +x marker_cols.py

     #run the script  
     ./marker_cols.py marker_file sheet_name output_file

     ``` 
  - Here's an example of already processed marker file: [marker_file](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/example_files/processed/Snoek_2019_markerUpdated.txt)

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
    parser = argparse.ArgumentParser(
      prog="Celegans-geno-processor",
      description="Convert genotype data to GN2/GEMMA-compatible format")
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
     python3 Celegans_genotype_processor.py -h ## more information, how to use
     python3 Celegans_genotype_processor.py geno_file marker_file output_file 
     ```
     - make it executable, then run it 
     ```sh
     #make it executable 
     chmod +x Celegans_genotype_processor.py

     #run the script  
     ./Celegans_genotype_processor.py geno_file marker_file output_file
     ``` 

### 02 Processing Experimental Phenotypes 
- On this part, we expect to have two files. One for the phenotypes and one for the corresponding descriptions

- processing the experimental phenotypes 
  - The could be provided in the supplementary files (as in Snoek), but sometimes, it is provided as a stand alone file 
  - At least the file should look like this: [phenotype_file]()

```python 
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

```

- HOW TO USE IT: 
  01. Download the script from this link: [celegans_pheno_processor](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/Scripts/Celegans_phenotype_processor.py)
  ```sh 
  git clone "paste here the link to the script above and run" 
  ```
  02. To run the script, either: 
     - use python directly;
     ```sh 
     python3 Celegans_phenotype_processor.py -h ## more information, how to use
     python3 Celegans_phenotype_processor.py pheno_file  output_file --sheet_name SHEET_NAME ## {run the --sheet_name option if your input file has .xlsx extension}

     ```
     - make it executable, then run it 
     ```sh
     #make it executable 
     chmod +x Celegans_phenotype_processor.py

     #run the tool 
     ./Celegans_phenotype_processor.py pheno_file output_file --sheet_name SHEET_NAME 

     ``` 
  - Example of phenotype raw file:[phenotype raw file](https://github.com/fetche-lab/GeneNetwork_24L/tree/main/Data_processing/C_elegans/example_files/Phenotype_Testfile.txt)
  - Examples of phenotype output files: [phenotype processed files](https://github.com/fetche-lab/GeneNetwork_24L/tree/main/Data_processing/C_elegans/example_files/processed/Test_pheno01.txt)

  
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










