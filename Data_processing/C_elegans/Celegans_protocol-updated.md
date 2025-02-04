## C elegans Data Processing and Entry Manual 

### This document contains many of the basic steps and hacks on how to process Celegans datasets, and prepare them in the format that is fit to be uploaded into the GeneNetwork2(GN2) production server.  

### In general, the following are file types GeneNetwork expects;  

01. Genotype file 
  - Raw file: Usually in `Variant Call Format (VCF)` file format, also could be customary to the data provider  
  - Processed file: in `*.geno and *BIMBAM.txt` format 

02. Classical Phenotype file 
  - Contain experimental phenotypes
  - Usually containing measurements on factors like; blood pressure, weight, height, growth rates, etc.  

03. Log2/Expression files (Non classical phenotypes) 
  - Contain log2 normalized gene expression values 

04. Annotation/Description files  
  - Contain descriptive information about the rest of the values in files mentioned above

### The following are general steps that can be used to process the files mentioned above to give out desired outputs for GN2 production server 

### 01 Processing genotype files 
- The script to process genotypes is written in Python  
- What is expected on the genotype input file(s); 
  - 2 files from which, genotype information is extracted 
    - `\*_map.txt`
     - A simple matrix with individual names as column headers and first column containing marker ids
     - An example file: [Snoek_map.txt](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/C_elegans/example_files/Snoek_map.txt){:target="_blank"}
    - `\*_marker.txt`
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
import pandas as pd

## load in the file 
marker_df = pd.read_excel("./your_markerfile.xlsx", sheet_name="your_sheet_name")

## save the file in txt format 
marker_df.to_csv("./your_processed_markerfile.txt", header = True, index = None, sep="\t") 

```
- N.B; make sure `openpyxl` and `pandas` are installed in your system. 
    - you can use `pip install xyz` to install them, where `xyz` represents the package(s) to be installed 

- Then you can run the following script to process the genotype file for GN2 
 
```python 

## import required python libraries 
import os 
import pandas as pd 

## load in the file(s) into a dataframe 
geno_df = pd.read_table('./file_map.txt')
marker_df = pd.read_table('./file_marker.txt')

## check for unique values 
unique_values = []
for values in geno_df.columns[1:]:
    unique_values.append(geno_df[values].unique())

#print(unique_values) ## 1, -1 are the only unique values 

## Encode the values into the gemma gn2 accepted format 
### -1 for A (females?) and 1 for B (males?)
for values in geno_df.columns:
    geno_df[values].replace({-1:'A', 1:'B', 0.5:'H'}, inplace=True)

## fill null values with NA 
geno_df.fillna('NA', inplace=True)

## add columns {chr, locus, cM, Mb} from the marker dataframe 
### calculate the cM and Mb 
marker_df['Mb'] = marker_df['Pos_WS258']/10**6
marker_df['cM'] = marker_df['Pos_WS258']/10**6 

### add the columns to the genotype dataframe 
marker_df.rename(columns={'Chromosome':'Chr', 'name':'Locus'}, inplace=True) 
marker_df.drop('Pos_WS258', axis=1, inplace=True)
marker_df = marker_df[['Chr', 'Locus', 'cM', 'Mb']] #rearrange the columns 
geno_df1 = geno_df.reset_index() #convert rownames to column index
geno_df1.drop('index', axis=1, inplace=True) #drop the index column
merged_df = pd.concat([marker_df, geno_df1], axis=1) #add cols in marker df into genotype df 

## save the dataframe into a file 
merged_df.to_csv('./processed/C_elegans-Genotypes.tsv', index = None, header = True, sep='\t')

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










