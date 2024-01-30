# Using Python to convert datasets into the R/qtl2 format 
## In this context, the dataset used was the one shared by Harm. The idea behind was to generate the genotype file(s), phenotype file(s), genetic/physical mapping file(s), all in the csv format, as well as the control file in yaml.
## The tabular datasets selected were; 
### table07 (genetic mapping file) 
### table08 (genotype file)
### table09 (phenotype file)

## Genetic mapping file 
```python
import pandas as pd  
## upload the datasets into dataframes 
# table07 
table07 = pd.read_excel('../../Projects_2024/data_engineering/data/R_qtl2/data/Harms_sup_data/supplementary_table_7-9/table-s7.xlsx')

## rename start col to position (cM)
table07.rename(columns = {'start':'position(cM)'}, inplace = True) 

## save table07 as a csv file 
table07.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/gmap_Harm_arabid.csv', index = None, header = True) 

```

## Column renaming function 
```python
## create a function to help rename strain column names with varied extensions corresponding to the type of experiments 
## create a function to rename the columns 
import pandas as pd

def rename_columns(dataframe, old_to_new_mapping):
    """
    Rename columns in a DataFrame based on the provided mapping.

    Parameters:
    - dataframe: The DataFrame to modify.
    - old_to_new_mapping: A dictionary where keys are old column names, and values are new column names.

    Returns:
    - Modified DataFrame.
    """
    dataframe.rename(columns=old_to_new_mapping, inplace=True)
    return dataframe


```

## Genotype file 
```python
import pandas as pd
## upload the dataset into dataframes 
# table08 
table08 = pd.read_excel('../../Projects_2024/data_engineering/data/R_qtl2/data/Harms_sup_data/supplementary_table_7-9/table-s8.xlsx')
# create a list of old names 
geno_col_old = [col for col in table08.columns]

# create a list of new names 
geno_col = [col.split('_') for col in table08.columns]
geno_col_new= [col[0] for col in geno_col]

## create a dictionary with old names as keys and new names as values 
geno_dic = dict(zip(geno_col_old, geno_col_new))

# Use the function to rename columns
geno_df_new = rename_columns(table08, geno_dic)
 
# save the resulting dataframe into a csv file 
geno_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/genotype_Harm_arabid.csv', index = None, header = True)


```

## Phenotype file 
### im experiment 
```python
#im experiment 
im_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_im_exp.tsv', sep = '\t') 

# create a list of old names 
im_col_old = [col for col in im_df.columns if '_im' in col]

# create a list of new names 
im_col = [col.split('_') for col in im_df.columns if '_im' in col]
im_col1_new= [col[0] for col in im_col]

## create a dictionary with old names as keys and new names as values 
im_dic = dict(zip(im_col_old, im_col1_new))

# Use the function to rename columns
im_df_new = rename_columns(im_df, im_dic)

# save the resulting dataframe into a csv file 
im_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_im_Harm_arabid.csv', index = None, header = True)


```
### pd experiment 
```python
#pd experiment 
pd_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_pd_exp.tsv', sep = '\t') 

# create a list of old names 
pd_col_old = [col for col in pd_df.columns if '_pd' in col]

# create a list of new names 
pd_col = [col.split('_') for col in pd_df.columns if '_pd' in col]
pd_col1_new= [col[0] for col in pd_col]

## create a dictionary with old names as keys and new names as values 
pd_dic = dict(zip(pd_col_old, pd_col1_new))

# Use the function to rename columns
pd_df_new = rename_columns(pd_df, pd_dic)

# save the resulting dataframe into a csv file 
pd_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_pd_Harm_arabid.csv', index = None, header = True)

```
### rp experiment 
```python
#rp experiment 
rp_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_rp_exp.tsv', sep = '\t') 

# create a list of old names 
rp_col_old = [col for col in rp_df.columns if '_rp' in col]

# create a list of new names 
rp_col = [col.split('_') for col in rp_df.columns if '_rp' in col]
rp_col1_new= [col[0] for col in rp_col]

## create a dictionary with old names as keys and new names as values 
rp_dic = dict(zip(rp_col_old, rp_col1_new))

# Use the function to rename columns
rp_df_new = rename_columns(rp_df, rp_dic)

# save the resulting dataframe into a csv file 
rp_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_rp_Harm_arabid.csv', index = None, header = True)

```
## Control file 
```markdown
# Data from Hartanto et al., 2020 G3 Genes|Genomes|Genetics, Volume 10, Issue 11, 1 November 2020, Pages 4215–4226
# Abstract of the paper at G3 Genes|Genomes|Genetics, https://doi.org/10.1534/g3.120.401477
# Available at https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Hartanto_et_al_2020/12844358

crosstype: RIL
geno: genotype_Harm_arabid.csv
pheno:
    - phenotype_im_Harm_arabid.csv 
    - phenotype_pd_Harm_arabid.csv 
    - phenotype_rp_Harm_arabid.csv 
gmap: gmap_Harm_arabid.csv

genotypes:
  AA: -1 
  BB: 1 
na.strings:
- '-'
- NA

```
