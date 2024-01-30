'''  
Program name: Arabidospsis_rqtl2_bundle.py 
Author: Felix Lisso 

Description: This script converts the selected tabular data shared by Harm in GeneNetwork data upload team, 
             to rqtl2 format and gets to be uploaded to GeneNetwork web service databases 

Input data: 
    - table07 as the genetic mapping file  
    - table08 as the genotype file 
    - table09 as the phenotype file (in three separate subsets) 
    
Output data: 
    - a zipped file containing three csv files and one control file in yaml format 

'''

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

## upload the datasets into dataframes 
# table07 
table07 = pd.read_excel('../../Projects_2024/data_engineering/data/R_qtl2/data/Harms_sup_data/supplementary_table_7-9/table-s7.xlsx')

#table08 
table08 = pd.read_excel('../../Projects_2024/data_engineering/data/R_qtl2/data/Harms_sup_data/supplementary_table_7-9/table-s8.xlsx')

#table09 
## already being processed into three separate datasets, more information on the readme file shared with this script 
#im experiment 
im_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_im_exp.tsv', sep = '\t') 

#pd experiment 
pd_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_pd_exp.tsv', sep = '\t') 

#rp experiment 
rp_df = pd.read_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/Harm_rp_exp.tsv', sep = '\t') 

## processing genetic mapping file (table07)

## rename start col to position (cM)
table07.rename(columns = {'start':'position(cM)'}, inplace = True) 

## save table07 as a csv file 
table07.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/gmap_Harm_arabid.csv', index = None, header = True) 

## processing the genotype file (table08) 
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


## processing the phenotype file(s)

# modify the columnnames to only remain with the actual strain names 
# create a list of old names 
rp_col_old = [col for col in rp_df.columns if '_rp' in col] #replace the _rp to _im for im experiment, and _pd for pd experiment 

# create a list of new names 
rp_col = [col.split('_') for col in rp_df.columns if '_rp' in col]
rp_col1_new= [col[0] for col in rp_col]

## create a dictionary with old names as keys and new names as values 
rp_dic = dict(zip(rp_col_old, rp_col1_new))

# Use the function to rename columns
rp_df_new = rename_columns(rp_df, rp_dic)

## the above section of code works for all three subsets of phenotype datasets, (im, pd, rp) 

# save the resulting dataframe into a csv file 
#rp
rp_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_rp_Harm_arabid.csv', index = None, header = True)
#pd
pd_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_pd_Harm_arabid.csv', index = None, header = True)
#im
im_df_new.to_csv('./data/R_qtl2/data/Harms_sup_data/rqlt2_data-H/phenotype_im_Harm_arabid.csv', index = None, header = True)

## control file in yaml format 
'''
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

'''
