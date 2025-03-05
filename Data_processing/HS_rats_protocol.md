# Heterogenous Stock (HS) Rats Protocol 
The following contains scripts used to process datasets on HS rats
These include: 
 - Expression data 
 - Experimental phenotypes 
 - Genotype data 
 - Phenotype Descriptions 

## 01 Processing genotypes 
- From the downloads, HS genotype files are in VCF format. So, there is a separate script to process vcf files that is also used for HS rats. Yet, there is still a need to optimize and automate this for future HS rats data. 

- Check the following link, for the genotype processing; [genotyping_vcf](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/genotype_smoothing.md)

 
## 02 Log2 dataset 
*(more information to be added here after several tests with a variety of data)*

```python 
# import required library/libraries 
import pandas as pd 

# load in the file into dataframe 
file_df = pd.read_table('file.rn7.expr.log2.bed.gz', compression='gzip', low_memory=False)

# drop the unwanted columns 
cols_to_drop = ['#chr', 'start', 'end']
file_df.drop(columns=cols_to_drop, axis=1, inplace=True)

# rename the gene_id to ID 
file_df.rename(columns={'gene_id':'ID'}, inplace=True)

# Round the values into 4 decimal places 
file_df_rounded = file_df.round(4)

# fill the empty cells with `NA` 
file_df_rounded.fillna('NA', inplace=True)

# save the processed file into tsv and xlsx format 
file_df_rounded.to_csv('processed/file.rn7.expr.log2.tsv', index=False, float_format='%.4f', sep='\t') #tsv 
file_df_rounded.to_excel('processed/file.rn7.expr.log2.xlsx', index=False) #xlsx

```

## 03 Experimental phenotypes 
*(more information to be added here after several tests with a variety of data)*

- The experimental dataset in this comes in xlsx sheets. The sheets of interest will be selected and processed

- Python will be used to process the files 

```python 
# import the required lib(s)
import pandas as pd 

# load in the excel file into a dataframe (select the needed sheets)
phenotype_df = pd.read_excel('./Supplementary_Tables_revise_25Aug.xlsx', sheet_name=1, skiprows=#ofrowstoskip) # selected the index of the sheet needed 

## filter some of the columns that are not regarded as phenotype perSec 
## in this case (Handler, Bleeder, Injector, Dissector)
## Modify the `Tissue.Harvest.Age` col into `Tissue.Harvest.Age(Weeks)`, then remove the weeks part in col values 
cols_to_drop = ['Handler', 'Bleeder', 'Injector', 'Dissector']
phenotype_df.drop(columns=cols_to_drop, axis=1, inplace=True)
phenotype_df.rename(columns={'Tissue.Harvest.Age':'Tissue.Harvest.Age(weeks)'}, inplace=True)
phenotype_df['Tissue.Harvest.Age(weeks)'] = phenotype_df['Tissue.Harvest.Age(weeks)'].str.extract('(\d+)')
phenotype_df['Tissue.Harvest.Age(weeks)'].unique()

# add the SE and N cols after each of the available columns, then add string x as the values for both added cols 


## Create a new DataFrame with the desired structure
new_cols = [] 

for cols in phenotype_df.columns: 
    new_cols.append(cols)
    new_cols.append(f"SE")
    new_cols.append(f"N") 

## Create a new DataFrame with the same number of rows and fill SE and N with 'x'
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
pheno_df_T = phenotype01_df.set_index('Rat ID').T 
pheno_df_T = pheno_df_T.reset_index()
pheno_df_T.rename(columns={'index':'Strain'}, inplace=True)

## drop the first two rows 'SE' and 'N' as they are not needed for the col headers
pheno_df_T = pheno_df_T.drop(index=[0,1])
pheno_df_T.reset_index(drop=True, inplace=True)

## Save your file 
pheno_df_T.to_csv('./Adipose_Liver.rn7.exp-phenotypes-updated.tsv', index=None, header=True, sep='\t')

```

## 04 Experimental phenotype descriptions 
*(more information to be added here after several tests with a variety of data)*

- This part contains the information describing the experimental phenotypes 

- load the summary file, process it to get the description information 

```python 
## import the required libraries 
import pandas as pd 

## load in the xlsx sheet of interest (in this case, index = 0)
desc_df = pd.read_excel('./Supplementary_Tables_revise_25Aug.xlsx', sheet_name=0, skiprows=[0,1,2,24,25,26,27])

## drop unwanted columns 
desc_df.drop(columns=['Unit', 'Number.of.rats', 'Transformation', 'Covariates', 'Heritability(SE)'], axis=1, inplace=True)

## process the missing columns use the following general format 
desc_df['#missing-col'] = '#values'

## process the publication description cols 
desc_updated_df['Post Publication Description'] = 'Cofactor, Metadata: ' + desc_updated_df['Post Publication Description'] + ' [' + desc_updated_df['Pre Publication Abbreviation'] + ']' # do same for pre publication description 

## rearrange the columns in the needed order 
#needed col order
cols = ['Pubmed ID','Pre Publication Description','Post Publication Description','Original Description','Pre Publication Abbreviation','Post Publication Abbreviation','Phenotype','Lab Code','Submitter','Owner','Authorized Users','Authors','Title','Abstract','Journal','Volume','Pages','Month','Year','Units']

#update col order 
desc_updated_df = desc_df[cols]

#save the file
desc_updated_df.to_csv('./Adipose_Liver.rn7.exp-phenotypes-descriptions.tsv', index=None, header=True, sep='\t')

```

## Processing the strain and strainXref 
*(This is yet to be fully automated as it highly depends on the backend GN team to decide how the strain and the strainXref files should be structured and entered into the database)*
