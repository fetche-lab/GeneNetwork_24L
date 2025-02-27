### 02 Processing classical phenotypes 

### In this case, we will use python/pandas to process the classical phenotype file 


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
