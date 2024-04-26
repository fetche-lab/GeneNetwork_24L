## FILTERING THE RECUNDANT VARIANTS FROM THE GENOTYPE DATA 
### The following contains programming scripts involved in processing medaka fish genotype files to remove redundant variants 
### The steps involved are as follows: 
#### 01 Filtering and generating genotype data from vcf files 
The link to this step is: https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/vcf_to_genotype.sh

#### 02 Encoding the genotypes to gemma format 
The link to this step is: https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/mikk_data_wrangling.py

#### 03 Phasing the vcf files 
In this step, `vcf_phase.py` tool was used to generate the phased version of the vcf file. For more information, here is the link: https://ppp.readthedocs.io/en/latest/PPP_pages/Functions/vcf_phase.html
```bash 
vcf_phase.py --vcf mikk_small_filtrd.vcf --phase-algorithm beagle --out mikk_small_filtrd_phased.vcf --out-format vcf 

```
#### 04 Collapsing the phased genotypes and generate a hapmap 
In this step, the file generated was to be used to inspect and filter the redundant variants by condensing the genotypes 
```R
## EXPERIMENTING WITH R TO GENERATE HAPLOTYPES FROM VCF FILES 
## SET THE WORKING DIRECTORY 
getwd()
setwd('/home/fetche/Desktop/gn_remote/mikk_exp_new/beagle_algorithm/')
## INSTALL AND LOAD THE REQUIRED LIBRARIES 
library(dplyr)
library(stringr)
require(reshape2)
require(data.table, quietly = T)
require(stringi)

hap_data <- fread('mikk_small_filtrd_phased.vcf')

processed_data <- data.table(gsub("|", "", do.call(paste0, hap_data[, -c(1:9)]), fixed = TRUE))

write.csv(processed_data, file = "mikk_small_filtrd_phased_hapmap.csv", row.names = T)

```

#### 05 Filtering the redundant variants and updating the genotype data 
In this step, a custom python script was used, as follows; 
```python
## import pandas 
import pandas as pd

## load the mikk_small hapmap into a dataframes 
mikk_df1 = pd.read_csv("mikk_small_filtrd_phased_hapmap.csv")

## begin the filtering process 
## a function to inspect redundant rows 
def same_value(row): 
    return len(set(row)) == 1 

# create the boolean filter
same_val_mask = mikk_df2['Haplotypes'].apply(same_value) 

## Select rows where all values are the same 
redundant_df = mikk_df2[same_val_mask] 

## Select rows where values are different 
unique_df = mikk_df2[~same_val_mask] 

## drop duplicates in the same value df and add the unique values to the diff val df 
## dropping the duplicates 
sorted_redundant_df = redundant_df.drop_duplicates(subset='Haplotypes') 

## add the sorted values into the diff val df 
final_df = pd.concat([unique_df, sorted_redundant_df])

## filter the redundant file from the actual file 
original_df = pd.read_csv('mikk_small_geno_marked2.csv')
filtered_df = original_df[original_df['Unnamed: 0'].isin(final_df.index)]

##save the new variant filtered genotype file 
filtered_df.to_csv('mikk_small_hapVariants.csv', index=None, header = True)

```
