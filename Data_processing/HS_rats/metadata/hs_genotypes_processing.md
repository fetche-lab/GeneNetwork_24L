## Processing and formatting hs rats vcf files 
> The following explains steps to extract genotypes from vcf files (for hs rats), and format the resulting file for genenetwork webserver 

### QC the file 
> filter low frequency markers, below `0.05` threshold
```sh 
bcftools view -i 'AF>0.05' file.vcf -o filtered_file.vcf 
```
> perform LD pruning, `r2=0.99, window=100` 
```sh 
 bcftools +prune -m 0.99 -w 100 filtered_file.vcf -Ov -o filtered_file_pruned.vcf 
```
### Phase genotypes  
> Use (https://faculty.washington.edu/browning/beagle/beagle.html)[beagle] for processing 
```sh
java -Xmx10g -jar ./beagle.27Feb25.75f.jar gt=./filtered_file_pruned.vcf out=./filtered_file_pruned_phased nthreads=10 
```
> Index the resulting phased vcf file 
```sh 
bcftools index filtered_file_pruned_phased.vcf 
```
### Generate hapmaps and genotypes 
> Use R packages to generate hapmap file from the phased vcf file 
```R 
#!/usr/bin/env Rscript 

## EXPERIMENTING WITH R TO GENERATE HAPLOTYPES FROM VCF FILES
## SET THE WORKING DIRECTORY
#getwd()
setwd("$HOME/phased_files/")
## INSTALL AND LOAD THE REQUIRED LIBRARIES
install.packages('R.utils')
library(dplyr)
library(stringr)
require(reshape2)
require(data.table, quietly = T)
require(stringi)

hap_data <- fread('filtered_file_pruned_phased.vcf') # sometimes add the quote="" to avoid warnings related to improper quotings

processed_data <- data.table(gsub("|", "", do.call(paste0, hap_data[, -c(1:9)]), fixed = TRUE))

write.csv(processed_data, file = "hapmap.csv", row.names = T)
```
> Extract genotypes from the phased vcf file 
```sh 
(echo -e "CHROM\tPOS\tID$(bcftools query -l filtered_file_pruned_phased.vcf | awk '{printf "\t%s", $0}')"; \
 bcftools query  -f '%CHROM\t%POS\t%ID[\t%GT]\n' filtered_file_pruned_phased.vcf) > genotypes.txt
```
> Filter the hapmap using custom python scripts 
```python 
## import packages 
import pandas as pd 
import numpy as np 

# load in the file 
hapmaps = pd.read_csv("hapmap.csv",low_memory=False)
genotypes = pd.read_csv("genotypes.txt", sep="\t", low_memory=False)

hapmaps.rename(columns={"Unnamed: 0":"ids"}, inplace=True)
hapmaps.insert(0, "markers", genotypes["ID"])

## a function to filter redundant/monomorphic snps 
def same_value(row): 
    return len(set(row)) == 1

## apply the function 
#same_val_mask = hapmaps['ids'].apply(same_value)
same_val_mask = hapmaps['ids'].astype(str).apply(same_value)
redundant_df = hapmaps[same_val_mask]
unique_df = hapmaps[~same_val_mask]

#filter snps using the unique hapmap df values 
filtered_genotypes = genotypes[genotypes["ID"].isin(unique_df["markers"])]

## process and format the resulting geno df 
# gemma encoding
genotypes_cols = filtered_genotypes.columns.difference(['CHROM', 'POS', 'ID'])
genotype_map = {
'0|0': 1,
'0|1': 2,
'1|0': 2,
'1|1': 3
}
filtered_genotypes[genotypes_cols] = filtered_genotypes[genotypes_cols].replace(genotype_map)

# process marker ids 
def format_marker(marker): 
    if isinstance(marker, str) and ":" in marker:
        chrom, pos = marker.split(":") 
        try: 
            return f'hsr{int(pos):08d}'
        except ValueError:
            return marker 
    return marker
filtered_genotypes["ID"] = filtered_genotypes["ID"].apply(format_marker)
#and Xsomes 
filtered_genotypes["Chr"] = filtered_genotypes["Chr"].str.replace('chr', '', regex=False)
filtered_genotypes.fillna("NA", inplace=True)

# headers and their proper order 
filtered_genotypes.rename(columns={"CHROM":"Chr", "POS":"Mb", "ID":"Locus"}, inplace=True)
filtered_genotypes["Mb"] = filtered_genotypes["Mb"]/1000000

new_order_cols = ["Chr", "Locus", "Mb"] + [col for col in filtered_genotypes.columns if col not in ["Chr", "Locus", "Mb"]]
filtered_genotypes01 = filtered_genotypes[new_order_cols]

## gn2 encoding 
geno_cols = filtered_genotypes01.columns.difference(["Chr", "Locus", "Mb"])
geno_map = {1:"A", 2:"H", 3:"B", "NA":"U"}
filtered_genotypes01[geno_cols] = filtered_genotypes01[geno_cols].replace(geno_map)

## add metadata 
hs_metadata = [
    '## This file has genotype data representing HS (Heterogeneous Stock) rats from Abe’s group (RatGTEx)',
    '## It represents 10 tissues exclusive of liver and adipose tissues, which have their own separate genotype file', 
    '## `ref` refers to the homozygous reference allele',
    '## `alt` refers to the homozygous alternative allele',
    '## `het` refers to the heterozygous alleles, containing one allele from ref and the other from alt',
    '## `unk` represents all unknown genotypes',
    '## the marker ids were formatted so that the include {`hsr` + `position value`} = `hsrxxxxx`', 
    '## `type` represents cross type, in this case its Heterogeneous Stock cross, derived from 8 inbred strains',
    ' ',
    '@name: HS-RATS-RatGTEx',
    '@type: Heterogeneous Stock (HS) cross',
    '@ref:A',
    '@alt:B',
    '@het:H',
    '@unk:U',
    ' '
]

hs_output = "hs_genotypes.geno"
with open(hs_output, 'w') as f:
    for line in hs_metadata: 
        f.write(line + '\n') 
    filtered_genotypes01.to_csv(f, sep='\t', index=False)

```
> The end 