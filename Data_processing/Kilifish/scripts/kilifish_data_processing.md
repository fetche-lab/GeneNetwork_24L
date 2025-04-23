### Kilifish Project 
The following contains information on all steps and actions taken to process kilifish datasets into the required GeneNetwork format before the uploading process 

#### 01 Genotype file 
More information on regarding the nature and structure of the genotype file; 

- The “Geno” list, contains the genotype of each (506) F2 fish. It is comprised of a matched array which contains a SNP marker location (roughly one marker each 100kb) which corresponds to the reference genome: MPIA_NFZ_2.0. This reference genome had been made available to you under the file name: `annotation.gtf`
  
- The reference genome from N. furzeri, counts 19 chromosomes (including the sex determining chromosome), chromosome 3 (CM025010.1) on this assembly. The identified SNP’s (500 to 900 per chromosome) were characterized by Danny Arends using the following pipeline:
>> 1. Filter0: Require 3 genotypes

>> 2. Filter1: Very coarse filter for 25,50,25

>> 3. Filter2: Remove duplicate markers

>> 4. Filter3: Sliding window, keep only the best marker in a 50kb window, we slide forward by 1/4th of a window

>> 5. Filter4: Squared correlation filter (remove bad markers < 0.05 squared correlation)

>> 6. Filter5: 10 iterations of Loess filter for getting rid of bad markers

- Each SNP marker location is matched with a genotype score of 1, 2 or 3:

>>> `1` indicates that at that location, the individual has homozygous N. furzeri SNP (the reference). 

>>> `2` indicates that at that location, the individual is heterozygous and carries one SNP allele from N. furzeri (the reference) and one allele from N. kadleci (the alternative allele). 

>>> `3` indicates that at that location, the individual has homozygous N. kadleci SNP (the alternative allele).

For processing, the following is needed for the file to be ready for upload into GN2. 
> 1. column with chromosome inforamtion (`chr`)

> 2. column with marker information (`locus`)

> 3. column with SNP position (in `Mb` and `cM`)

> 4. letter encoded snp values (in this case, `1=A`, `2=H`, `3=B`, `NaN=U`)

> 5. metadata to describe the genotype file at the top of the actual values 

Now, let's process and format the genotype file; 

```python 
## import pandas 
import pandas as pd 

## load in the file into dataframe 
genotype_df = pd.read_csv("killi.v0.1.cross", low_memory=False, header=None) 

## inspect the file, filter unwanted rows 
index_to_drop = [0, 2, 3, 4, 5]
genotype_df.drop(index_to_drop, inplace=True)
genotype_df.reset_index(inplace=True) #restore the index order 
genotype_df.drop("index", axis=1, inplace=True)

## convert row 01 to be column header 
genotype_df.columns = genotype_df.iloc[0]
genotype_df.drop(0, inplace=True)
genotype_df.reset_index(inplace=True)
genotype_df.drop("index", axis=1, inplace=True)

## rename the nan values on columns 02 and 03 
rename_dict = {genotype_df.columns[1]:"Chr", genotype_df.columns[2]:"Mb"} #dict with renaming conversions
genotype_df.rename(columns=rename_dict, inplace=True)

## Convert the chromosome ids into numbers as per gn2 requirements 
genotype_df["CHR"] = genotype_df["Chr"].str.extract(r'(\d{2})\.1$').astype(int) - 7
genotype_df.drop("Chr", axis=1, inplace=True) #drop the old Xsome names 
genotype_df.insert(0, 'Chr', genotype_df['CHR']) 
genotype_df.drop("CHR", axis=1, inplace=True)

## Encode the genotype values 
mapping_criteria = {"1":"A", "2":"H", "3":"B", "NA":"U"}
geno_cols = genotype_df.columns[3:]
genotype_df[geno_cols] = genotype_df[geno_cols].replace(mapping_criteria)

## Work on the marker IDs 
genotype_df["Locus"] = 'rsk' + genotype_df['ID'].astype(str).str.zfill(5) #adds `rsk` string and buffers the id numbers to reflect the maximum id value, 5 digits in length 

### add metadata before saving to file 
## define metadata 
metadata = [
    '@name:F2-KILIFISH-UTHS',
    '@type:F2 cross',
    '@ref:A',
    '@alt:B',
    '@het:H',
    '@unk:U'
]

# Save to text file with metadata followed by DataFrame
output_file = "Kilifish_Genotypes.geno"

with open(output_file, 'w') as f:
    # Write metadata lines
    for line in metadata:
        f.write(line + '\n')
    # Append DataFrame as tab-separated values
    genotype_df.to_csv(f, sep='\t', index=False)

```

#### 02 Phenotype file(s)
`TODO....,`  
