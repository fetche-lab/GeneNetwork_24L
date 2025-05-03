## Kilifish Project 
The following contains information on all steps and actions taken to process kilifish datasets into the required GeneNetwork format before the uploading process 

### 01 Genotype file 
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
    '## LIA stands for "Leiniz Institute of Aging" from which this data was generated', 
    '## `ref` refers to the homozygous reference allele',
    '## `alt` refers to the homozygous alternative allele',
    '## `het` refers to the heterozygous alleles, containing one allele from ref and the other from alt',
    '## `unk` represents all unknown genotypes', 
    '## the crossing used in this data is of second generation from the parent generation',
    ' ', 
    '@name:F2-KILIFISH-UTHS',
    '@type:F2 cross',
    '@ref:A',
    '@alt:B',
    '@het:H',
    '@unk:U', 
    ' ' 
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

### 02 Phenotype file(s)
#### a) CaseAttributes 
> CaseAttributes provide context useful in data comparison during mapping. 

> This dataset had `sex`, and `Age` range in days as case attributes 

> The processing was as follows; 
```python 
## import pandas 
import pandas as pd 

## load in the file with attributes 
 attr_df = pd.read_csv("./2025-04-08RQTL_mapping_phenotype_master_table_Felix.csv", low_memory=False)

## inspect the file, then select the columns of interest 
### filter unwanted rows 
attr_df = attr_df[attr_df['Generation'] == 'F2'] #selects all rows with F2 as their value in Generation column
### select cols of interest 
cols_to_pick = ["ID", "DOB", "DOS", "sex"]
attr_df = attr_df[cols_to_pick]

## adjust DOB and DOS into time range in days 
# Convert DOB and DOS columns to datetime format, handling different formats
attr_df['DOB'] = pd.to_datetime(attr_df['DOB'], format='%d.%m.%y', errors='coerce')
attr_df['DOS'] = pd.to_datetime(attr_df['DOS'], format='%d.%m.%y', errors='coerce')

# Calculate the difference in days between DOS and DOB
attr_df['Age(days)'] = (attr_df['DOS'] - attr_df['DOB']).dt.days
#attr_df['Age(months)'] = ((attr_df['DOS'] - attr_df['DOB']).dt.days / 30.44).round(2)

## drop DOB and DOS
cols_to_drop = ["DOB", "DOS"]
attr_df.drop(columns=[cols_to_drop], axis=1, inplace=True)

## save the df into a csv file 
attr_df.to_csv("./Kilifish_CaseAttributes.csv", index=False)

```

#### b) Experimental phenotypes 
> These are non expression measurements that are also useful in mapping and association experiments between genotypes and phenotypes 

> The processing is as follows; 
```python

## import pandas 
import pandas as pd 

## load in the file with attributes 
 exp_df = pd.read_csv("./2025-04-08RQTL_mapping_phenotype_master_table_Felix.csv", low_memory=False)

### filter unwanted rows 
exp_df = exp_df[exp_df['Generation'] == 'F2']

### select wanted columns 
cols_to_pick = ["ID", "Generation", "DOB", "DOS", "Brain.weight", "sex", "Fish.length.no.tail"]
exp_df = exp_df[cols_to_pick]

## adjust DOB and DOS into time range in days 
# Convert DOB and DOS columns to datetime format, handling different formats
exp_df['DOB'] = pd.to_datetime(exp_df['DOB'], format='%d.%m.%y', errors='coerce')
exp_df['DOS'] = pd.to_datetime(exp_df['DOS'], format='%d.%m.%y', errors='coerce')

# Calculate the difference in days between DOS and DOB
exp_df['Age(days)'] = (exp_df['DOS'] - exp_df['DOB']).dt.days
#exp_df['Age(months)'] = ((exp_df['DOS'] - exp_df['DOB']).dt.days / 30.44).round(2)

## drop DOB and DOS
cols_to_drop = ["DOB", "DOS"]
exp_df.drop(columns=[cols_to_drop], axis=1, inplace=True)

## add the SE and N cols after each of the available columns, then add string x as the values for both added cols 
#Create a new DataFrame with the desired structure
new_cols = []

for cols in exp_pheno_df.columns: 
    new_cols.append(cols) 
    new_cols.append(f"SE") 
    new_cols.append(f"N")

#Create a new DataFrame with the same number of rows and fill SE and N with 'x'
exp_pheno_df1 = pd.DataFrame(columns=new_cols)

for cols in exp_pheno_df.columns:

    # Assign values from original DataFrame
    exp_pheno_df1[cols] = exp_pheno_df[cols] 

    # Assign 'x' to SE and N columns
    exp_pheno_df1[f"SE"] = 'x' 
    exp_pheno_df1[f"N"] = 'x'

## Transpose the df and maintain the indexing 
exp_pheno_df2 = exp_pheno_df1.set_index("ID").T
exp_pheno_df2 = exp_pheno_df2.reset_index()
exp_pheno_df2.rename(columns={'index':'Strain'}, inplace = True)ce=True)

## drop the first two rows 'SE' and 'N' as they are not needed for the col headers
exp_pheno_df2 = exp_pheno_df2.drop(index=[0, 1])
exp_pheno_df2.reset_index(drop=True, inplace=True)

## save the file 
exp_pheno_df2.to_csv("../processed/Kilifish_Experimental_Phenotypes.csv", index=False)

```
#### c) Expression phenotypes (transcriptomics and proteomics)

> The dataset includes transcriptomics and proteomics datasets 

> Both were generated from samples of either of the Kilifish brain's hemisphere 

> all columns with "XM.." are transcriptomic data

> all columns with "A.." are proteomics data 

> Normalization techniques used: 
>> Transcriptomics dataset was `rpm (reads per million)` normalized 
>> Proteomics dataset was `spike-in control` normalized 

> Processing the expression dataset is as follows: 
```python 
## import packages 
import pandas as pd 

## load in data 
expression_df = pd.read_csv("2025-04-08RQTL_mapping_phenotype_master_table_Felix.csv", low_memory=False) 

## inspect the df 
### in this case; select all rows with F2 as their generation 
expression_df = expression_df[expression_df["Generation"] == "F2"]

## drop columns not needed in this process 
cols_to_drop = ["Generation", "DOB", "DOS", "Brain.weight", "sex", "Fish.length.no.tail"]

expression_df.drop(columns=cols_to_drop, axis=1, inplace=True)

## now extract the transcriptomic data 
xm_cols = [cols for cols in expression_df.columns if cols.startswith("XM_")]#selects columns starting with "XM_"
xm_df = expression_df[xm_cols]#creates a new df with the xm columns 
xm_df.insert(0,"ID", expression_df['ID'])# adds the individual IDs 

## now extract proteomics data 
proteomics_cols = [cols for cols in expression_df.columns if not cols.startswith("XM_")]# extracts all columns not starting with "XM_" 
proteomics_df = expression_df[proteomics_cols]#creates new df with proteomics columns 

## standardize to 3 decimal places for both 
##xm_df
xm_df.iloc[:, 1:] = xm_df.iloc[:,1:].round(3)
xm_df.fillna("NA", inplace=True)
##proteomics_df 
proteomics_df.iloc[:, 1:] = proteomics_df.iloc[:, 1:].round(3)
proteomics_df.fillna("NA", inplace=True)

## save the dfs into files (as csv and text)
xm_df.to_csv("../processed/Kilifish_transcriptomic_data.csv", index=False)
xm_df.to_csv("../processed/Kilifish_transcriptomic_data.txt", index=False, header=None, sep="\t")
proteomics_df.to_csv("../processed/Kilifish_proteomic_data.csv", index=False)
proteomics_df.to_csv("../processed/Kilifish_proteomic_data.txt", index=False, header=None, sep="\t")

```

### 03 Metadata 
#### a) Experimental Phenotype descriptions 

> This involves decriptive information on the phenotypes selected to be used in this project 

> The process involves updating the column headers reflecting the following information: 

>>`['Pubmed ID', 'Pre Publication Description','Post Publication Description', 'Original Description','Pre Publication Abbreviation', 'Post Publication Abbreviation','Lab Code', 'Submitter', 'Owner', 'Autorized Users', 'Authors', 'Title','Abstract', 'Journal', 'Volume', 'Pages', 'Month', 'Year', 'Units']`

> A simple excel manipulation is enough to get this done. It is important however, to adhere to the gn2 descriptions guidelines, more information on this link: [GN2_guidelines]("https://info.genenetwork.org/sop/data_submission/")

b) Annotation information??
`TODO` 
