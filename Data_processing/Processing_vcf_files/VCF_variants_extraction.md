### EXTRACTING VARIANTS FROM VCF FILES FOR GENENETWORK2 
### The following information provides an overview of how to process vcf files to extract variants/genotypes and use them in gn2-related data analysis 
### The following are the main steps/procedures to use to succeed in getting the variants/genotypes 
### 01 Extracting genotypes from vcf files 
In this step, `bcftools` is used to first; filter variants with less than 5% allelic frequency, and second; extract genotypes with their associated positions using special options.
```bash
#!/bin/bash 

## The following commands were used to filter and retrieve genotype information from two vcf files 

## Variant filtering (variants with less than 5% allelic frequency [AF])

## file.vcf 
bcftools view -i 'AF>0.05' file.vcf -o filtered_file.vcf 

## Extracting genotype information from the filtered vcf files 

bcftools query -l filtered_file.vcf | tr '\n' ',' | sed 's/,$//' | awk '{print "CHROM,POS,ID," $0}' > sample01_genotypes.csv
bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' filtered_file.vcf | awk -F'\t' 'BEGIN {OFS=","} {print}' >> sample01_genotypes.csv
```
An example of a genotype file would look like this: 

```bash
   MIKK4_1	MIKK4_2	MIKK5_1	MIKK7_2	MIKK8_2	MIKK10_1	...
0	 1/1	    1/1	    1/1	    ./.	    ./.	    1/1		...
1	 1/1	    1/1	    1/1	    1/1	    ./.	    ./.		...
2	 1/1	    1/1	    1/1	    1/1	    ./.	    ./.		...
3	 ./.	    1/1	    ./.	    ./.	    ./.	    1/1	  ...
4	 1/1	    1/1	    1/1	    1/1	    ./.	    1/1		...

```
### 02 Encoding the genotypes into gemma format 
`gemma (Genome-wide Efficient Mixed-Model Association)` is a bioinformatics tool used in GWAS studies, more information at; https://github.com/genetics-statistics/GEMMA/tree/master. In this step, the genotype values are converted into gemma format, to facilitate easy experimentation with gemma. A custom Python script was used to achieve this step; 
```python
## import pandas 
import pandas as pd 

## import the geno file as a dataframe 
file_df = pd.read_csv('./filtrd_file_genotypes.txt', sep='\t') 

## read the information in the file with headers, then extract them as a list
headers_lst = []
with open('headers.txt', 'r') as file1: 
    headers = file1.readlines()
    
for line in headers: 
    headers_lst.append(line[:-1]) # this removes the last `\n` character 
print(len(headers_lst))

## remove the column with empty values 
file_df1 = file_df.drop(columns=['Unnamed: 79'])

## append the headers corresponding to the genotype file 
file_df1.columns = headers_lst

## check for the unique values 
unique_lst=[]
for names in file_df1.columns: 
    unique_names = file_df1[names].unique()
    unique_lst.append(unique_names)
for data in unique_lst:
    print(data)
    
## replace values in the geno file with the gemma required standard
#df.replace({'old_value1': 'new_value1', 'old_value2': 'new_value2'}, inplace=True)
file_df1.replace({'0/3':'0.0','2/6': '1.0', '1/2':'1', '1/4':'1', '3/3':'1.0',  
                  '4/4':'1.0', '0/4':'0.0', '0/5':'0.0', '2/4':'1.0', '6/6':'1.0', 
                  '5/5':'1.0', '2/3':'0.0', '0/6':'0.0', '4/6':'1.0', '3/4':'1.0', 
                  '2/5':'1.0', '1/6':'1.0', '3/5':'1.0', '1/5':'1.0', '3/6':'1.0', 
                  '4/5':'1.0', '2/0':'1', '1/3':'1.0', '0/0':'0', '0/1': '0.5',
                  '1/1':'1', '0/1':'0.5', '1/2':'1','5/6':'1.0'}, inplace=True)

## save the dataframe into a csv file format 
file_df1.to_csv('file_genotype_gemma_encoded.csv', index = None, header = True)
```
The output should look something like this; 
```bash
	MIKK4_1	MIKK4_2	MIKK5_1	MIKK7_2	MIKK8_2	MIKK10_1	MIKK11_1	MIKK11_2	MIKK13_2	MIKK14_1	...	MIKK132_4_1	MIKK132_5	MIKK134_1	MIKK135_1	MIKK135_2	MIKK137_4	MIKK138_1	MIKK139_4	MIKK140_1	MIKK140_3
Ol223467v1_chr1_12989	1	1	1	NA	NA	1	1	NA	1	1	...	NA	NA	1	1	1	1	NA	1	1	1
Ol223467v1_chr1_13217	1	1	1	1	NA	NA	NA	1	1	1	...	1	1	NA	NA	1	NA	1	NA	1	1
Ol223467v1_chr1_13234	1	1	1	1	NA	NA	NA	1	1	1	...	1	1	NA	NA	1	NA	1	NA	1	1
Ol223467v1_chr1_13314	NA	1	NA	NA	NA	1	1	1	1	1	...	1	0	NA	NA	NA	NA	NA	1	1	NA
Ol223467v1_chr1_13359	

```

### In certain cases, the vcf files contain redundant variants. So, it is important to filter all non-informative variants, to ensure a more accurate and precise downstream analysis using the variants/genotypes. The following steps include key processes involved in achieving this. 

### 01 Vcf preprocessing and genotype phasing 
Use `bcftools` to filter out low-quality variants as in `step 01` above. The next step is to convert the unphased genotypes `e.g 0/1, 1/1, etc` into phased genotypes `e.g 0|1, 1|1, etc`. There are available bioinformatics tools for this task. However, in this situation, `py-popgen/vcf_phase.py`: https://ppp.readthedocs.io/en/latest/PPP_pages/Functions/vcf_phase.html, was used to process the genotypes. The phased genotypes can then be processed to generate hapmaps and identify and filter out the redundant variants. 
```bash
vcf_phase.py --vcf file_filtrd.vcf --phase-algorithm beagle --out file_filtrd_phased.vcf --out-format vcf 

```
NB; sometimes, the vcf file is too large for these tools to handle. The best approach is to split the file into chunks, process each chunk, and merge the files back to the original state. Examples of the scripts used in this practice are as follows; 
A) Splitting the vcf file 
```python
import os 

def split_vcf(input_file, output_dir, chunk_size):
    """
    Splits a large vcf file into smaller files while maintaining the metadata and headers. 

    Args:
        input_file (str): Path to the input vcf file.
        output_dir (str): Path to the output directory where the smaller files will be saved. 
        chunk_size (int): Number of variants to include in each smaller file. 
    
    """ 

    os.makedirs(output_dir, exist_ok=True)

    with open(input_file, 'r') as f: 
        #Read the metadata and header lines 
        metadata_and_header = [] 
        line = f.readline()
        while line.startswith('##') or line.startswith('#'):
            metadata_and_header.append(line)
            line = f.readline()

        # split the vcf file into smaller files 
        file_count = 1 
        chunk = [] 
        while line: 
            chunk.append(line)
            if len(chunk) == chunk_size: 
                output_file = os.path.join(output_dir, f'chunk_{file_count}.vcf')
                with open(output_file, 'w') as out_f:
                    out_f.writelines(metadata_and_header)
                    out_f.writelines(chunk)
                chunk = []
                file_count += 1 
            line = f.readline()

        #write the last chunk if any 
        if chunk: 
            output_file = os.path.join(output_dir, f'chunk_{file_count}.vcf')
            with open(output_file, 'w') as out_f:
                out_f.writelines(metadata_and_header)
                out_f.writelines(chunk)

# run the function 
split_vcf('file_filtrd.vcf','./vcf_chunks_500k', 500000)

```
B) Phasing the chunks 
```bash
#!/bin/bash

#set the input and output directories
input_dir="/path to the input files (chunks)"
output_dir="path to the phased_output"

#iterate through the files in the input directory 
for file in "$input_dir"*.vcf; do 

    #check if the entry is a regular file 
    if [ -f "$file" ]; then 

        #get the filename without extension 
        filename=$(basename "$file" | cut -d"." -f1)

        #process the file and save the output 
        vcf_phase.py --vcf "$file" --phase-algorithm beagle --out-prefix "$output_dir/${filename}_phased" --out-format vcf
        echo "Processed file: $file" 
    fi 
done 

echo "Processing complete." 
```
c) Merge back the processed vcf chunks 
```bash
## Compress the phased vcf chunks
 for file in chunk_*_phased.vcf; do bgzip -c $file > $file.gz; done

## Index the gzipped vcf chunks using bcftools 
for file in chunk_*_phased.vcf.gz; do bcftools index $file; done

## Merge the chunks using bcftools
bcftools concat $(ls chunk_*_phased.vcf.gz | sort -V) -Ov -o file_merged_phased.vcf 

```

### 02 Generating hapmaps from the phased vcf files 
In this case, a set of R packages and custom scripts were used to generate hapmap from the phased genotypes. This script can handle large files too. 
```R
## EXPERIMENTING WITH R TO GENERATE HAPLOTYPES FROM VCF FILES 
## SET THE WORKING DIRECTORY 
getwd()
setwd('/path to the working directory/')
## INSTALL AND LOAD THE REQUIRED LIBRARIES 
library(dplyr)
library(stringr)
require(reshape2)
require(data.table, quietly = T)
require(stringi)

hap_data <- fread('file_merged_phased.vcf')

processed_data <- data.table(gsub("|", "", do.call(paste0, hap_data[, -c(1:9)]), fixed = TRUE))

write.csv(processed_data, file = "file_merged_phased_hapmap.csv", row.names = T)

```
### 03 Filter the redundant variants 
This step involved processing the hapmap and filtering out the repetitive rows that represent the redundant variants. A custom Python script is used to achieve this. 
```python

## load the hapmap vcf file into a dataframe
file_df = pd.read_csv("file_filtrd_phased_hapmap.csv")

## begin the filtering process 
## a function to inspect redundant rows 
def same_value(row): 
    return len(set(row)) == 1 

# create the boolean filter
same_val_mask = file_df['Haplotypes'].apply(same_value) 

## Select rows where all values are the same 
redundant_df = file_df[same_val_mask] 

## Select rows where values are different 
unique_df = file_df[~same_val_mask] 

## drop duplicates in the same value df and add the unique values to the diff val df 
## dropping the duplicates 
sorted_redundant_df = redundant_df.drop_duplicates(subset='Haplotypes') 

## add the sorted values into the diff val df 
final_df = pd.concat([unique_df, sorted_redundant_df])

## filter the redundant file from the actual file 
original_df = pd.read_csv('file_geno_marked2.csv')
filtered_df = original_df[original_df['Unnamed: 0'].isin(final_df.index)]

##save the new variant filtered genotype file 
filtered_df.to_csv('file_hapVariants.csv', index=None, header = True)

```
The following script can be used when processing relatively large file 
```python
### Proceed from the `final_df` in the script above; 
## Define the file path and column names (if applicable)
file_path = './file_geno_marked2.csv' 
chunk_size = 100000 # Adjust this according to your system's memory capacity

# Initialize an empty list to store processed chunks
processed_chunks = []

# Read the file in chunks
for chunk in pd.read_csv(file_path, chunksize=chunk_size):
    
    #fill empty spaces with NA
    chunk.fillna('NA', inplace=True)
    
    #rename the marker header 
    chunk.reset_index()
    chunk=chunk.rename(columns={'Unnamed: 0':'Marker_Ids'})
    
    #filter non informative genotypes 
    chunk=chunk[chunk['Marker_Ids'].isin(final_clean_df.index)]
    
     #append the processed chunks into a list
    processed_chunks.append(chunk)
    
#save the list as a dataframe 
final_processed_df = pd.concat(processed_chunks)

#save it into a csv file 
final_processed_df.to_csv('file_geno_hapVariants.csv', index = None, Header = True)

```
The last part of the script above is to save the merged dataframe into a file. If it fails to work due to memory limits, the following R script can help solve the problem 

```R
getwd()
setwd('/setting path to the working directory')

# Load the required library
library(dplyr)

# Read the files into data frames
file1 <- read.csv("file_geno_marked2.csv")
file2 <- read.csv("file_hapvarMaps.csv")

# Filter file1 based on the first column of file2
filtered_file1 <- file1 %>% 
  filter(X %in% file2$Marker_ids)

# Write the filtered data to a new file
write.csv(filtered_file1, "file_hapVariants.csv", row.names = FALSE)

```
This workflow is still under progressive development. Contributions to its improvement will greatly be appreciate.
