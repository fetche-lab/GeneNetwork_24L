## Genotype smoothing for VCF file formats 

## Step 01: Filter the vcf file 

- bcftools is used to the filtering 

### 01 Filter out Allele Frequencies (AF) less than 0.05, quality (QUAL) less than 20 and depth coverage (INFO/DP) less than 10  

```bash 
bcftools view -i 'AF>0.05 && QUAL>20 && INFO/DP > 10' input_file.vcf -o output_file_filtered.vcf 

```
### 02 Filter out redundant variants 

- `bcftools +prune` plugin is used to filter the markers with very high LD scores
  - r2 = 0.99 
  - pruning window = 50bp; important to consider heterozygosity and homozygosity
  - NB; the values above are subjective to change depending on the context of the study 
- example code;

```bash 
 bcftools +prune -m 0.99 -w 50 output_file_filtered.vcf -Ov -o output_file_filtered_pruned.vcf 

```
## Step 02: Impute and phase the vcf 

- In this step, we use beagle.

### 01 For files less than 1 Gb 

```bash 
java -Xmx8g -jar ./beagle/beagle.06Aug24.a91.jar gt=./mikk-germline/chunks/chunk_16.vcf out=./mikk-germline/chunk_16_phased nthreads=10

```
### 02 For files larger than 1 Gb 

- Break the large file into smaller chunks of at least 1 Gb each 

```python
import os
import argparse

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
        # Read the metadata and header lines
        metadata_and_header = []
        line = f.readline()
        while line.startswith('##') or line.startswith('#'):
            metadata_and_header.append(line)
            line = f.readline()

        # Split the vcf file into smaller files
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

        # Write the last chunk if any
        if chunk:
            output_file = os.path.join(output_dir, f'chunk_{file_count}.vcf')
            with open(output_file, 'w') as out_f:
                out_f.writelines(metadata_and_header)
                out_f.writelines(chunk)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split a VCF file into smaller chunks.')
    parser.add_argument('--input_file', required=True, help='Path to the input VCF file.')
    parser.add_argument('--output_dir', required=True, help='Directory to save the output chunks.')
    parser.add_argument('--chunk_size', type=int, required=True, help='Number of variants per chunk.')

    args = parser.parse_args()

    split_vcf(args.input_file, args.output_dir, args.chunk_size)


```
- On the commandline, you will run the above script as follows; 

```bash 
## an example of a run 
python3 -m split_vcf --input_file ./joint_germline-filtrd-prnd.vcf --output_dir chunks/ --chunk_size 250000

## to get more information about the script; 
python3 -m split_vcf --help 

```
- Iterate the phasing command on the chunked files 

```bash 
#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <beagle_jar_path>"
    echo ""
    echo "Arguments:"
    echo "  <input_directory>    The directory containing the input VCF files."
    echo "  <output_directory>   The directory where the phased output files will be saved."
    echo "  <beagle_jar_path>    The full path to the Beagle JAR file."
    echo ""
    echo "Example:"
    echo "  $0 /path/to/input_directory/ /path/to/output_directory/ /path/to/beagle/beagle.06Aug24.a91.jar"
    exit 1
fi

# Set the input and output directories from command line arguments
input_dir="$1/"
output_dir="$2/"
beagle_jar="$3"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate through the files in the input directory
for file in "$input_dir"*.vcf; do 

    # Check if the entry is a regular file 
    if [ -f "$file" ]; then 

        # Get the filename without extension 
        filename=$(basename "$file" | cut -d"." -f1)

        # Process the file and save the output using Beagle
        java -Xmx8g -jar "$beagle_jar" gt="$file" out="${output_dir}${filename}_phased" nthreads=10
        echo "Processed file: $file" 
    fi 
done 

echo "Processing complete."

```
- On the commndline, run the above script as follows; 

```bash 
## an example of the run 
./beagle_phase.sh /path/to/input_directory/ /path/to/output_directory/ /path/to/beagle/beagle.06Aug24.a91.jar

## for more help, run; 
./beagle_phase.sh -h 
```
- Next, index the chunks using `bcftools index` 

```bash 
for file in chunk_*_phased.vcf.gz; do bcftools index $file; done
##OR
bcftools index phased.vcf.gz ## for smaller files 
```

- Then merge back the chunks into one whole file. (for larger files)

```bash 
bcftools concat $(ls chunk_*_phased.vcf.gz | sort -V) -Ov -o file_merged_phased.vcf 
```
- NB; For more information regarding the beagle tool, how to install and use it, visit the following link [beagle]('https://faculty.washington.edu/browning/beagle/beagle.html')

## Step 03: Extract haplomaps from the phased vcf file 

### 01 Create a hapmap file 

- R packages are used 

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

hap_data <- fread('file_merged_phased.vcf') # sometimes add the quote="" to avoid warnings related to improper quotings 

processed_data <- data.table(gsub("|", "", do.call(paste0, hap_data[, -c(1:9)]), fixed = TRUE))

write.csv(processed_data, file = "file_merged_phased_hapmap.csv", row.names = T)

```

### 02 Smoothening the hapmap file by filtering non informative rows 

- Add the position column from the vcf file to the hapmap. This will help to track the positions that are filtered out when extracting genotypes from the phased vcf file 

- 
```python 
## load files into dataframe 
hapmap_df = pd.read_csv('./joint-germline_phased_hapmap.csv') #hapmap file 
geno_df = pd.read_table('./joint_germline_headers.txt', low_memory=False, header=None) ## contains position column

## rename headers 
geno_df.rename(columns={1:'position', 0:'Xsome'}, inplace=True)

## insert the position column in the hapmap dataframe 
hapmap_df.insert(1, 'marker', geno_df['position'])
hapmap_df.rename(columns={'V1':'Haplotypes'}, inplace = True)
hapmap_df.drop('Unnamed: 0', axis=1, inplace=True)

## perform smooothing 
mask = []
    
for haplotype in hapmap_df['Haplotypes']:
    # Check if haplotype contains only 0s or only 1s
    if set(haplotype) == {'0'} or set(haplotype) == {'1'}:
        mask.append(False)
        continue
    
    # Count the number of 0s and 1s
    count_0 = haplotype.count('0')
    count_1 = haplotype.count('1')
    
    # Calculate total length of the haplotype string
    total_count = len(haplotype)
    
    # Calculate ratio of 0s
    ratio_0 = count_0 / total_count

    # Calculate ratio of 1s 
    ratio_1 = count_1 / total_count
    
    # Filter out if ratio of 0s is more than 60%
    if ratio_0 > 0.4:
        mask.append(False)
    elif ratio_1 > 0.5: 
        mask.append(False)
    else:
        mask.append(True)

# Apply mask to filter DataFrame
filtered_df = hapmap_df[mask]

```
- Save the filtered df into a file and use it to filter the phased vcf genotypes after extracting them as in the next step 

## Step 04: Extract genotypes and related information from the vcf file 

### 01 Extraction 

- use bash custom commands and bcftools 

```bash 
## Extract sample names using bcftools 
bcftools query -l output.vcf.gz | tr '\n' '\t' > sample_headers.txt

## Create a variable to hold the fields names to be combined with  sample names 
header="CHROM\tPOS\tID\tREF\tALT\t$(paste -sd '\t' sample_headers.txt)"

## Extract the genotypes using bcftools and append the header variables to it 
(echo -e "$header" && bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[\t%GT]\n' output.vcf.gz) > output_genotypes.txt #for if the %ID has empty values, do not select it 

```
### 02 Encoding the genotypes into GEMMA format 

- Check for the unique values first. Then encode based on what is available 

```python 
## NB; Beware of column displacement when loading the file into a dataframe 

#Load the files into a dataframe 
genotype_df = pd.read_table('./joint_germline_genotypes.txt', low_memory=False)
hapmap_filtrd_df = pd.read_csv('./joint-germline_phased-filtrd_hapmap.csv')

## check for the unique values 
unique_lst=[]
for names in genotype_df.columns: #genotype_df.columns[5:] for the case of MIKK 
    unique_names = file_df1[names].unique()
    unique_lst.append(unique_names)
for data in unique_lst:
    print(data)

genotype_df.replace({'0|3':'0.0','2|6': '1.0', '1|2':'1', '1|4':'1', '3|3':'1.0', '4|4':'1.0', '0|4':'0.0', '0|5':'0.0', '2|4':'1.0', '6|6':'1.0', '5|5':'1.0', '2|3':'0.0', '0|6':'0.0', '4|6':'1.0', '3|4':'1.0', '2|5':'1.0', '1|6':'1.0', '3|5':'1.0', '1|5':'1.0', '3|6':'1.0','4|5':'1.0', '2|0':'1', '1|3':'1.0', '0|0':'0', '0|1': '0.5', '1|1':'1', '0|1':'0.5', '1|2':'1','5|6':'1.0','6|5':'1', '0|2':'0', '2|2':'1', '1|0':'1', '3|0':'1', '2|1':'1', '6|1':'1', '3|1':'1', '3|2':'1','5|0':'1', '4|3':'1', '4|2':'1', '4|1':'1', '5|2':'1', '4|0':'1', '6|0':'1', '5|1':'1', '6|3':'1', '5|4':'1', '5|3':'1', '6|2':'1', '6|4':'1'}, inplace=True) ## Note: this is subjective to the type of variant calls present in the file. "checking for unique values" part of this script will determine that 

## filter the non-informative variants from the coded df using the hapmap marker column 

filtered_genotype_df = genotype_df[genotype_df['POS'].isin(hapmap_filtrd_df['marker'])]

## save the file 
## TODO: HOW TO SAVE THE FILE IN GENO FORMAT TOO 

filtered_genotype_df.to_csv('./joint_germline_genotypes-BIMBAM.tsv', index=None, header = True, sep='\t')
```





