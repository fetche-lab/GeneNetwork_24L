### EXTRACTING VARIANTS FROM VCF FILES FOR GENENETWORK2 
### The following information provides an overview on how to process vcf files inorder to extract variants/genotypes and use them in gn2 related data analysis 
### The following are the main steps/procedures to use in order to succeed getting the variants/genotypes 
### 01 Extracting genotypes from vcf files 
In this step, `bcftools` is used to first; filter variants with less than 5% allelic frequency, and second; extract genotypes with their associated positions using special options.
```bash
#!/bin/bash 

## The following commands were used to filter and retrieve genotype information from two vcf files 

## Variant filtering (variants with less than 5% allelic frequency [AF])

## file.vcf 
bcftools view -i 'AF>0.05' file.vcf -o filtered_file.vcf 

## Extracting genotype information from the filtered vcf files 

# extract sample headers 
bcftools query -l filtered_file.vcf | tr '\n' '\t' > sample_headers.txt  

# add the field headers to a new file 
echo -e 'CHROM\tPOS\tREF" > fields.txt  
echo -e "\n" >> fields.txt

# concatinate the field and sample headers 
cat fields.txt sample_headers.txt > headers.txt

# extract the genotypes from the files 
bcftools query -f '%CHROM\t%POS\t%REF\t[\t%GT]\n' filtered_file.vcf | cat headers.txt -> filtered_file_genotype.txt 

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
`gemma (Genome-wide Efficient Mixed-Model Association)` is a bioinformatics tool used in gwas studies, more information in; https://github.com/genetics-statistics/GEMMA/tree/master. In this step, the genotypes values are converted into gemma format, so as to facilitate easy experimentation with gemma. A custom python script was used to achieve this step; 
```python

```
