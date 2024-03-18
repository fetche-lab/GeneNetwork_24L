#!/bin/bash 

## The following commands were used to filter and retrieve genotype information from two vcf files 

## Variant filtering (variants with less than 5% allelic frequency [AF])

## mikk_small.vcf 
bcftools view -i 'AF>0.05' mikk_small.vcf -o filtered_mikk_small.vcf 

## mikk_clean.vcf 
bcftools view -i 'AF>0.05' mikk_clean.vcf -o filtered_mikk_clean.vcf

## Extracting genotype information from the filtered vcf files 

# extract sample headers 
bcftools query -l filtered_mikk_small.vcf | tr '\n' '\t' > sample_headers.txt # same sample headers for mikk_clean.vcf 

# add the field headers to a new file 
echo -e 'CHROM\tPOS\tREF" > fields.txt  
echo -e "\n" >> fields.txt

# concatinate the field and sample headers 
cat fields.txt sample_headers.txt > headers.txt

# extract the genotypes from the files 
bcftools query -f '%CHROM\t%POS\t%REF\t[\t%GT]\n' filtered_mikk_small.vcf | cat headers.txt -> filtered_mikk_small_genotype.txt 
bcftools query -f '%CHROM\t%POS\t%REF\t[\t%GT]\n' filtered_mikk_clean.vcf | cat headers.txt -> filtered_mikk_clean_genotype.txt 
