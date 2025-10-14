## Genotype smoothing for Heterogeneous Stock (HS) Rats 

* The following contains a quick roadmap and important things to note, when dealing with multiparent genotypes. 

* There is still room for improvement on the current state of the art of genotyping, and this is one of the major steps to be taken, that is, the importance of genotype smoothing. 

* Genotyping smoothing, is simply an approach used to separate signal from noise in genotypes, so that one ends up having a clean and accurate set of genotypes for further downstream analyses. 

### More context to genotyping, phasing, imputation, and haplotyping in multiparent species
 
### Improving HS rats genotyping and smoothing 

* Define the rat genome assembly explicitly in metadata to clarify coordinate space and support reproducibility.

* Perform QC by counting genotype calls per type and per rat—deviations in expected ratios can flag data quality issues.

* Use pairs of markers to define haplotype blocks more reliably; this improves recombination detection and genotype confidence.

* Consider the expected recombination density (~2000 events per animal at Gen70) to assess whether marker density is sufficient.

* Evaluate genotype data quality and sparsity empirically, e.g., by analyzing the detection power for cis-eQTLs and checking for haplotype "runs."

### Genotyping accuracy in HS rats 

* Consider raising the minor allele frequency (MAF) threshold to 10–20%, given the 8 inbred progenitors and expected low haplotype diversity per locus in HS rats.

* Avoid using markers spaced <1000 bp apart, as drastic genotype shifts over short distances may indicate technical artifacts rather than real recombination.

* Be cautious with interpreting telomeric regions—they can be harder to genotype reliably and shouldn't be generalized from limited data.

### Haplotype Refinement and Smoothing: Key Takeaways

* Defining proximal/distal SNPs per haplotype block and applying MAF > 20% greatly improved genotype clarity, enabling reliable identification of long haplotype runs.

* Spurious recombinations often result from either mispositioned markers or poor-quality SNPs—these should be flagged for exclusion or correction.

* Recombination smoothing is now the priority: count and analyze per-rat and per-marker recombination flips to identify poor-quality samples (e.g., erratic genotypes like in animal 57241-5) and overly disruptive markers.

### Integrating HS Founders to Improve Genotype Smoothing

* The HS rat population was derived from eight fully inbred progenitor strains, seven of which have been sequenced, enabling precise identification of most haplotype origins.

* Each HS rat has undergone approximately 1,440 meiotic recombinations over 80 generations, leading to expected haplotype block sizes of ~1.5–2.0 Mb.

* Genotype files currently reflect both the original founder SNP variants and the new recombination-derived haplotypes—these layers need to be disentangled.

* Separating founder variants from HS-specific recombinations will reduce mapping noise and enhance power in QTL analyses.

* Including the 8 founder genotypes (as reference columns) in the current genotype matrix will allow alignment of observed haplotypes with founder origins.

* Mapping should first use clean HS haplotypes for detection, and then refine with founder SNPs for fine-mapping resolution.

### Why Mapping HS Rats to Founders Matters

* Each HS rat genome is a mosaic of haplotype blocks inherited from eight fully inbred founder strains, and can be (in principle) traced back to founder origins.

* Founder-specific haplotypes, especially from genetically distinct strains like BN, are easier to identify due to unique SNP signatures.

* Dense SNP data would ideally allow “color-coding” each chromosome segment by founder origin, improving interpretability and mapping resolution.

* Fully phased genotypes (assigning parent-of-origin to each haplotype) provide the highest power for trait mapping—but are often unavailable in HS data.

* Founder haplotypes carry complex structural variants (e.g., mobile elements) that SNPs alone may miss, enabling richer variant discovery beyond SNP-level variation.

## Practical steps to clean, infer and generate haplotypes for HS rats 

### Step 1: QC/QA on the raw data 
* There are many tools used for checking the quality of genotypes before any downstream analysis. Deciding which tool to pick, depends on the file format (vcf, bed, etc), type of downstream analysis, and computational efficiency. 
* In this case, [BCFtools](https://www.htslib.org/doc/1.0/bcftools.html) was selected for its simplicity, efficiency, and multiple options for analysis 
1. Filtering MAF (minor allele frequency) less than 20% 
```sh 
bcftools filter -e 'MAF < 0.2' -Oz -o output.vcf.gz input.vcf.gz 
```
2. Filtering low quality snps 
```sh 
bcftools view -i 'QUAL>20 && DP>10' input.vcf.gz -o filtered.vcf.gz 
```
 * Many other options to explore, however this will depend with one's objectives to filtering the vcf file 

### Step 2: Processing founders 
* The original vcf founder file contains 42 founder strains, which are homozygous. But we only need to select 8 of the 42. 
* The 8 founders include: 
 - ACI/N
 - BN/SsN
 - BUF/N
 - F344/N
 - M520/N
 - MR/N
 - WKY/N
 - WN/N
*  We still use bcftools to select these 8 founders
```sh 
## create a file containing list of the 8 founder strain names 
## use it to filter the strains to 8 founders 
bcftools view -S samples.txt input.vcf -o output.vcf 
```
* Another important thing is to sync the ref genome used to generate founder vcf and the current outbred vcf. Usually, the founder vcf needs to be lifted over to match the upto date ref genome used on the outbred vcf 
* [CrossMap](https://crossmap.sourceforge.net/) is a program used to convert genome coordinates between different assemblies thus, snycing the two vcf references 

* What is needed as input for the uplifting? 
    - assembly.chain/liftover file (containing coordinate conversions) 
    - vcf file to be processed (in this case, the founder file) 
    - up-to-date reference genome 
* The following is run to generate an liftover founder vcf that matches the latest reference genome, in terms of genome coordinates
```sh 
CrossMap vcf rn5ToRn7.over.chain ./HS_rats_8-founders_snps_homozygous.filtered.vcf.gz GCF_015227675.2_mRatBN7.2_genomic.chr.fa ./HS_rats_8-founders_snps_homozygous.filtered_lifted.vcf.gz
``` 
 * [More on liftover concept](https://genome.ucsc.edu/cgi-bin/hgLiftOver).  

### Step 3: Phasing founders and outbred vcfs 
* This step involved using beagle, to perform statistical phasing to estimate parents of origin for each allele type and/or snp per chromosome 

* This does not end here, as in the next steps, we need to design a way to infer and separate founder confoundings on the outbred snps. 
```sh 
java -Xmx4g -jar ./beagle.27Feb25.75f.jar gt=./input{founder/outbred}.vcf.gz out=./output_phased impute=true nthreads=10
``` 

### Step 4: Dealing with founder confounding on the outbred HS rats snps 
* The HS outbred snps are a genetic mosaic of the 8 founder strains, after 80 generations of crossing. So for mapping reasons, we need to separate the founder layers in the outbred snps, as these present confoundings that are technically noise on one aspect, though useful when fine mapping is involved 

* The scripts contain two parts; 
> * The first part is the core part that actually generates snps independent of the founder's confounding effect 
```sh 
(base) $ ./distinctive_snps.py -h 
usage: distinctive_snps.py [-h] --hs-vcf HS_VCF --founder-vcf FOUNDER_VCF --chrom
                           CHROM [--output-dir OUTPUT_DIR]
                           [--min-diff-pct MIN_DIFF_PCT]
                           [--snp-call-rate SNP_CALL_RATE] [--skip-snp-filtering]

Identify distinctive SNPs in HS vs Founder VCFs

options:
  -h, --help            show this help message and exit
  --hs-vcf HS_VCF       VCF file for HS samples (gzipped)
  --founder-vcf FOUNDER_VCF
                        VCF file for Founder samples (gzipped)
  --chrom CHROM         Chromosome to process
  --output-dir OUTPUT_DIR
                        Directory for outputs
  --min-diff-pct MIN_DIFF_PCT
                        Min percent founders differing to call distinctive SNP
  --snp-call-rate SNP_CALL_RATE
                        Min SNP call rate (0-1). Use 0 to disable
  --skip-snp-filtering  Skip SNP-level call rate filtering


```
> * Source code is found in this link:  [distinctive_snps.py](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/HS_rats/codes/distinctive_snps.py)

> * The second part is the one that splits the founders and hs vcf files to be processed chromosome-wise, then merged at the end of the process; 
```sh 
(base) $ ./inferring_pipeline.py -h 
usage: inferring_pipeline.py [-h] --hs-vcf HS_VCF --founder-vcf FOUNDER_VCF
                             --chromosomes CHROMOSOMES [CHROMOSOMES ...]
                             --output-dir OUTPUT_DIR [--hs-prefix HS_PREFIX]
                             [--founder-prefix FOUNDER_PREFIX]

Run SNP filtering pipeline with runtime/memory monitoring.

options:
  -h, --help            show this help message and exit
  --hs-vcf HS_VCF       Input HS VCF file (gzipped).
  --founder-vcf FOUNDER_VCF
                        Input Founder VCF file (gzipped).
  --chromosomes CHROMOSOMES [CHROMOSOMES ...]
                        List of chromosomes (e.g 1 2 3).
  --output-dir OUTPUT_DIR
                        Output directory.
  --hs-prefix HS_PREFIX
                        Prefix for HS chromosomes (default: 'chr').
  --founder-prefix FOUNDER_PREFIX
                        Prefix for Founder chromosomes (default:'').

```
> * Source code is found in this link: [inferring_pipeline.py](https://github.com/fetche-lab/GeneNetwork_24L/blob/main/Data_processing/HS_rats/codes/inferring_pipeline.py)
 
### Step 5: Generating Haplotype blocks 
* With PLINK, it was possible to generate marginal haplotype blocks for the new outbred vcf file free of founders confounding 
* The scripts used:
 
```sh 
plink --vcf filtered_hs.vcf.gz --blocks no-pheno-req --blocks-max-kb 500 --blocks-strong-lowci 0.60 --blocks-strong-highci 0.95 --out HS_haplo_blocks
 
```

* Filter the resulting haploblocks:

```python 
#!/usr/bin/env python3 

import pandas as pd 

blocks = pd.read_csv("10_4-10-20.blocks.det", sep=r"\s+", engine="python")

#inspect 
print(blocks.head())

#basic qc summaries 
print("Number of blocks:", len(blocks))
print("Median block size (kb):", blocks["KB"].median())
print("Max block size (KB):", blocks["KB"].max())
print("Median number of SNPs per block:", blocks["NSNPS"].median())

# Filtering criteria 
qc_blocks = blocks[ 
    (blocks["KB"] <= 2000) & # no more than 2mb  
    (blocks["NSNPS"] >= 3)  # at least 3 snps 
]

print("Blocks kept after QC:", len(qc_blocks))

# save the filtered list 
qc_blocks.to_csv("./filtered_10_4-10-20_blocks.tsv", sep="\t", index=False)

```
* Generate haplotype edges only from the filtered  `blocks.det` file  

```python 

#!/usr/bin/env python3

import pandas as pd

# Load your blocks.det file
blocks = pd.read_csv("/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/10_4-all.blocks.txt", sep="\t")

# Combine BP1 and BP2 into one series of edges
#edges = pd.concat([blocks["BP1"], blocks["BP2"]]).reset_index(drop=True)

edges_list = [] 
for _, row in blocks.iterrows():
    edges_list.append({"CHR":row["CHR"], "haploEdges": row["BP1"]}) 
    edges_list.append({"CHR":row["CHR"], "haploEdges": row["BP2"]})

# convert to dataframe 
#edges_df = pd.DataFrame({"haploEdges": edges})
edges_df = pd.DataFrame(edges_list) 

# Inspect
#print(edges_df.head())

edges_df = edges_df.reset_index(drop=True) 

# Save to file
edges_df.to_csv("/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/10_4-all.edges.txt", sep="\t", index=False)

```


### Step 6: Generate bimbam files and run gemma 
* [Gemma](https://github.com/genetics-statistics/GEMMA) runs of either plink or bimbam format 
* In this case, we opt to tranform our files into bimbam format. 
1. Bimbam genotype file and bimbam snp file

```python 
#!/usr/bin/env python3

import pandas as pd

# === Input/Output Paths ===
GENO_FILE = "/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/10_4-all.haplotypes.txt"
BLOCKS_FILE = "/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/10_4-all.edges.txt"

OUT_GENOTYPES = "/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/genotypes/bimbam/10_4-all.bimbam.haplotypes1.txt"
OUT_SNPS = "/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/genotypes/bimbam/10_4-all.bimbam.snps1.txt"

# === Load Haplo Edges ===
blocks = pd.read_csv(BLOCKS_FILE, sep="\t")
haplo_edges = set(blocks["haploEdges"].astype(int))  # convert to set for faster lookup

# === Define genotype code mapping ===
geno_map = {"0|0": 1, "0|1": 2, "1|0": 2, "1|1": 3}

# === Helper: format SNP IDs ===
def format_marker(marker):
    if isinstance(marker, str) and ":" in marker:
        chrom, pos = marker.split(":")
        try:
            return f"hsr{int(pos):09d}"
        except ValueError:
            return marker
    return marker

# === Process in Chunks ===
chunksize = 10000  # adjust depending on memory (100k rows at a time)

# Prepare output files (overwrite if they exist)
open(OUT_GENOTYPES, "w").close()
open(OUT_SNPS, "w").close()

for i, chunk in enumerate(pd.read_csv(GENO_FILE, sep="\t", chunksize=chunksize)):
    print(f"[INFO] Processing chunk {i+1}...")

    # Filter rows by haploEdges
    chunk = chunk[chunk["POS"].isin(haplo_edges)]
    if chunk.empty:
        continue

    # # Replace genotype encodings
    geno_cols = chunk.columns.difference(["CHROM", "POS", "ALT", "REF", "ID"])
    # chunk[geno_cols] = chunk[geno_cols].replace(geno_map)

   
   # Replace genotype encodings safely
    chunk.loc[:, geno_cols] = (
        chunk[geno_cols].replace(geno_map).infer_objects(copy=False)
    )

    # Format IDs safely
    chunk.loc[:, "ID"] = chunk["ID"].apply(format_marker)


    # # Reformat SNP IDs
    # chunk["ID"] = chunk["ID"].apply(format_marker)

    # === BIMBAM genotypes ===
    col_order = ["ID", "ALT", "REF"] + [col for col in chunk.columns if col not in ["CHROM", "POS", "ALT", "REF", "ID"]]
    bimbam_genotypes = chunk[col_order]

    # Append to file (no header, no index)
    bimbam_genotypes.to_csv(OUT_GENOTYPES, sep="\t", header=False, index=False, mode="a")

    # === BIMBAM SNPs ===
    bimbam_snps = chunk[["ID", "POS", "CHROM"]].copy()
    bimbam_snps["CHROM"] = bimbam_snps["CHROM"].str.replace("chr", "", regex=False)

    bimbam_snps.to_csv(OUT_SNPS, sep="\t", header=False, index=False, mode="a")

print("[INFO] Finished processing all chunks!")


```

* Now we run gemma 

> * Calculate kinship 
```sh 
./gemma-0.98.5-linux-static-AMD64 -g 10_4-all.bimbam.genotypes.txt -p hs_phenotypes_bimbam.txt -gk -o 10_4.all.kinship 

```
> * Calculate LMM (univariate linear mix model)
```sh 
./gemma-0.98.5-linux-static-AMD64 -g 10_4-all.bimbam.haplotypes1.txt -p hs_phenotypes_bimbam.txt -a 10_4-all.bimbam.snps1.txt -k output/10_4.all.kinship.cXX.txt -lmm -o 10_4.all.kinship.lmm 

```

* Now we generate plots 

```R

#!/usr/bin/env Rscript

# -------------------------
# Setup
# -------------------------
setwd("/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/genotypes/bimbam/output/")

if(!require(qqman)) install.packages("qqman", repos='https://cloud.r-project.org')
if(!require(dplyr)) install.packages("dplyr", repos='https://cloud.r-project.org')
library(qqman)
library(dplyr)

# -------------------------
# Function to load + process GWAS results
# -------------------------
process_gwas <- function(file, label) {
  df <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  df <- df %>%
    mutate(
      CHR = as.numeric(gsub("chr", "", chr)),
      BP = as.numeric(ps),
      SNP = rs,
      P = as.numeric(p_wald)
    ) %>%
    filter(!is.na(P) & P > 0) %>%
    arrange(CHR, BP)
  
  df$FDR <- p.adjust(df$P, method = "fdr")
  df$logp <- -log10(df$P)
  df$fdr_thr <- max(df$P[df$FDR < 0.05], na.rm = TRUE)
  df$dataset <- label
  return(df)
}

# -------------------------
# Load both datasets
# -------------------------
smoothed <- process_gwas("10_4.all.smoothed1.lmm.assoc.txt", "Smoothed haplotypes")
unsmoothed <- process_gwas("10_4.all.unsmoothed.lmm.assoc.txt", "Unsmoothed genotypes")

# -------------------------
# Shared y-limit and threshold
# -------------------------
combined <- bind_rows(smoothed, unsmoothed)
y_max <- ceiling(max(combined$logp, na.rm = TRUE))
shared_fdr_thr <- -log10(min(unique(c(smoothed$fdr_thr, unsmoothed$fdr_thr)), na.rm = TRUE))

cat("Shared FDR line at -log10(p) =", round(shared_fdr_thr, 2), "\n")
cat("Shared Y-axis limit =", y_max, "\n")

# -------------------------
# Plot smoothed
# -------------------------
#pdf("compared_manhattan/HS_smoothed_FDR_qqman.pdf", width=12, height=6)
png("compared_manhattan/HS_smoothed_FDR_qqman.png", 
    width = 3000, height = 1500, res = 300)  # High-resolution

manhattan(
  smoothed,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  main = "HS smoothed haplotypes (FDR-adjusted)",
  col = c("#4C72B0", "#55A868"),
  cex = 0.3, cex.axis = 0.8,
  genomewideline = shared_fdr_thr, # same for both
  ylim = c(0, y_max),
  suggestiveline = FALSE
)
dev.off()

# -------------------------
# Plot unsmoothed
# -------------------------
#pdf("compared_manhattan/HS_unsmoothed_FDR_qqman.pdf", width=12, height=6)
png("compared_manhattan/HS_unsmoothed_FDR_qqman.png", 
    width = 3000, height = 1500, res = 300)  # High-resolution
manhattan(
  unsmoothed,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  main = "HS unsmoothed genotypes (FDR-adjusted)",
  col = c("#4C72B0", "#55A868"),
  cex = 0.3, cex.axis = 0.8,
  genomewideline = shared_fdr_thr,
  ylim = c(0, y_max),
  suggestiveline = FALSE
)
dev.off()

cat("\nBoth plots saved with identical y-axis scaling and FDR threshold line.\n")

```

### Step 7: Generate genotypes files for gn2 format 
* Generate a simple genotype file from the outbred vcf with distinctive snps 
```sh 
## Generate genotype file from the vcf 
(echo -e "CHROM\tPOS\tID$(bcftools query -l filtered_hs_chr1.vcf.gz | awk '{printf "\t%s", $0}')";  bcftools query  -f '%CHROM\t%POS\t%ID[\t%GT]\n' filtered_hs_chr1.vcf.gz) > hs_genotypes.txt
```

* Generate final genotype file using the haplo edges as inclusive markers 

```python
#!/usr/bin/env python3

import pandas as pd

# Input files
geno = "/path_to_genotype_file/genotypes.txt"
edges = "/path_to_edges_file/edge_markers.txt"
hs_output = "/path_to_the_output_file/genotypes.geno"

# Load haplotype block positions into sets for fast lookup
blocks = set(pd.read_csv(edges, sep="\t")["haploEdges"].values)

# Mapping genotypes
geno_map = {'0|0': 'A', '0|1': 'H', '1|0': 'H', '1|1': 'B'}

# Metadata header
hs_metadata = [
    '## This file has genotype data representing HS (Heterogeneous Stock) rats from Abe’s group (RatGTEx)',
    '## Two version of the file. round10_4 represents rats in liver and adipose cohorts. round10_5 represents rats in Basolateral amygdala, Brain hemisphere, Eye, Infralimbic cortex,  Lateral habenula, Nucleus accumbens core, Orbitofrontal cortex, Prelimbic cortex, Posterior ventral tegmental area, and Rostromedial tegmental nucleus cohorts', 
    '## There are chromosomes 1 to 20',
    '## Rat assembly used is rn7 (mRatBN7.2) genome assembly',
    '## `ref` = homozygous reference allele',
    '## `alt` = homozygous alternative allele',
    '## `het` = heterozygous alleles (ref + alt)',
    '## `unk` = unknown genotype',
    '## Marker IDs formatted as {hsr + position}',
    '@name: HS-RATS-RatGTEx',
    '@type: Heterogeneous Stock (HS) cross',
    '@ref:A',
    '@alt:B',
    '@het:H',
    '@unk:U',
    ''
]

# Write metadata first
with open(hs_output, "w") as f:
    for line in hs_metadata:
        f.write(line + "\n")

# Helper function to process a chunk
def process_chunk(chunk, keep_positions):
    # Filter rows
    chunk = chunk[chunk["POS"].isin(keep_positions)]

    if chunk.empty:
        return None

    # Replace genotypes
    geno_cols = chunk.columns.difference(['CHROM', 'POS', 'ID'])
    chunk.loc[:, geno_cols] = chunk[geno_cols].replace(geno_map)

    # Reformat marker IDs
    def format_marker(marker):
        if isinstance(marker, str) and ":" in marker:
            chrom, pos = marker.split(":")
            try:
                return f'hsr{int(pos):09d}'
            except ValueError:
                return marker
        return marker
    chunk.loc[:, "ID"] = chunk["ID"].apply(format_marker)

    # Clean CHROM and scale POS
    chunk.loc[:, "CHROM"] = chunk["CHROM"].str.replace("chr", "", regex=False)
    chunk.loc[:, "POS"] = chunk["POS"] / 1_000_000

    # Rename columns
    chunk = chunk.rename(columns={"CHROM": "Chr", "POS": "Mb", "ID": "Locus"})
    chunk.insert(2, "cM", chunk["Mb"])

    # Reorder columns
    new_order = ["Chr", "Locus", "Mb", "cM"] + [c for c in chunk.columns if c not in ["Chr", "Locus", "Mb", "cM"]]
    return chunk[new_order]

# Process both genotype files in chunks
for geno_file, keep_positions in (geno, blocks):
    for chunk in pd.read_csv(geno_file, sep="\t", chunksize=50000):
        processed = process_chunk(chunk, keep_positions)
        if processed is not None:
            processed.to_csv(hs_output, sep="\t", index=False, mode="a", header=True)

```

