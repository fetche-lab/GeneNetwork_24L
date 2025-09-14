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

* Spurious recombinations (e.g., at Marker 25 and 41) often result from either mispositioned markers or poor-quality SNPs—these should be flagged for exclusion or correction.

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

* So, the following were steps taken to infer genotypes to their founders and estimate how unique or similar they were to the founders. This helped us separate snps independent of founders' confounding 
    - Select common markers from both files, sorts them, and deals with any missing snp  
    - Extracts genotypes from the phased HS  vcf into maternal and paternal alleles (e.g; 0|1 => C|T)  
    - For HS founders vcf, assuming their homozygosity, only picks one allele per SNP (e.g; 0|0 => C) 
    - Compares each HS rat's maternal and paternal alleles to the 8 founders alleles at each SNP 
    - Identifies SNPs where an HS haplotype differs from at least 80% of founders (6 out of 8, set by min_diff_pct=80) 
    - Generates an HS vcf file that only contains distinctive snps from founders (at least 80%) 

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
> * Source code is found in this link: 

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
> * Source code is found in this link: 

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

```python 
#!/usr/bin/env python3 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data (assuming 'chr' and 'ps' columns exist)
df = pd.read_csv(
    "/home/fetche-lab/Desktop/gn_remote/HS_rats/data/processed/hs_qc_qa/phased/simple_processing/data/processed/final_vcfs_10_4/genotypes/bimbam/output/10_4.all.smoothed1.lmm.assoc.txt",
    sep="\t"
)

# Convert p-values to -log10 scale
df["-log10p"] = -np.log10(df["p_wald"])

# Bonferroni genome-wide threshold
bonf_threshold = -np.log10(0.05 / len(df))
suggestive_threshold = -np.log10(1e-5)

# Identify the top SNP
top_snp = df.loc[df["p_wald"].idxmin()]

# Sort by chromosome and position
df["chr"] = df["chr"].astype(str).str.replace("chr", "")  # cleanup if needed
df["chr"] = df["chr"].astype(int)  # ensure numeric sorting
df = df.sort_values(["chr", "ps"])

# Create a cumulative position for plotting
chr_offsets = {}
cumulative_bp = 0
xticks = []
xlabels = []

for chrom in sorted(df["chr"].unique()):
    chr_min = df.loc[df["chr"] == chrom, "ps"].min()
    chr_max = df.loc[df["chr"] == chrom, "ps"].max()
    
    chr_offsets[chrom] = cumulative_bp - chr_min
    cumulative_bp += chr_max
    
    # midpoint for x-axis labeling
    xticks.append((df.loc[df["chr"] == chrom, "ps"].median() + chr_offsets[chrom]))
    xlabels.append(f"{chrom}")

df["pos_cum"] = df.apply(lambda row: row["ps"] + chr_offsets[row["chr"]], axis=1)

# Plot Manhattan
plt.figure(figsize=(18, 9))

colors = ["#4C72B0", "#55A868"]  # alternating colors
for i, chrom in enumerate(sorted(df["chr"].unique())):
    subset = df[df["chr"] == chrom]
    plt.scatter(
        subset["pos_cum"],
        subset["-log10p"],
        c=colors[i % len(colors)],
        s=10,
        alpha=0.7
        #label=f"chr{chrom}" if i < 2 else None  # avoid cluttering legend
    )

# Highlight top SNP
plt.scatter(
    top_snp["ps"] + chr_offsets[top_snp["chr"]],
    top_snp["-log10p"],
    c="red",
    s=50,
    edgecolor="black",
    label=f"Top SNP ({top_snp['rs']})"
)

# Threshold lines
plt.axhline(y=bonf_threshold, color="red", linestyle="--", label="Bonferroni")
plt.axhline(y=suggestive_threshold, color="orange", linestyle="--", label="Suggestive")

# Labels & formatting
plt.xticks(xticks, xlabels, rotation=0)
plt.xlabel("Chromosome")
plt.ylabel("-log10(p-value)")
plt.title("HS original genotypes GWAS Manhattan Plot (20 chromosomes)")
plt.legend()
plt.tight_layout()
plt.show()

```

