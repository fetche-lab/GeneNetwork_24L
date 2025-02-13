## **CELEGANS (*Caenorhabditis elegans*)**
### MAGIC (Multiparent Advanced Generation Intercross) Lines 

Metadata on Celegans dataset in GeneNetwork by Snoek et al. 

Publication Title: "A multi-parent recombinant inbred line
population of C. elegans allows identificationof novel QTLs for complex life history traits."

[Publication and supplementary access](https://doi.org/10.1186/s12915-019-0642-8)

#### Summary Overview 

The nematode *Caenorhabditis elegans* is widely used to study genetic variation in complex traits through QTL analysis of recombinant inbred lines (RILs). To investigate the impact of local genetic diversity, researchers developed a new multi-parental RIL population from four wild isolates collected in France, identifying nearly 9000 SNPs and doubling mapping resolution. This approach revealed novel QTLs for traits like lifespan and pathogen response, demonstrating that multi-parent RILs improve trait mapping and enhance understanding of local adaptation.

#### Key Findings

Multi-parent RIL populations provide more informative SNP markers than traditional two-parent RILs, improving genetic mapping across organisms. Our study shows that using mpRILs increases the number of detectable QTLs and enhances the resolution for identifying candidate causal SNPs. This approach improves the genetic characterization of complex traits with greater precision.

#### Experimental Design/Methods 

C. elegans strains were cultured on OP50 at 20 °C. A multi-parental RIL was constructed using a crossing scheme (Fig. 1, Table S1) involving four parental isolates (JU1511, JU1941, JU1926, JU1931). To allow recombination and generate extra crossovers, reciprocal inter-crosses were performed, and males were induced by heat stress. F2 worms were then single-worm inbred for 6 generations to create homozygous genotypes. From 383 lines, 200 mpRILs were randomly selected for mRNA sequencing to assess genetic variation.

RNA isolation: 
- C. elegans RNA was isolated from mpRILs and parental strains grown on OP50 at 16°C, then bleached and raised at 24°C for 48 hours8. RNA extraction was performed using the Maxwell® 16 LEV simplyRNA Tissue Kit with a modified lysis step involving proteinase K incubation8. This method, building upon guanidinium thiocyanate-phenol-chloroform extraction14, allowed for high-quality RNA isolation suitable for downstream applications like sequencing

Platform: 
- Illumina HiSeq™ 2000

Access to Raw sequencing data:
- (SRA; www.ncbi.nlm.nih.gov/bioproject/ PRJNA495983/) with ID 
PRJNA495983

#### Datasets 

Genotypes: 
- SNP calling: Untrimmed paired-end reads were mapped to the N2 reference genome (WS220) using Tophat, allowing mismatches and edit distance. SNPs were then called using SAMtools mpileup with bcftools and vcfutils, requiring a minimum read depth, followed by quality filtering to minimize false positives.
- Genetic Map: To construct a C. elegans genetic map from RNA-seq data SNPs, the method from Serin et al. was adapted. SNPs present in parental lines with high quality (quality > 199) and a minimum presence in mpRILs were selected, alongside those correlating with neighboring SNPs. Further filtering based on heterozygosity was conducted, and the remaining SNPs (Table S2) were used to create both a direct SNP map and a parental origin genetic map (Tables S3, S4), determining parental origin via stretches of ten SNPs. [Access to figures and tables mentioned](https://doi.org/10.1186/s12915-019-0642-8)

Phenotypes: 
- The 200 mpRILs and parental strains were phenotyped for a diverse set of life-history traits, reflecting natural conditions. Measurements were taken for length, width, length-to-width ratio, volume, lifespan (both standard and under dietary restriction), survival under heat shock and oxidative stress, the occurrence of males, developmental speed (assayed on E. coli OP50 and Erwinia rhapontici), and population growth (on E. coli OP50, Erwinia rhapontici, Sphingomonas sp., Bacillus thuringiensis DSM-350E, and pathogenic B. thuringiensis NRRL B-18247). Assays were conducted under varied bacterial food and abiotic conditions. 

#### Analyses 

QTL Mapping: 
- QTL analysis was performed to identify genetic loci associated with the observed phenotypic variation. Broad peaks were further refined by testing for associations between phenotypic variation and the parental origin at a given location in the mpRIL genomes to identify smaller, better-resolved QTLs. They did not account for epistatic interactions among the loci.

Statistical Analysis: 
- A standard interval mapping approach was employed using a Hidden Markov Model (HMM) to calculate the probability of each mpRIL inheriting each of the four parental alleles at 1 cM intervals across the genome. The association between the posterior probabilities and the phenotypic data was then tested using a likelihood ratio test (LRT). The significance thresholds were determined using a permutation approach. 

#### Acknowledgements 

We thank Marie-Anne Félix for providing the parental lines used to make the
mpRILs. We thank the compilers of WormBase for making it a versatile and
important resource for C. elegans.

Funding: 
- We acknowledge financial support from the Deutsche Forschungsgemeinschaft
to HS, grant number SCHU 1415/11 and project Al within the CRC 1182.
Furthermore, financial support from the NWO-ALW (project 855.01.151) to RJMV




