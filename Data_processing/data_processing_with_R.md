## ANNOTATING AFFYMETRIX MICROARRAY DATASETS USING AVAILABLE R ANNOTATION PACKAGES 
### The following contains scripts prepared in R language with annotation packages that were used in processing few of the selected microarray datasets taken from NCBI GEO repositories 
### For experimentation purposes, `C elegans` with accession number `GSE163253`, and `A thaliana` with accession number `GSE71188` were used 

### R script for C elegans 
```R
# annotating c elegans microarray dataset 
# case study dataset - GSE163253

# experiment the packages 
library(BiocManager)
BiocManager::install('celegans.db')
BiocManager::install('celeganscdf')
BiocManager::install('celegansprobe')
BiocManager::install('pd.celegans', force = TRUE)
BiocManager::install('affy')
BiocManager::install('GEOquery')
BiocManager::install(c('gcrma', 'dplyr','tibble'))

# load the packages  
library(celegans.db) # gene level annotation package 
library(celeganscdf) # chip description file package 
library(affy)
library(GEOquery)
library(celegansprobe) # probe level annotation package 
library(pd.celegans)
library(gcrma)
library(dplyr)
library(tibble)

#start the experimentation analysis, to normalize, transform, and annotate data 
#set working environment 
setwd('/home/fetche/Desktop/Software_dev/Pjotr_software_team/Projects_2024/data_engineering/data/C_elegans_data/GSE163253_dataset/')

#unpack the CEL files 
untar("GSE163253_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

#set the working directory to data folder 
setwd('/home/fetche/Desktop/Software_dev/Pjotr_software_team/Projects_2024/data_engineering/data/C_elegans_data/GSE163253_dataset/data/')

#read the dataset using the chip description file for celegans microarray data 
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="celegans")

#perform rma normalization 
data.rma.norm <- gcrma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma.data.curated=exprs(data.rma.norm)

#perform log2 transformation 
rma.data.log2 = log2(rma.data.curated)

#Format values to 6 decimal places
rma.data.curated1=format(rma.data.log2, digits=7)
print(rma.data.curated1)

#Map probe sets to gene symbols or other annotations
#check the available mapping options in the annotation packages 
ls('package:celegans.db')
columns(celegans.db)

#extract the rownames as probes 
Probe_IDs = row.names(rma.data.curated1)

#extract important annotation information from the package(s)
Gene_Symbols <- unlist(mget(probes, celegansSYMBOL, ifnotfound = NA))
Gene_Names <- unlist(mget(probes, celegansGENENAME, ifnotfound = NA))
Entrez_IDs <- unlist(mget(probes, celegansENTREZID, ifnotfound = NA))
Ensembl_IDs<- unlist(mget(probes, celegansENSEMBL, ifnotfound = NA))
Refseq_IDs <- unlist(mget(probes, celegansREFSEQ, ifnotfound = NA))
CHR <- unlist(mget(probes, celegansCHR, ifnotfound = NA))
CHR_start <- unlist(mget(probes, celegansCHRLOC, ifnotfound = NA))
CHR_end <- unlist(mget(probes, celegansCHRLOCEND, ifnotfound = NA))
CHR_length <- unlist(mget(probes, celegansCHRLENGTHS, ifnotfound = NA))
PubMed_IDs <- unlist(mget(probes, celegansPMID, ifnotfound = NA))
UniProt_IDs <- unlist(mget(probes, celegansUNIPROT, ifnotfound = NA))



#Combine gene annotations with raw data
rma=cbind(probes,Gene_Symbols,Gene_Names,Entrez_IDs,Ensembl_IDs,
          Refseq_IDs,CHR,CHR_start,CHR_end,PubMed_IDs,UniProt_IDs,
          celegansprobe$sequence,rma.data.curated1)

#Write RMA-normalized, mapped data to file
write.table(rma, file = "GSE163253_celegans_rma_annot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) # annotation file
write.csv(rma.data.curated2, file = "GSE163253_log2_data.csv", quote = FALSE,row.names = FALSE, col.names = TRUE) # dataset file 

```

### R script for A thaliana 
```R
# annotating A thaliana microarray dataset 
# case study dataset - GSE71188
# genechip - Affymetrix Arabidopsis ATH1 Genome Array
BiocManager::install('ath1121501.db', force = TRUE)

# load packages 
library(BiocManager)
library(ath1121501.db) # gene level annotation package 
library(ath1121501probe) # probe level annotation package 
library(ath1121501cdf) # chip description file package 
library(org.At.tair.db) # organism level annotation package 
library(affy)
library(GEOquery)
library(gcrma)
library(dplyr)
library(tibble)


# Search for available packages related to 'arabidopsis'
available_packages <- BiocManager::available()
arabidopsis_packages <- grep("arabidopsis", available_packages, ignore.case = TRUE, value = TRUE)
print(arabidopsis_packages)

# start the experimentation analysis to normalize, transform, and annotate the data 
#set working environment 
setwd('/home/fetche/Desktop/Software_dev/Pjotr_software_team/Projects_2024/data_engineering/data/Arabidopsis_data/GSE71188_dataset/')

#unpack the CEL files 
untar("GSE71188_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

#set the working directory to data folder 
setwd('/home/fetche/Desktop/Software_dev/Pjotr_software_team/Projects_2024/data_engineering/data/Arabidopsis_data/GSE71188_dataset/data/')

#read the dataset using the chip description file for arabidopsis microarray data 
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="ath1121501")

#perform rma normalization 
data.rma.norm <- gcrma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma.data.expr=exprs(data.rma.norm)


#perform log2 transformation 
rma.data.log2 = log2(rma.data.expr)

#Format values to 6 decimal places
rma.data.curated1=format(rma.data.log2, digits=7)

#print(rma.data.curated1)

#extract the rownames as probes 
Probe_IDs = row.names(rma.data.curated1)

#extract important annotation information from the package(s)
Gene_Symbols <- unlist(mget(Probe_IDs, ath1121501SYMBOL, ifnotfound = NA))
Gene_Names <- unlist(mget(Probe_IDs, ath1121501GENENAME, ifnotfound = NA))
Entrez_IDs <- unlist(mget(Probe_IDs, ath1121501ENTREZID, ifnotfound = NA))
CHR <- unlist(mget(Probe_IDs, ath1121501CHR, ifnotfound = NA))
CHR_start <- unlist(mget(Probe_IDs, ath1121501CHRLOC, ifnotfound = NA))
CHR_end <- unlist(mget(Probe_IDs, ath1121501CHRLOCEND, ifnotfound = NA))
PubMed_IDs <- unlist(mget(Probe_IDs, ath1121501PMID, ifnotfound = NA))

#Combine gene annotations with raw data
rma1=cbind(Probe_IDs,Gene_Symbols,Gene_Names,Entrez_IDs,CHR,CHR_start,CHR_end,PubMed_IDs,
           ath1121501probe$sequence,rma.data.curated1)


#Write RMA-normalized, mapped data to file
write.table(rma1, file = "GSE71188_arabid_annot.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) # annotation file 
write.csv(rma.data.curated1, file = "GSE71188_arabid_log2_data.csv", quote = FALSE, col.names = TRUE) # dataset file 

```
