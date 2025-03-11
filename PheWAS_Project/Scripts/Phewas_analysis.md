## Phewas Analysis on Mice/Rats data in GeneNetwork 
The following contains steps involved in carrying out phewas (phenome-wide association study) analysis on mice and rats data from GN2. 

### A bit of context 
- The analysis involves extracting phenotype/trait information associated with a particular genetic variant or set of variants of interest, then study how strongly associated they are. So, its reverse GWAS (Genome-Wide Association Study)

### Steps involved in this analysis 
01. Extracting relevant records from genenetwork using gene names provided by source, then quering them on global search, specifying the species of interest (in this case, mice and rats)

02. Filtering traits that are associated with the studied condition of interest, in this case, `restless leg syndrome`. 

  - The script used; 

  ```R 
    library(dplyr)

    input_dir <- "/home/fetche-lab/Desktop/phd_bioinformatics/PheWAS_Project/Data/Raw/phewas01/BTBD9_MICE"
    output_dir <- "/home/fetche-lab/Desktop/phd_bioinformatics/PheWAS_Project/Data/Processed"
    dir.create(output_dir, showWarnings=FALSE)

    csv_files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)

    process_file <- function(file_path) {
    data <- read.csv(file_path, skip=5, header=TRUE, stringsAsFactors=FALSE)
    data$Record.ID <- as.character(data$Record.ID)
    data$Locus_at_Peak <- as.character(data$Locus_at_Peak)
    
    phewas_data <- data %>%
        filter(Symbol == "Btbd9" & grepl("restless leg syndrome", Description, ignore.case=TRUE)) %>%
        select(Record.ID, Locus_at_Peak, P_value_of_MAX)
    
    output_file <- file.path(output_dir, paste0("phewas_", tools::file_path_sans_ext(basename(file_path)), ".csv"))
    write.csv(phewas_data, output_file, row.names=FALSE)
    
    message(paste("Processed", basename(file_path), ":", nrow(phewas_data), "rows"))
    return(phewas_data)
    }

    results <- lapply(csv_files, process_file)
    combined_data <- bind_rows(results)
    if (nrow(combined_data) > 0) {
    write.csv(combined_data, file.path(output_dir, "phewas_btbd9_combined.csv"), row.names=FALSE)
    print(head(combined_data))
    }

    # Filter significant associations (p < 0.05)
    sig_hits <- combined_data %>%
    filter(P_value_of_MAX < 0.05)

    # save the file with filtered hits
    write.csv(sig_hits, "phewas_btbd9_combined-filtered.csv", row.names=FALSE)


  ```
  - In summary, the R script above did the following; 
    - Each Genename (from source) had a set of groups for mice, and in each group there were multiple files containing trait records and their hits associated with the genetic variants representing that gene. 
    - The script extracts traits, variants and significant hits (as pvalues).
    - It then filters out traits whose pvalue is more than 0.5 and retains the ones that are less than 0.5 as these reflect significant association to the variants. 
    - This will then be applied to the rest of the groups. 

