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

write.csv(sig_hits, "phewas_btbd9_combined-filtered.csv", row.names=FALSE)


# PheWAS plot for combined data (Optional)
library(ggplot2)
ggplot(sig_hits, aes(x=Locus_at_Peak, y=-log10(P_value_of_MAX), color=Record.ID)) +
  geom_point(size=2) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  labs(title="PheWAS for Btbd9 Across Datasets", x="QTL Locus", y="-log10(P-value)")
ggsave(file.path(output_dir, "phewas_btbd9_combined.pdf"), width=12, height=8)
