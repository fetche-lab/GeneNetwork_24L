#!/usr/bin/env Rscript

# ===============================================================
# Ensembl Annotation for Rat Gene IDs (BXD / HS compatible)
#  - Handles mixed gene names, LOC IDs, and RGD IDs
#  - One-to-one mapping, preserves all input rows
# ===============================================================

setwd("/home/fetche-lab/Desktop/gn_remote/HS_rats/phenotypes_hs_arthur/expression/")

# --- Setup packages ---
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt", repos="https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos="https://cloud.r-project.org")

library(biomaRt)
library(dplyr)
library(readr)

# --- Input file ---
gene_file <- "gene_ids_liver.txt"
genes <- read_tsv(gene_file, col_names = "gene_id", show_col_types = FALSE)

cat("\n Loaded", nrow(genes), "genes from file:", gene_file, "\n")

# --- Connect to Ensembl (Rat) ---
cat("Connecting to Ensembl biomart (rat)...\n")
ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")

# --- Define the annotation fields we want ---
attribs <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "external_synonym",
  "rgd_symbol",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand",
  "description"
)

# --- Split input by ID type ---
rgd_ids <- unique(genes$gene_id[grepl("^RGD", genes$gene_id, ignore.case = TRUE)])
non_rgd_ids <- unique(genes$gene_id[!genes$gene_id %in% rgd_ids])

cat("Non-RGD IDs :", length(non_rgd_ids),
    "\nRGD IDs     :", length(rgd_ids), "\n")

# --- Query Ensembl: external_gene_name ---
cat("Fetching annotations by gene name / LOC...\n")
annot_gene <- getBM(
  attributes = attribs,
  filters = "external_gene_name",
  values = non_rgd_ids,
  mart = ensembl
)

# --- Query Ensembl: RGD IDs ---
annot_rgd <- data.frame()
if (length(rgd_ids) > 0) {
  cat("Fetching annotations by RGD symbol...\n")
  annot_rgd <- getBM(
    attributes = attribs,
    filters = "rgd_symbol",
    values = rgd_ids,
    mart = ensembl
  )
}

# --- Combine and clean ---
#annot_all <- bind_rows(annot_gene, annot_rgd) %>%
 # mutate(
  #  strand = ifelse(strand == 1, "+", ifelse(strand == -1, "-", NA))
  #) %>%
  #distinct()

# --- Harmonize column types before combining ---
common_cols <- intersect(names(annot_gene), names(annot_rgd))

annot_gene <- annot_gene %>%
  mutate(across(all_of(common_cols), as.character))

annot_rgd <- annot_rgd %>%
  mutate(across(all_of(common_cols), as.character))

# --- Combine and clean ---
annot_all <- bind_rows(annot_gene, annot_rgd) %>%
  mutate(
    strand = ifelse(strand == "1", "+", ifelse(strand == "-1", "-", NA))
  ) %>%
  distinct()


# --- Reduce to one-to-one (unique gene symbol) ---
annot_1to1 <- annot_all %>%
  group_by(external_gene_name) %>%
  slice_head(n = 1) %>%
  ungroup()

# --- Merge with input list (preserve all genes) ---
annotated_genes <- genes %>%
  left_join(annot_1to1, by = c("gene_id" = "external_gene_name")) %>%
  rename_with(~ gsub("\\.x$|\\.y$", "", .x)) %>%
  distinct()

# --- Summary ---
total_genes <- nrow(genes)
annotated_count <- annotated_genes %>%
  filter(!is.na(ensembl_gene_id)) %>%
  distinct(gene_id) %>%
  nrow()
missing_count <- total_genes - annotated_count
percent_annotated <- round(annotated_count / total_genes * 100, 2)

cat("\nAnnotation Summary:")
cat("\n   • Total input genes      :", total_genes)
cat("\n   • Successfully annotated  :", annotated_count, sprintf("(%s%%)", percent_annotated))
cat("\n   • Missing / not found     :", missing_count, "\n")

# --- Save results ---
out_file <- sub("\\.txt$", "_annotated_rat_1to1_v4.csv", gene_file)
write_csv(annotated_genes, out_file)

cat("\nAnnotation complete!")
cat("\n   → Output file:", out_file, "\n\n")

