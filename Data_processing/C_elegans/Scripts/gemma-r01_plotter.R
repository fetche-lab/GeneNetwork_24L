# Install required packages if not already installed
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load libraries
library(data.table)
library(ggplot2)

# Read the data (assuming it's saved as a tab-delimited file)
# Replace 'your_file.txt' with your actual file name
data <- fread("celegans_kinship_lmm01.assoc.txt", sep = "\t")


# Prepare data for plotting
# Extract chromosome number and convert X to 23
data[, chr_numeric := as.numeric(gsub("X_", "", sub("_.*", "", chr)))]
data[is.na(chr_numeric), chr_numeric := 23]  # X chromosome as 23

# Convert p-values to -log10 scale
data[, neg_log_p := -log10(p_wald)]

# Sort by chromosome and position
setorder(data, chr_numeric, ps)

# Calculate cumulative position
data[, cumulative_pos := 0]
current_pos <- 0
for (chrom in unique(data$chr_numeric)) {
  chrom_data <- data[chr_numeric == chrom]
  data[chr_numeric == chrom, cumulative_pos := ps + current_pos]
  current_pos <- max(data[chr_numeric == chrom]$cumulative_pos)
}

# Calculate chromosome midpoints for x-axis labels
chrom_midpoints <- data[, .(midpoint = (min(cumulative_pos) + max(cumulative_pos)) / 2), by = chr_numeric]
chrom_midpoints[, chr_label := ifelse(chr_numeric == 23, "X", as.character(chr_numeric))]

# Determine a more appropriate threshold
# Check the maximum -log10(p-value) in your data to set a reasonable threshold
max_neg_log_p <- max(data$neg_log_p, na.rm = TRUE)
cat("Maximum -log10(p-value) in data:", max_neg_log_p, "\n")

# Set threshold 
threshold <- 2  # This corresponds to p = 0.01; adjust as needed

# Create the Manhattan plot
p <- ggplot(data, aes(x = cumulative_pos, y = neg_log_p, color = factor(chr_numeric))) +
  geom_point(alpha = 0.65, size = 1) +
  scale_color_manual(values = rep(c("#1f77b4", "#ff7f0e"), length.out = length(unique(data$chr_numeric)))) +
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed", alpha = 0.7) +
  scale_x_continuous(breaks = chrom_midpoints$midpoint, labels = chrom_midpoints$chr_label) +
  labs(x = "Chromosome", y = "-log10(p-value)", title = "Celegans Manhattan Plot") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center the title
  )

# Display the plot
print(p)

# Save the plot as PDF
ggsave("manhattan_plotR.pdf", p, width = 15, height = 6, device = "pdf")