# Install required packages if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load libraries
library(ggplot2)
library(dplyr)

# Step 1: Read the data
# Replace 'your_file.txt' with the actual file path
data <- read.table("", header = TRUE, sep = "\t")

# Step 2: Prepare the data
# Convert p-values to -log10 scale, handle very small p-values to avoid Inf
data <- data %>%
  mutate(
    neg_log10_p = -log10(pmax(p_wald, 1e-300)),  # Avoid -Inf by setting a floor
    position = as.numeric(ps),                   # Ensure position is numeric
    chromosome = factor(chr, levels = unique(chr))  # Factorize chromosome
  )

# For a single chromosome (X), no need for cumulative position, but included for generality
# If multi-chromosome data were present, we'd calculate cumulative positions
data <- data %>%
  arrange(chromosome, position)

# Step 3: Define significance threshold (e.g., p = 0.05 or Bonferroni-corrected)
threshold <- -log10(0.05)  # Simple threshold
# For Bonferroni correction: -log10(0.05 / nrow(data))

# Step 4: Create the Manhattan plot
manhattan_plot <- ggplot(data, aes(x = position / 1e6, y = neg_log10_p)) +
  # Scatter plot with points
  geom_point(aes(color = chromosome), size = 1.5, alpha = 0.6) +
  # Add significance threshold line
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed", size = 0.5) +
  # Customize axes
  labs(
    x = "Position on Chromosome X (Mb)",
    y = "-log10(p-value)",
    title = "Manhattan Plot of GWAS Results (Chromosome X)"
  ) +
  # Customize theme to resemble GeneNetwork2 style
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    legend.position = "none"  # Hide legend since it's only X; adjust if multi-chromosome
  ) +
  # Set x-axis to Mb scale
  scale_x_continuous(breaks = seq(0, max(data$position / 1e6), by = 2)) +
  # Color scale (for multi-chromosome, alternate colors; here just one color)
  scale_color_manual(values = c("X" = "darkblue"))

# Step 5: Display the plot
print(manhattan_plot)

# Optional: Save the plot to a file
ggsave("manhattan_plot.png", manhattan_plot, width = 10, height = 6, dpi = 300)