## The following are scripts used to process information from GeneNetwork on mice, for Phewas analysis 
### STEP01: understanding my data 
install.packages("dplyr")
library(dplyr)

data <- read.csv("/home/fetche-lab/Desktop/phd_bioinformatics/PheWAS_Project/Data/Raw/phewas01/BTBD9_MICE/AKXD_traits.csv", header=TRUE, skip=5, stringsAsFactors=FALSE)
colnames(data) <- col_names  # Apply manually if 62 columns match

unique(data$Symbol)  # Should be "Btbd9"T
unique(data$Description)  # Check RLS mentions
unique(data$Locus_at_Peak)  # Unique variants
summary(data$P_value_of_MAX)  # P-value range

### Step 02: Extract key columns for Phewas analysis 
#### The columns include; Record_ID, Locus_at_Peak, and P_value_of_MAX) 
phewas_data <- data %>%
  select(Record.ID, Locus_at_Peak, P_value_of_MAX)

### Step 03: Filter significant hits 
# Filter significant associations (p < 0.05)
sig_hits <- phewas_data %>%
  filter(P_value_of_MAX < 0.05)

### lets plot our results and see 
library(ggplot2)

# Basic PheWAS plot
ggplot(phewas_data, aes(x=Locus_at_Peak, y=-log10(P_value_of_MAX), color=Record.ID)) +
  geom_point(size=2) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="PheWAS for Btbd9 (AKXD)", x="QTL Locus", y="-log10(P-value)")

ggsave("phewas_btbd9_akxd.pdf", width=10, height=6)
