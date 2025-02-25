# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


library(tidyverse)
library(cowplot)
library(ggrepel)

# Load Snakemake variables
counts <- snakemake@input[["counts"]]

# Load count data
data <- read.delim(counts)

# Remove barcodes where sum of all numeric columns is zero
data <- data %>%
  filter(rowSums(select_if(., is.numeric)) > 0)

# Identify numeric column with largest sum
max.col <- which.max(colSums(select_if(data, is.numeric)))
max.sum <- colSums(data[names(max.col)])

# Normalise to largest sample:
# divide max.sum by each numeric column sum to get correction factor
# then multiple each numeric column by that correction factor
for (column in names(data)[sapply(data, is.numeric)]) {
  sum_ <- colSums(data[column])
  factor <- max.sum / sum_
  data[column] <- data[column] * factor
}

# Reshape data to get a sample and bin column
data_long <- data %>%
  pivot_longer(cols = -c(barcode_id, orf_id, gene),
               names_to = "variable",
               values_to = "count") %>%
  separate(variable, into = c("sample", "bin"), sep = "_")

# Prepare data for PCA:
# wide format, samples as rows, variables as columns)
# Use pivot_wider, but this time, we create unique sample identifiers
data_wide <- data_long %>%
  unite("sample_id", c(barcode_id, bin), sep = "_", remove = FALSE) %>% #create a unique sample id.
  pivot_wider(id_cols = sample_id, names_from = sample, values_from = count) %>%
  select(-sample_id) #Remove now the sample id

# Remove rows if all values across all conditions are zero
data_wide <- data_wide %>%
  filter(rowSums(select_if(., is.numeric)) > 0)

# Perform PCA
pca_result <- prcomp(data_wide, scale. = TRUE, center = TRUE)
summary(pca_result)

df <- pca_result$rotation %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")

# Create plot
colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
             "#FFFF33", "#A65628", "#F781BF", "#999999")
p <- ggplot(df, aes(x = PC1, 
                    y = PC2, 
                    color = Sample)) +
  geom_point(size = 8) +
  theme_cowplot(18) +
  labs(title = "PCA of barcode profiles",
       x = paste0("PC1 (", round(pca_result$sdev[1] / sum(pca_result$sdev) * 100, 2), "%)"),
       y = paste0("PC2 (", round(pca_result$sdev[2] / sum(pca_result$sdev) * 100, 2), "%)")) +
  scale_color_manual(values = colours[1:length(df$Sample)]) +
  geom_label_repel(data = df, 
                  aes(label = Sample), 
                  size = 5,) +
  theme(legend.position = "none")

# Save plot
ggsave(snakemake@output[["pdf"]], p, width = 6, height = 4)

# Close log file
sink(log, type = "output")
sink(log, type = "message")