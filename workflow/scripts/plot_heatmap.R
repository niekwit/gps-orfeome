# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Load gene summary data of all conditions
data <- read_csv(snakemake@input[["csv"]], show_col_types = FALSE)

# Remove rows (ORFs) with all columns starting with
# abs(delta_PSI_mean) < threshold
threshold <- snakemake@wildcards[["ht"]]
print(paste0(
  "Filtering data with absolute delta_PSI_mean threshold: ",
  threshold
))
data <- data %>%
  filter(if_any(starts_with("delta_PSI_mean_"), ~ abs(.) >= threshold))

# Combine orf_id and gene into gene_id
data <- data %>%
  unite("gene_id", c("orf_id", "gene"), sep = "_", remove = TRUE)

# Remove columns that start with z_score_corr_
data <- data %>%
  select(-starts_with("z_score_corr_"))

# Convert gene_id to rownames
data <- data %>% column_to_rownames("gene_id")

# Convert NA values to 0
data[is.na(data)] <- 0

# Prepare colours
max_abs_value <- max(abs(data), na.rm = TRUE)
colour_min <- -max_abs_value
colour_max <- max_abs_value
num_colours_in_palette <- 100
colour_palette <- colorRampPalette(c("navy", "white", "red3"))(
  num_colours_in_palette
)

# Define breaks for the colour palette
custom_breaks <- seq(
  colour_min,
  colour_max,
  length.out = num_colours_in_palette + 1
)

# Remove delta_PSI_mean_ substring from column names
colnames(data) <- gsub("delta_PSI_mean_", "", colnames(data))

# Create heatmap to get clusters
show_row_names <- snakemake@params[["rownames"]]
results <- pheatmap(
  data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  breaks = as.numeric(custom_breaks),
  filename = NA,
  silent = TRUE
)

# Get gene clusters
k <- snakemake@params[["clusters"]]
clusters <- cutree(results$tree_row, k = k)
df <- data.frame(gene_id = rownames(data), cluster = clusters) %>%
  arrange(cluster)


# Re-draw heatmap with clusters annotation
annotation_row_df <- data.frame(
  Cluster = factor(clusters) # Convert cluster assignments to a factor
)

cluster_colours <- brewer.pal(k, "Set2")
names(cluster_colours) <- levels(annotation_row_df$Cluster)

results <- pheatmap(
  data,
  angle_col = 45,
  cluster_rows = TRUE,
  treeheight_row = 40,
  cluster_cols = FALSE,
  fontsize = 8,
  border_color = NA,
  color = colour_palette,
  breaks = as.numeric(custom_breaks),
  show_rownames = show_row_names,
  fontsize_row = 4,
  filename = snakemake@output[["pdf"]],
  width = snakemake@params[["width"]],
  height = snakemake@params[["height"]],
  main = paste0("Heatmap of dPSI\n(abs. threshold: ", threshold, ")"),
  annotation_row = annotation_row_df,
  annotation_colors = list(
    Cluster = cluster_colours
  )
)


# Add dPSI and z-scores
data$gene_id <- rownames(data)
df <- df %>%
  left_join(data, by = "gene_id") %>%
  # convert zero values back to NA
  mutate(across(starts_with("delta_PSI_mean_"), ~ ifelse(. == 0, NA, .))) %>%
  mutate(across(starts_with("z_score_corr_"), ~ ifelse(. == 0, NA, .)))

# Write output
write_csv(df, snakemake@output[["csv"]])
