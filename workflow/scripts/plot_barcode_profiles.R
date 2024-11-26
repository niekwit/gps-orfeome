# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(reshape2)
library(cowplot)


data.file <- snakemake@input[["csv"]]
comparison <- snakemake@wildcards[["comparison"]]
test.sample <- str_split(comparison, "_vs_")[[1]][1]
ref.sample <- str_split(comparison, "_vs_")[[1]][2]
outdir <- snakemake@output[["outdir"]]

# Create output directory
dir.create(outdir, showWarnings = FALSE)

# Load data and calculate proportion of reads in each bin
data <- read_csv(data.file, show_col_types = FALSE) %>%
  group_by(barcode) %>%
  rowwise() %>%
  mutate(across(starts_with(paste0(ref.sample, "_")), ~ . / sum(c_across(starts_with(paste0(ref.sample, "_")))))) %>%
  mutate(across(starts_with(paste0(test.sample, "_")), ~ . / sum(c_across(starts_with(paste0(test.sample, "_")))))) 


# Write genes to text file
#genes <- tmp %>%
#  select(gene.id) %>%
#  distinct() %>%
#  pull() %>%
#  sort()
#write.table(genes, file.path(dir, "genes.txt"), 
#            row.names = FALSE, 
#            col.names = FALSE,
#            quote = FALSE)

# Get columns that contain hit information
info.columns <- data %>%
  as.data.frame() %>%
  select(c("orf", "gene"),
         starts_with("stabilised_in_"), 
         starts_with("destabilised_in_")) %>%
  colnames()

# Plot each gene and its proportion of reads in bins
for (column in info.columns[3:length(info.columns)]) {
  # Create sub directory for each column
  dir <- file.path(outdir, column)
  dir.create(dir, showWarnings = FALSE)
  
  # Subset genes/ORF IDs
  df <- data %>%
    filter(get(column) == TRUE)
  
  # Get gene IDs
  genes <- df %>%
    select(gene) %>%
    distinct() %>%
    pull()
  
  
  print(paste0("Plotting ORF ID: ", gene))
  # Get gene data
  df <- tmp[tmp$gene.id == gene,] %>%
    select(gene.id, starts_with("Un."), starts_with("Hyp.")) %>%
    melt() %>%
    separate(variable, into = c("condition", "bin"), sep = "\\.") %>%
    group_by(bin, condition) %>%
    # add replicate number to condition, e.g. Un_1
    mutate(condition = paste0(condition, "_", row_number()))
  
  # Define colours
  samples <- length(unique(df$condition))
  grey.colours <- scales::seq_gradient_pal("black", "grey", "Lab")(seq(0,1,length.out=samples / 2))
  red.colours <- scales::seq_gradient_pal("red4", "red", "Lab")(seq(0,1,length.out=samples / 2))
  
  # Create plot
  p <- ggplot(df, aes(x = bin, 
                      y = value, 
                      group = condition,
                      colour = condition)) +
    geom_point(size = 4) +
    geom_line(linewidth = 1) +
    labs(title = gene,
         y = "Proportion of reads",
         x = "Bin") +
    theme_cowplot(18) + 
    scale_colour_manual(values = c(red.colours, grey.colours)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  dir.create(file.path(dir, "plots"), showWarnings = FALSE)
  ggsave(file.path(dir, "plots", paste0(gene, ".png")), 
         p, 
         width = 10, 
         height = 6,
         bg = "white")
}






# Close log file
sink(log, type = "output")
sink(log, type = "message")
