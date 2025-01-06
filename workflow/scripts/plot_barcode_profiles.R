# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print(Sys.time())

library(tidyverse)
library(reshape2)
library(cowplot)
library(doParallel)
library(foreach)

# Load Snakemake variables
data.file <- snakemake@input[["csv"]]
comparison <- snakemake@wildcards[["comparison"]]
test.sample <- str_split(comparison, "_vs_")[[1]][1]
ref.sample <- str_split(comparison, "_vs_")[[1]][2]
outdir <- snakemake@params[["outdir"]]
threads <- snakemake@threads

# Set up parallel backend
print(paste0("Setting up parallel backend with ", threads, " threads"))
cl <- makeCluster(threads)
registerDoParallel(cl)

# Create output directory
print(paste0("Creating output directory: ", outdir))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load data and calculate proportion of reads in each bin
print("Calculating proportion of reads in each bin")
data <- read_csv(data.file, show_col_types = FALSE) %>%
  group_by(barcode_id) %>%
  mutate(
    # Pre-compute row sums for ref.sample and test.sample columns
    ref_sum = rowSums(across(starts_with(paste0(ref.sample, "_")), ~ ., .names = "ref_{col}"), na.rm = TRUE),
    test_sum = rowSums(across(starts_with(paste0(test.sample, "_")), ~ ., .names = "test_{col}"), na.rm = TRUE)
  ) %>%
  # Normalize ref.sample and test.sample columns
  mutate(
    across(starts_with(paste0(ref.sample, "_")), ~ . / ref_sum, .names = "{col}"),
    across(starts_with(paste0(test.sample, "_")), ~ . / test_sum, .names = "{col}")
  ) %>%
  # Drop intermediate columns if needed
  select(-ref_sum, -test_sum)

# Get columns that contain hit information
info.columns <- data %>%
  as.data.frame() %>%
  select(c("orf_id", "gene"),
         starts_with("stabilised_in_"), 
         starts_with("destabilised_in_")) %>%
  colnames()

# Plot each gene and its proportion of reads in bins
print("Plotting gene profiles")
for (column in info.columns[3:length(info.columns)]) {
  # Create sub directory for each column
  dir <- file.path(outdir, column)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Subset genes/ORF IDs
  tmp <- data %>%
    filter(get(column) == TRUE)

  # Get ORF IDs
  tmp <- tmp %>%
    mutate(gene.id = paste0(gene, "_", orf_id))

  foreach(id = tmp$gene.id, .packages = c("tidyverse", 
                                          "reshape2", 
                                          "cowplot", 
                                          "scales")) %dopar% {

    df <- tmp[tmp$gene.id == id, ] %>%
      select(gene.id, starts_with(ref.sample), starts_with(test.sample)) %>%
      select(-ends_with("distance")) %>%
      reshape2::melt() %>%
      separate(variable, into = c("barcode", "bin"), sep = "\\_") %>%
      group_by(bin, barcode) %>%
      mutate(barcode = paste0(barcode, "_", row_number()))
    
    # Define colour gradient for lines
    barcode.count <- length(unique(df$barcode)) / 2
    grey.colours <- scales::seq_gradient_pal("black", "grey", "Lab")(seq(0, 1, length.out = barcode.count))
    red.colours <- scales::seq_gradient_pal("red4", "red", "Lab")(seq(0, 1, length.out = barcode.count))
    
    # Create plot
    p <- ggplot(df, aes(x = bin,
                        y = value,
                        group = barcode,
                        colour = barcode)) +
      geom_point(size = 4) +
      geom_line(linewidth = 1) +
      labs(title = id,
           y = "Proportion of reads",
           x = "Bin") +
      theme_cowplot(18) + 
      scale_colour_manual(values = c(red.colours, grey.colours)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Save plot
    ggsave(file.path(dir, paste0(id, ".png")), 
            p, 
            width = 10, 
            height = 6,
            bg = "white")
  }
}

# Stop the parallel backend
stopCluster(cl)

# Create an empty file to signal that the script has completed
file.create(snakemake@output[["flag"]])

print("Plotting complete")
print(Sys.time())

# Close log file
sink(log, type = "output")
sink(log, type = "message")