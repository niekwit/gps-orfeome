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

# Check if barcodes with twin peaks have been removed
# (assumes that at least one barcode with twin peaks in the entire data set exists)
dpeak.values <- unique(data[["twin_peaks"]])
if (length(dpeak.values) == 1) {
  print("Barcodes with twin peaks detection was skipped...")
  dpeaks.shapes.removed <- FALSE
} else {
  print("Barcodes with twin peaks detection was performed...")
  dpeaks.shapes.removed <- TRUE
}

# Plot each gene and its proportion of reads in bins
for (column in info.columns[3:length(info.columns)]) {
  # Create sub directory for each column
  dir <- file.path(outdir, column)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Subset genes/ORF IDs
  tmp <- data %>%
    filter(get(column) == TRUE)

  # Create ID (gene + ORF) for plots
  tmp <- tmp %>%
    mutate(gene.id = paste0(gene, "_", orf_id))

  foreach(id = tmp$gene.id, .packages = c("tidyverse", 
                                          "reshape2", 
                                          "cowplot", 
                                          "scales")) %dopar% {
    print(paste0("Plotting ORF ID: ", id))
    df <- tmp[tmp$gene.id == id, ]
    deltaPSI.mean <- round(unique(df$delta_PSI_mean), 2)
    deltaPSI.sd <- round(unique(df$delta_PSI_SD), 2)
    df <- df %>%
      dplyr::select(gene.id, starts_with(ref.sample), starts_with(test.sample), twin_peaks) %>%
      reshape2::melt() %>%
      separate(variable, into = c("barcode", "bin"), sep = "\\_") %>%
      group_by(bin, barcode) %>%
      # add replicate number to condition, e.g. Un_1
      mutate(barcode = paste0(barcode, "_", row_number()))
    
    # Define colour gradient for lines
    barcode.count <- length(unique(df$barcode)) / 2
    grey.colours <- scales::seq_gradient_pal("black", "grey", "Lab")(seq(0, 1, length.out = barcode.count))
    red.colours <- scales::seq_gradient_pal("red4", "red", "Lab")(seq(0, 1, length.out = barcode.count))
    
    # Create plot
    if (dpeaks.shapes.removed == FALSE){
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
      
    } else {
      p <- ggplot(df, aes(x = bin,
                          y = value,
                          group = barcode,
                          colour = barcode,
                          shape = twin_peaks,
                          alpha = twin_peaks)) +
        geom_point(size = 6) +
        geom_line(linewidth = 1) +
        labs(title = id,
             y = "Proportion of reads",
             x = "Bin") +
        theme_cowplot(18) + 
        scale_colour_manual(values = c(red.colours, grey.colours)) +
        scale_alpha_manual(values = c(1, 0.1)) +
        scale_shape_manual(values = c(16, 15)) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    # Add label for delta PSI +- SD inside top left plot
    p <- p +
      geom_text(aes(x = 1.5,
                    y = max(df$value, na.rm = TRUE) * 0.95,
                    label = paste0("dPSI: ", deltaPSI.mean, " ± ", deltaPSI.sd)),
                colour = "black",
                size = 6,
                hjust = 0.5,
                vjust = 0.5)
    
    # Save plot
    ggsave(file.path(dir, paste0(id, ".pdf")),
           p,
           width = 10,
           height = 6)
  }
}

# Stop the parallel backend
stopCluster(cl)

print("Plotting complete")
print(Sys.time())

# Close log file
sink(log, type = "output")
sink(log, type = "message")