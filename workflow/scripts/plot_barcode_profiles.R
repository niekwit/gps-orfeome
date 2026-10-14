# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print(Sys.time())

library(tidyverse)
library(reshape2)
library(doParallel)
library(foreach)
source(file.path(snakemake@scriptdir, "theme_gpsw.R"))
theme_path <- file.path(snakemake@scriptdir, "theme_gpsw.R")


# Load Snakemake variables
bin_number <- snakemake@params[["bin_number"]]
rank_file <- snakemake@input[["ranked"]]
prop_file <- snakemake@input[["proportions"]]
comparison <- snakemake@wildcards[["comparison"]]
test_sample <- str_split(comparison, "_vs_")[[1]][1]
ref_sample <- str_split(comparison, "_vs_")[[1]][2]
outdir <- snakemake@params[["outdir"]]
threads <- snakemake@threads

# Set up parallel backend
print(paste0("Setting up parallel backend with ", threads, " threads"))
cl <- makeCluster(threads)
registerDoParallel(cl)

# Create output directory
print(paste0("Creating output directory: ", outdir))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load proportion data
data <- read_csv(prop_file) %>%
  select(
    barcode_id,
    orf_id,
    gene,
    starts_with(test_sample),
    starts_with(ref_sample),
    twin_peaks,
    delta_PSI_mean,
    delta_PSI_SD
  )

# Get columns that contain hit information
data_ranked <- read_csv(rank_file) %>%
  select(orf_id, stabilised, destabilised, good_barcodes, agreeing_barcodes)
data <- data %>%
  left_join(data_ranked, by = "orf_id") %>%
  # Only require the plotted bin proportions/twin_peaks/hit
  # columns to be complete. delta_PSI_SD is legitimately NA for genes
  # with a single good barcode (sd() of one value), so it must not
  # be part of this completeness check, or such genes would be dropped
  # entirely before they ever reach the plotting loop below.
  filter(complete.cases(select(
    .,
    -delta_PSI_mean,
    -delta_PSI_SD
  )))

# Ranks are NA for non-hits, so they are joined after the complete.cases filter
data_ranks <- read_csv(rank_file) %>%
  select(orf_id, stabilised_rank, destabilised_rank)
data <- data %>%
  left_join(data_ranks, by = "orf_id")

# Check if barcodes with twin peaks have been removed
# (assumes that at least one barcode with twin
# peaks in the entire data set exists)
dpeak_values <- unique(data[["twin_peaks"]])
if (length(dpeak_values) == 1) {
  print("Barcode analysis with twin peaks detection was skipped...")
  dpeaks_shapes_removed <- FALSE
} else {
  print("Barcodes with twin peaks detection was performed...")
  dpeaks_shapes_removed <- TRUE
}

# Plot each gene and its proportion of reads in bins
columns <- c("stabilised", "destabilised")
for (column in columns) {
  # Create sub directory for each column
  dir <- file.path(outdir, column)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  # Subset genes/ORF IDs
  tmp <- data %>%
    filter(get(column) == TRUE)

  # Create ID (gene + ORF) for plots
  tmp <- tmp %>%
    mutate(
      gene.id = paste0(gene, "_", orf_id),
      rank = .data[[paste0(column, "_rank")]]
    )

  foreach(
    id = tmp$gene.id,
    .packages = c("tidyverse", "reshape2", "cowplot", "scales")
  ) %dopar%
    {
      source(theme_path)
      df <- tmp[tmp$gene.id == id, ]
      hit_rank <- unique(df$rank)
      dpsi_mean <- round(unique(df$delta_PSI_mean), 2)
      dpsi_sd <- round(unique(df$delta_PSI_SD), 2)
      n_good <- unique(df$good_barcodes)
      n_agree <- unique(df$agreeing_barcodes)
      df <- df %>%
        dplyr::select(
          gene.id,
          starts_with(ref_sample),
          starts_with(test_sample),
          twin_peaks
        ) %>%
        reshape2::melt() %>%
        separate(variable, into = c("barcode", "bin"), sep = "\\_") %>%
        group_by(bin, barcode) %>%
        # add replicate number to condition, e.g. Un_1
        mutate(barcode = paste0(barcode, "_", row_number()))

      # Define colour gradient for lines. The colour scale below is a named
      # vector keyed by the actual test/ref barcode labels (not just
      # concatenated in order), so test is always red and control is
      # always black/grey
      barcode_count <- length(unique(df$barcode)) / 2
      grey_colours <- scales::seq_gradient_pal("black", "grey", "Lab")(seq(
        0,
        1,
        length.out = barcode_count
      ))
      red_colours <- scales::seq_gradient_pal("red4", "red", "Lab")(seq(
        0,
        1,
        length.out = barcode_count
      ))
      colour_values <- setNames(
        c(red_colours, grey_colours),
        c(
          paste0(test_sample, "_", seq_len(barcode_count)),
          paste0(ref_sample, "_", seq_len(barcode_count))
        )
      )

      # Create plot
      if (dpeaks_shapes_removed == FALSE) {
        p <- ggplot(
          df,
          aes(x = bin, y = value, group = barcode, colour = barcode)
        ) +
          geom_point(size = 4) +
          geom_line(linewidth = 1) +
          labs(title = id, y = "Proportion of reads", x = "Bin") +
          theme_gpsw() +
          scale_colour_manual(values = colour_values) +
          theme(plot.title = element_text(hjust = 0.5))
      } else {
        p <- ggplot(
          df,
          aes(
            x = bin,
            y = value,
            group = barcode,
            colour = barcode,
            shape = twin_peaks,
            alpha = twin_peaks
          )
        ) +
          geom_point(size = 6) +
          geom_line(linewidth = 1) +
          labs(title = id, y = "Proportion of reads", x = "Bin") +
          theme_gpsw() +
          scale_colour_manual(values = colour_values) +
          scale_alpha_manual(values = c(1, 0.1)) +
          scale_shape_manual(values = c(16, 15)) +
          theme(plot.title = element_text(hjust = 0.5))
      }

      # Add label for delta PSI +- SD and number of good barcodes that agree
      # in direction with the mean delta PSI inside top left plot
      p <- p +
        geom_text(
          aes(
            x = 1.5,
            y = max(df$value, na.rm = TRUE) * 0.95,
            label = paste0(
              "dPSI: ",
              dpsi_mean,
              " ± ",
              dpsi_sd,
              "\n",
              n_agree,
              "/",
              n_good,
              " barcodes agree"
            )
          ),
          colour = "black",
          size = 6,
          hjust = 0.5,
          vjust = 0.5
        )

      # Save plot
      ggsave(
        file.path(dir, paste0(hit_rank, "_", id, ".pdf")),
        p,
        width = 10,
        height = 6
      )
    }
}

# Stop the parallel backend
stopCluster(cl)

print("Plotting complete")
print(Sys.time())
