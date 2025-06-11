# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print(Sys.time())

library(tidyverse)
library(cowplot)
library(doParallel)
library(foreach)

# Load Snakemake variables
bin.number <- snakemake@params[["bin_number"]]
rank.files <- unique(snakemake@input[["ranked"]])
prop.files <- unique(snakemake@input[["proportions"]])
outdir <- snakemake@params[["outdir"]]
threads <- snakemake@threads

# Set up parallel backend
print(paste0("Setting up parallel backend with ", threads, " threads"))
cl <- makeCluster(threads)
registerDoParallel(cl)

# Create output directory
print(paste0("Creating output directory: ", outdir))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load all proportion data
data.list <- lapply(prop.files, read_csv) %>%
  lapply(function(x) {
    select(x, 1:(4 + 2 * bin.number), twin_peaks, delta_PSI_mean, delta_PSI_SD)
  })

# Combine orf_id and gene into gene.id and remove orf_id/gene columns
create.geneid <- function(x) {
  x <- x %>%
    unite("gene.id", c("orf_id", "gene"), sep = "_", remove = TRUE)
  return(x)
}
data.list <- lapply(data.list, create.geneid)

# Load all rank data
rank.list <- lapply(rank.files, read_csv) %>%
  lapply(function(x) {
    select(
      x,
      c("orf_id", "gene"),
      starts_with("stabilised_in_"),
      starts_with("destabilised_in_")
    )
  })

# Combine orf_id and gene into gene.id and remove orf_id/gene columns
rank.list <- lapply(rank.list, create.geneid)

# Join data
data.list <- mapply(
  function(x, y) {
    left_join(x, y, by = "gene.id") %>%
      filter(complete.cases(.))
  },
  data.list,
  rank.list,
  SIMPLIFY = FALSE
)

# Remove lines where twin_peaks are TRUE
data.list <- lapply(data.list, function(x) {
  dpeak.values <- unique(x[["twin_peaks"]])
  if (length(dpeak.values) == 1) {
    print("Barcodes with twin peaks detection was skipped...")
    dpeaks.shapes.removed <- FALSE
  } else {
    print("Barcodes with twin peaks detection was performed...")
    dpeaks.shapes.removed <- TRUE
    x <- x %>%
      filter(twin_peaks == FALSE)
  }
  return(x)
})

calculate_means <- function(df) {
  comparison <- unique(df$Comparison)
  ref_prefix <- str_split(comparison, "_vs_")[[1]][2]
  test_prefix <- str_split(comparison, "_vs_")[[1]][1]

  # Identify ref and test columns based on prefix
  ref_cols <- grep(paste0("^", ref_prefix), names(df), value = TRUE)
  test_cols <- grep(paste0("^", test_prefix), names(df), value = TRUE)

  # Combine ref and test columns for easier processing
  all_cols <- c(ref_cols, test_cols)

  # Group by gene.id and calculate the mean for all relevant columns
  result_df <- df %>%
    group_by(gene.id) %>%
    summarise(across(
      all_of(all_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}"
    )) %>%
    ungroup()

  return(result_df)
}

# Calculate means for each data frame
mean.list <- lapply(data.list, calculate_means)

# Check if there is overlap in the ref samples
col.names <- lapply(mean.list, names)
ref.overlap <- intersect(col.names[[1]], col.names[[2]])[
  2:length(intersect(col.names[[1]], col.names[[2]]))
]

# Remove overlapping ref samples in second data frame, if any
if (length(ref.overlap) > 0) {
  mean.list[[2]] <- mean.list[[2]] %>%
    select(-all_of(ref.overlap))
}

# Join data
# only keep orfs were both conditions have data
df <- inner_join(mean.list[[1]], mean.list[[2]], by = "gene.id") %>%
  pivot_longer(
    cols = -gene.id,
    names_to = "sample",
    values_to = "proportion"
  ) %>%
  separate(sample, into = c("sample", "bin"), sep = "\\_") %>%
  group_by(bin, sample) %>%
  ungroup()

# Overlapping orf_ids
ids <- unique(df$gene.id)

# Plot overlapping orf_ids
colours <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#FFFF33",
  "#A65628",
  "#F781BF",
  "#999999"
)
foreach(
  id = ids,
  .packages = c("tidyverse", "cowplot", "RColorBrewer")
) %dopar%
  {
    tmp <- df[df$gene.id == id, ]

    # Define colours for lines
    length <- length(unique(tmp$sample))
    colours.plot <- colours[1:length]

    # Create plot
    p <- ggplot(
      tmp,
      aes(x = bin, y = proportion, group = sample, colour = sample)
    ) +
      geom_point(size = 4) +
      geom_line(linewidth = 1) +
      labs(title = id, y = "Mean proportion of reads", x = "Bin") +
      theme_cowplot(18) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_manual(values = colours.plot)

    ggsave(file.path(outdir, paste0(id, ".pdf")), p, width = 10, height = 6)
  }

# Stop the parallel backend
stopCluster(cl)

print("Plotting complete")
print(Sys.time())
