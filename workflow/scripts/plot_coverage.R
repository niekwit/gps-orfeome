# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

fasta <- snakemake@input[["fasta"]]
bin.number <- snakemake@params[["bin_number"]]

# Load count table
if (bin.number == 1) {
  data <- read.delim(snakemake@input[["tsv"]]) %>%
    dplyr::select(-all_of(c("sgRNA", "gene")))
} else {
  data <- read.delim(snakemake@input[[1]]) %>%
    dplyr::select(-all_of(c("barcode_id", "orf_id", "gene")))
}

# Create df to store coverage
# Remove prepended X from samples names (only happens when they start with a number)
df <- as.data.frame(matrix(data = NA, ncol = 2, nrow = ncol(data))) %>%
  rename(sample = 1, coverage = 2) %>%
  mutate(sample = str_remove(names(data), "^X"))

# Get number of sgRNAs from fasta file (library size)
sgrnas <- length(readLines(fasta)) / 2

# Calculate sequence coverage for each sample
for (n in seq(names(data))) {
  tmp <- data %>%
    dplyr::select(n)

  sum.reads <- sum(tmp)
  coverage <- sum.reads / sgrnas

  df[n, "coverage"] <- coverage
}

# Plot sequence coverage and save
p <- ggplot(df, aes(x = sample, y = coverage)) +
  geom_bar(stat = "identity", fill = "#419179", colour = "black") +
  theme_cowplot(16) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ylab("Fold sequence coverage") +
  xlab(NULL)

ggsave(snakemake@output[[1]], p)


# Close log file
sink(log, type = "output")
sink(log, type = "message")
