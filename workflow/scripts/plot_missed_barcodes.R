# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(reshape2)
library(cowplot)

bin.number <- snakemake@params[["bin_number"]]

# Load count table
if (bin.number == 1) {
  counts <- read.delim(snakemake@input[[1]]) %>%
    dplyr::select(-all_of(c("sgRNA", "gene")))
  
  # Get number of sgRNAs with zero counts
  df <- colSums(counts == 0) %>%
    melt() %>%
    rownames_to_column(var = "sample")
} else {
  counts <- read.delim(snakemake@input[[1]]) %>%
    dplyr::select(-all_of(c("barcode_id", "orf_id", "gene")))
  
  # Get unique conditions (no bin number)
  conditions <- unique(gsub("_[0-9]$", "", names(counts)))
  
  # Sum of all the columns that contain a condition
  for (c in conditions) {
    counts[[c]] <- rowSums(counts %>% select(contains(c)))
  }
  # Only keep the summed columns
  counts <- counts %>% select(all_of(conditions))
  
  # Get number of sgRNAs with zero counts for each condition
  zeros.count <- colSums(counts ==0) 
    
  # Count rows where all values are zero
  zeros.count["Missing in\nall conditions"] <- sum(rowSums(counts == 0) == ncol(counts))
  
  # Convert to data frame for plotting
  df <- data.frame(
    sample = names(zeros.count),
    value = as.numeric(zeros.count)
  )
  
  # Change order factor levels of sample column so that "Missing in all conditions" is last
  df$sample <- factor(df$sample, levels = c(conditions, "Missing in\nall conditions"))
}

# Create bar graph
p <- ggplot(data = df, aes(x = sample, y = value)) + 
  geom_bar(stat = "identity",
           fill = "#419179",
           colour = "black") +
  theme_cowplot(16) +
  xlab(NULL) +
  ylab("Missed barcodes") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(guide = guide_axis(angle = 45))

# Save plot
ggsave(snakemake@output[[1]], p)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")