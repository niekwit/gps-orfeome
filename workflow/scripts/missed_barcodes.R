# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(reshape2)
library(cowplot)

# Load count table
counts <- read.delim(snakemake@input[[1]]) %>%
  dplyr::select(-c("sgRNA", "gene"))

# Get number of sgRNAs with zero counts
df <- colSums(counts == 0) %>%
  melt() %>%
  rownames_to_column(var = "sample") 

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




