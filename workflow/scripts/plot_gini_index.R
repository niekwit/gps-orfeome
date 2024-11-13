# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(ineq)
library(cowplot)

# Load count table
counts <- read.delim(snakemake@input[[1]])

# Create df to store Gini indices
df <- as.data.frame(matrix(data = NA,
                           ncol = 2,
                           nrow = ncol(counts) - 2))
names(df) <- c("sample","Gini_index")
df$sample <- names(counts)[3:ncol(counts)]

# Calculate Gini index for each sample
for (i in seq_len(nrow(df))){
  df[i,"Gini_index"] <- Gini(as.vector(counts[i + 2])[[1]])
}

# Plot Gini index bar graph and save
p <- ggplot(data = df, aes(x = sample, 
                           y = Gini_index)) + 
  geom_bar(stat = "identity",
           fill = "#419179",
           colour = "black") +
  theme_cowplot(16) +
  ylab("Gini index") +
  xlab(NULL) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggtitle("Evenness of barcodes") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave(snakemake@output[[1]], p)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

