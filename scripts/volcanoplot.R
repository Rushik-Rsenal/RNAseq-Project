se <- readRDS("data/GSE96870_se.rds")
se <- se[rowSums(assay(se, "counts")) > 5, ]

# Create DESeq Data set
dds <- DESeq2::DESeqDataSet(se,
                            design = ~ sex + time)

# THE FIX: Assign the result back to dds
dds <- DESeq(dds)

# Now this will work perfectly!
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))

# Volcano Plot
#1. Convert to a clean data frame (removing NAs)
res_df <- as.data.frame(resTime)
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# 2. Add significance coloring thresholds (FDR < 0.05 and Log2FC > 1)
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Upregulated"
res_df$Significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated"

# 3. THE CRUCIAL STEP: Create a specific label column using your 'genes' variable
# We set everything to NA first, then only fill in the names that match your top 10
res_df$Label <- NA
res_df$Label[rownames(res_df) %in% genes] <- rownames(res_df)[rownames(res_df) %in% genes]

# 4. Load the plotting libraries
library(ggplot2)
library(ggrepel)

# 5. Draw the Volcano Plot with labels!
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Not Significant" = "grey80", 
                                "Upregulated" = "red")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  
  # The magic function that adds non-overlapping text labels
  geom_text_repel(aes(label = Label), 
                  color = "black",       # Forces text to be black so it's readable
                  box.padding = 0.8,     # Adds space around the text
                  point.padding = 0.3,   # Adds space between text and the point
                  max.overlaps = Inf,    # Forces it to draw all 10 labels no matter what
                  show.legend = FALSE) + # Stops a weird 'a' from appearing in your legend
  
  labs(title = "Volcano Plot: Day 8 vs Day 0",
       subtitle = "Top 10 most significant genes labeled",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Gene Status")

ggsave("volcano_plot_Day8_vs_Day0.png", 
       width = 8,        # Width in inches
       height = 6,       # Height in inches
       dpi = 300)        # dpi = 300 ensures it is high-resolution/publication quality