library(SummarizedExperiment)
library(DESeq2)
library(vsn)
library(ggplot2)

library(RColorBrewer)
library(hexbin)

# library(iSEE)
library(ComplexHeatmap)
library(cowplot)
library(ComplexHeatmap)
library(apeglm)

se <- readRDS("data/GSE96870_se.rds")
se <- se[rowSums(assay(se, "counts")) > 5, ]

# Create DESeq Data set
dds <- DESeq2::DESeqDataSet(se,
                            design = ~ sex + time)

# THE FIX: Assign the result back to dds
dds <- DESeq(dds)

# Now this will work perfectly!
resTime <- results(dds, contrast = c("time", "Day8", "Day0"))
summary(resTime)
head(resTime[order(resTime$pvalue), ])



vsd <- vst(dds, blind = TRUE)

# Get top DE genes
genes <- resTime[order(resTime$pvalue), ] |>
  head(10) |>
  rownames()
heatmapData <- assay(vsd)[genes, ]

# Scale counts for visualization
heatmapData <- t(scale(t(heatmapData)))

# Add annotation
heatmapColAnnot <- data.frame(colData(vsd)[, c("time", "sex")])
heatmapColAnnot <- HeatmapAnnotation(df = heatmapColAnnot)

# Plot as heatmap
geneExpressionHeatmap <- ComplexHeatmap::Heatmap(heatmapData,
                        top_annotation = heatmapColAnnot,
                        cluster_rows = TRUE, cluster_columns = FALSE)

# Save the heatmap
# 1. Open the PDF file (adjust width and height in inches as needed)
pdf("gene_expression_heatmap.pdf", width = 8, height = 6)

# 2. Explicitly draw your heatmap object
ComplexHeatmap::draw(geneExpressionHeatmap)

# 3. Close the graphic device to finish saving the file
dev.off()


