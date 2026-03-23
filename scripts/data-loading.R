
# Load the data
# Counts has each sample as a column and Each gene expression as a row
# row.names  = 1 because column 1 has the identifiers for each gene
counts <- read.csv("data/GSE96870_counts_cerebellum.csv", row.names = 1)
View(head(counts, 2))

# Load the Sample annotations
# This is the metadata for each sample
# Each sample is column in counts, so we will call this coldata
# row.names  = 1 because column 1 has the unique identifiers for each sample

coldata <- read.csv("data/GSE96870_coldata_cerebellum.csv",
                    row.names = 1)
View(head(coldata, 2))

# Load the Gene Annotations
# We use row.name = 5 because in rowranges.tsv, the 5th column is the gene name
rowranges <- read.delim("data/GSE96870_rowranges.tsv",
                        sep = "\t",
                        colClasses = c(ENTREZID = "character"),
                        header = TRUE,
                        quote = "",
                        row.name = 5)
View(head(rowranges, 2))

# Combine all these tables into one object called SummarisedExperiment
# count object will be saved in assays slot
# coldata object will be stored in (sample metadata)
# rowranges object describing the genes will be stored in the rowRanges slot
# (features metadata)


# We must match the column names of count with rownames of coldata
all.equal(colnames(counts), rownames(coldata))

# We must match row names of counts with rows names of the rowRanges 
# gene
all.equal(rownames(counts), rownames(rowranges))

# Create our summarised experiment
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(rowranges, "GRanges"),
  colData = coldata
)

# Access the counts
head(assay(se, "counts"))
head(colData(se))
head(rowData(se))

# Make better sample IDs that show sex, time and mouse ID
se$Label <- paste(se$sex, se$time, se$mouse, sep = "_")
se$Label

# set the column names
colnames(se) <- se$Label


# Ordering the samples based on sex and time
se$Group <- paste(se$sex,
                  se$time,
                  sep='_')

# The way we want to order the objects goes like this:
#"Female_Day0"
#"Male_Day0"

#"Female_Day4"
#"Male_Day4"

#"Female_Day8"
#"Male_Day8"

se$Group <- factor(se$Group, levels = c(
  "Female_Day0", "Male_Day0",
  "Female_Day4", "Male_Day4",
  "Female_Day8", "Male_Day8"
))

se <- se[, order(se$Group)]
colData(se)
se$Label <- factor(se$Label, levels = se$Label)

saveRDS(se, "data/GSE96870_se.rds")
rm(se) # remove the object!
se <- readRDS("data/GSE96870_se.rds")

### End of Data Loading



