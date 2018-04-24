library("tximport")
library("readr")
library("ensembldb")
library("DESeq2")
library("vsn")
library("dplyr")
library("ggplot2")
library("calibrate")


# Ensembledb loading of gene annotation
# points to gtf file
gtffile <- "quantification/gene_annotation/Arabidopsis_thaliana.TAIR10.39.gtf.gz"
# Generate SQLite database file
DB <- ensDbFromGtf(gtf = gtffile)
# Load DB file 
EDB <- EnsDb(DB)
# Convert DB file to data frame containing transcript info
tx2gene <- transcripts(EDB, return.type = "DataFrame")

dir <- "quantification"
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
rownames(samples) <- samples$Name
files <- file.path(dir, "quants", samples$Name, "quant.sf")
names(files) <- paste0(samples$Name)
all(file.exists(files))

txi.tx <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi.tx)
head(txi.tx$counts)

dds <- DESeqDataSetFromTximport(txi.tx,
  colData = samples,
  design = ~ condition
)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)


dds <- DESeq(dds, fitType = "mean")
res <- results(dds)

plotMA(res, ylim = c(-5, 5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topGene, pos = 2, col = "dodgerblue")
})

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue),
  pch = 20, main = "Volcano plot",
  xlim = c(-10, 10), ylim = c(0, 3)
))

# Colors
with(subset(res, padj < .05), points(log2FoldChange, -log10(pvalue),
  pch = 20,
  col = "red"
))
with(
  subset(res, abs(log2FoldChange) > 1),
  points(log2FoldChange, -log10(pvalue), pch = 20, col = "orange")
)
with(
  subset(res, padj < .05 & abs(log2FoldChange) > 1),
  points(log2FoldChange, -log10(pvalue), pch = 20, col = "green")
)

# Label points with the textxy function from the calibrate plot
# with(subset(res, padj<.05 & abs(log2FoldChange)>1),
#     textxy(log2FoldChange, -log10(pvalue), cex=.8))

# reorder res object to show top log2 fold changes
res <- res[order(res$log2FoldChange), ]
tail(res, 10)

# To find Gene, search transcript ID in est.fa.gz,
# then NCBI blast search for gene.
