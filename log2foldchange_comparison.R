library("tximport")
library("readr")
library("DESeq2")
library("vsn")
library("dplyr")
library("ggplot2")
library("calibrate")

# tximport pipeline -------------------------------------------------------

# get tx2gene file using ensembledb
source("lib/gene_annotation.R")

# point tximport to location of quant files
dir <- "quantification"
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
rownames(samples) <- samples$Name
files <- file.path(dir, "quants", samples$Name, "quant.sf")
names(files) <- paste0(samples$Name)
# all(file.exists(files))

# Run tximport function on quant files
txi_tx <- tximport(files, type = "salmon", tx2gene = tx2gene,
                   txIdCol = "tx_id", geneIdCol = "gene_id")
# names(txi_tx)
# head(txi_tx$counts)


# Diff expr ---------------------------------------------------------------

dds <- DESeqDataSetFromTximport(txi_tx,
  colData = samples,
  design = ~ condition
)
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)


dds <- DESeq(dds, fitType = "mean")
res <- results(dds)


# plotting ----------------------------------------------------------------


plotMA(res, ylim = c(-5, 5))
top_gene <- rownames(res)[which.min(res$padj)]
with(res[top_gene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, top_gene, pos = 2, col = "dodgerblue")
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
