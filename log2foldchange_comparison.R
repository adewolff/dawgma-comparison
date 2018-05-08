library("ensembldb")
library("tximport")
library("readr")
library("DESeq2")
library("vsn")
library("dplyr")
library("ggplot2")
library("calibrate")

# tximport pipeline -------------------------------------------------------

# get tx2gene file using ensembledb
if (length(Sys.glob("gene_annotation/*.sqlite")) > 0){
  # load sqlite file
  EDB <- EnsDb(Sys.glob("gene_annotation/*.sqlite"))
  # Convert DB file to data frame containing transcript info
  tx2gene <- data.frame(transcripts(EDB, return.type = "DataFrame"))
  tx2gene <- dplyr::select(tx2gene, tx_id, gene_id, tx_biotype, tx_seq_start,
                           tx_seq_end, tx_cds_seq_start, tx_cds_seq_end,
                           tx_name)
} else{
  # if sqlite file not present in root, source("lib/gene_annotation.R") instead
  source("lib/gene_annotation.R")
}

# point tximport to location of quant files
samples <- read.table("samples.txt", header = TRUE)
rownames(samples) <- samples$Name
files <- file.path("quants", samples$Name, "quant.sf")
names(files) <- paste0(samples$Name)
all(file.exists(files))

# Run tximport function on quant files
txi_tx <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi_tx)
head(txi_tx$counts)


# Diff expr ---------------------------------------------------------------

dds <- DESeqDataSetFromTximport(txi_tx,
  colData = samples,
  design = ~ condition
)

# Filter counts with low reads
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
nrow(dds)

dds <- DESeq(dds, fitType = "mean")
res <- results(dds)

# shrink LFC estimates
resLFC <- lfcShrink(dds, coef = paste0(resultsNames(dds)[2]), type = "apeglm")

# order results by smallest p-val
res_ordered <-res[order(res$pvalue), ]

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
