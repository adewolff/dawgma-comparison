library("ensembldb")
library("tximport")
library("readr")
library("DESeq2")
library("vsn")
library("dplyr")
library("ggplot2")
library("calibrate")
library('plotly')
library("IHW")

# Set Working Directory to source file location.
# Only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

# Convert tximport data to DESeqDataSet
dds <- DESeqDataSetFromTximport(txi_tx,
  colData = samples,
  design = ~ condition
)

# Filter counts with low reads
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10, ]
nrow(dds)

# Run DESeq analysis
dds <- DESeq(dds, fitType = "parametric")
res <- results(dds)

# Shrink LFC estimates
resLFC <- lfcShrink(dds, coef = paste0(resultsNames(dds)[2]), type = "apeglm")

# Order results by smallest p-val
res_ordered <- res[order(res$pvalue), ]

# Perform independent hypothesis weighting
resIHW <- results(dds, filterFun = ihw)

# Put resLFC in form of normal dataframe for ggplot2
resnorm <- data.frame(gene_id = rownames(resLFC),
                      base_mean = resLFC$baseMean,
                      log2_fold_change = resLFC$log2FoldChange,
                      lfc_se = resLFC$lfcSE,
                      pvalue = resLFC$pvalue,
                      padj = resLFC$padj)


# plotting ----------------------------------------------------------------

# Create MA plot of shrunken results
plotMA(resLFC, ylim = c(-12, 12))
top_gene <- rownames(res)[which.max(resnorm$log2_fold_change)]
with(res[top_gene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, top_gene, pos = 2, col = "dodgerblue")
})

# Create volcano plot
ggplot(data = resnorm, aes(x = log2_fold_change, y = -log2(pvalue))) +
  geom_point(mapping = aes(color = base_mean)) +
  scale_color_gradient(limits = c(10, 1088)) +
  xlim(-6, 6) +
  xlab("log2FoldChange") +
  ylab("-log2(p-value)")

# Create Plotly interactive map
plot_ly(
  resnorm, x = ~log2_fold_change, y = ~-log2(pvalue),
  # Hover text:
  text = ~paste("gene: ", gene_id),
  color = ~base_mean
) %>% 
  layout(xaxis = list(range = c(-6,6)))

# reorder res object to show top log2 fold changes
resnorm <- resnorm[order(resnorm$log2_fold_change), ]
View(tail(resnorm, 10))


# To find Gene, search transcript ID in est.fa.gz,
# then NCBI blast search for gene.
