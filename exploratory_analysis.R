# This script relies on data obtained from "importer.R". Please run this script 
# before attempting to run this one.

library("pheatmap")
library("RColorBrewer")

# Removes rows that have no/nearly no information about amt of gene expression.
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Log transformation to stabilize variance
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Determine similarity between samples
sampleDists <- dist(t(assay(rld)))
sampleDists

# Create heatmap of sample similarity
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)