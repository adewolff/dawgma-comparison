library("ensembldb")
library("readr")
library("dplyr")


# Ensembledb loading of gene annotation
# points to gtf file
gtffile <- "quantification/gene_annotation/Arabidopsis_thaliana.TAIR10.39.gtf.gz"
# Generate SQLite database file
DB <- ensDbFromGtf(gtf = gtffile)
# Load DB file 
EDB <- EnsDb(DB)
# Convert DB file to data frame containing transcript info
tx2gene <- data.frame(transcripts(EDB, return.type = "DataFrame"))
tx2gene <- dplyr::select(tx2gene, tx_name, gene_id)