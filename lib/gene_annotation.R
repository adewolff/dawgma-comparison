# gtf file location: 
# ftp://ftp.ensemblgenomes.org/pub/plants/release-39/gtf/arabidopsis_thaliana/

# Ensembledb loading of gene annotation
# points to gtf file
if(length(Sys.glob("gene_annotation/*.gtf.gz")) > 0){
  gtffile <- Sys.glob("gene_annotation/*.gtf.gz")
} else if(length(Sys.glob("gene_annotation/*.gtf")) > 0){
  gtffile <- Sys.glob("gene_annotation/*.gtf")
}
# Generate SQLite database file
DB <- ensDbFromGtf(gtf = gtffile, path = "gene_annotation")
# Load DB file 
EDB <- EnsDb(DB)
# Convert DB file to data frame containing transcript info
tx2gene <- data.frame(transcripts(EDB, return.type = "DataFrame"))
tx2gene <- dplyr::select(tx2gene, tx_name, gene_id)