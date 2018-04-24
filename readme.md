# DAWGMA PRJDB2508 analysis

RNA-SEQ analysis of [European Nucleotide Archive's PRJDB2508 study](https://www.ebi.ac.uk/ena/data/view/PRJDB2508).

##How to use R code for analyzing rna seq data:
1. Put quantification files in quantification/quants
2. Put gene annotation gtf file in quantification/gene_annotation and/or ensure
sqlite database file is present in root.
3. Put samples.txt file in quantification directory.
4. Open log2foldchange_comparison.R
    - On Line 28 modify file.path to correspond to naming convention of quant
  files.
    - Ensure tx2gene vector in log2log2foldchange_comparison.R(line 18) and
  gene_annotation.R(line 17) both select proper column names for tx and gene id
    - Ensure line 27-29 are referencing correct column of samples.txt
  (Name by default)


##quantification directory
This directory contains all files related to quantification.

- index.sh indexes the *Arabidopsis thaliana* reference transcript
- sequence.sh sequences the reads obtained from the study website.
- Arabidopsis_thaliana.....fa.gz is the reference transcript used to create the Salmon index
- quants directory contains results from Salmon quantification of transcripts obtained from the study website
