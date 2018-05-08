# THIS BRANCH UTILIZES PARALLELIZATION WHEN AVAILABLE. USABILITY ON ALL SYSTEMS NOT GUARANTEED. IF YOU DON'T KNOW WHAT THIS MEANS PLEASE USE MASTER BRANCH. 

# DAWGMA RNA-seq Differential expression analysis
This program allows DAWGMA members to run a Differential expression analysis of
their rna seq results.

## How to use R code for analyzing rna seq data:
1. Put quantification files in quants directory
2. Put gene annotation .gtf file or sqlite file in gene_annotation directory
3. Put samples.txt file in root directory.
4. Open log2foldchange_comparison.R
    - Set working directory to source file
    (Session>Set working directory>To Source File Location)
    - On Line 27 modify file.path to correspond to naming convention of quant
  files and their directories (or uphold xx_quant/quant.sf structure).
    - Ensure tx2gene vector in log2log2foldchange_comparison.R(line 18) and
  gene_annotation.R(line 17) both select proper column names for tx and gene id
  (found using `names(tx2gene)`)
    - Ensure line 27-29 are referencing correct column of samples.txt
  (Name by default)
5. Source log2foldchange_comparison.R


## FAQ
#### Where do I obtain a "gene annotation .gf file?"
```
Go here: ftp://ftp.ensemblgenomes.org/pub/, and select a release. Which
release you pick doesn't matter. However, keep note of which release you use, as
this is important for experiment reproducibility. If you're unsure of which
release to pick, pick the latest. From here, navigate to the kingdom in which
your species is. For example, Saccharomyces cerevisiae would be under fungi.
Next, go to the "gtf" directory, and from here select your species. Finally,
download the gtf.gz file. DO NOT download the abinitio.gtf.gz file.   
```

#### What should my samples.txt look like?
```
Check out the samples.txt included in this repository.
(reference edition coming soonTM)
```
