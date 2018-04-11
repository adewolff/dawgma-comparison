# DAWGMA PRJDB2508 analysis

RNA-SEQ analysis of [European Nucleotide Archive's PRJDB2508 study](https://www.ebi.ac.uk/ena/data/view/PRJDB2508).

##quantification directory
This directory contains all files related to quantification.

- index.sh indexes the *Arabidopsis thaliana* reference transcript
- sequence.sh sequences the reads obtained from the study website.
- Arabidopsis_thaliana.....fa.gz is the reference transcript used to create the Salmon index
- quants directory contains results from Salmon quantification of transcripts obtained from the study website
