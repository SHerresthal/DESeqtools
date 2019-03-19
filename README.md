#  DESeqtools

A package for analysis of RNASeq-data based on the package DESeq2. (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) 

## Package installation

```r
install.packages("devtools")
library(devtools)
devtools::install_github("SHerresthal/RNA-DESeq2")
```

## Example Analysis

This tutorial contains an example analysis which can be used as a template for standard bulk RNA-seq analysis using a dataset of human MDSC-like cells (published here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92852). 
See https://sherresthal.github.io for the example workflow. It contains: 

1.  Installation of the DESeqtools-package
2.  Project Information
3.  Obligatory Data Structure
4. Data Import
5. Build DESeq Data Set
6. Exploratory Data Analysis
7. Batch Effect Removal (treat with caution!)
8. Differential expression analysis
9. Export of the results
10. Save image and session info
