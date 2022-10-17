# biorunR

An R package for running bioinformatic analysis. Includes wrapper functions for performing differential expression (DESeq2) and geneset enrichment (GSEA, ssGSEA, topGO), and functions for plotting data.

## Dependencies
 - DESeq2 (& ashr)
 - dplyr
 - fgsea
 - ggplot2
 - GSVA
 - magrittr
 - msigdbr
 - org.Hs.eg.db
 - RColorBrewer
 - topGO

## Installation
```R
if(!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("nbay13/biorunR")
```
## Usage
...