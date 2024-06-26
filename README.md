# biorunR <img src="imgs/biorunr logo.png" align="right" width="150" height="150" />

An R package for running bioinformatic analyses. Includes wrapper functions for performing differential expression (DESeq2), geneset enrichment (GSEA, GSVA/ssGSEA, topGO), and data visualization (ggplot2).

## Dependencies
 - classInt
 - DESeq2 (& ashr)
 - dplyr
 - fgsea
 - ggplot2
 - GSVA
 - Hmisc
 - magrittr
 - msigdbr
 - org.Hs.eg.db
 - RColorBrewer
 - tibble
 - topGO

## Installation
```R
if(!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("nbay13/biorunr")
```
## Usage
...
