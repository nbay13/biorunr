# biorunR

An R package for running bioinformatic analysis. Includes wrapper functions for performing differential expression (DESeq2), geneset enrichment (GSEA, GSVA/ssGSEA, topGO), and data vizualization (ggplot2).

## Dependencies
 - classInt
 - DESeq2 (& ashr)
 - dplyr
 - fgsea
 - ggplot2
 - GSVA
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