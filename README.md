# DOTSeq
DOTSeq is an R package for identifying differentially translated open reading frames (ORFs) from ribosome profiling and matched RNA-seq datasets. Unlike most existing tools that operate at gene level, DOTSeq performs analysis at the ORF level, enabling the detection of differential translation efficiency between conditions, and occupancy shifts of ribosomes on ORFs within a single gene. DOTSeq integrates both ribosome profiling and RNA-seq read counts into quasi-binomial generalised linear model (GLM) using a modified design matrix and model fitting formula inspired by Riborex and satuRn. At present, DOTSeq accepts count data generated with featureCounts. Development is underway to extend support for additional quantification tools, including mmquant and HTseq. The package also provides functions for visualisation and exploration of results.

## DEPENDENCIES
* R (>= 4.5.0)
* biomaRt (>=2.65.0)
* DEXSeq (>= 1.55.1)
* IRanges (>= 2.43.0)
* GenomicRanges (>= 1.61.1)
* SummarizedExperiment (>= 1.39.1)
* rtracklayer (>= 1.69.1)
* satuRn (>= 1.17.0)
* locfdr (>= 1.1-8)

## INSTALLATION
Please ensure the dependencies listed above are installed using the following steps before installing DOTSeq:
```r
# Install BiocManager if not already available
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Initialise usage of Bioconductor devel version
BiocManager::install(version = "devel")

# Install required Bioconductor packages
BiocManager::install(c(
  "biomaRt@2.65.0",
  "IRanges@2.43.0",
  "GenomicRanges@1.61.1",
  "SummarizedExperiment@1.39.1",
  "rtracklayer@1.69.1",
  "DEXSeq@1.55.1",
  "satuRn@1.17.0",
  "locfdr@1.1-8"
))

# When prompted with "Update all/some/none? [a/s/n]:", enter 'n' to skip updates.

# Install devtools if not already available
install.packages("devtools")
library(devtools)
options(unzip = "internal")

# Install DOTSeq from GitHub
devtools::install_github("compgenom/DOTSeq")
```

## DOCUMENTATION
Please refer to [vignettes](https://github.com/compgenom/DOTSeq/tree/main/vignettes) for how to use DOTSeq.

## CONTACTS AND BUG REPORTS
- Chun Shen Lim: chunshen [dot] lim [at] otago [dot] ac [dot] nz
- Gabrielle Chieng: gabrielle [dot] chieng [at] postgrad [dot] otago [dot] ac [dot] nz


