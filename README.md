# DOTSeq

`DOTSeq` is an R package for identifying **differentially translated open reading frames (ORFs)** from ribosome profiling (Ribo-seq) and matched RNA-seq datasets. Unlike traditional gene-level approaches, `DOTSeq` performs analysis at the **ORF level**, enabling detection of:

- **Differential ORF Usage (DOU)** — changes in ORF usage within the same gene.
- **Differential Translation Efficiency (DTE)** — changes in ribosome loading relative to RNA level across conditions.

`DOTSeq` models Ribo-seq and RNA-seq read counts using a **beta-binomial generalised linear model (GLM)** implemented via [`glmmTMB`](https://cran.r-project.org/web/packages/glmmTMB/index.html). It supports experimental designs with multiple conditions, and uses an interaction term (`condition:strategy`) to isolate translation-specific effects.

Post hoc contrasts are computed using [`emmeans`](https://cran.r-project.org/web/packages/emmeans/index.html), and empirical Bayes shrinkage is applied via [`ashr`](https://cran.r-project.org/web/packages/ashr/index.html).


## DEPENDENCIES
* R (>= 4.5.0)
* biomaRt (>=2.65.0)
* DEXSeq (>= 1.55.1)
* IRanges (>= 2.43.0)
* GenomicRanges (>= 1.61.1)
* SummarizedExperiment (>= 1.39.1)
* rtracklayer (>= 1.69.1)

## INSTALLATION
Please ensure the dependencies listed above are installed using the following steps before installing `DOTSeq`:
```r
# Create a directory for R packages if not already
dir.create(file.path(Sys.getenv("HOME"), "R/4.5"), showWarnings = TRUE, recursive = TRUE)

# Install BiocManager if not already available
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = file.path(Sys.getenv("HOME"), "R/4.5"))

# Initialise usage of Bioconductor devel version
BiocManager::install(version = "devel")

# Install required Bioconductor packages
BiocManager::install(c(
  "biomaRt",
  "IRanges",
  "GenomicRanges",
  "SummarizedExperiment",
  "rtracklayer",
  "DEXSeq"), lib = file.path(Sys.getenv("HOME"), "R/4.5"))

# Install devtools if not already available
install.packages("devtools")
library(devtools)
options(unzip = "internal")

# Install DOTSeq from GitHub
devtools::install_github("compgenom/DOTSeq")
```

## DOCUMENTATION
Please refer to [vignettes](https://github.com/compgenom/DOTSeq/tree/main/vignettes) for how to use `DOTSeq`.

## CONTACTS AND BUG REPORTS
- Chun Shen Lim: chunshen [dot] lim [at] otago [dot] ac [dot] nz
- Gabrielle Chieng: gabrielle [dot] chieng [at] postgrad [dot] otago [dot] ac [dot] nz


