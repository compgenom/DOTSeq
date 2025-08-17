# DOTSeq
DOTSeq is an R package for identifying differentially translated open reading frames (ORFs) from ribosome profiling and matched RNA-seq datasets. Unlike most existing tools that operate at gene level, DOTSeq performs analysis at the ORF level, enabling the detection of differential translation efficiency between conditions, and occupancy shifts of ribosomes on ORFs within a single gene. DOTSeq integrates both ribosome profiling and RNA-seq read counts into quasi-binomial generalised linear model (GLM) using a modified design matrix and model fitting formula inspired by Riborex and satuRn. At present, DOTSeq accepts count data generated with featureCounts. Development is underway to extend support for additional quantification tools, including mmquant and HTseq. The package also provides functions for visualisation and exploration of results.

## DEPENDENCIES
* R (>= 4.5.0)
* DEXSeq (>= 1.55.1)
* IRanges (>= 2.43.0)
* GenomicRanges (>= 1.61.1)
* SummarizedExperiment (>= 1.39.1)
* rtracklayer (>= 1.69.1)
* satuRn (>= 1.17.0)
* locfdr (>= 1.1-8)

## INSTALLATION
Please make sure you have the above dependencies installed prior to installing DOTSeq.
Then start R and enter:
```r
  install.packages('devtools')
  library(devtools)
  options(unzip='internal')
  devtools::install_github('compgenom/DOTSeq')
```

## DOCUMENTATION
Please refer to vignettes/DOTSeq.html for how to use DOTSeq.

## CONTACTS AND BUG REPORTS
Chun Shen Lim
chunshen.lim@otago.ac.nz

Gabrielle Chieng
gabrielle.chieng@postgrad.otago.ac.nz

>>>>>>> 4b8498600bd71a8e56750000ec3f18fe43a42059
