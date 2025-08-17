# DOTSeq
DOTSeq is a R package for identifying differentially translated open reading frames (ORFs) from paired Ribo- and RNAseq data. Most existing tools operate at gene level, DOTSeq performs analysis at the ORF level, enabling the detection of differential translation between conditions, and occupancy shifts of ribosomes on ORFs within a single gene.

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

