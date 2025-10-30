#' DOTSeq: Differential Analysis of Translation Efficiency and Usage of Open Reading Frames
#'
#' DOTSeq provides a flexible beta-binomial generalized linear model (GLM) 
#' framework for modeling the expected proportion of ribosome profiling 
#' (Ribo-seq) to RNA-seq counts for individual open reading frames (ORFs)
#' relative to other ORFs within the same gene. It also includes a negative 
#' binomial GLM framework for detecting changes in translation efficiency 
#' across experimental conditions.
#'
#' @docType package
#' @name DOTSeq
#' @useDynLib DOTSeq, .registration = TRUE
#' @import GenomicRanges GenomicAlignments GenomeInfoDb
#' @importFrom Rcpp sourceCpp
#' @importFrom methods is
#' @importFrom stats model.matrix
#' @importFrom utils head tail
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
"_PACKAGE"
