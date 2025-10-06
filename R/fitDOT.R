#' Fit DOTSeq Differential ORF Translation Models
#'
#' This function performs the complete DOTSeq analysis pipeline for differential ORF translation:
#' loading count data, aligning with sample metadata, filtering ORFs, normalizing counts,
#' calculating translational efficiency (TE), and fitting beta-binomial and negative binomial
#' generalized linear models using \code{DOTSeq::fitDOU} and \code{DESeq2::DESeq}.
#'
#' @param count_table Path to a count table file or a data frame. Must contain columns:
#'   \code{Geneid}, \code{Chr}, \code{Start}, \code{End}, \code{Strand}, \code{Length},
#'   plus one column per sample.
#' @param condition_table Path to a sample metadata file or a data frame. Must include columns:
#'   \code{run}, \code{strategy}, \code{condition}, \code{replicate}.
#' @param flattened_gtf Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param formula A formula object specifying the design, e.g., \code{~ condition * strategy}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term.
#'   Contrasted against the baseline condition (default: \code{NULL}).
#' @param baseline Character string specifying the desired reference level (default: \code{NULL}).
#' @param min_count Minimum count threshold for filtering ORFs (default: \code{1}).
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#'   \describe{
#'     \item{\code{TRUE}}{Keep ORFs where all replicates in at least one condition pass \code{min_count}.}
#'     \item{\code{FALSE}}{Keep ORFs where all replicates in at least one condition-strategy group pass \code{min_count}.}
#'     \item{\code{NULL}}{Keep ORFs where total counts across replicates pass \code{min_count}.}
#'   }
#' @param dispersion_modeling String specifying the dispersion modeling approach.
#'   Options include \code{"auto"}, \code{"shared"}, or \code{"custom"} (default: \code{"auto"}).
#' @param dispformula Optional formula object for custom dispersion modeling.
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test comparing full vs reduced models
#'   to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics including tests for overdispersion,
#'   zero inflation, and residual properties (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization.
#'   If \code{NULL}, parallelism is disabled. Default: \code{list(n = 4L, autopar = TRUE)}.
#'   Parallelization provides a noticeable speed-up only when \code{optimizers = TRUE}. 
#'   When using only the default optimizer (\code{nlminb}), parallelism has little impact on speed.
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers
#'   in \code{glmmTMB}: \code{nlminb}, \code{bobyqa}, and \code{optim} (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#'
#' @return A named \code{list} containing:
#' \describe{
#'   \item{sumExp}{\code{SummarizedExperiment} object containing pre-filtered normalized counts and sample metadata.}
#'   \item{dds}{\code{DESeqDataSet} object used for modeling differential ORF translation.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- fitDOT(
#'   count_table = "counts.txt",
#'   condition_table = "samples.txt",
#'   flattened_gtf = "flattened.gff",
#'   bed = "orfs.bed",
#'   min_count = 1,
#'   stringent = TRUE,
#'   verbose = TRUE
#' )
#' head(result$sumExp)
#' }
#'
#' @importFrom DESeq2 DESeq results design
#' @importFrom SummarizedExperiment rowRanges assay rowData rowData<-
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#' @importFrom utils capture.output 
#' 
#' @export
#' 
fitDOT <- function(count_table, 
                   condition_table, 
                   flattened_gtf, 
                   bed, 
                   formula = ~ condition * strategy,
                   target = NULL,
                   baseline = NULL,
                   min_count = 1, 
                   stringent = TRUE, 
                   dispersion_modeling = "auto",
                   dispformula = NULL,
                   lrt = FALSE,
                   diagnostic = FALSE,
                   parallel = list(n=4L, autopar=TRUE),
                   optimizers = FALSE,
                   verbose = TRUE) {
  
  
  dot <- DOTSeq::DOTSeqDataSet(count_table, 
                               condition_table, 
                               flattened_gtf, 
                               bed, 
                               formula = formula,
                               target = target,
                               baseline = baseline,
                               min_count = min_count, 
                               stringent = stringent, 
                               verbose = verbose)

  sumExp_filtered <- dot$sumExp_raw[rowRanges(dot$sumExp_raw)$is_kept == TRUE, ]
  
  if (isTRUE(verbose)) {
    message(" - Fit beta-binomial generalised linear models")
    message(" - Use ", parallel$n, " threads")
    start_dou <- Sys.time()
  }
  
  sumExp <- fitDOU(object = sumExp_filtered,
                           formula = metadata(sumExp_filtered)$formula,
                           target = target,
                           dispersion_modeling= dispersion_modeling,
                           dispformula = dispformula,
                           lrt = lrt,
                           diagnostic = diagnostic,
                           parallel = parallel,
                           optimizers = optimizers,
                           verbose = verbose)
  
  # L <- contrastMatrix(sumExp, fmla, baseline = baseline)
  
  if (verbose) {
    end_dou <- Sys.time()
    
    elapsed_dou <- runtime(end_dou, start_dou)
    
    if (!is.null(elapsed_dou$mins)) {
      message(sprintf(" - DOU runtime: %d mins %.3f secs", elapsed_dou$mins, elapsed_dou$secs))
    } else {
      message(sprintf(" - DOU runtime: %.3f secs", elapsed_dou$secs))
    }
    
    message(" - Start differential translation efficiency analysis")
    start_dte <- Sys.time()
  }
  
  
  # Run the DESeq2 analysis
  dds <- DESeq(dot$dds)
  
  if (verbose) {
    end_dte <- Sys.time()
    
    elapsed_dte <- runtime(end_dte, start_dte)
    
    if (!is.null(elapsed_dte$mins)) {
      message(sprintf(" - DTE runtime: %d mins %.3f secs", elapsed_dte$mins, elapsed_dte$secs))
    } else {
      message(sprintf(" - DTE runtime: %.3f secs", elapsed_dte$secs))
    }
    
    # DOU summary
    dou_res <- extract_results(sumExp)
    msg <- capture.output(print(table(dou_res$model_type)))
    msg <- msg[2:length(msg)]
    message(paste(" - DOU model fitting summary:\n  ", paste(msg, collapse = "\n")))
    
    # DTE summary
    dte_res <- results(dds)
    dte_res$model_type <- ifelse(is.na(dte_res$padj), "NA", "nbinom")
    
    msg_dte <- capture.output(print(table(dte_res$model_type)))
    msg_dte <- msg_dte[2:length(msg_dte)]
    message(paste(" - DTE model fitting summary:\n", paste(msg_dte, collapse = "\n")))
  }
  
  return(list(sumExp = sumExp,
              dds = dds))
}



