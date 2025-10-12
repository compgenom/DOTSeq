#' Perform Differential ORF Translation Analysis with DOTSeq
#'
#' DOTSeq is a unified framework for modeling differential ORF translation using ribosome profiling and RNA-seq data.
#' It includes a novel beta-binomial modeling approach for Differential ORF Usage (DOU), and integrates DESeq2-based
#' modeling for Differential Translation Efficiency (DTE). This wrapper function runs the full DOTSeq workflow,
#' including data loading, ORF filtering, normalization, and model fitting for one or both modules.
#' 
#' @param dotseq_dataset A named \code{list} containing pre-constructed DOTSeq input objects.
#' This list must include:
#' \describe{
#'   \item{\code{sumExp}}{A \code{RangedSummarizedExperiment} object containing pre-filtered raw counts, 
#'   sample metadata, and ORF-level annotations. This object is used for modeling Differential ORF Usage (DOU) 
#'   within the DOTSeq framework.}
#'   \item{\code{dds}}{A \code{DESeqDataSet} object used for modeling Differential Translation Efficiency (DTE) 
#'   within the DOTSeq framework via DESeq2.}
#' }
#' If \code{dotseq_dataset} is provided, the function will skip raw input parsing and use these objects directly. 
#' Otherwise, all of \code{count_table}, \code{condition_table}, \code{flattened_gtf}, and \code{bed} must be supplied 
#' to construct the dataset from scratch.
#' @param count_table Path to a count table file or a data frame. Must contain columns:
#'   \code{Geneid}, \code{Chr}, \code{Start}, \code{End}, \code{Strand}, \code{Length},
#'   plus one column per sample.
#' @param condition_table Path to a sample metadata file or a data frame. Must include columns:
#'   \code{run}, \code{strategy}, \code{condition}, \code{replicate}.
#' @param flattened_gtf Optional path to a flattened GFF/GTF file containing exon definitions.
#' @param bed Path to a BED file with ORF annotations.
#' @param formula A formula object specifying the design, e.g., \code{~ condition * strategy}.
#' @param modules Character vector specifying which DOTSeq modules to run.
#'   Options include \code{"DOU"} (Differential ORF Usage) and \code{"DTE"} (Differential Translation Efficiency).
#'   Both are components of the DOTSeq framework. Default is \code{c("DOU", "DTE")}, which runs both.
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
#' @param dispersion_modeling String specifying the dispersion modeling approach for DOU.
#'   Options include \code{"auto"}, \code{"shared"}, or \code{"custom"} (default: \code{"auto"}).
#' @param dispformula Optional formula object for custom dispersion modeling in DOU.
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test comparing full vs reduced models
#'   to assess translation-specific effects in DOU (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics in DOU, including tests for overdispersion,
#'   zero inflation, and residual properties (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization in DOU.
#'   If \code{NULL}, parallelism is disabled. Default: \code{list(n = 4L, autopar = TRUE)}.
#'   Parallelization provides a noticeable speed-up only when \code{optimizers = TRUE}. 
#'   When using only the default optimizer (\code{nlminb}), parallelism has little impact on speed.
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers
#'   in \code{glmmTMB}: \code{nlminb}, \code{bobyqa}, and \code{optim} (default: \code{FALSE}).
#' @param nullweight Numeric. Prior weight on the null hypothesis for empirical Bayes shrinkage in DOU.
#'   Higher values yield more conservative lfsr estimates. Default is \code{500}.
#' @param contrasts_method Character string specifying the method for post hoc contrasts in DOU.
#'   Default is \code{"pairwise"}; other methods supported by \code{emmeans} may be used.
#' @param seed Optional integer to set the random seed for reproducibility in model fitting (default: \code{NULL}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages and runtime summaries (default: \code{TRUE}).
#'
#' @return A named \code{list} containing:
#' \describe{
#'   \item{sumExp}{A \code{SummarizedExperiment} object containing filtered raw counts, sample metadata, and DOU results. 
#'   Used for modeling Differential ORF Usage within the DOTSeq framework.}
#'   \item{dds}{A \code{DESeqDataSet} object containing normalized counts and metadata, 
#'   used for modeling Differential Translation Efficiency within the DOTSeq framework via DESeq2.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- DOTSeq(
#'   count_table = "counts.txt",
#'   condition_table = "samples.txt",
#'   flattened_gtf = "flattened.gff",
#'   bed = "orfs.bed",
#'   min_count = 1,
#'   stringent = TRUE,
#'   seed = 42,
#'   verbose = TRUE
#' )
#' head(result$sumExp)
#' }
#'
#' @importFrom DESeq2 DESeq design resultsNames lfcShrink
#' @importFrom SummarizedExperiment rowRanges assay rowData rowData<-
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#' 
#' @export
#' 
DOTSeq <- function(
    dotseq_dataset = NULL,
    count_table = NULL, 
    condition_table = NULL, 
    flattened_gtf = NULL, 
    bed = NULL, 
    formula = ~ condition * strategy,
    modules = c("DOU", "DTE"),
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
    nullweight = 500,
    contrasts_method = "pairwise",
    seed = NULL,
    verbose = TRUE
    ) {
  
  if (!("DOU" %in% modules) && !("DTE" %in% modules)) {
    stop("Please specify at least one module: 'DOU' and/or 'DTE'.")
  }
  
  if (verbose) {
    start_dou <- Sys.time()
  }
  
  if (!is.null(dotseq_dataset)) {
    if (!is.list(dotseq_dataset) || !"sumExp" %in% names(dotseq_dataset)) {
      stop("'dotseq_dataset' must be a list containing at least a 'sumExp' RangedSummarizedExperiment.")
    }
    dot <- dotseq_dataset
  } else {
    # Check that all required raw inputs are provided
    if (any(sapply(list(count_table, condition_table, flattened_gtf, bed), is.null))) {
      stop("Either provide a 'dotseq_dataset' object or all of 'count_table', 'condition_table', 'flattened_gtf', and 'bed'.")
    }
    dot <- DOTSeqDataSet(
      count_table = count_table, 
      condition_table = condition_table, 
      flattened_gtf = flattened_gtf, 
      bed = bed, 
      formula = formula,
      target = target,
      baseline = baseline,
      min_count = min_count, 
      stringent = stringent, 
      verbose = verbose
    )
  }
  
  if ("DOU" %in% modules) {
    if ("DOUResults" %in% names(rowData(dot$sumExp)) || "interaction_results" %in% names(metadata(dot$sumExp))) {
      message("skipping Differential ORF Usage (DOU) analysis; model fitting has already been performed on this object. To re-run DOU, please provide a fresh DOTSeqDataSet without fitted results.")
    } else {
      # Total input ORFs
      nrow_input <- nrow(rowData(dot$sumExp))
      
      dot$sumExp <- dot$sumExp[rowRanges(dot$sumExp)$is_kept == TRUE, ]
      
      # Kept ORFs
      nrow_kept <- nrow(rowData(dot$sumExp))
      
      if (verbose) {
        message("starting Differential ORF Usage (DOU) analysis")
      }
      
      rowData(dot$sumExp)[["DOUResults"]] <- fitDOU(
        countData = assay(dot$sumExp),
        orf2gene = rowData(dot$sumExp), 
        anno = colData(dot$sumExp),
        formula = metadata(dot$sumExp)$formula,
        emm_specs = metadata(dot$sumExp)$emm_specs,
        target = target,
        dispersion_modeling = dispersion_modeling, 
        dispformula = dispformula,
        lrt = lrt,
        diagnostic = diagnostic,
        parallel = parallel,
        optimizers = optimizers,
        seed = seed,
        verbose = verbose
      )
      
      dot$sumExp <- testDOU(
        dot$sumExp,
        emm_specs = metadata(dot$sumExp)$emm_specs,
        contrasts_method = contrasts_method,
        nullweight = nullweight,
        verbose = verbose
      )
      
      interaction_results <- metadata(dot$sumExp)$interaction_results
      all_contrasts <- unique(interaction_results$contrast)
      
      if (verbose) {
        end_dou <- Sys.time()
        
        elapsed_dou <- runtime(end_dou, start_dou)
        
        if (!is.null(elapsed_dou$mins)) {
          message(sprintf("DOU runtime: %d mins %.3f secs", elapsed_dou$mins, elapsed_dou$secs))
        } else {
          message(sprintf("DOU runtime: %.3f secs", elapsed_dou$secs))
        }
        
        # DOU summary
        for (c_name in all_contrasts) {
          nrow_fitted <- nrow(interaction_results[interaction_results$contrast == c_name, ])
          message("DOU model fitting summary: ", c_name)
          message("  models fitted: ", paste0(nrow_fitted))
          message("  fit errors (non-convergence or invalid Hessian): ", paste0(nrow_kept - nrow_fitted))
          message("  not fitted (filtered out): ", paste0(nrow_input - nrow_kept))
        }
      }
    }
    
  } else {
    if (verbose) {
      message("DOU module not selected; returning SummarizedExperiment object without model fitting. This object can still be used as input for DOTSeq() to run the DOU module later.")
    }
  }
  
  if ("DTE" %in% modules) {
    if ("interaction_results" %in% names(metadata(dot$dds))) {
      message("skipping Differential Translation Efficiency (DTE) analysis: model fitting has already been performed on this object. To re-run DTE, please provide a fresh DOTSeqDataSet without fitted results")
    } else {
      if (verbose) {
        message("starting Differential Translation Efficiency (DTE) analysis")
        start_dte <- Sys.time()
      }
      
      # Run the DESeq2 analysis
      dot$dds <- DESeq(dot$dds)
      
      terms <- resultsNames(dot$dds)
      matched_term <- terms[grepl("\\.", terms)]
      
      if (verbose) {
        message("starting post hoc analysis")
      }

      contrast_vectors_list <- contrast_vectors(dot$dds)

      contrast_results <- list()
      for (c_name in names(contrast_vectors_list)) {
        if (verbose) {
          message("performing empirical Bayesian shrinkage on the effect size for ", c_name)
        }

        contrast_results_df <- lfcShrink(dot$dds, contrast = contrast_vectors_list[[c_name]], type = "ashr", quiet = TRUE)
        contrast_results_df$orf_id <- rownames(contrast_results_df)
        rownames(contrast_results_df) <- NULL
        contrast_results_df$contrast <- c_name
        contrast_results <- c(contrast_results, contrast_results_df)
      }
      contrast_results <- do.call(rbind, contrast_results)

      metadata(dot$dds)$interaction_results <- contrast_results
      
      all_contrasts <- unique(contrast_results$contrast)
      
      if (verbose) {
        end_dte <- Sys.time()

        elapsed_dte <- runtime(end_dte, start_dte)

        if (!is.null(elapsed_dte$mins)) {
          message(sprintf("DTE runtime: %d mins %.3f secs", elapsed_dte$mins, elapsed_dte$secs))
        } else {
          message(sprintf("DTE runtime: %.3f secs", elapsed_dte$secs))
        }

        # DTE summary
        for (c_name in all_contrasts) {
          nrow_input <- nrow(contrast_results[contrast_results$contrast == c_name, ])
          nrow_fitted <- nrow(contrast_results[!is.na(contrast_results$padj) & contrast_results$contrast == c_name, ])
          
          message("DTE model fitting summary: ", c_name)
          message("  models fitted: ", paste0(nrow_fitted))
          message("  padj = NA (filtered or flagged): ", paste0(nrow_input - nrow_fitted))
        }
      }
    }
    
  } else {
    if (verbose) {
      message("DTE module not selected; returning DESeqDataSet object without model fitting. This object can still be used as input for DOTSeq() to run the DTE module later.")
    }
  }

  return(dot)
}


