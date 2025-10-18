#' Perform Differential ORF Translation Analysis with DOTSeq
#'
#' @description
#' DOTSeq is a statistical framework for modeling Differential ORF 
#' Translation using ribosome profiling (Ribo-seq) and RNA-seq data. It 
#' includes a novel beta-binomial modeling approach for detecting 
#' Differential ORF Usage (DOU), and integrates DESeq2-based modeling for 
#' Differential Translation Efficiency (DTE). This wrapper function 
#' optionally executes one or both modules, including data loading, 
#' ORF-level filtering, model fitting, post hoc contrasts, and adaptive 
#' shrinkage of effect sizes. Plotting and downstream analysis are handled 
#' separately via the \code{\link{plotDOT}} function. 
#' 
#' @seealso \code{\link{DOTSeqDataSet}}, \code{\link{fitDOU}}, 
#' \code{\link{testDOU}}, \code{\link{plotDOT}}
#'
#' @param dotseq_dataset A named \code{list} containing pre-constructed DOTSeq
#'     input objects. This list must include:
#'     \describe{
#'         \item{\code{sumExp}}{
#'             A \code{RangedSummarizedExperiment} object containing filtered
#'             raw counts, sample metadata, and ORF-level annotations.
#'         }
#'         \item{\code{dds}}{
#'             A \code{DESeqDataSet} object used for modeling DTE via DESeq2.
#'         }
#'     }
#'     If \code{dotseq_dataset} is provided, the function skips raw input
#'     parsing and uses these objects directly.
#'
#' @param count_table Path to a count table file or a data frame. Must contain
#'     columns: \code{Geneid}, \code{Chr}, \code{Start}, \code{End},
#'     \code{Strand}, \code{Length}, plus one column per sample.
#'
#' @param condition_table Path to a sample metadata file or a data frame.
#'     Must include columns: \code{run}, \code{strategy}, \code{condition},
#'     \code{replicate}.
#'
#' @param flattened_gtf Optional path to a flattened GFF/GTF file.
#'
#' @param bed Path to a BED file with ORF annotations.
#'
#' @param formula A formula object specifying the design.
#'     Default is \code{~ condition * strategy}.
#'
#' @param modules Character vector specifying which DOTSeq modules to run.
#'     Options include \code{"DOU"} and \code{"DTE"}.
#'     Default is \code{c("DOU", "DTE")}.
#'
#' @param target Character string specifying the non-reference condition level
#'     to extract the corresponding interaction term.
#'     Default is \code{NULL}.
#'
#' @param baseline Character string specifying the desired reference level.
#'     Default is \code{NULL}.
#'
#' @param min_count Minimum count threshold for filtering ORFs.
#'     Default is \code{1}.
#'
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#'     \describe{
#'         \item{\code{TRUE}}{
#'             Keep ORFs where all replicates in at least one condition pass
#'             \code{min_count}.
#'         }
#'         \item{\code{FALSE}}{
#'             Keep ORFs where all replicates in at least one condition-strategy
#'             group pass \code{min_count}.
#'         }
#'         \item{\code{NULL}}{
#'             Keep ORFs where total counts across replicates pass
#'             \code{min_count}.
#'         }
#'     }
#'
#' @param dispersion_modeling String specifying the dispersion modeling
#'     approach for DOU. Options include \code{"auto"}, \code{"shared"},
#'     or \code{"custom"}. Default is \code{"auto"}.
#'
#' @param dispformula Optional formula object for custom dispersion modeling.
#'
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test
#'     comparing full vs reduced models in DOU.
#'     Default is \code{FALSE}.
#'
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics in DOU.
#'     Default is \code{FALSE}.
#'
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel
#'     optimization in DOU. Default is \code{list(n = 4L, autopar = TRUE)}.
#'
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization
#'     using multiple optimizers in \code{glmmTMB}.
#'     Default is \code{FALSE}.
#'
#' @param nullweight Numeric. Prior weight on the null hypothesis for empirical
#'     Bayes shrinkage in DOU. Default is \code{500}.
#'
#' @param contrasts_method Character string specifying the method for post hoc
#'     contrasts in DOU. Default is \code{"revpairwise"}.
#'
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'     Default is \code{TRUE}.
#'
#' @return A named \code{list} containing:
#'     \describe{
#'         \item{\code{sumExp}}{
#'             A \code{SummarizedExperiment} object with DOU results.
#'         }
#'         \item{\code{dds}}{
#'             A \code{DESeqDataSet} object with DTE results.
#'         }
#'     }
#'
#' @importFrom DESeq2 DESeq design resultsNames lfcShrink
#' @importFrom SummarizedExperiment rowRanges assay rowData rowData<-
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#'
#' @export
#'
#' @examples
#' # Load SummarizedExperiment
#' library(SummarizedExperiment)
#'
#' dir <- system.file("extdata", package = "DOTSeq")
#'
#' cnt <- read.table(
#'     file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
#'     header = TRUE,
#'     comment.char = "#"
#' )
#' names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
#'
#' flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#'
#' m <- DOTSeq(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     bed = bed,
#'     modules = "DTE"
#' )
#'
#' names(m)
#' metadata(m$dds)
#'
#' set.seed(42)
#' random_rows <- sample(seq_len(nrow(m$sumExp)), size = 100)
#' m$sumExp <- m$sumExp[random_rows, ]
#'
#' m <- DOTSeq(dotseq_dataset = m)
#' metadata(m$sumExp)
#'
#' m <- DOTSeqDataSet(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     bed = bed
#' )
#'
#' m <- DOTSeq(dotseq_dataset = m, modules = "DTE")
#' metadata(m$dds)
#'
#' @references
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W.,
#' Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017).
#' glmmTMB balances speed and flexibility among packages for zero-inflated
#' generalized linear mixed modeling. The R Journal, 378–400.
#' \url{https://doi.org/10.32614/RJ-2017-066}
#'
#' Love, M. I., Huber, W., Anders, S. (2014).
#' Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#' Genome Biology, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#'
#' Lenth, R., Piaskowski, J. (2025).
#' emmeans: Estimated Marginal Means, aka Least-Squares Means.
#' R package version 2.0.0. \url{https://rvlenth.github.io/emmeans/}
#'
#' Stephens, M. (2016).
#' False discovery rates: a new deal. Biostatistics, 18(2).
#' \url{https://doi.org/10.1093/biostatistics/kxw041}
#'
#' Hartig, F. (2025).
#' DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed)
#' Regression Models. R package version 0.4.7.
#' \url{https://github.com/florianhartig/dharma}
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
    parallel = list(n = 4L, autopar = TRUE),
    optimizers = FALSE,
    nullweight = 500,
    contrasts_method = "revpairwise",
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
            stop(
                "'dotseq_dataset' must be a list containing ", 
                "at least a 'sumExp' RangedSummarizedExperiment."
            )
        }
        dot <- dotseq_dataset
    } else {
        # Check that all required raw inputs are provided
        if (any(vapply(list(count_table, condition_table, flattened_gtf, bed), is.null, logical(1)))) {
            stop(
                "Either provide a 'dotseq_dataset' object ", 
                "or all of 'count_table', 'condition_table', ", 
                "'flattened_gtf', and 'bed'."
            )
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
            message(
                "skipping Differential ORF Usage (DOU) analysis; ", 
                "model fitting has already been performed on this object. ", 
                "To re-run DOU, please provide a fresh DOTSeqDataSet ", 
                "without fitted results."
            )
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
                count_table = assay(dot$sumExp),
                rowdata = rowData(dot$sumExp),
                anno = colData(dot$sumExp),
                formula = metadata(dot$sumExp)$formula,
                emm_specs = metadata(dot$sumExp)$emm_specs,
                dispersion_modeling = dispersion_modeling,
                dispformula = dispformula,
                lrt = lrt,
                diagnostic = diagnostic,
                parallel = parallel,
                optimizers = optimizers,
                verbose = verbose
            )

            dot$sumExp <- testDOU(
                dot$sumExp,
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
                    message("  models fitted: ", nrow_fitted)
                    message("  non-convergence or invalid Hessian: ", nrow_kept - nrow_fitted)
                    message("  not fitted (filtered out): ", nrow_input - nrow_kept)
                }
            }
        }
    } else {
        if (verbose) {
            message(
                "DOU module not selected; ", 
                "returning SummarizedExperiment object without model fitting. ", 
                "This object can still be used as input for DOTSeq() ", 
                "to run the DOU module later."
            )
        }
    }

    if ("DTE" %in% modules) {
        if ("interaction_results" %in% names(metadata(dot$dds))) {
            message(
                "skipping Differential Translation Efficiency (DTE) analysis: ", 
                "model fitting has already been performed on this object. ", 
                "To re-run DTE, please provide a fresh DOTSeqDataSet ", 
                "without fitted results"
            )
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
                    message("  models fitted: ", nrow_fitted)
                    message("  padj = NA (filtered or flagged): ", nrow_input - nrow_fitted)
                }
            }
        }
    } else {
        if (verbose) {
            message(
                "DTE module not selected; returning DESeqDataSet object without ",
                "model fitting. This object can still be used as input for DOTSeq() ",
                "to run the DTE module later."
            )
        }
    }

    return(dot)
}
