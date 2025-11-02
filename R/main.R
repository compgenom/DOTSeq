#' Perform Differential ORF Translation Analysis with DOTSeq
#'
#' @description
#' \pkg{DOTSeq} is a statistical framework for modeling Differential ORF 
#' Translation using ribosome profiling (Ribo-seq) and RNA-seq data. It 
#' supports two modes of input:
#' 
#' - A named list of raw input components: \code{count_table}, 
#' \code{condition_table}, \code{flattened_gtf}, and \code{bed}.
#' - A pre-constructed \code{DOTSeqDataSets} object containing processed data.
#' 
#' The function automatically detects the input type and proceeds with the 
#' appropriate workflow. It performs ORF-level filtering, model fitting, 
#' post hoc contrasts, and adaptive shrinkage of effect sizes. Plotting and 
#' downstream analysis are handled separately via the \code{\link{plotDOT}} function.
#' 
#' @seealso \code{\link{DOTSeqDataSets}}, \code{\link{fitDOU}}, 
#' \code{\link{testDOU}}, \code{\link{plotDOT}}
#'
#' @param datasets Either:
#' \describe{
#'     \item{\code{DOTSeqDataSets} object}{
#'     A pre-constructed \code{\link{DOTSeqDataSets-class}} object created using 
#'     \code{\link{DOTSeqDataSets}}. It must include:
#'         \describe{
#'             \item{\code{DOU}}{
#'                 A \code{\link{DOUData-class}} object containing filtered
#'                 raw counts, sample metadata, and ORF-level annotations.
#'             }
#'             \item{\code{DTE}}{
#'                 A \code{\link[DESeq2]{DESeqDataSet-class}} object used for 
#'                 modeling DTE via DESeq2.
#'             }
#'         }
#'     If a \code{DOTSeqDataSets} object is provided, the function skips raw input
#'     parsing and uses these objects directly.
#'     }
#' }
#'
#' @param formula A formula object specifying the design.
#' Default is \code{~ condition * strategy}.
#'
#' @param modules Character vector specifying which \pkg{DOTSeq} modules 
#' to run. Options include \code{"DOU"} and \code{"DTE"}.
#' Default is \code{c("DOU", "DTE")}.
#'
#' @param target Character string specifying the non-reference condition 
#' level to extract the corresponding interaction term.
#' Default is \code{NULL}.
#'
#' @param baseline Character string specifying the desired reference level.
#' Default is \code{NULL}.
#'
#' @param min_count Minimum count threshold for filtering ORFs.
#' Default is \code{1}.
#'
#' @param stringent Logical or \code{NULL}; determines the filtering strategy:
#' \describe{
#'     \item{\code{TRUE}}{
#'         Keep ORFs where all replicates in at least one condition pass
#'         \code{min_count}.
#'     }
#'     \item{\code{FALSE}}{
#'         Keep ORFs where all replicates in at least one condition-strategy
#'         group pass \code{min_count}.
#'     }
#'     \item{\code{NULL}}{
#'         Keep ORFs where total counts across replicates pass
#'         \code{min_count}.
#'     }
#' }
#'
#' @param dispersion_modeling String specifying the dispersion modeling
#' approach for DOU. Options include \code{"auto"}, \code{"shared"},
#' or \code{"custom"}. Default is \code{"auto"}.
#'
#' @param dispformula Optional formula object for custom dispersion modeling.
#'
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test
#' comparing full vs reduced models in DOU.
#' Default is \code{FALSE}.
#'
#' @param diagnostic Logical; if \code{TRUE}, enables model diagnostics in DOU.
#' Default is \code{FALSE}.
#'
#' @param parallel A list passed to \code{\link[glmmTMB]{glmmTMBControl}} to configure parallel
#' optimization in DOU. Default is \code{list(n = 4L, autopar = TRUE)}.
#'
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization
#' using multiple optimizers in \code{\link[glmmTMB]{glmmTMBControl}}.
#' Default is \code{FALSE}.
#'
#' @param nullweight Numeric. Prior weight on the null hypothesis for empirical
#' Bayes shrinkage in DOU. Default is \code{500}.
#'
#' @param contrasts_method Character string specifying the method for post hoc
#' contrasts in DOU. Default is \code{"revpairwise"}.
#'
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' Default is \code{TRUE}.
#'
#' @return A \code{DOTSeqDataSets} object containing:
#' \describe{
#'     \item{\code{DOU}}{
#'         A \code{\link{DOUData-class}} object with DOU results.
#'     }
#'     \item{\code{DTE}}{
#'         A \code{\link{DTEData-class}} object with DTE results.
#'     }
#' }
#'
#' @importFrom DESeq2 DESeq design resultsNames lfcShrink
#' @importFrom SummarizedExperiment rowRanges assay rowData rowData<-
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#'
#' @export
#'
#' @examples
#' # Load test data
#' dir <- system.file("extdata", package = "DOTSeq")
#'
#' cnt <- read.table(
#'     file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
#'     header = TRUE,
#'     comment.char = "#"
#' )
#' names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
#'
#' gtf <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#'
#' # Use raw input list
#' raw <- list(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' 
#' d <- DOTSeq(datasets = raw, modules = "DTE")
#'
#' show(d)
#' 
#' # Get the DOUData object
#' dou <- getDOU(d)
#' 
#' # Subset DOUData and edit DOTSeqDataSets in place
#' set.seed(42)
#' random_rows <- sample(seq_len(nrow(dou)), size = 50)
#' getDOU(d) <- dou[random_rows, ]
#' 
#' # Run the DOU module
#' d <- DOTSeq(datasets = d)
#' 
#' # Create a DOTSeqDataSets object and use it as input
#' d <- DOTSeqDataSets(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' d <- DOTSeq(datasets = d, modules = "DTE")
#'
#' @references
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, 
#' C. W., Nielsen, A., Skaug, H. J., Mächler, M. and Bolker, B. M. (2017).
#' glmmTMB balances speed and flexibility among packages for zero-inflated
#' generalized linear mixed modeling. The R Journal, 378–400.
#' DOI: 10.32614/RJ-2017-066
#'
#' Love, M. I., Huber, W., Anders, S. (2014).
#' Moderated estimation of fold change and dispersion for RNA-seq data with 
#' DESeq2. Genome Biology, 15:550. DOI: 10.1186/s13059-014-0550-8
#'
#' Lenth, R., Piaskowski, J. (2025).
#' emmeans: Estimated Marginal Means, aka Least-Squares Means.
#' R package version 2.0.0. \url{https://rvlenth.github.io/emmeans/}
#'
#' Stephens, M. (2016).
#' False discovery rates: a new deal. Biostatistics, 18(2).
#' DOI: 10.1093/biostatistics/kxw041
#'
#' Hartig, F. (2025).
#' DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed)
#' Regression Models. R package version 0.4.7.
#' \url{https://github.com/florianhartig/dharma}
#'
DOTSeq <- function(
        datasets,
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
    
    if (diagnostic && !requireNamespace("DHARMa", quietly = TRUE)) {
        stop(
            "Model diagnostics require the 'DHARMa' package. ", 
            "Please install it by running: install.packages('DHARMa')"
        )
    }
    
    if (!("DOU" %in% modules) && !("DTE" %in% modules)) {
        stop("Please specify at least one module: 'DOU' and/or 'DTE'.")
    }

    if (verbose) {
        start_dou <- Sys.time()
    }
    
    if (inherits(datasets, "DOTSeqDataSets")) {
        dot <- datasets
        if (!is(dot, "DOTSeqDataSets")) {
            stop("'datasets' must be a DOTSeqDataSets object.")
        }
        if (!is(getDOU(dot), "DOUData")) {
            stop("The 'DOU' slot must be a DOUData object.")
        }
        if (!inherits(getDTE(dot), "DTEData") && !inherits(getDTE(dot), "DESeqDataSet")) {
            stop("The 'DTE' slot must be a DTEData or DESeqDataSet object.")
        }

    } else {
        stop("'datasets' must be either a DOTSeqDataSets object.")
    }
    
    dou <- getDOU(dot)

    if ("DOU" %in% modules) {
        if ("DOUResults" %in% names(rowData(dou))) {
            message(
                "skipping Differential ORF Usage (DOU) analysis; ", 
                "model fitting has already been performed on this object. ", 
                "To re-run DOU, please provide a fresh DOTSeqDataSets ", 
                "without fitted results."
            )
        } else {
            # Total input ORFs
            nrow_input <- nrow(rowData(dou))

            dou <- dou[rowRanges(dou)$is_kept == TRUE, ]

            # Kept ORFs
            nrow_kept <- nrow(rowData(dou))
            
            if (nrow_kept == 0) {
                stop("No data after filtering. Please check the count table.")
            }

            if (verbose) {
                message("starting Differential ORF Usage (DOU) analysis")
            }

            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = fmla(dou),
                specs = specs(dou),
                dispersion_modeling = dispersion_modeling,
                dispformula = dispformula,
                lrt = lrt,
                diagnostic = diagnostic,
                parallel = parallel,
                optimizers = optimizers,
                verbose = verbose
            )

            dou <- testDOU(
                data = dou,
                contrasts_method = contrasts_method,
                nullweight = nullweight,
                verbose = verbose
            )

            interaction_results <- getContrasts(dou, type = "interaction")
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
    
    dte <- getDTE(dot)

    if ("DTE" %in% modules) {
        if (inherits(getDTE(dot), "DTEData")) {
            message(
                "skipping Differential Translation Efficiency (DTE) analysis: ", 
                "model fitting has already been performed on this object. ", 
                "To re-run DTE, please provide a fresh DOTSeqDataSets ", 
                "without fitted results"
            )
        } else {
            if (verbose) {
                message("starting Differential Translation Efficiency (DTE) analysis")
                start_dte <- Sys.time()
            }

            # Run the DESeq2 analysis
            dte <- DESeq(dte)

            terms <- resultsNames(dte)
            matched_term <- terms[grepl("\\.", terms)]

            if (verbose) {
                message("starting post hoc analysis")
            }

            contrast_vectors_list <- contrast_vectors(dte)

            contrast_results <- list()
            for (c_name in names(contrast_vectors_list)) {
                # if (verbose) {
                #     message("performing empirical Bayesian shrinkage on the effect size for ", c_name)
                # }

                contrast_results_df <- lfcShrink(dte, contrast = contrast_vectors_list[[c_name]], type = "ashr", quiet = TRUE)
                contrast_results_df$orf_id <- rownames(contrast_results_df)
                rownames(contrast_results_df) <- NULL
                contrast_results_df$contrast <- c_name
                contrast_results <- c(contrast_results, contrast_results_df)
            }
            contrast_results <- do.call(rbind, contrast_results)

            dte <- new(
                "DTEData",
                dte,
                formula = design(dte),
                specs = metadata(dte)$specs,
                interaction = contrast_results
            )

            metadata(dte)$specs <- NULL

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

    return(new("DOTSeqDataSets", DOU = dou, DTE = dte))
}
