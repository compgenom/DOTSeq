#' Compute Differential ORF Usage (DOU) Contrasts Using Estimated Marginal 
#' Means
#'
#' @description
#' Performs Differential ORF Usage (DOU) analysis by computing contrasts
#' between ribosome profiling and RNA-seq modalities using estimated 
#' marginal means (EMMs) from fitted models (via 
#' \code{\link[emmeans]{emmeans}}). Supports both interaction-specific 
#' and strategy-specific contrasts, and applies empirical Bayes shrinkage 
#' via \code{\link[ashr]{ash}} to stabilize effect size estimates.
#'
#' @seealso \code{\link{DOTSeq}}, \code{\link{DOTSeqDataSet}}, 
#' \code{\link{fitDOU}}, \code{\link{plotDOT}}
#' 
#' @param sumExp A SummarizedExperiment object containing fitted model objects,
#'     typically stored in \code{rowData(sumExp)[['DOUResults']]}.
#'
#' @param contrasts_method Character string specifying the method for computing
#'     contrasts. Default is \code{"revpairwise"}.
#'
#' @param nullweight Numeric. Prior weight on the null hypothesis for empirical
#'     Bayes shrinkage. Higher values yield more conservative lfsr estimates.
#'     Default is \code{500}.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'     Default is \code{TRUE}.
#'
#' @return A \code{SummarizedExperiment} object (the input \code{sumExp} or
#'     \code{m$sumExp}) with two new \code{S4Vectors::DataFrame} objects 
#'     stored in its \code{metadata} slot. These tables contain the 
#'     long-format results for all computed contrasts:
#'     \describe{
#'         \item{\code{interaction_results}}{
#'             A long-format \code{S4Vectors::DataFrame} containing the
#'             Differential ORF Usage (DOU) effect sizes (Ribo-seq contrast
#'             minus RNA-seq contrast) for all computed interaction contrasts.
#'             Columns include \code{contrast}, and the full set of shrunken
#'             and unshrunken metrics (e.g., \code{betahat}, \code{sebetahat},
#'             \code{WaldPadj}, \code{PosteriorMean}, \code{lfsr}).
#'         }
#'         \item{\code{strategy_results}}{
#'             A long-format \code{S4Vectors::DataFrame} containing the
#'             strategy-specific effect sizes (e.g., estimates for Ribo-seq
#'             only) for all computed contrasts. Columns include 
#'             \code{strategy} (e.g., "ribo", "rna"), \code{contrast}, and 
#'             the full set of shrunken and unshrunken metrics.
#'         }
#'     }
#'
#' @details
#' The results for post hoc contrasts are stored in long format via explicit
#' \code{contrast} and/or \code{strategy} columns. Non-converged models are
#' omitted.
#'
#' @importFrom emmeans emmeans contrast
#' @importFrom ashr ash get_pm get_qvalue get_lfdr get_lfsr
#' @importFrom methods is
#' @importFrom S4Vectors Rle DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom pbapply pblapply
#' @importFrom stats AIC aggregate anova as.dendrogram as.formula 
#' @importFrom stats p.adjust pnorm complete.cases
#'
#' @export
#' @examples
#' # Load SummarizedExperiment to enable subsetting and access to 
#' # components like rowRanges and rowData
#' library(SummarizedExperiment)
#'
#' # Read in count matrix, condition table, and annotation files
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
#' # extract only samples processed using cyclohexamide
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL # remove the treatment column
#'
#' # Create SummarizedExperiment objects.These objects can be used as input
#' # for DOTSeq and fitDOU
#' m <- DOTSeqDataSet(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     bed = bed
#' )
#'
#' # Keep ORFs where all replicates in at least one condition pass min_count
#' # Single-ORF genes are removed
#' m$sumExp <- m$sumExp[rowRanges(m$sumExp)$is_kept == TRUE, ]
#'
#' # Randomly sample 100 ORFs for fitDOU
#' n <- 100
#' set.seed(42)
#' random_rows <- sample(seq_len(nrow(m$sumExp)), size = n)
#'
#' # Subset the SummarizedExperiment object
#' m$sumExp <- m$sumExp[random_rows, ]
#'
#' # Model fitting using fitDOU
#' rowData(m$sumExp)[["DOUResults"]] <- fitDOU(
#'     count_table = assay(m$sumExp),
#'     rowdata = rowData(m$sumExp),
#'     anno = colData(m$sumExp),
#'     formula = ~ condition * strategy,
#'     emm_specs = ~ condition * strategy,
#'     dispersion_modeling = "auto",
#'     lrt = FALSE,
#'     optimizers = FALSE,
#'     diagnostic = FALSE,
#'     parallel = list(n = 4L, autopar = TRUE),
#'     verbose = TRUE
#' )
#'
#' # Run post hoc contrasts, Wald tests, and effect size shrinkage
#' m$sumExp <- testDOU(m$sumExp, verbose = TRUE)
#'
#' @references
#' Lenth R, Piaskowski J (2025). emmeans: Estimated Marginal Means, aka
#' Least-Squares Means. R package version 2.0.0,
#' \url{https://rvlenth.github.io/emmeans/}
#'
#' Stephens, M. (2016) False discovery rates: a new deal.
#' Biostatistics, 18:2. \doi{10.1093/biostatistics/kxw041}
#'
#'
testDOU <- function(
    sumExp,
    contrasts_method = "revpairwise",
    nullweight = 500,
    verbose = TRUE
) {
    if (verbose) {
        message("starting post hoc analysis")
    }

    result_list <- rowData(sumExp)[["DOUResults"]]

    valid_indices <- which(vapply(result_list, function(obj) {
        !is.null(posthoc(obj)) && inherits(posthoc(obj), "emmGrid")
    }, logical(1)))
    
    emm_list <- lapply(valid_indices, function(i) posthoc(result_list[[i]]))

    all_results <- list()
    all_results[["interaction_specific"]] <- list()

    ### Manual Interaction contrasts
    first_emm <- emm_list[[1]]
    contrast_df_template <- as.data.frame(
        summary(
            contrast(
                first_emm, 
                method = contrasts_method, 
                by = "strategy", 
                enhance.levels = FALSE, 
                adjust = "none"
            )
        )
    )
    all_contrast_names <- unique(as.character(contrast_df_template$contrast))

    for (c_name in all_contrast_names) {
        # if (verbose) {
        #     message("calculating the effect size and se for ", c_name)
        # }

        # Use lapply to iterate through all genes and extract the betas and SEs
        all_contrasts <- pblapply(emm_list, function(emm) {
            contrast_df <- summary(
                contrast(
                    emm,
                    method = contrasts_method, 
                    by = "strategy", 
                    enhance.levels = FALSE, 
                    adjust = "none"
                )
            )

            # Set baseline
            strategy_vals <- as.character(contrast_df$strategy)

            # Define a flexible regex pattern: matches "rna", "RNA", "RNA-seq", and "0"
            rna_pattern <- "(?i)rna|^0$"

            # Find matching values to RNA-seq
            rna_like <- unique(strategy_vals[grepl(rna_pattern, strategy_vals)])

            # Handle multiple matches
            if (length(rna_like) > 1) {
                warning("Multiple RNA-seq levels found. Using the first match as reference and the last as .")
                rna_like <- rna_like[1]
            } else if (length(rna_like) == 0) {
                stop("No RNA-seq level found. Please use rna, RNA, RNA-seq, or 0 to represent RNA-seq.")
            }

            # Find the Ribo-seq level
            ribo_like <- setdiff(unique(strategy_vals), rna_like)

            if (length(ribo_like) > 1) {
                warning("Multiple Ribo-seq levels found. Use ", ribo_like[1], " as target.")
                ribo_like <- ribo_like[1]
            } else if (length(ribo_like) == 0) {
                stop("No Ribo-seq level found. Please check your data.")
            }

            # Use as.character() for robust subsetting
            beta_ribo <- contrast_df$estimate[(as.character(contrast_df$strategy) == ribo_like) & (as.character(contrast_df$contrast) == c_name)]
            se_ribo <- contrast_df$SE[(as.character(contrast_df$strategy) == ribo_like) & (as.character(contrast_df$contrast) == c_name)]

            beta_rna <- contrast_df$estimate[(as.character(contrast_df$strategy) == rna_like) & (as.character(contrast_df$contrast) == c_name)]
            se_rna <- contrast_df$SE[(as.character(contrast_df$strategy) == rna_like) & (as.character(contrast_df$contrast) == c_name)]

            beta_dou <- beta_ribo - beta_rna
            se_dou <- sqrt(se_ribo^2 + se_rna^2)

            return(list(beta = beta_dou, se = se_dou))
        })

        # Extract the betas and SEs for the current contrast name
        betas_for_ashr <- vapply(all_contrasts, `[[`, numeric(1), "beta")
        ses_for_ashr <- vapply(all_contrasts, `[[`, numeric(1), "se")
        pvalues <- 2 * (1 - pnorm(abs(betas_for_ashr / ses_for_ashr)))

        # Run ashr on the full set of data for this contrast
        if (any(!is.na(betas_for_ashr))) {
            # if (verbose) {
            #     message("performing empirical Bayesian shrinkage on the effect size for ", c_name)
            # }

            ash_result <- ash(
                betas_for_ashr, 
                ses_for_ashr, 
                nullweight = nullweight, 
                pointmass = TRUE
            )

            # Collect results in a data frame
            res_df <- data.frame(
                betahat = betas_for_ashr,
                sebetahat = ses_for_ashr,
                WaldP = pvalues,
                WaldPadj = p.adjust(pvalues, method = "BH"),
                PosteriorMean = get_pm(ash_result),
                qvalue = get_qvalue(ash_result),
                lfdr = get_lfdr(ash_result),
                lfsr = get_lfsr(ash_result)
            )

            all_results[["interaction_specific"]][[c_name]] <- res_df
        } else {
            warning("No valid betas for manual contrast:", c_name)
            all_results[["interaction_specific"]][[c_name]] <- NULL
        }
    }

    ### Strategy-specific contrasts
    all_results[["strategy_specific"]] <- list()
    strategy_combos <- unique(summary(contrast_df_template)[, c("contrast", "strategy")])

    for (i in seq_len(nrow(strategy_combos))) {
        c_name <- as.character(strategy_combos$contrast[i])
        by_name <- as.character(strategy_combos$strategy[i])

        # if (verbose) {
        #     message("calculating the effect size and se for contrast: ", c_name, ", strategy: ", by_name)
        # }

        # Check if the list for the current contrast name exists
        if (is.null(all_results[["strategy_specific"]][[c_name]])) {
            all_results[["strategy_specific"]][[c_name]] <- list()
        }

        betas <- vapply(emm_list, function(emm) {
            res <- summary(
                contrast(
                    emm, 
                    method = contrasts_method, 
                    by = "strategy", 
                    enhance.levels = FALSE, 
                    adjust = "none"
                )
            )
            res$estimate[res$contrast == c_name & res$strategy == by_name]
        }, FUN.VALUE = numeric(1))

        ses <- vapply(emm_list, function(emm) {
            res <- summary(
                contrast(
                    emm, 
                    method = contrasts_method, 
                    by = "strategy", 
                    enhance.levels = FALSE, 
                    adjust = "none"
                )
            )
            res$SE[res$contrast == c_name & res$strategy == by_name]
        }, FUN.VALUE = numeric(1))

        pvalues <- 2 * (1 - pnorm(abs(betas / ses)))

        # if (verbose) {
        #     message("performing empirical Bayesian shrinkage on the effect size for ", c_name)
        # }

        if (any(!is.na(betas))) {
            ash_result <- ash(betas, ses)
            res_df <- data.frame(
                betahat = betas,
                sebetahat = ses,
                WaldP = pvalues,
                WaldPadj = p.adjust(pvalues, method = "BH"),
                PosteriorMean = get_pm(ash_result),
                qvalue = get_qvalue(ash_result),
                lfdr = get_lfdr(ash_result),
                lfsr = get_lfsr(ash_result)
            )
        } else {
            res_df <- data.frame(
                beta = NA, 
                se = NA, 
                PosteriorMean = NA, 
                lfsr = NA, 
                qvalue = NA
            )
        }

        all_results[["strategy_specific"]][[c_name]][[by_name]] <- res_df
    }

    # Flatten interaction_specific
    interaction_list <- lapply(names(all_results$interaction_specific), function(c_name) {
        df <- all_results$interaction_specific[[c_name]]
        df$orf_id <- rownames(df)
        rownames(df) <- NULL
        df$contrast <- c_name
        return(df)
    })
    interaction_df <- do.call(rbind, interaction_list)
    interaction_df <- DataFrame(interaction_df)

    # Apply Rle compression
    interaction_df$contrast <- Rle(interaction_df$contrast)

    # Flatten strategy_specific
    strategy_list <- lapply(names(all_results$strategy_specific), function(c_name) {
        lapply(names(all_results$strategy_specific[[c_name]]), function(by_name) {
            df <- all_results$strategy_specific[[c_name]][[by_name]]
            df$orf_id <- rownames(df)
            rownames(df) <- NULL
            df$contrast <- c_name
            df$strategy <- by_name
            return(df)
        })
    })
    strategy_df <- do.call(rbind, unlist(strategy_list, recursive = FALSE))
    strategy_df <- DataFrame(strategy_df)

    strategy_df$contrast <- Rle(strategy_df$contrast)
    strategy_df$strategy <- Rle(strategy_df$strategy)

    metadata(sumExp)$interaction_results <- interaction_df
    metadata(sumExp)$strategy_results <- strategy_df
    
    if (verbose) {
        message("effect sizes calculated and shrunk successfully")
    }

    return(sumExp)
}


#' Generate Contrast Vectors for Pairwise and Baseline Comparisons
#'
#' @description
#' Constructs a named list of contrast vectors for differential analysis using
#' a DESeq2 \code{DESeqDataSet} object. It identifies interaction terms (e.g.,
#' \code{condition:strategy}) from the model design and builds pairwise
#' contrasts between them, as well as contrasts against a specified baseline
#' condition.
#'
#' This is useful for extracting custom contrasts not directly available via
#' \code{resultsNames(dds)} and for post hoc comparisons in complex designs
#' involving interaction terms.
#'
#' @param dds A \code{DESeqDataSet} object containing count data and sample
#'     annotations.
#'
#' @param formula Optional. A model formula used to generate the design matrix
#'     (e.g., \code{~ condition * strategy}). If not provided, the function
#'     will attempt to extract it from \code{metadata(dds)$formula}.
#'
#' @param baseline Optional. A character string specifying the baseline
#'     condition for comparisons. If \code{NULL}, the first level of
#'     \code{condition} in \code{colData(dds)} is used.
#'
#' @param delim A character string used to identify interaction terms in
#'     coefficient names. Default is \code{"."}.
#'
#' @return A named list of numeric contrast vectors, where each vector
#'     corresponds to a pairwise or baseline comparison. These vectors can be
#'     used with \code{results(dds, contrast = ...)} for custom differential
#'     testing.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom stats model.matrix
#' @importFrom utils combn
#' @importFrom DESeq2 resultsNames
#'
#' @keywords internal
#' @examples
#' \dontrun{
#'     # Generate contrast vectors from a DESeqDataSet
#'     contrast_list <- contrast_vectors(
#'         dds,
#'         formula = ~ condition * strategy
#'     )
#'
#'     # Use a contrast vector with DESeq2
#'     res <- lfcShrink(
#'         dot$dds,
#'         contrast = contrast_list[[1]],
#'         type = "ashr"
#'     )
#' }
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550.
#' \doi{10.1186/s13059-014-0550-8}
#'
contrast_vectors <- function(
    dds, 
    formula = NULL, 
    baseline = NULL, 
    delim = "."
) {
    
    # Build design matrix
    anno <- colData(dds)
    if (is.null(baseline)) {
        baseline <- levels(anno$condition)[1]
    }

    strategy_levels <- levels(anno$strategy)
    if (length(strategy_levels) > 2) {
        stop("Expect two strategy levels got: ", paste(strategy_levels, collapse = ", "))
    } else if (length(strategy_levels) == 1) {
        stop("Expect two strategy levels got: ", strategy_levels)
    } else {
        ribo_level <- strategy_levels[2]
    }

    # Get the names of all coefficients in the model matrix
    all_terms <- resultsNames(dds)
    num_terms <- length(all_terms)
    contrast_factors <- grep(paste0("\\", delim), all_terms, value = TRUE)

    if (!is.null(formula)) {
        design <- model.matrix(formula, data = anno)
    } else if ("formula" %in% names(metadata(dds))) {
        design <- model.matrix(metadata(dds)$formula, data = anno)
    } else {
        stop("Please provide design formula.")
    }

    attrs <- names(attr(design, "contrasts"))

    contrast_vectors_list <- list()

    if (length(contrast_factors) > 1) {
        comparison_pairs <- combn(contrast_factors, 2, simplify = FALSE)

        for (pairs in comparison_pairs) {
            contrast_vector <- numeric(num_terms)
            index1 <- which(all_terms == pairs[1])
            index2 <- which(all_terms == pairs[2])
            contrast_vector[index1] <- -1
            contrast_vector[index2] <- 1

            target1 <- sub(attrs[1], "", pairs[1])
            target1 <- sub(paste0(delim, attrs[2], ribo_level), "", target1)
            target2 <- sub(attrs[1], "", pairs[2])
            target2 <- sub(paste0(delim, attrs[2], ribo_level), "", target2)

            contrast_name <- paste0(target2, " - ", target1)
            contrast_vectors_list[[contrast_name]] <- contrast_vector
        }

        for (term in all_terms) {
            if (any(grep(paste0("\\", delim), term))) {
                contrast_vector <- numeric(num_terms)
                idx <- which(all_terms == term)
                contrast_vector[idx] <- 1

                target <- sub(attrs[1], "", term)
                target <- sub(paste0(delim, attrs[2], ribo_level), "", target)
                contrast_name <- paste0(target, " - ", baseline)

                contrast_vectors_list[[contrast_name]] <- contrast_vector
            }
        }
    } else {
        for (term in all_terms) {
            if (any(grep(paste0("\\", delim), term))) {
                contrast_vector <- numeric(num_terms)
                idx <- which(all_terms == term)
                contrast_vector[idx] <- 1

                target <- sub(attrs[1], "", term)
                target <- sub(paste0(delim, attrs[2], ribo_level), "", target)
                contrast_name <- paste0(target, " - ", baseline)

                contrast_vectors_list[[contrast_name]] <- contrast_vector
            }
        }
    }

    return(contrast_vectors_list)
}
