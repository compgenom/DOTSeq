#' Compute Differential ORF Usage (DOU) Contrasts Using Estimated Marginal Means
#'
#' @description
#' Performs Differential ORF Usage (DOU) analysis by computing contrasts between ribosome profiling 
#' and RNA-seq modalities using estimated marginal means (EMMs) from fitted GLMM models. 
#' Supports both interaction-specific and strategy-specific contrasts, and applies empirical Bayes 
#' shrinkage via the `ashr` package to stabilize effect size estimates.
#'
#' @param sumExp A SummarizedExperiment object containing fitted model objects,
#'   typically stored in `rowData(sumExp)[['fitDOUModels']]`.
#' @param emm_specs A formula specifying the structure of the estimated marginal means.
#'   Default is \code{~condition * strategy}.
#' @param contrasts_method Character string specifying the method for computing contrasts.
#'   Default is \code{"pairwise"}.
#' @param workers Integer. Number of parallel workers to use. Default is \code{1}.
#' @param nullweight Numeric. Prior weight on the null hypothesis for empirical Bayes shrinkage.
#'   Higher values yield more conservative lfsr estimates. Default is \code{500}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Default is \code{TRUE}.
#'
#' @return A **SummarizedExperiment** object (the input \code{sumExp} or \code{m$sumExp}) 
#'   with two new \code{S4Vectors::DataFrame} objects stored in its \code{metadata} slot.
#'   These tables contain the long-format results for all computed contrasts:
#' \describe{
#'    \item{\code{interaction_results}}{A long-format \code{S4Vectors::DataFrame} containing the 
#'      **Differential ORF Usage (DOU) effect sizes** (Ribo-seq contrast minus RNA-seq contrast) 
#'      for all computed interaction contrasts. Columns include \code{contrast}, 
#'      the full set of shrunken and unshrunken metrics (e.g., \code{betahat}, \code{sebetahat}, 
#'      \code{WaldPadj}, \code{PosteriorMean}, \code{lfsr}).}
#'    \item{\code{strategy_results}}{A long-format \code{S4Vectors::DataFrame} containing the 
#'      **strategy-specific effect sizes** (e.g., estimates for Ribo-seq only) for all computed contrasts.
#'      Columns include \code{strategy} (e.g., "ribo", "rna"), \code{contrast},  
#'      and the full set of shrunken and unshrunken metrics.}
#' }
#' 
#' @details
#' The results for post hoc contrasts are stored in a long format via explicit \code{contrast} and/or 
#' \code{strategy} columns. Non-converged models are omitted.
#'
#' @importFrom emmeans emmeans contrast
#' @importFrom ashr ash get_pm get_qvalue get_lfdr get_lfsr
#' @importFrom methods is
#' @importFrom S4Vectors Rle DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom pbapply pblapply
#' @importFrom stats AIC aggregate anova as.dendrogram as.formula complete.cases
#' @importFrom stats p.adjust pnorm
#' 
#' @export
#' 
testDOU <- function(
    sumExp, 
    emm_specs = ~condition * strategy, 
    contrasts_method = "pairwise", 
    nullweight = 500,
    BPPARAM = bpparam(), 
    verbose = TRUE
    ) {
  
  if (verbose) {
    message(" - Starting post hoc analysis")
  }
  
  result_list <- rowData(sumExp)[['DOUResults']]
  
  valid_indices <- which(vapply(result_list, function(obj) {
    !is.null(obj@posthoc) && inherits(obj@posthoc, "emmGrid")
  }, logical(1)))
  
  emm_list <- lapply(valid_indices, function(i) result_list[[i]]@posthoc)
  
  all_results <- list()
  all_results[["interaction_specific"]] <- list()
  
  ### Manual Interaction contrasts
  first_emm <- emm_list[[1]]
  contrast_df_template <- as.data.frame(summary(contrast(first_emm, method = contrasts_method, by = "strategy", adjust = "none")))
  all_contrast_names <- unique(as.character(contrast_df_template$contrast))
  
  for (c_name in all_contrast_names) {
    if (verbose) {
      message(" - Calculating the effect size and standard error for ", c_name)
    }
    
    # Use lapply to iterate through all genes and extract the betas and SEs
    all_contrasts <- pblapply(emm_list, function(emm) {
      contrast_df <- summary(contrast(emm, method = contrasts_method, by = "strategy", adjust = "none"))
      
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
    betas_for_ashr <- sapply(all_contrasts, `[[`, "beta")
    ses_for_ashr <- sapply(all_contrasts, `[[`, "se")
    pvalues <- 2 * (1 - pnorm(abs(betas_for_ashr / ses_for_ashr)))
    
    # Run ashr on the full set of data for this contrast
    if (any(!is.na(betas_for_ashr))) {
      if (verbose) {
        message(" - Performing empirical Bayesian shrinkage on the effect size for ", c_name)
      }
      
      ash_result <- ash(betas_for_ashr, ses_for_ashr, nullweight = nullweight, pointmass = TRUE)
      
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
      warning(paste("No valid betas for manual contrast:", c_name))
      all_results[["interaction_specific"]][[c_name]] <- NULL
    }
  }
  
  ### Strategy-specific contrasts
  all_results[["strategy_specific"]] <- list()
  strategy_combos <- unique(summary(contrast_df_template)[, c("contrast", "strategy")])
  
  for (i in seq_len(nrow(strategy_combos))) {
    c_name <- as.character(strategy_combos$contrast[i])
    by_name <- as.character(strategy_combos$strategy[i])
    
    if (verbose) {
      message(" - Calculating the effect size and standard error for contrast: ", c_name, ", strategy: ", by_name)
    }
    
    # Check if the list for the current contrast name exists
    if (is.null(all_results[["strategy_specific"]][[c_name]])) {
      all_results[["strategy_specific"]][[c_name]] <- list()
    }
    
    betas <- vapply(emm_list, function(emm) {
      res <- summary(contrast(emm, method = contrasts_method, by = "strategy", adjust = "none"))
      res$estimate[res$contrast == c_name & res$strategy == by_name]
    }, FUN.VALUE = numeric(1))
    
    ses <- vapply(emm_list, function(emm) {
      res <- summary(contrast(emm, method = contrasts_method, by = "strategy", adjust = "none"))
      res$SE[res$contrast == c_name & res$strategy == by_name]
    }, FUN.VALUE = numeric(1))
    
    pvalues <- 2 * (1 - pnorm(abs(betas / ses)))
    
    if (verbose) {
      message(" - Performing empirical Bayesian shrinkage on the effect size for ", c_name)
    }
    
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
      res_df <- data.frame(beta = NA, se = NA, PosteriorMean = NA, lfsr = NA, qvalue = NA)
    }
    
    all_results[["strategy_specific"]][[c_name]][[by_name]] <- res_df
  }
  
  # Flatten interaction_specific
  interaction_list <- lapply(names(all_results$interaction_specific), function(c_name) {
    df <- all_results$interaction_specific[[c_name]]
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
  
  return(sumExp)
}



