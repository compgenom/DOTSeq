#' Test for Differential ORF Usage (DOU) Using Estimated Marginal Means
#'
#' @description
#' Performs differential ORF usage (DOU) analysis using estimated marginal means (EMMs)
#' from a fitted model object. Computes contrast-specific effect sizes and applies
#' empirical Bayes shrinkage using the `ashr` package. Supports both interaction-specific
#' and strategy-specific contrasts.
#'
#' @param m A named list of fitted model objects (e.g., from `glmmTMB`), one per ORF.
#' @param emm_specs A formula specifying the EMM structure (default: \code{~condition * strategy}).
#' @param contrasts_method Character. Method for computing contrasts (default: \code{"pairwise"}).
#' @param workers Integer. Number of parallel workers to use (default: \code{1}).
#' @param BPPARAM A BiocParallel parameter object for parallel execution (default: \code{BiocParallel::bpparam()}).
#' @param verbose Logical. If \code{TRUE}, print progress messages (default: \code{FALSE}).
#'
#' @return A named list with two components:
#' \describe{
#'   \item{interaction_specific}{A list of contrast-specific data frames containing effect sizes, standard errors, shrunk estimates, local false sign rates (lfsr), and q-values.}
#'   \item{strategy_specific}{A nested list of contrast-strategy combinations with the same structure as above.}
#' }
#'
#' @details
#' The function first filters out invalid or NULL EMM objects. For each valid contrast,
#' it computes the difference in effect sizes between ribosome profiling and RNA-seq
#' (i.e., DOU), estimates the standard error, and applies shrinkage using `ashr::ash`.
#'
#' Strategy-specific contrasts are also computed separately for each modality.
#'
#' @import emmeans
#' @import ashr
#' @importFrom BiocParallel register
#' @importFrom BiocParallel bplapply
#' @importFrom pbapply pblapply
#' @export
testDOU <- function(m, emm_specs = ~condition * strategy, 
                    contrasts_method = "pairwise", 
                    workers = 1,
                    BPPARAM = BiocParallel::bpparam(), 
                    verbose = FALSE) {
  
  fitted_models_list <- rowData(m$sumExp)[['fitDOUModels']]
  
  if (isTRUE(verbose)) {
    message(" - Perform contrasts")
    if (workers > 1) {
      BPPARAM <- BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
      BiocParallel::register(BPPARAM)
      emm_list <- BiocParallel::bplapply(fitted_models_list, function(m) {
        if (is(m, "StatModel") && m@type == "glmmTMB") {
          tryCatch(emmeans(m@model, specs = emm_specs), error = function(e) NULL)
        } else {
          NULL
        }
      }, BPPARAM = BPPARAM)
    } else {
      emm_list <- pbapply::pblapply(fitted_models_list, function(m) {
        if (is(m, "StatModel") && m@type == "glmmTMB") {
          tryCatch(emmeans(m@model, specs = emm_specs), error = function(e) NULL)
        } else {
          NULL
        }
      })
    }
  } else {
    # No progress bar, default backend
    emm_list <- lapply(fitted_models_list, function(m) {
      if (is(m, "StatModel") && m@type == "glmmTMB") {
        tryCatch(emmeans(m@model, specs = emm_specs), error = function(e) NULL)
      } else {
        NULL
      }
    })
  }
  
  valid_genes <- !sapply(emm_list, is.null)
  emm_list <- emm_list[valid_genes]
  gene_ids <- names(emm_list)
  
  if (length(emm_list) == 0) {
    warning("No valid emmeans objects could be generated.")
    return(NULL)
  }
  
  all_results <- list()
  all_results[["interaction_specific"]] <- list()
  
  ### Manual Interaction contrasts
  first_emm <- emm_list[[1]]
  contrast_df_template <- as.data.frame(summary(contrast(first_emm, method = contrasts_method, by = "strategy", adjust = "none")))
  all_contrast_names <- unique(as.character(contrast_df_template$contrast))
  
  for (c_name in all_contrast_names) {
    if (verbose) {
      message(" - Calculate the effect size and standard error for ", c_name)
    }
    
    # Use lapply to iterate through all genes and extract the betas and SEs
    all_contrasts <- lapply(emm_list, function(emm) {
      contrast_df <- summary(contrast(emm, method = contrasts_method, by = "strategy", adjust = "none"))
      
      # Use as.character() for robust subsetting
      beta_ribo <- contrast_df$estimate[(as.character(contrast_df$strategy) == "1") & (as.character(contrast_df$contrast) == c_name)]
      se_ribo <- contrast_df$SE[(as.character(contrast_df$strategy) == "1") & (as.character(contrast_df$contrast) == c_name)]
      
      beta_rna <- contrast_df$estimate[(as.character(contrast_df$strategy) == "0") & (as.character(contrast_df$contrast) == c_name)]
      se_rna <- contrast_df$SE[(as.character(contrast_df$strategy) == "0") & (as.character(contrast_df$contrast) == c_name)]
      
      beta_dou <- beta_ribo - beta_rna
      se_dou <- sqrt(se_ribo^2 + se_rna^2)
      
      return(list(beta = beta_dou, se = se_dou))
    })
    
    # Extract the betas and SEs for the current contrast name
    betas_for_ashr <- sapply(all_contrasts, `[[`, "beta")
    ses_for_ashr <- sapply(all_contrasts, `[[`, "se")
    
    # Run ashr on the full set of data for this contrast
    if (any(!is.na(betas_for_ashr))) {
      if (verbose) {
        message(" - Perform empirical Bayesian shrinkage on the effect size for ", c_name)
      }
      
      ash_result <- ashr::ash(betas_for_ashr, ses_for_ashr)
      
      # Collect results in a data frame
      res_df <- data.frame(
        beta = betas_for_ashr,
        se = ses_for_ashr,
        shrunkBeta = ashr::get_pm(ash_result),
        lfsr = ashr::get_lfsr(ash_result),
        qvalue = ashr::get_qvalue(ash_result)
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
      message(" - Calculate the effect size and standard error for contrast: ", c_name, ", strategy: ", by_name)
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
    
    if (verbose) {
      message(" - Perform empirical Bayesian shrinkage on the effect size for ", c_name)
    }
    
    if (any(!is.na(betas))) {
      ash_result <- ashr::ash(betas, ses)
      res_df <- data.frame(
        beta = betas,
        se = ses,
        shrunkBeta = ashr::get_pm(ash_result),
        lfsr = ashr::get_lfsr(ash_result),
        qvalue = ashr::get_qvalue(ash_result)
      )
    } else {
      res_df <- data.frame(beta = NA, se = NA, shrunkBeta = NA, lfsr = NA, qvalue = NA)
    }
    
    all_results[["strategy_specific"]][[c_name]][[by_name]] <- res_df
  }
  return(all_results)
}


#' Extract Model Results into a Unified Data Frame
#'
#' @description
#' Extracts scalar results from a named list of model objects (typically from `glmmTMB` fits)
#' into a tidy data frame. Handles nested lists, missing or `NULL` models, and ensures
#' consistent columns across all entries.
#'
#' @param models_list A named list of model objects, each representing an ORF.
#'   Each object should contain a `@type` slot and optionally a `@results` slot.
#' @param verbose Logical. If \code{TRUE}, messages will be printed for skipped or failed models.
#'
#' @return A data frame where each row corresponds to an ORF and its associated model results.
#'   Columns include scalar parameters extracted from the model, \code{ORF_ID}, and \code{modelType}.
#'
#' @details
#' The function uses a recursive helper (`flatten_scalars`) to extract scalar values
#' from nested lists. It supports models of type \code{"glmmTMB"} and \code{"glmmTMB_joint"},
#' and gracefully handles \code{NULL} or unsupported model types by returning minimal rows.
#'
#' Missing columns across models are filled with \code{NA} to ensure a consistent structure.
#'
#' @import SummarizedExperiment
#' @export
extract_results <- function(models_list, verbose = TRUE) {
  
  # Helper to recursively flatten and extract scalar values
  flatten_scalars <- function(x, prefix = NULL) {
    if (is.atomic(x) && length(x) == 1) {
      name <- if (is.null(prefix)) "value" else prefix
      return(setNames(list(x), name))
    } else if (is.list(x)) {
      result <- list()
      # Use `names(x)` and `seq_along` to handle unnamed elements gracefully
      list_names <- names(x)
      if (is.null(list_names)) {
        list_names <- paste0("elem", seq_along(x))
      }
      for (i in seq_along(x)) {
        n <- list_names[i]
        sub_prefix <- if (is.null(prefix)) n else paste0(prefix, ".", n)
        result <- c(result, flatten_scalars(x[[i]], sub_prefix))
      }
      return(result)
    } else {
      return(list()) # Return an empty list for non-atomic, non-list objects
    }
  }
  
  # Process a single model
  process_model <- function(model_obj, name) {
    # --- ADDED: Check for NULL objects here ---
    if (is.null(model_obj)) {
      if (verbose) message("Skipping ORF_ID ", name, " due to NULL model object.")
      return(data.frame(
        ORF_ID = name,
        modelType = "NULL", # Or "Failed", etc.
        stringsAsFactors = FALSE
      ))
    }
    
    model_type <- model_obj@type
    
    if ((model_type == "glmmTMB_joint") | (model_type == "glmmTMB")) {
      params <- model_obj@results
      param_values <- flatten_scalars(params)
      
      param_values$ORF_ID <- name
      param_values$modelType <- model_type
      
      # --- CORRECTED: Create data.frame from the single list to guarantee one row ---
      df <- as.data.frame(param_values, stringsAsFactors = FALSE)
      
    } else {
      # Return minimal row with NA for failed or single-ORF models
      df <- data.frame(
        ORF_ID = name,
        modelType = model_type,
        stringsAsFactors = FALSE
      )
    }
    
    return(df)
  }
  
  # Apply to all models
  results_df_list <- lapply(names(models_list), function(name) {
    process_model(models_list[[name]], name)
  })
  
  # Get all unique column names
  all_cols <- unique(unlist(lapply(results_df_list, names)))
  
  # Fill missing columns with NA
  results_df_list_filled <- lapply(results_df_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]
  })
  
  # Combine all rows
  results_df <- do.call(rbind, results_df_list_filled)
  rownames(results_df) <- NULL
  
  return(results_df)
}
