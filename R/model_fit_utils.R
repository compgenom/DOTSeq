#' @title Calculate Moment-Based Dispersion Estimate
#'
#' @description
#' This helper function computes a moment-based estimate of dispersion (intra-class correlation)
#' for a given set of ORF counts and total gene counts across samples. It avoids division by zero
#' and clamps negative dispersion estimates to zero.
#'
#' @param counts A numeric vector of ORF counts across samples.
#' @param totals A numeric vector of total gene counts across the same samples.
#'
#' @return A \code{list} containing:
#' \describe{
#'   \item{mean_prop}{Mean proportion of ORF counts relative to total gene counts.}
#'   \item{rawRho}{Unclamped dispersion estimate (can be negative).}
#'   \item{raw_intra_class_corr}{Clamped dispersion estimate (non-negative).}
#' }
#' 
#' @importFrom stats var
#' 
#' @keywords internal
#' 
.calculate_mean_dispersion <- function(counts, totals) {
  # Avoid division by zero
  if (length(counts) < 2) {
    return(list(mean_prop = NA_real_, raw_rho = NA_real_, raw_intra_class_corr = NA_real_))
  }
  
  # The proportions
  proportions <- counts / totals
  
  # The mean proportion for the group
  mu_hat <- mean(proportions)
  
  # The sample variance of the proportions
  var_hat <- var(proportions)
  
  # The mean of the total counts for the group
  N_bar <- mean(totals)
  
  # Moment-based estimate for rho
  # Based on: var_hat = (mu_hat * (1 - mu_hat) / N_bar) * (1 + (N_bar - 1) * rho)
  # Solve for rho_hat
  rho_hat <- (var_hat * N_bar) / (mu_hat * (1 - mu_hat)) - 1
  rho_hat <- rho_hat / (N_bar - 1)
  
  # Ensure rho is non-negative
  return(list(mean_prop = mu_hat, 
              rawRho = rho_hat,
              raw_intra_class_corr = max(0, rho_hat)))
}



#' @title Fit a Beta-Binomial GLMM with Multiple Optimizers
#'
#' @description
#' This internal function fits a beta-binomial generalized linear mixed model (GLMM) using \code{glmmTMB}, 
#' with support for multiple optimizers to improve convergence. It is used as a backend for fitting 
#' models to ORF-level count data in DOTSeq.
#'
#' @param formula A formula object specifying the fixed effects structure of the model.
#' @param dispformula A formula object specifying the dispersion model (e.g., \code{~strategy}).
#' @param data A data frame containing the model data, including counts and covariates.
#' @param optimizers A character vector of optimizer names to try sequentially. Supported values include 
#'   \code{"nlminb"}, \code{"bobyqa"}, and \code{"optim"} (default: \code{c("nlminb", "bobyqa", "optim")}).
#' @param max_iter An integer specifying the maximum number of iterations for the optimizer (default: \code{1000}).
#'
#' @return A fitted \code{glmmTMB} model object if convergence is successful. If all optimizers fail, 
#'   the last attempted model (with convergence failure) is returned.
#'
#' @details
#' The function attempts to fit the model using each optimizer in the specified order. If a model converges 
#' successfully (i.e., \code{model$fit$convergence == 0} and \code{model$sdr$pdHess == TRUE}), it is returned immediately. 
#' Otherwise, the function continues trying the next optimizer.
#'
#' @importFrom glmmTMB glmmTMB glmmTMBControl betabinomial
#'
#' @keywords internal
#' 
fit_glmm <- function(formula, dispformula, data, family = betabinomial(), parallel = list(n=4L, autopar=TRUE),
                     optimizers = c("nlminb", "bobyqa", "optim"), max_iter = 1000) {
  for (opt in optimizers) {
    opt_args <- if (opt == "optim") list(method = "BFGS") else list()
    model <- try(glmmTMB(formula, dispformula = dispformula, family = family,
                         data = data,
                         control = glmmTMBControl(optimizer = opt,
                                                  optArgs = opt_args,
                                                  optCtrl = list(iter.max = max_iter, eval.max = max_iter),
                                                  parallel = parallel)),
                 silent = TRUE)
    
    if (inherits(model, "glmmTMB") && model$fit$convergence == 0 && isTRUE(model$sdr$pdHess)) {
      return(model)
    }
  }
  
  # If none of the optimizers worked
  return(model)
}


#' @title Fit Beta-Binomial Models for a Single Gene
#'
#' @description
#' This internal worker function fits beta-binomial generalized linear models (GLMs) or generalized linear mixed models (GLMMs)
#' for all ORFs within a single gene using \code{glmmTMB}. It supports multiple dispersion
#' modeling strategies and optional model diagnostics.
#'
#' @param ribo_mat A numeric matrix of Ribo-seq counts for a single gene (rows = ORFs, columns = samples).
#' @param rna_mat A numeric matrix of RNA-seq counts for the same gene (same dimensions as \code{ribo_mat}).
#' @param anno A data frame containing sample annotations. Must include columns such as \code{condition}, \code{strategy}, and \code{replicate}.
#' @param formula A formula object specifying the model design, e.g., \code{~ condition * strategy}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param dispersion_modeling Character string specifying the dispersion modeling strategy.
#'   Options are:
#'   \describe{
#'     \item{\code{"auto"}}{Fit both strategy-dependent and shared dispersion models, and select the best via likelihood ratio test.}
#'     \item{\code{"shared"}}{Assume constant dispersion across all predictor levels.}
#'     \item{\code{"custom"}}{Use a user-specified dispersion formula via \code{dispformula}.}
#'   }
#' @param dispformula Optional formula object specifying a custom dispersion model (used when \code{dispersion_modeling = "custom"}).
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, runs DHARMa diagnostics to assess model fit (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization,
#'   e.g., \code{list(parallel = TRUE, ncpus = 4)}. Default: \code{list(n = 4L, autopar = TRUE)}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{FALSE}).
#'
#' @return A named \code{list} of \code{StatModel} objects, one for each ORF in the gene.
#' 
#' @importFrom stats AIC aggregate anova as.formula
#' @importFrom stats order.dendrogram p.adjust pnorm predict rbinom relevel rgamma
#' 
#' @importFrom DHARMa simulateResiduals testDispersion testZeroInflation testResiduals
#' @importFrom glmmTMB glmmTMB glmmTMBControl fixef betabinomial
#' 
#' @keywords internal
#' 
.fitBetaBinomial <- function(ribo_mat, 
                             rna_mat, 
                             anno, 
                             formula = ~ condition * strategy, 
                             target = NULL,
                             dispersion_modeling = c("auto", "shared", "custom"), 
                             dispformula = NULL, 
                             lrt = FALSE,
                             diagnostic = FALSE, 
                             seed = NULL, 
                             optimizers = FALSE,
                             parallel = list(n=4L, autopar=TRUE), 
                             verbose = FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dispersion_modeling <- match.arg(dispersion_modeling)
  
  # Check if gene has only one ORF
  num_orfs <- nrow(ribo_mat)
  # Handle single-ORF genes
  if(nrow(ribo_mat) == 1 && nrow(rna_mat) == 1) {
    .out <- .StatModel(
      type = "singleORF",
      results = list(),
      model = NA
    )
    names(.out) <- rownames(ribo_mat)
    return(list(.out))
  }
  
  # Get a list of ORFs to process
  orf_names <- rownames(ribo_mat)
  num_orfs <- length(orf_names)
  num_samples <- ncol(ribo_mat)
  
  # Extract the strategy levels from anno
  strategy_levels <- levels(anno$strategy)
  
  # Define a flexible regex pattern: matches "rna", "RNA", "RNA-seq", and "0"
  rna_pattern <- "(?i)rna|^0$"
  
  # Find matching values
  if (isTRUE(grepl(rna_pattern, strategy_levels[1]))) {
    # Prepare the data for a single gene dynamically based on the design matrix
    long_data <- data.frame(
      counts = c(as.vector(ribo_mat), as.vector(rna_mat)), # Ribo first
      ORF = factor(rep(orf_names, times = num_samples * 2)),
      row.names = NULL
    )
  } else {
    long_data <- data.frame(
      counts = c(as.vector(rna_mat), as.vector(ribo_mat)),
      ORF = factor(rep(orf_names, times = num_samples * 2)),
      row.names = NULL
    )
  }
  
  # Replicate the design variables for each ORF
  long_design <- as.data.frame(anno[rep(seq_len(nrow(anno)), each = num_orfs), ])
  
  # Bind the replicated design matrix to the long data frame
  long_data <- cbind(long_data, long_design)

  # Check for negative counts
  if (isTRUE(any(long_data$counts < 0))) {
    print(long_data)
    stop("Encounter negative counts")
  }
  
  # Check for zero-count ORFs and remove them
  total_counts_by_orf <- aggregate(counts ~ ORF, data = long_data, sum)
  orfs_to_keep <- as.character(total_counts_by_orf$ORF[total_counts_by_orf$counts > 0])
  
  if (length(orfs_to_keep) < 2) {
    return(list(.StatModel(
      type = "fitError",
      results = list(),
      model = NA
    )))
  }
  
  # Fit the model for each non-zero ORF's proportion against the others
  models_gene <- lapply(orfs_to_keep, function(current_orf) {
    
    # Get the data for the current ORF
    model_data_this_orf <- long_data[long_data$ORF == current_orf, ]
    
    # Get the total gene counts for each sample
    if (isTRUE(grepl(rna_pattern, strategy_levels[1]))) {
      total_gene_counts <- c(colSums(ribo_mat), colSums(rna_mat))
    } else {
      total_gene_counts <- c(colSums(rna_mat), colSums(ribo_mat))
    }
    
    model_data_this_orf$success <- model_data_this_orf$counts
    model_data_this_orf$failure <- total_gene_counts - model_data_this_orf$success
    
    # # Construct random effects string
    # fixed_part <- as.character(formula)[2]
    # if (is.null(randomEffects)) {
    #   random_part <- ""
    # } else {
    #   random_part <- paste0(" + ", paste0(paste0("(1|", randomEffects, ")"), collapse = " + "))
    # }
    # 
    # full_formula_str <- paste0("cbind(success, failure) ~", fixed_part, random_part)
    # null_formula_str <- paste0("cbind(success, failure) ~", fixed_part, " - condition:strategy", random_part)
    # 
    # full_formula <- as.formula(full_formula_str)
    # null_formula <- as.formula(null_formula_str)
    
    full_formula <- as.formula(paste("cbind(success, failure) ~", as.character(formula)[2]))
    
    if (isTRUE(lrt)) {
      null_formula_str <- paste0(as.character(formula)[2], " - condition:strategy")
      null_formula <- as.formula(paste("cbind(success, failure) ~", null_formula_str)) 
    }
    
    # Use a single, robust tryCatch block to handle all potential errors
    result <- tryCatch({
      
      # Initialize all potential models and metrics to NULL
      model_strategy <- NULL
      model_shared <- NULL
      model_null <- NULL
      model_null_shared <- NULL
      
      # Initialize the results list with the new structure
      if (isTRUE(lrt)) {
        results <- list(model_fit = list(), tests = list(), estimates = list(), dispersion = list(), diagnostics = list())
      } else {
        results <- list(model_fit = list(), estimates = list(), dispersion = list(), diagnostics = list())
      }
      
      
      if (isTRUE(optimizers)) {
        # Fit Models based on dispersion_modeling
        if (dispersion_modeling %in% c("auto")) { #, "strategy"
          model_strategy <- fit_glmm(formula = full_formula, dispformula = ~strategy, data = model_data_this_orf)
          if ((is.null(model_strategy) | model_strategy$fit$convergence == 1 | isFALSE(model_strategy$sdr$pdHess))) {
            model_shared <- fit_glmm(full_formula, dispformula = ~1, family = betabinomial, data = model_data_this_orf)
            if (isTRUE(lrt)) {
              model_null_shared <- tryCatch({
                fit_glmm(formula = null_formula, dispformula = ~1, data = model_data_this_orf)
              }, error = function(e) {
                return(NA)
              })
            }
          }
          if (isTRUE(lrt)) {
            model_null <- tryCatch({
              fit_glmm(formula = null_formula, dispformula = ~ strategy, data = model_data_this_orf)
            }, error = function(e) {
              return(NA)
            })
          }
        }
        
        if (dispersion_modeling %in% c("shared")) { # "auto", 
          model_shared <- fit_glmm(formula = full_formula, dispformula = ~1, data = model_data_this_orf)
          if (isTRUE(lrt)) {
            model_null_shared <- tryCatch({
              fit_glmm(formula = null_formula, dispformula = ~1, data = model_data_this_orf)
            }, error = function(e) {
              return(NA)
            })
          }
        }
        
        if (dispersion_modeling == "custom") {
          if  (!is.null(dispformula)) {
            model_custom <- fit_glmm(formula = full_formula, dispformula = dispformula, data = model_data_this_orf)
            if (isTRUE(lrt)) {
              model_null <- tryCatch({
                fit_glmm(formula = null_formula, dispformula = dispformula, data = model_data_this_orf)
              }, error = function(e) {
                return(NA)
              })
            }
          } else {
            stop("Please provide dispformula.")
          }
        }
      } else {
        if (dispersion_modeling %in% c("auto")) { #, "strategy"
          model_strategy <- glmmTMB(full_formula, dispformula = ~strategy, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          if ((is.null(model_strategy) | model_strategy$fit$convergence == 1 | isFALSE(model_strategy$sdr$pdHess))) {
            model_shared <- glmmTMB(full_formula, dispformula = ~1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            if (isTRUE(lrt)) {
              model_null_shared <- tryCatch({
                glmmTMB(null_formula, dispformula = ~ 1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
              }, error = function(e) {
                return(NA)
              })
            }
          }
          if (isTRUE(lrt)) {
            model_null <- tryCatch({
              glmmTMB(null_formula, dispformula = ~ strategy, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            }, error = function(e) {
              return(NA)
            })
          }
        }
        
        if (dispersion_modeling %in% c("shared")) { #"auto", 
          model_shared <- glmmTMB(full_formula, dispformula = ~1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          if (isTRUE(lrt)) {
            model_null_shared <- tryCatch({
              glmmTMB(null_formula, dispformula = ~ 1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            }, error = function(e) {
              return(NA)
            })
          }
        }
        
        if (dispersion_modeling == "custom") {
          if  (!is.null(dispformula)) {
            model_custom <- glmmTMB(full_formula, dispformula = dispformula, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            if (isTRUE(lrt)) {
              model_null <- tryCatch({
                glmmTMB(null_formula, dispformula = dispformula, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
              }, error = function(e) {
                return(NA)
              })
            }
          } else {
            stop("Please provide dispformula.")
          }
        } 
      }
      # Build return object based on chosen strategy
      model_to_return <- NULL
      
      if (dispersion_modeling == "auto") {
        # Full comparison mode
        
        # Add strategy-specific metrics
        if (!is.null(model_strategy) && model_strategy$fit$convergence == 0 && isTRUE(model_strategy$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_strategy)
          if (isTRUE(lrt) && !is.null(model_null) && model_null$fit$convergence == 0 && isTRUE(model_null$sdr$pdHess)) {
            results$tests$pvalue <- anova(model_strategy, model_null)$P[2]
          }
          if (isTRUE(diagnostic)) {
            sim_out <- simulateResiduals(fittedModel = model_strategy, plot = FALSE)
            results$diagnostics$diagnostics_strategy <- list(overdispersion = testDispersion(sim_out, plot = FALSE)$p.value, 
                                                             zeroInflation = testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                             uniformity = testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                             residualsDispersion = testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                             outliers = testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
        
        # Add shared-specific metrics
        if (!is.null(model_shared) && model_shared$fit$convergence == 0 && isTRUE(model_shared$sdr$pdHess)) {
          results$model_fit$aic_shared_disp <- AIC(model_shared)
          
          # results$tests$lrt_shared_disp <- anova(model_strategy, model_shared)$P[2]
          if (isTRUE(lrt) && !is.null(model_null_shared) && model_null_shared$fit$convergence == 0 && isTRUE(model_null_shared$sdr$pdHess)) {
            results$tests$lrt_shared <- anova(model_shared, model_null_shared)$P[2]
          }
          
          if (isTRUE(diagnostic)) {
            sim_out <- simulateResiduals(fittedModel = model_shared, plot = FALSE)
            results$diagnostics$diagnostics_shared <- list(overdispersion = testDispersion(sim_out, plot = FALSE)$p.value, 
                                                           zeroInflation = testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                           uniformity = testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                           residualsDispersion = testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                           outliers = testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
        
        valid_strategy <- !is.null(model_strategy) && model_strategy$fit$convergence == 0 && isTRUE(model_strategy$sdr$pdHess)
        valid_shared <- !is.null(model_shared) && model_shared$fit$convergence == 0 && isTRUE(model_shared$sdr$pdHess)
        
        if (valid_strategy) { # isTRUE(lrt) && !is.na(results$tests$lrt_shared_disp) && results$tests$lrt_shared_disp < 0.05 && 
          model_to_return <- model_strategy
        } else if (valid_shared) {
          model_to_return <- model_shared
        # } else if (valid_strategy) {
        #   model_to_return <- model_strategy
        } else {
          model_to_return <- NULL
        }
        
      } else if (dispersion_modeling == "shared") {
        # Only fit shared model
        model_to_return <- model_shared
        if (!is.null(model_to_return) && model_to_return$fit$convergence == 0 && isTRUE(model_to_return$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_to_return)
          if (isTRUE(lrt) && !is.null(model_null_shared) && model_null_shared$fit$convergence == 0 && isTRUE(model_null_shared$sdr$pdHess)) {
            results$tests$pvalue <- anova(model_to_return, model_null_shared)$P[2]
            # results$tests$pvalue_best <- results$tests$pvalue
          }
          if (isTRUE(diagnostic)) {
            sim_out <- simulateResiduals(fittedModel = model_to_return, plot = FALSE)
            results$diagnostics$diagnostics_shared <- list(overdispersion = testDispersion(sim_out, plot = FALSE)$p.value, 
                                                           zeroInflation = testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                           uniformity = testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                           residualsDispersion = testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                           outliers = testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
      } else if (dispersion_modeling == "custom") {
        # Only fit custom model
        model_to_return <- model_custom
        if (!is.null(model_to_return) && model_to_return$fit$convergence == 0 && isTRUE(model_to_return$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_to_return)
          if (isTRUE(lrt) && !is.null(model_null) && model_null$fit$convergence == 0 && isTRUE(model_null$sdr$pdHess)) {
            results$tests$pvalue <- anova(model_to_return, model_null)$P[2]
            # results$tests$pvalue_best <- results$tests$pvalue
          }
          if (isTRUE(diagnostic)) {
            sim_out <- simulateResiduals(fittedModel = model_to_return, plot = FALSE)
            results$diagnostics$diagnostics_custom <- list(overdispersion = testDispersion(sim_out, plot = FALSE)$p.value, 
                                                             zeroInflation = testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                             uniformity = testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                             residualsDispersion = testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                             outliers = testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
      }
      
      # Check if a valid model was chosen and fitted
      if (is.null(model_to_return) || model_to_return$fit$convergence != 0 || !isTRUE(model_to_return$sdr$pdHess)) {
        warning(paste("Model fit failed for ORF", current_orf, "due to convergence or non-positive-definite Hessian."))
        return(.StatModel(type = "fitError", 
                          results = results, 
                          model = NA))
      }
      
      raw_rho_rna <- .calculate_mean_dispersion(model_data_this_orf[model_data_this_orf$strategy == 0, ]$counts, model_data_this_orf[model_data_this_orf$strategy == 0, ]$success + model_data_this_orf[model_data_this_orf$strategy == 0, ]$failure)
      raw_rho_ribo <- .calculate_mean_dispersion(model_data_this_orf[model_data_this_orf$strategy == 1, ]$counts, model_data_this_orf[model_data_this_orf$strategy == 1, ]$success + model_data_this_orf[model_data_this_orf$strategy == 1, ]$failure)
      
      results$estimates$mean_prop_rna <- raw_rho_rna$mean_prop
      results$estimates$mean_prop_ribo <- raw_rho_ribo$mean_prop
      results$dispersion$raw_rho_rna <- raw_rho_rna$raw_rho
      results$dispersion$raw_rho_ribo <- raw_rho_ribo$raw_rho
      
      disp_coefs <- fixef(model_to_return)$disp
      names(disp_coefs) <- paste0("disp~", names(disp_coefs))
      
      disp_intercept <- disp_coefs["disp~(Intercept)"]
      results$dispersion$fitted_disp_rna <- exp(disp_intercept)
      
      if ((dispersion_modeling == "auto" && !is.null(model_strategy))) { # dispersion_modeling == "strategy" || 
        disp_ribo_effect <- disp_coefs["disp~strategy1"]
        results$dispersion$fitted_disp_ribo <- exp(disp_intercept + disp_ribo_effect)
      } else {
        results$dispersion$fitted_disp_ribo <- results$dispersion$fitted_disp_rna
      }
      
      results$dispersion$fittedVsRawDispDiff <- abs(results$dispersion$fitted_disp_rna - results$dispersion$raw_rho_rna) + abs(results$dispersion$fitted_disp_ribo - results$dispersion$raw_rho_ribo)
      
      return(.StatModel(type = "glmmTMB", 
                        results = results, 
                        model = model_to_return))
      
    },
    error = function(e) {
      warning(paste("Model fit failed for ORF", current_orf, "with error:", e$message))
      return(.StatModel(
        type = "fitError",
        results = list(),
        model = NA
      ))
    })
    
    return(result)
  })
  
  names(models_gene) <- orfs_to_keep
  return(models_gene)
}


#' @title Fit Differential ORF Usage Models
#'
#' @description
#' This internal function fits beta-binomial models for differential ORF usage (DOU)
#' across all genes. It supports multiple dispersion modeling strategies and optional
#' diagnostics using DHARMa. This function is adapted from the \code{satuRn} package
#' to support beta-binomial GLM/GLMMs via \code{glmmTMB}.
#'
#' @param countData A numeric matrix of ORF-level counts (rows = ORFs, columns = samples).
#' @param orf2gene A data frame mapping ORF IDs to gene IDs. Must contain columns \code{orf_id} and \code{gene_id}.
#' @param anno A data frame containing sample annotations. Must include columns such as \code{condition}, \code{strategy}, and \code{replicate}.
#' @param formula A formula object specifying the model design, e.g., \code{~ condition * strategy}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param dispformula Optional formula object specifying a custom dispersion model (used when \code{dispersion_modeling = "custom"}).
#' @param dispersion_modeling Character string specifying the dispersion modeling strategy.
#'   Options are:
#'   \describe{
#'     \item{\code{"auto"}}{Fit both strategy-dependent and shared dispersion models, and select the best via likelihood ratio test.}
#'     \item{\code{"strategy"}}{Model dispersion as a function of sequencing strategy.}
#'     \item{\code{"shared"}}{Assume constant dispersion across all predictor levels.}
#'     \item{\code{"custom"}}{Use a user-specified dispersion formula via \code{dispformula}.}
#'   }
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, runs DHARMa diagnostics to assess model fit (default: \code{FALSE}).
#' @param parallel Logical; if \code{TRUE}, enables parallel processing using \code{BiocParallel} (default: \code{FALSE}).
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{FALSE}).
#'
#' @return A named \code{list} of \code{StatModel} objects, one per ORF.
#' 
#' @importFrom S4Vectors mcols mcols<-
#' 
#' @keywords internal
#' 
.fitDOU_internal <- function(countData,
                             orf2gene,
                             anno,
                             formula,
                             target,
                             dispformula,
                             dispersion_modeling,
                             lrt,
                             diagnostic, 
                             parallel,
                             optimizers,
                             seed,
                             verbose) {
  
  stopifnot(class(countData)[1] %in% c("matrix", "data.frame", "dgCMatrix", "DelayedMatrix"))
  
  single_orf <- !mcols(orf2gene)$gene_id %in% mcols(orf2gene)$gene_id[duplicated(mcols(orf2gene)$gene_id)]
  if (any(single_orf)) {
    message("Genes with only one type of transcript/isoform detected. Results of such transcripts will be set to NA and flagged as lonelyTranscript.")
  }
  
  geneForEachOrf <- orf2gene$gene_id[match(rownames(countData), rownames(orf2gene))]
  geneForEachOrf <- as.character(geneForEachOrf)
  stopifnot(length(geneForEachOrf) == nrow(countData))
  
  # countData_ribo <- countData[, grep("ribo", colnames(countData))]
  # countData_rna <- countData[, grep("rna", colnames(countData))]
  
  # Get the logical indices from the 'strategy' column of the sample annotation (anno)
  strategy_levels <- levels(anno$strategy)
  if (length(strategy_levels) > 2) {
    warning("Expected two strategy levels got: ", strategy_levels)
    warning("Use ", strategy_levels[1], " for RNA-seq and ", strategy_levels[2], " for Ribo-seq")
  }
  ribo_idx <- anno$strategy == strategy_levels[2]
  rna_idx <- anno$strategy == strategy_levels[1]
  
  # Subset the countData matrix using these indices
  countData_ribo <- countData[, ribo_idx]
  countData_rna <- countData[, rna_idx]
  
  gene_ids <- levels(factor(geneForEachOrf))
  
  fit_model_for_gene <- function(gene) {
    idx <- which(geneForEachOrf == gene)
    ribo_mat <- countData_ribo[idx, , drop = FALSE]
    rna_mat <- countData_rna[idx, , drop = FALSE]
    
    tryCatch({
      models_gene <- .fitBetaBinomial(ribo_mat = ribo_mat, rna_mat = rna_mat, anno = anno, 
                                      formula = formula, target = target,
                                      dispersion_modeling = dispersion_modeling, dispformula = dispformula, lrt = lrt,
                                      diagnostic = diagnostic, seed = seed, 
                                      parallel = parallel, optimizers = optimizers,
                                      verbose = verbose)
      return(models_gene)
    },
    error = function(e) {
      warning(paste("Model fit failed for gene", gene, "with a critical error:", e$message))
      # Return a named list with a single fitError model
      fit_error_model <- setNames(
        list(.StatModel(type = "fitError", results = list(), model = NA)),
        gene
      )
      return(list(fit_error_model))
    })
    
    return(models_gene)
  }
  
  # if (parallel) {
  #   models <- BiocParallel::bplapply(gene_ids, fit_model_for_gene, BPPARAM = BPPARAM)
  # } else if (verbose) {
  if (verbose) {
    models <- pbapply::pblapply(gene_ids, fit_model_for_gene)
  } else {
    models <- lapply(gene_ids, fit_model_for_gene)
  }
  
  # Filter out any NULL elements that might have slipped through
  all_models <- models[!sapply(models, is.null)]
  
  # Unlist to get a single named list of models for all ORFs
  all_models <- unlist(all_models, recursive = FALSE)
  
  # Create a final list with all ORFs from the original countData
  final_models <- setNames(
    lapply(rownames(countData), function(x) {
      # Look up the model from the list of fitted models
      model_obj <- all_models[[x]]
      if (is.null(model_obj)) {
        # If no model was found, it's a zero-count ORF, create a fitError model
        return(.StatModel(type = "fitError", results = list(), model = NA))
      }
      return(model_obj)
    }),
    rownames(countData)
  )
  
  return(final_models)
}


#' @title Fit Differential ORF Usage Models
#'
#' @description
#' This function fits beta-binomial generalized linear models (GLMs) or generalized linear mixed models (GLMMs) for
#' differential ORF usage (DOU) across all genes in a \code{SummarizedExperiment} object.
#' It supports multiple dispersion modeling strategies and optional DHARMa diagnostics.
#'
#' This function is adapted from the \code{satuRn} package to support beta-binomial GLM/GLMMs via \code{glmmTMB}.
#'
#' @param object A \code{SummarizedExperiment}, \code{RangedSummarizedExperiment}, or \code{SingleCellExperiment} object.
#'   The assay slot must contain transcript-level expression counts as a \code{matrix}, \code{DataFrame},
#'   \code{sparseMatrix}, or \code{DelayedMatrix}. The \code{rowData} must include columns \code{orf_id}
#'   and \code{gene_id}, and \code{colData} must contain sample annotations. The model formula should be
#'   specified in the metadata or passed directly.
#' @param formula A formula object specifying the model design, e.g., \code{~ 0 + group}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param dispersion_modeling Character string specifying the dispersion modeling strategy.
#'   Options: \code{"auto"}, \code{"strategy"}, \code{"shared"}, or \code{"custom"}.
#' @param dispformula Optional formula object for custom dispersion modeling.
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, enables DHARMa diagnostics (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization,
#'   e.g., \code{list(n = 4L, autopar = TRUE)}.
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#' @param ... A placeholder for further enhancements.
#' 
#' @return An updated \code{SummarizedExperiment} object with a new list of models
#'   stored in \code{rowData(object)[["fitDOUModels"]]}.
#'
#'
#' @rdname fitDOU
#' @author Jeroen Gilis (original), Chun Shen Lim (modifications)
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom DHARMa simulateResiduals testDispersion testZeroInflation testResiduals
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats model.matrix
#' @importFrom pbapply pblapply
#' @importFrom BiocParallel bplapply bpparam
#' 
# setGeneric("fitDOU", function(object,
#                               formula,
#                               target,
#                               dispersion_modeling,
#                               dispformula,
#                               lrt = FALSE,
#                               diagnostic = FALSE,
#                               parallel = list(n=4L, autopar=TRUE),
#                               optimizers = FALSE,
#                               seed = NULL,
#                               verbose = TRUE) {
#   standardGeneric("fitDOU")
# })
setGeneric("fitDOU", function(object, ...) standardGeneric("fitDOU"))

# Wrapper function
setMethod(
  f = "fitDOU",
  signature = "SummarizedExperiment",
  function(object,
           formula,
           target,
           dispersion_modeling,
           dispformula,
           lrt = FALSE,
           diagnostic = FALSE,
           parallel = list(n=4L, autopar=TRUE),
           optimizers = FALSE,
           seed = NULL,
           verbose = TRUE) {
    if (ncol(colData(object)) == 0) 
      stop("colData is empty")
    
    if (!"gene_id" %in% colnames(rowData(object)) && all(grepl("O", names(rowData(object))))) {
      stop("rowData does not contain columns gene_id and orf_id")
    }
    if (!all(rownames(object) == rownames(rowData(object)))) {
      stop("not all row names of the expression matrix match 
                 the orf_id column of the object's rowData")
    }
    
    rowData(object)[["fitDOUModels"]] <- .fitDOU_internal(
      countData = assay(object),
      orf2gene = rowData(object), 
      anno = colData(object),
      formula = formula,
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
    
    return(object)
  }
)
