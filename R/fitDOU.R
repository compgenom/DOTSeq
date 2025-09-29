#' @title The StatModel class for DOTSeq
#'
#' @description
#' The \code{StatModel} class represents a statistical model fitted to ORF-level data
#' using a beta-binomial GLMM via \code{glmmTMB}. It is used to store model results,
#' diagnostics, and metadata for each ORF.
#'
#' Objects of this class are typically created by the user-level function \code{fitDOU()},
#' or manually using the \code{StatModel()} constructor. In the DOTSeq pipeline,
#' each ORF is assigned a \code{StatModel} object, which is stored in a \code{DataFrame}
#' and embedded in the \code{rowData} slot of a \code{SummarizedExperiment}.
#'
#' @slot type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#'   If the model is successfully fitted, the type is typically \code{"glmmTMB"}.
#' @slot results A \code{list} containing model results, parameters, and test statistics.
#' @slot model An object of class \code{ANY} storing the fitted model object (e.g., from \code{glmmTMB}).
#'
#' @export
#'
#' @examples
#' ## Create a dummy StatModel object
#' myModel <- StatModel(
#'   type = "glmmTMB",
#'   results = list(
#'     model_fit = list(aic = 96.03),
#'     tests = list(pvalueBest = 0.355)
#'   )
#' )
#' myModel
.StatModel <- setClass("StatModel",
                       slots = c(
                         type = "character",
                         results = "list",
                         model = "ANY"
                       )
)

#' @title Construct a StatModel Object
#'
#' @description
#' This function constructs a new \code{StatModel} object, which stores the results
#' of a statistical model fitted to ORF-level data. It is typically used internally
#' by the DOTSeq pipeline or manually for testing and diagnostics.
#'
#' @param type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#' @param results A \code{list} containing model results, parameters, and test statistics.
#' @param model An optional fitted model object (e.g., from \code{glmmTMB}). Default is \code{NA_real_}.
#'
#' @return A \code{StatModel} S4 object.
#'
#' @examples
#' ## Create a dummy StatModel object
#' myModel <- StatModel(
#'   type = "glmmTMB",
#'   results = list(x = 3, y = 7, b = 4)
#' )
#' myModel
#'
#' @importFrom methods new
#' @export
StatModel <- function(type = "fitError",
                      results = list(),
                      model = NA_real_) {
  out <- new("StatModel")
  out@type <- type
  out@results <- results
  out@model <- model
  return(out)
}


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
#'   \item{meanProp}{Mean proportion of ORF counts relative to total gene counts.}
#'   \item{rawRho}{Unclamped dispersion estimate (can be negative).}
#'   \item{rawIntraClassCorr}{Clamped dispersion estimate (non-negative).}
#' }
.calculate_mean_dispersion <- function(counts, totals) {
  # Avoid division by zero
  if (length(counts) < 2) {
    return(list(meanProp = NA_real_, rawDisp = NA_real_, clampedDisp = NA_real_))
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
  return(list(meanProp = mu_hat, 
              rawRho = rho_hat,
              rawIntraClassCorr = max(0, rho_hat)))
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
#' @param dispersionStrategy Character string specifying the dispersion modeling strategy.
#'   Options are:
#'   \describe{
#'     \item{\code{"auto"}}{Fit both strategy-dependent and shared dispersion models, and select the best via likelihood ratio test.}
#'     \item{\code{"strategy"}}{Model dispersion as a function of sequencing strategy.}
#'     \item{\code{"shared"}}{Assume constant dispersion across all predictor levels.}
#'     \item{\code{"custom"}}{Use a user-specified dispersion formula via \code{dispformula}.}
#'   }
#' @param dispformula Optional formula object specifying a custom dispersion model (used when \code{dispersionStrategy = "custom"}).
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, runs DHARMa diagnostics to assess model fit (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization,
#'   e.g., \code{list(parallel = TRUE, ncpus = 4)}. Default: \code{list(n = 8L, autopar = TRUE)}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{FALSE}).
#'
#' @return A named \code{list} of \code{StatModel} objects, one for each ORF in the gene.
.fitBetaBinomial <- function(ribo_mat, 
                             rna_mat, 
                             anno, 
                             formula = ~ condition * strategy, 
                             target = NULL,
                             dispersionStrategy = c("auto", "strategy", "shared", "custom"), 
                             dispformula = NULL, 
                             lrt = FALSE,
                             diagnostic = FALSE, 
                             seed = NULL, 
                             optimizers = FALSE,
                             parallel = list(n=8L, autopar=TRUE), 
                             verbose = FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  dispersionStrategy <- match.arg(dispersionStrategy)
  
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
  
  # Prepare the data for a single gene dynamically based on the design matrix
  long_data <- data.frame(
    counts = c(as.vector(rna_mat), as.vector(ribo_mat)),
    ORF = factor(rep(orf_names, times = num_samples * 2)),
    row.names = NULL
  )
  
  # Replicate the design variables for each ORF
  long_design <- as.data.frame(anno[rep(seq_len(nrow(anno)), each = num_orfs), ])
  
  # Bind the replicated design matrix to the long data frame
  long_data <- cbind(long_data, long_design)

  # long_data$gene <- sapply(strsplit(as.character(long_data$ORF), ":"), `[`, 1)
  # long_data$gene <- factor(long_data$gene)
  
  # Corrected check for negative counts
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
    total_gene_counts <- c(colSums(rna_mat), colSums(ribo_mat))
    
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
      
      fit_glmm <- function(formula, dispformula, data,
                           optimizers = c("nlminb", "bobyqa", "optim"), max_iter = 1000) {
        for (opt in optimizers) {
          opt_args <- if (opt == "optim") list(method = "BFGS") else list()
          model <- try(glmmTMB(formula, dispformula = dispformula, family = betabinomial(),
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
      
      if (isTRUE(optimizers)) {
        # Fit Models based on dispersionStrategy
        if (dispersionStrategy %in% c("auto", "strategy")) {
          model_strategy <- fit_glmm(formula = full_formula, dispformula = ~strategy, data = model_data_this_orf)
          if (isTRUE(lrt)) {
            model_null <- fit_glmm(formula = null_formula, dispformula = ~ strategy, data = model_data_this_orf)
          }
        }
        
        if (dispersionStrategy %in% c("auto", "shared")) {
          model_shared <- fit_glmm(formula = full_formula, dispformula = ~1, data = model_data_this_orf)
          if (isTRUE(lrt)) {
            model_null_shared <- fit_glmm(formula = null_formula, dispformula = ~1, data = model_data_this_orf)
          }
        }
        
        if (dispersionStrategy == "custom") {
          if  (!is.null(dispformula)) {
            model_custom <- fit_glmm(formula = full_formula, dispformula = dispformula, data = model_data_this_orf)
            if (isTRUE(lrt)) {
              model_null <- fit_glmm(formula = null_formula, dispformula = dispformula, data = model_data_this_orf)
            }
          } else {
            stop("Please provide dispformula.")
          }
        }
      } else {
        if (dispersionStrategy %in% c("auto", "strategy")) {
          model_strategy <- glmmTMB(full_formula, dispformula = ~strategy, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          if (isTRUE(lrt)) {
            model_null <- glmmTMB(null_formula, dispformula = ~ strategy, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          }
        }
        
        if (dispersionStrategy %in% c("auto", "shared")) {
          model_shared <- glmmTMB(full_formula, dispformula = ~1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          if (isTRUE(lrt)) {
            model_null_shared <- glmmTMB(null_formula, dispformula = ~1, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
          }
        }
        
        if (dispersionStrategy == "custom") {
          if  (!is.null(dispformula)) {
            model_custom <- glmmTMB(full_formula, dispformula = dispformula, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            if (isTRUE(lrt)) {
              model_null <- glmmTMB(null_formula, dispformula = dispformula, family = betabinomial, data = model_data_this_orf, control = glmmTMBControl(parallel = parallel))
            }
          } else {
            stop("Please provide dispformula.")
          }
        } 
      }
      # Build return object based on chosen strategy
      model_to_return <- NULL
      
      if (dispersionStrategy == "auto") {
        # Full comparison mode
        
        # Add strategy-specific metrics
        if (!is.null(model_strategy) && model_strategy$fit$convergence == 0 && isTRUE(model_strategy$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_strategy)
          if (isTRUE(lrt)) {
            results$tests$pvalue <- anova(model_strategy, model_null)$P[2]
          }
          if (isTRUE(diagnostic)) {
            sim_out <- DHARMa::simulateResiduals(fittedModel = model_strategy, plot = FALSE)
            results$diagnostics$diagnostics_strategy <- list(overdispersion = DHARMa::testDispersion(sim_out, plot = FALSE)$p.value, 
                                                             zeroInflation = DHARMa::testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                             uniformity = DHARMa::testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                             residualsDispersion = DHARMa::testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                             outliers = DHARMa::testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
        
        # Add shared-specific metrics
        if (!is.null(model_shared) && model_shared$fit$convergence == 0 && isTRUE(model_shared$sdr$pdHess)) {
          results$model_fit$aicSharedDisp <- AIC(model_shared)
          
          if (isTRUE(lrt)) {
            results$tests$lrtSharedDisp <- anova(model_strategy, model_shared)$P[2]
            results$tests$lrtShared <- anova(model_shared, model_null_shared)$P[2]
          }
          if (isTRUE(diagnostic)) {
            sim_out <- DHARMa::simulateResiduals(fittedModel = model_shared, plot = FALSE)
            results$diagnostics$diagnostics_shared <- list(overdispersion = DHARMa::testDispersion(sim_out, plot = FALSE)$p.value, 
                                                           zeroInflation = DHARMa::testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                           uniformity = DHARMa::testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                           residualsDispersion = DHARMa::testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                           outliers = DHARMa::testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
        
        # Determine the final p-value based on the gene-by-gene selection logic
        if (isTRUE(lrt) && !is.null(results$tests$lrtSharedDisp) && !is.na(results$tests$lrtSharedDisp) && results$tests$lrtSharedDisp < 0.05) {
          if (!is.null(results$tests$pvalue) && !is.na(results$tests$pvalue)) {
            results$tests$pvalueBest <- results$tests$pvalue
          } else {
            results$tests$pvalueBest <- results$tests$lrtShared
          }
        } else {
          results$tests$pvalueBest <- results$tests$lrtShared
        }
        
        # The 'auto' use the strategy model as the default.
        model_to_return <- model_strategy
        
      } else if (dispersionStrategy == "strategy") {
        # Only fit strategy model
        model_to_return <- model_strategy
        if (!is.null(model_to_return) && model_to_return$fit$convergence == 0 && isTRUE(model_to_return$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_to_return)
          if (isTRUE(lrt)) {
            results$tests$pvalue <- anova(model_to_return, model_null)$P[2]
            results$tests$pvalueBest <- results$tests$pvalue
          }
          if (isTRUE(diagnostic)) {
            sim_out <- DHARMa::simulateResiduals(fittedModel = model_to_return, plot = FALSE)
            results$diagnostics$diagnostics_strategy <- list(overdispersion = DHARMa::testDispersion(sim_out, plot = FALSE)$p.value, 
                                                             zeroInflation = DHARMa::testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                             uniformity = DHARMa::testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                             residualsDispersion = DHARMa::testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                             outliers = DHARMa::testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
        
      } else if (dispersionStrategy == "shared") {
        # Only fit shared model
        model_to_return <- model_shared
        if (!is.null(model_to_return) && model_to_return$fit$convergence == 0 && isTRUE(model_to_return$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_to_return)
          if (isTRUE(lrt)) {
            results$tests$pvalue <- anova(model_to_return, model_null_shared)$P[2]
            results$tests$pvalueBest <- results$tests$pvalue
          }
          if (isTRUE(diagnostic)) {
            sim_out <- DHARMa::simulateResiduals(fittedModel = model_to_return, plot = FALSE)
            results$diagnostics$diagnostics_shared <- list(overdispersion = DHARMa::testDispersion(sim_out, plot = FALSE)$p.value, 
                                                           zeroInflation = DHARMa::testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                           uniformity = DHARMa::testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                           residualsDispersion = DHARMa::testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                           outliers = DHARMa::testResiduals(sim_out, plot = FALSE)$outliers$p.value)
          }
        }
      } else if (dispersionStrategy == "custom") {
        # Only fit custom model
        model_to_return <- model_custom
        if (!is.null(model_to_return) && model_to_return$fit$convergence == 0 && isTRUE(model_to_return$sdr$pdHess)) {
          results$model_fit$aic <- AIC(model_to_return)
          if (isTRUE(lrt)) {
            results$tests$pvalue <- anova(model_to_return, model_null)$P[2]
            results$tests$pvalueBest <- results$tests$pvalue
          }
          if (isTRUE(diagnostic)) {
            sim_out <- DHARMa::simulateResiduals(fittedModel = model_to_return, plot = FALSE)
            results$diagnostics$diagnostics_strategy <- list(overdispersion = DHARMa::testDispersion(sim_out, plot = FALSE)$p.value, 
                                                             zeroInflation = DHARMa::testZeroInflation(sim_out, plot = FALSE)$p.value,
                                                             uniformity = DHARMa::testResiduals(sim_out, plot = FALSE)$uniformity$p.value,
                                                             residualsDispersion = DHARMa::testResiduals(sim_out, plot = FALSE)$dispersion$p.value,
                                                             outliers = DHARMa::testResiduals(sim_out, plot = FALSE)$outliers$p.value)
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
      
      if (isTRUE(lrt)) {
        # Extract and calculate dispersion-related metrics for the chosen model
        coeffs <- summary(model_to_return)$coefficients$cond
        interaction_term <- grep("condition.*?:strategy.*?", rownames(coeffs), value = TRUE)
        
        selected_term <- NULL
        
        if (length(interaction_term) == 1) {
          selected_term <- interaction_term[1]
          
        } else if (length(interaction_term) > 1 && !is.null(target)) {
          target_match <- grep(paste0("condition", target, ":strategy"), interaction_term, value = TRUE)
          if (length(target_match) == 1) {
            selected_term <- target_match
          } else {
            warning("Target condition not found. Using first available interaction term.")
            selected_term <- interaction_term[1]
          }
          
        } else {
          if (isTRUE(verbose)) {
            message("Detected interaction terms: ", paste(interaction_term, collapse = ", "), ". Using first available term.")
          }
          selected_term <- interaction_term[1]
        }
        
        # Final assignment
        results$shrinkage$estimate <- coeffs[selected_term, "Estimate"]
        results$shrinkage$se <- coeffs[selected_term, "Std. Error"]
      }
      
      raw_rho_rna <- .calculate_mean_dispersion(model_data_this_orf[model_data_this_orf$strategy == 0, ]$counts, model_data_this_orf[model_data_this_orf$strategy == 0, ]$success + model_data_this_orf[model_data_this_orf$strategy == 0, ]$failure)
      raw_rho_ribo <- .calculate_mean_dispersion(model_data_this_orf[model_data_this_orf$strategy == 1, ]$counts, model_data_this_orf[model_data_this_orf$strategy == 1, ]$success + model_data_this_orf[model_data_this_orf$strategy == 1, ]$failure)
      
      results$estimates$meanPropRna <- raw_rho_rna$meanProp
      results$estimates$meanPropRibo <- raw_rho_ribo$meanProp
      results$dispersion$rawDispRna <- raw_rho_rna$rawDisp
      results$dispersion$rawDispRibo <- raw_rho_ribo$rawDisp
      
      disp_coefs <- fixef(model_to_return)$disp
      names(disp_coefs) <- paste0("disp~", names(disp_coefs))
      
      disp_intercept <- disp_coefs["disp~(Intercept)"]
      results$dispersion$fittedDispRna <- exp(disp_intercept)
      
      if (dispersionStrategy == "strategy" || (dispersionStrategy == "auto" && !is.null(model_strategy))) {
        disp_ribo_effect <- disp_coefs["disp~strategy1"]
        results$dispersion$fittedDispRibo <- exp(disp_intercept + disp_ribo_effect)
      } else {
        results$dispersion$fittedDispRibo <- results$dispersion$fittedDispRna
      }
      
      results$dispersion$fittedVsRawDispDiff <- abs(results$dispersion$fittedDispRna - results$dispersion$rawDispRna) + abs(results$dispersion$fittedDispRibo - results$dispersion$rawDispRibo)
      
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
#' @param dispformula Optional formula object specifying a custom dispersion model (used when \code{dispersionStrategy = "custom"}).
#' @param dispersionStrategy Character string specifying the dispersion modeling strategy.
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
#' @param BPPARAM An object specifying the parallel back-end, e.g., \code{MulticoreParam()} or \code{SnowParam()}.
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{FALSE}).
#'
#' @return A named \code{list} of \code{StatModel} objects, one per ORF.
#'
#' @importFrom ashr ash
.fitDOU_internal <- function(countData,
                             orf2gene,
                             anno,
                             formula,
                             target,
                             dispformula,
                             dispersionStrategy,
                             lrt,
                             diagnostic, 
                             parallel,
                             optimizers,
                             seed,
                             verbose) {
  
  stopifnot(class(countData)[1] %in% c("matrix", "data.frame", "dgCMatrix", "DelayedMatrix"))
  
  single_orf <- !orf2gene$gene_id %in% orf2gene$gene_id[duplicated(orf2gene$gene_id)]
  if (any(single_orf)) {
    message("Genes with only one type of transcript/isoform detected. Results of such transcripts will be set to NA and flagged as lonelyTranscript.")
  }
  
  geneForEachOrf <- orf2gene$gene_id[match(rownames(countData), orf2gene$isoform_id)]
  geneForEachOrf <- as.character(geneForEachOrf)
  stopifnot(length(geneForEachOrf) == nrow(countData))
  
  countData_ribo <- countData[, grep("ribo", colnames(countData))]
  countData_rna <- countData[, grep("rna", colnames(countData))]
  gene_ids <- levels(factor(geneForEachOrf))
  
  fit_model_for_gene <- function(gene) {
    idx <- which(geneForEachOrf == gene)
    ribo_mat <- countData_ribo[idx, , drop = FALSE]
    rna_mat <- countData_rna[idx, , drop = FALSE]
    
    tryCatch({
      models_gene <- .fitBetaBinomial(ribo_mat = ribo_mat, rna_mat = rna_mat, anno = anno, 
                                      formula = formula, target = target,
                                      dispersionStrategy = dispersionStrategy, dispformula = dispformula, lrt = lrt,
                                      diagnostic = diagnostic, seed = seed, 
                                      parallel = parallel, optimizers = optimizers,
                                      verbose = verbose)
      return(models_gene)
    },
    error = function(e) {
      warning(paste("Model fit failed for gene", gene, "with a critical error:", e$message))
      # Return a named list with a single fitError model
      fit_error_model <- .StatModel(type = "fitError", results = list(), model = NA)
      names(fit_error_model) <- gene
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
#'   \code{sparseMatrix}, or \code{DelayedMatrix}. The \code{rowData} must include columns \code{isoform_id}
#'   and \code{gene_id}, and \code{colData} must contain sample annotations. The model formula should be
#'   specified in the metadata or passed directly.
#' @param formula A formula object specifying the model design, e.g., \code{~ 0 + group}.
#' @param target Character string specifying the non-reference condition level to extract the corresponding interaction term from the model. 
#' This is contrasted against the baseline condition (default: \code{NULL}).
#' @param dispersionStrategy Character string specifying the dispersion modeling strategy.
#'   Options: \code{"auto"}, \code{"strategy"}, \code{"shared"}, or \code{"custom"}.
#' @param dispformula Optional formula object for custom dispersion modeling.
#' @param lrt Logical; if \code{TRUE}, performs a likelihood ratio test to compare the full model (with interaction) against a reduced model 
#' (without interaction) to assess translation-specific effects (default: \code{FALSE}).
#' @param diagnostic Logical; if \code{TRUE}, enables DHARMa diagnostics (default: \code{FALSE}).
#' @param parallel A list passed to \code{glmmTMBControl} to configure parallel optimization,
#'   e.g., \code{list(n = 8L, autopar = TRUE)}.
#' @param optimizers Logical; if \code{TRUE}, enables brute-force optimization using multiple optimizers in \code{glmmTMB} (default: \code{FALSE}).
#' @param seed Optional integer to set the random seed for reproducibility (default: \code{NULL}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#'
#' @return An updated \code{SummarizedExperiment} object with a new list of models
#'   stored in \code{rowData(object)[["fitDOUModels"]]}.
#'
#' @examples
#' data(sumExp_example, package = "satuRn")
#' sumExp <- fitDOU(
#'   object = sumExp_example,
#'   formula = ~ 0 + group,
#'   dispersionStrategy = "auto",
#'   parallel = list(n = 4L, autopar = TRUE),
#'   verbose = TRUE
#' )
#'
#' @rdname fitDOU
#' @author Jeroen Gilis (original), Chun Shen Lim (modifications)
#'
#' @import SummarizedExperiment
#' @import DHARMa
#' @import glmmTMB
#' @importFrom stats model.matrix
#' @importFrom pbapply pblapply
#' @importFrom BiocParallel bplapply bpparam
#' @export
setGeneric("fitDOU", function(object, ...) standardGeneric("fitDOU"))

# Wrapper function
setMethod(
  f = "fitDOU", # Changed from fitDTU
  signature = "SummarizedExperiment",
  function(object,
           formula,
           target,
           dispersionStrategy,
           dispformula,
           lrt = FALSE,
           diagnostic = FALSE,
           parallel = list(n=8L, autopar=TRUE),
           optimizers = FALSE,
           seed = NULL,
           verbose = TRUE) {
    if (ncol(SummarizedExperiment::colData(object)) == 0) 
      stop("colData is empty")
    
    if (!"gene_id" %in% colnames(rowData(object)) | 
        !"isoform_id" %in% colnames(rowData(object))) {
      stop("rowData does not contain columns gene_id and isoform_id")
    }
    if (!all(rownames(object) == rowData(object)[, "isoform_id"])) {
      stop("not all row names of the expression matrix match 
                 the isoform_id column of the object's rowData")
    }
    
    rowData(object)[["fitDOUModels"]] <- .fitDOU_internal(
      countData = assay(object),
      orf2gene = rowData(object)[, c("isoform_id", "gene_id")],
      anno = colData(object), # Changed from design
      formula = formula,
      target = target,
      dispersionStrategy = dispersionStrategy, 
      dispformula = dispformula,
      lrt = lrt,
      diagnostic = diagnostic,
      parallel = parallel,
      optimizers = optimizers,
      seed = seed, # New
      verbose = verbose
    )
    
    return(object)
  }
)
