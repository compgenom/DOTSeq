#' Generate Contrast Matrix for Differential Testing
#'
#' This function constructs a contrast matrix for quasi-binomial or GLM-based
#' differential testing. It uses the sample annotation (`colData`) from a
#' `SummarizedExperiment` object and a model formula to determine pairwise
#' contrasts between experimental conditions.
#'
#' @param sumExp A `SummarizedExperiment` object containing count data and
#'   sample annotation.
#' @param formula A formula object specifying the design, e.g., ~ 0 + replicate + condition * strategy.
#' @param baseline Optional character specifying the baseline condition for contrasts.
#'   If NULL, the function attempts to infer it automatically (default: NULL).
#' @param verbose Logical; if TRUE, prints progress messages. Default is FALSE.
#'
#' @return A numeric contrast matrix for differential testing in `testDOT`.
#'
#' @examples
#' # Construct contrast matrix from a SummarizedExperiment object
#' L <- contrastMatrix(sumExp, ~ 0 + replicate + condition * strategy, baseline = "control", verbose = TRUE)
#'
#' @export
contrastMatrix <- function(sumExp, formula, baseline = NULL, verbose = FALSE) {
  # Build design matrix
  anno <- colData(sumExp)
  design <- model.matrix(formula, data = anno)
  
  contrastFactors <- grep(":strategy1", colnames(design), value = TRUE)
  
  L_all <- NULL
  valid_rows <- colnames(design)
  
  for (contrastFactor in contrastFactors) {
    if (verbose) cat("\nProcessing:", contrastFactor, "\n")
    
    if (!is.factor(design[, contrastFactor])) {
      design[, contrastFactor] <- factor(design[, contrastFactor])
    }
    
    conditions <- unique(sapply(strsplit(names(design[, contrastFactor]), "\\."), `[`, 2))
    nconditions <- length(conditions)
    
    if (verbose) {
      cat(" - Number of valid conditions:", nconditions, "\n")
      cat(" - Valid conditions:", paste(conditions, collapse = ", "), "\n")
    }
    
    # Infer baseline if not provided
    if (is.null(baseline)) {
      # valid_levels <- paste0(contrastFactor, conditions)
      # missing_levels <- setdiff(valid_levels, colnames(design))
      
      valid_levels <- sapply(strsplit(contrastFactors, ":"), `[`, 1)
      valid_levels <- sub("^condition", "", valid_levels)
      missing_levels <- setdiff(conditions, valid_levels)
      
      if (length(missing_levels) == 1) {
        inferred_baseline <- sub(contrastFactor, "", missing_levels)
        local_baseline <- inferred_baseline
        if (verbose) cat(" - Inferred baseline:", local_baseline, "\n")
      } else {
        local_baseline <- conditions[1]
        if (verbose) cat(" - Could not infer baseline uniquely. Using:", local_baseline, "\n")
      }
    } else {
      if (!(baseline %in% conditions)) {
        stop(paste(" - Baseline", baseline, "not found among valid conditions:", paste(conditions, collapse = ", ")))
      }
      local_baseline <- baseline
      if (verbose) cat(" - Using user-specified baseline:", local_baseline, "\n")
    }
    
    if (nconditions > 1) {
      combinations <- combn(conditions, 2, simplify = FALSE)
      L <- matrix(0, nrow = length(valid_rows), ncol = length(combinations),
                  dimnames = list(valid_rows, rep("", length(combinations))))
      
      counter <- 1
      for (combo in combinations) {
        cond1 <- as.character(combo[[1]])
        cond2 <- as.character(combo[[2]])
        
        # row1 <- paste0(contrastFactor, cond1)
        # row2 <- paste0(contrastFactor, cond2)
        row1 <- grep(paste0(cond1, ".*:strategy1$"), colnames(design), value = TRUE)
        row2 <- grep(paste0(cond2, ".*:strategy1$"), colnames(design), value = TRUE)
        
        if (!is.null(local_baseline) && local_baseline %in% combo) {
          non_baseline <- ifelse(cond1 == local_baseline, cond2, cond1)
          contrast_name <- paste0(non_baseline, "_vs_", local_baseline)
          
          # Find the correct row for the non-baseline condition
          row <- grep(paste0(non_baseline, ".*:strategy1$"), colnames(design), value = TRUE)
          
          if (length(row) == 1 && row %in% valid_rows) {
            L[row, counter] <- 1
            if (verbose) cat(" - Assigned 1 to", row, "for", contrast_name, "\n")
          } else if (verbose) {
            cat(" - Row not found or ambiguous for non-baseline condition:", non_baseline, "\n")
          }
          
          colnames(L)[counter] <- contrast_name
        
        } else {
          # fallback to pairwise contrast
          contrast_name <- paste0(cond1, "_vs_", cond2)
          if (row1 %in% valid_rows) L[row1, counter] <- 1
          if (row2 %in% valid_rows) L[row2, counter] <- -1
          colnames(L)[counter] <- contrast_name
          
          if (verbose) {
            cat(" - Assigned",
                if (row1 %in% valid_rows) paste("+1 to", row1) else "",
                if (row2 %in% valid_rows) paste("-1 to", row2) else "",
                "for", contrast_name, "\n")
          }
        }
        
        counter <- counter + 1
        
      }
    } else {
      if (verbose) cat(" - Skipping contrast generation: not enough valid conditions.\n")
      next
    }
    
  }
  
  # # Check if the contrast matrix has a row named "(Intercept)"
  # if ("(Intercept)" %in% rownames(L)) {
  #   disp_intercept <- rep(0, ncol(L))
  #   L <- rbind(L, disp_intercept)
  #   rownames(L)[nrow(L)] <- "disp~(Intercept)"
  # }
  # 
  # if (modelType == "betabinomial") {
  #   disp_coefs <- rep(0, ncol(L))
  #   L <- rbind(L, disp_coefs)
  #   rownames(L)[nrow(L)] <- "disp~strategy1"
  # }
  
  if (verbose) {
    cat("Contrast matrix L:\n")
    print(L)
  }
  
  return(contrastMatrix = L)
}
