#' Generate Contrast Matrix for Differential Testing
#'
#' This function constructs a contrast matrix for quasi-binomial or GLM-based
#' differential testing. It uses the sample annotation (`colData`) from a
#' `SummarizedExperiment` object and a model formula to determine pairwise
#' contrasts between experimental conditions.
#'
#' @param sumExp A `SummarizedExperiment` object containing count data and
#'   sample annotation.
#' @param fmla A formula object specifying the design, e.g., ~ 0 + condition:effect.
#' @param baseline Optional character specifying the baseline condition for contrasts.
#'   If NULL, the function attempts to infer it automatically.
#' @param verbose Logical; if TRUE, prints progress messages. Default is FALSE.
#'
#' @return A numeric contrast matrix for differential testing in `testDOT`.
#'
#' @examples
#' # Construct contrast matrix from a SummarizedExperiment object
#' L <- contrastMatrix(sumExp, ~ 0 + effect1 + effect2, baseline = "control", verbose = TRUE)
#'
#' @export
contrastMatrix <- function(sumExp, fmla, baseline = NULL, verbose = FALSE) {
  # Build design matrix
  anno <- colData(sumExp)
  design <- model.matrix(fmla, data = anno)
  
  contrastFactors <- grep("effect", colnames(anno), value = TRUE)
  contrastFactors <- grep("effect2", contrastFactors, value = TRUE, invert = TRUE)
  
  L_all <- NULL
  valid_rows <- colnames(design)
  
  for (contrastFactor in contrastFactors) {
    if (verbose) cat("\nProcessing:", contrastFactor, "\n")
    
    conditions <- levels(anno[[contrastFactor]])
    conditions <- grep("none", conditions, value = TRUE, invert = TRUE)
    nconditions <- length(conditions)
    
    if (verbose) {
      cat(" - Number of valid conditions (excluding 'none'):", nconditions, "\n")
      cat(" - Valid conditions:", paste(conditions, collapse = ", "), "\n")
    }
    
    # Infer baseline if not provided
    if (is.null(baseline)) {
      valid_levels <- paste0(contrastFactor, conditions)
      missing_levels <- setdiff(valid_levels, colnames(design))
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
        
        row1 <- paste0(contrastFactor, cond1)
        row2 <- paste0(contrastFactor, cond2)
        
        if (!is.null(local_baseline) && local_baseline %in% combo) {
          non_baseline <- ifelse(cond1 == local_baseline, cond2, cond1)
          row <- paste0(contrastFactor, non_baseline)
          contrast_name <- paste0(non_baseline, "_vs_", local_baseline)
          
          if (row %in% valid_rows) {
            L[row, counter] <- 1
            if (verbose) cat(" - Assigned 1 to", row, "for", contrast_name, "\n")
          } else if (verbose) {
            cat(" - Row not found in design:", row, "\n")
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
    
    if (verbose) {
      cat("Contrast matrix L:\n")
      print(L)
    }
    
    if (!is.null(L)) {
      if (is.null(L_all)) {
        L_all <- L
      } else {
        all_rows <- union(rownames(L_all), rownames(L))
        all_rows <- intersect(all_rows, valid_rows)
        
        L_all_expanded <- matrix(0, nrow = length(all_rows), ncol = ncol(L_all),
                                 dimnames = list(all_rows, colnames(L_all)))
        L_expanded <- matrix(0, nrow = length(all_rows), ncol = ncol(L),
                             dimnames = list(all_rows, colnames(L)))
        
        common_rows_L_all <- intersect(rownames(L_all_expanded), rownames(L_all))
        common_cols_L_all <- intersect(colnames(L_all_expanded), colnames(L_all))
        if (length(common_rows_L_all) > 0 && length(common_cols_L_all) > 0) {
          L_all_expanded[common_rows_L_all, common_cols_L_all] <- L_all[common_rows_L_all, common_cols_L_all]
        }
        
        common_rows_L <- intersect(rownames(L_expanded), rownames(L))
        common_cols_L <- intersect(colnames(L_expanded), colnames(L))
        if (length(common_rows_L) > 0 && length(common_cols_L) > 0) {
          L_expanded[common_rows_L, common_cols_L] <- L[common_rows_L, common_cols_L]
        }
        
        L_all <- cbind(L_all_expanded, L_expanded)
      }
    }
  }
  
  if (verbose) {
    cat("\nFinal contrast matrix:\n")
    print(L_all)
  }
  
  return(contrastMatrix = L_all)
}
