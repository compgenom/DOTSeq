#' Generate Coefficients for Simulated Differential ORF Translation
#'
#' @description
#' Simulates log-fold change coefficients for differential ORF translation (DOT) analysis.
#' This function assigns baseline coefficients to all ORFs and applies specific log-fold
#' changes to multi-cistronic ORFs based on a selected regulatory scenario.
#'
#' @param orfs A data frame containing ORF annotations. Row names must align with the
#'   filtered count matrices. Must include columns `groupID` and `labels`.
#' @param scenario A character string specifying the regulatory scenario to simulate.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"uORF_up_mORF_down"}
#'     \item \code{"dORF_up_mORF_up"}
#'     \item \code{"uORF_up_mORF_flat"}
#'     \item \code{"uORF_down_mORF_flat"}
#'     \item \code{"uORF_down_mORF_down"}
#'     \item \code{"uORF_down_mORF_up"}
#'     \item \code{"uORF_flat_mORF_up"}
#'     \item \code{"uORF_flat_mORF_down"}
#'     \item \code{"dORF_up_mORF_down"}
#'     \item \code{"dORF_up_mORF_flat"}
#'     \item \code{"dORF_down_mORF_flat"}
#'     \item \code{"dORF_down_mORF_down"}
#'     \item \code{"dORF_down_mORF_up"}
#'     \item \code{"dORF_flat_mORF_up"}
#'     \item \code{"dORF_flat_mORF_down"}
#'   }
#' @param gcoeff Numeric. The log-fold change magnitude to apply to regulated ORFs.
#' @param shape Numeric. Shape parameter for the gamma distribution used to simulate baseline coefficients.
#' @param scale Numeric. Scale parameter for the gamma distribution used to simulate baseline coefficients.
#' @param seed Optional integer. If provided, sets the random seed for reproducibility.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{gcoeffs_ribo}{Named numeric vector of log-fold change coefficients for ribosome profiling data.}
#'   \item{gcoeffs_rna}{Named numeric vector of baseline log-fold change coefficients for RNA-seq data.}
#'   \item{labels}{Named binary vector indicating truly regulated ORFs (\code{1}) vs. non-regulated (\code{0}).}
#' }
#'
#' @examples
#' # Example usage:
#' # generate_coefficients(orfs_df, scenario = "uORF_up_mORF_down", seed = 123)
#'
#' @export
generate_coefficients <- function(orfs, 
                                  scenario = "uORF_up_mORF_down", 
                                  gcoeff = 1.5, 
                                  shape = 0.6, 
                                  scale = 0.5, 
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  df <- orfs
  
  # Identify multi-cistronic genes
  gene_counts <- table(df$groupID)
  multi_cistronic_genes <- names(gene_counts[gene_counts > 1])
  
  # Define the log-fold change for each ORF type per scenario
  scenario_lfc <- list(
    "uORF_up_mORF_up" = c(uORF = gcoeff, mORF = gcoeff, dORF = 0),
    "uORF_up_mORF_down" = c(uORF = gcoeff, mORF = -gcoeff, dORF = 0),
    "uORF_up_mORF_flat" = c(uORF = gcoeff, mORF = 0, dORF = 0),
    "uORF_down_mORF_flat" = c(uORF = -gcoeff, mORF = 0, dORF = 0),
    "uORF_down_mORF_down" = c(uORF = -gcoeff, mORF = -gcoeff, dORF = 0),
    "uORF_down_mORF_up" = c(uORF = -gcoeff, mORF = gcoeff, dORF = 0),
    "uORF_flat_mORF_up" = c(uORF = 0, mORF = gcoeff, dORF = 0),
    "uORF_flat_mORF_down" = c(uORF = 0, mORF = -gcoeff, dORF = 0),
    "dORF_up_mORF_up" = c(uORF = 0, mORF = gcoeff, dORF = gcoeff),
    "dORF_up_mORF_down" = c(uORF = 0, mORF = -gcoeff, dORF = gcoeff),
    "dORF_up_mORF_flat" = c(uORF = 0, mORF = 0, dORF = gcoeff),
    "dORF_down_mORF_flat" = c(uORF = 0, mORF = 0, dORF = -gcoeff),
    "dORF_down_mORF_down" = c(uORF = 0, mORF = -gcoeff, dORF = -gcoeff),
    "dORF_down_mORF_up" = c(uORF = 0, mORF = gcoeff, dORF = -gcoeff),
    "dORF_flat_mORF_up" = c(uORF = 0, mORF = gcoeff, dORF = 0),
    "dORF_flat_mORF_down" = c(uORF = 0, mORF = -gcoeff, dORF = 0)
  )
  
  if (!scenario %in% names(scenario_lfc)) {
    stop("Invalid scenario specified. Please choose from available scenarios.")
  }
  
  # Generate a baseline of coefficients for all ORFs
  gcoeffs <- rgamma(nrow(df), shape = shape, scale = scale) * sample(c(-1, 1), nrow(df), replace = TRUE)
  names(gcoeffs) <- rownames(df)
  gcoeffs_ribo <- gcoeffs
  
  # Initialize labels and a temporary list to track regulated ORFs
  labels <- rep(0, nrow(df))
  names(labels) <- rownames(df)
  regulated_orfs <- c()
  
  # Apply effects based on the scenario rules
  current_scenario_lfc <- scenario_lfc[[scenario]]
  for (orf_type in names(current_scenario_lfc)) {
    lfc <- current_scenario_lfc[orf_type]
    
    if (lfc != 0) {
      orf_ids <- rownames(df)[df$groupID %in% multi_cistronic_genes & df$labels == orf_type]
      gcoeffs_ribo[orf_ids] <- gcoeffs_ribo[orf_ids] + lfc
      regulated_orfs <- c(regulated_orfs, orf_ids)
    }
  }
  
  labels[regulated_orfs] <- 1
  
  return(list(
    gcoeffs_ribo = gcoeffs_ribo,
    gcoeffs_rna = gcoeffs,
    labels = labels
  ))
}


#' Simulate Differential ORF Translation (DOT)
#'
#' @description
#' Simulates ribosome profiling and matched RNA-seq count matrices with specified
#' differential ORF translation (DOT) effects. The simulation can include batch effects
#' and supports multiple experimental conditions and replicates.
#'
#' @param ribo A matrix or data frame of ribosome profiling counts (genes x samples).
#' @param rna A matrix or data frame of RNA-seq counts (genes x samples).
#' @param regulation_type Character. Specifies the type of DOT effect to simulate.
#'   Passed to the `scenario` argument of `generate_coefficients`.
#' @param orfs A data frame of ORF annotations. Must include `groupID` and `labels` columns.
#' @param te_genes Numeric. Percentage of genes to be assigned as differentially translated (default: 10).
#' @param bgenes Numeric. Percentage of genes to carry a batch effect (default: 10).
#' @param num_samples Integer. Number of biological replicates per condition (default: 2).
#' @param conditions Integer. Number of experimental conditions (default: 2).
#' @param gcoeff Numeric. Magnitude of log-fold change for DOT effects (default: 1.5).
#' @param bcoeff Numeric. Magnitude of batch effect coefficient (default: 0.9).
#' @param num_batches Integer. Number of batches (default: 2).
#' @param shape Numeric. Shape parameter for gamma distribution used to simulate baseline coefficients (default: 0.6).
#' @param scale Numeric. Scale parameter for gamma distribution used to simulate baseline coefficients (default: 0.5).
#' @param batch_scenario Character. Specifies the batch effect design. Must be one of:
#'   \itemize{
#'     \item \code{"balanced"}
#'     \item \code{"confounded"}
#'     \item \code{"random"}
#'     \item \code{"unbalanced"}
#'     \item \code{"nested"}
#'     \item \code{"modality_specific"}
#'   }
#' @param diagplot_ribo Logical. If \code{TRUE}, generate diagnostic plots for ribo data (default: \code{FALSE}).
#' @param diagplot_rna Logical. If \code{TRUE}, generate diagnostic plots for RNA data (default: \code{FALSE}).
#' @param seed Optional integer. Sets the random seed for reproducibility (default: \code{NULL}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{simData}{Combined simulated count matrix (genes x samples).}
#'   \item{colData}{Data frame containing sample-level metadata (condition, replicate, strategy, batch).}
#'   \item{labels}{Named binary vector indicating true positive (1) and negative (0) ORFs for DOT.}
#'   \item{logFC}{Named vector of true log-fold changes for the simulated DOT effect.}
#' }
#'
#' @import DESeq2
#' @export
simDOT <- function(
    ribo,
    rna,
    regulation_type = NULL,
    orfs = NULL,
    te_genes = 10,
    bgenes = 10,
    num_samples = 2,
    conditions = 2,
    gcoeff = 1.5,
    bcoeff = 0.9,
    num_batches = 2,
    shape = 0.6, 
    scale = 0.5,
    batch_scenario = "balanced",
    diagplot_ribo = FALSE,
    diagplot_rna = FALSE,
    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Ensure all inputs are matrices and have consistent gene names
  ribo <- as.matrix(ribo)
  rna <- as.matrix(rna)
  
  # Find genes that exist in ALL three input objects
  original_genes <- intersect(rownames(ribo), rownames(rna))
  if (!is.null(regulation_type) & !is.null(orfs)) {
    common_genes_all <- intersect(original_genes, rownames(orfs))
  } else {
    common_genes_all <- original_genes
  }
  if(length(common_genes_all) == 0) {
    stop("No common genes found across all three input datasets (ribo, rna, orfs).")
  }
  
  # Filter all inputs to this common set of genes, keeping them as matrices
  ribo_subset <- ribo[common_genes_all, , drop = FALSE]
  rna_subset <- rna[common_genes_all, , drop = FALSE]
  if (!is.null(regulation_type) & !is.null(orfs)) {
    message(" - Simulate differential ORF usage")
    orfs_filtered <- orfs[common_genes_all, , drop = FALSE]
  } else {
    message(" - Simulate differential translation efficiency")
  }
  
  # Perform the filtering based on count thresholds
  ribo_pass <- apply(ribo_subset, 1, function(x) length(x[x > 5]) >= 2)
  rna_pass <- apply(rna_subset, 1, function(x) length(x[x > 5]) >= 2)
  
  # Find the names of genes that passed the filter
  common_filtered_genes <- intersect(names(which(ribo_pass)), names(which(rna_pass)))
  
  if(length(common_filtered_genes) == 0) {
    stop("No genes passed the filtering criteria in both ribo and rna datasets.")
  }
  
  # Final alignment of all data to the filtered gene list
  counts_ribo_filtered <- ribo_subset[common_filtered_genes, , drop = FALSE]
  counts_rna_filtered <- rna_subset[common_filtered_genes, , drop = FALSE]
  if (!is.null(regulation_type) & !is.null(orfs)) {
    orfs_filtered <- orfs_filtered[common_filtered_genes, , drop = FALSE]
  } else {
  }
  
  total_samples <- num_samples * conditions * num_batches
  group <- rep(rep(seq(0, conditions - 1), each = num_samples), num_batches)
  
  # Explicitly convert group to factor for correct model.matrix behavior
  group <- as.factor(group)
  
  # Batch assignment logic
  if (batch_scenario == "balanced") {
    message(" - Use batch_scenario: ", batch_scenario)
    batch <- rep(seq(1, num_batches), each = num_samples * conditions)
    batch <- as.factor(batch)
    if (num_batches == 1) {
      mod <- model.matrix(~ -1 + group)
    } else {
      mod <- model.matrix(~ -1 + batch + group)
    }
    
  } else if (batch_scenario == "confounded") {
    message(" - Use batch_scenario: ", batch_scenario)
    if (num_batches > 1) {
      # In a confounded design, batch and group are the same factor
      batch <- group
      mod <- model.matrix(~ -1 + group)
    } else {
      stop("Use num_batches > 1 for confounded to generate simulated data")
    }
    
  } else if (batch_scenario == "random") {
    if (num_batches > 1) {
      message(" - Use batch_scenario: ", batch_scenario)
      batch <- sample(seq(1, num_batches), total_samples, replace = TRUE)
      batch <- as.factor(batch)
      mod <- model.matrix(~ -1 + batch + group)
    } else {
      stop("Use num_batches > 1 for random to generate simulated data")
    }
    
  } else if (batch_scenario == "unbalanced") {
    if (num_batches > 1) {
      message(" - Use batch_scenario: ", batch_scenario)
      major <- rep(1, round(0.6 * total_samples))
      minor <- sample(seq(2, num_batches), total_samples - length(major), replace = TRUE)
      batch <- c(major, minor)
      batch <- as.factor(batch)
      mod <- model.matrix(~ -1 + batch + group)
    } else {
      stop("Use num_batches > 1 for unbalanced to generate simulated data")
    }
    
  } else if (batch_scenario == "nested") {
    message(" - Use batch_scenario: ", batch_scenario)
    # Create a single factor representing the nested structure
    nested_factor <- as.factor(paste0("C", rep(seq(0, conditions - 1), each = num_samples * num_batches),
                                      "_B", rep(rep(seq(1, num_batches), each = num_samples), conditions)))
    mod <- model.matrix(~ -1 + nested_factor)
    batch <- nested_factor
    # Group and replicate info are for colData only in this scenario
    group <- as.factor(rep(rep(seq(0, conditions - 1), each = num_samples), num_batches))
    
  } else if (batch_scenario == "modality_specific") {
    if (num_batches > 1) {
      message(" - Use batch_scenario: ", batch_scenario)
      # Create a consistent batch vector for both modalities
      batch <- rep(seq_len(num_batches), each = num_samples * conditions)
      batch <- as.factor(batch)
      mod_ribo <- model.matrix(~ -1 + batch + group)
      mod_rna <- model.matrix(~ -1 + group)
    } else {
      stop("num_batches must be > 1 for modality_specific scenario.")
    }
  }
  
  params_ribo <- polyester::get_params(counts_ribo_filtered)
  params_rna <- polyester::get_params(counts_rna_filtered)
  
  if (!is.null(regulation_type) & !is.null(orfs)) {
    coeffs_list <- generate_coefficients(
      orfs = orfs_filtered,
      scenario = regulation_type,
      gcoeff = gcoeff,
      shape = shape, 
      scale = scale, 
      seed = seed
    )
    gcoeffs_ribo <- coeffs_list$gcoeffs_ribo
    gcoeffs_rna <- coeffs_list$gcoeffs_rna
    labels <- coeffs_list$labels
    
    bgenes <- round(bgenes * nrow(counts_ribo_filtered) / 100)
    bselect <- if(bgenes > 0) sample(1:nrow(counts_ribo_filtered), bgenes) else c()
    bcoeffs <- rep(0, nrow(counts_ribo_filtered))
    if(length(bselect) > 0) bcoeffs[bselect] <- bcoeff * sample(c(1, -1), bgenes, replace = TRUE)
    
  } else {
    dTE_genes <- round(te_genes * nrow(counts_ribo_filtered) / 100)
    bgenes <- round(bgenes * nrow(counts_ribo_filtered) / 100)
    
    gcoeffs_rna <- rgamma(nrow(counts_ribo_filtered), shape = shape, scale = scale) * sample(c(-1, 1), nrow(counts_ribo_filtered), replace = TRUE)
    select <- sample(1:nrow(counts_ribo_filtered), dTE_genes)
    
    bcoeffs <- rep(0, nrow(counts_ribo_filtered))
    bselect <- if (bgenes <= length(select)) sample(select, bgenes) else c(select, sample(setdiff(1:nrow(counts_ribo_filtered), select), bgenes - length(select)))
    bcoeffs[bselect] <- bcoeff * sample(c(1, -1), bgenes, replace = TRUE)
    
    gcoeffs_ribo <- gcoeffs_rna
    gcoeffs_ribo[select] <- gcoeffs_rna[select] + sample(c(gcoeff, -gcoeff), length(select), replace = TRUE)
    
    labels <- rep(0, nrow(counts_ribo_filtered))
    names(labels) <- rownames(counts_ribo_filtered)
    labels[select] <- 1
  }
  
  if (batch_scenario == "confounded" && length(bselect) > 0) {
    # The batch effect and condition effect are the same.
    gcoeffs_ribo[bselect] <- gcoeffs_ribo[bselect] + bcoeffs[bselect]
    gcoeffs_rna[bselect] <- gcoeffs_rna[bselect] + bcoeffs[bselect]
  }
  
  # Build the beta matrices using the list output and bcoeffs
  if (batch_scenario == "modality_specific") {
    # Create beta matrices with correct dimensions for each modality
    coeffs_ribo <- matrix(0, nrow = nrow(counts_ribo_filtered), ncol = ncol(mod_ribo))
    colnames(coeffs_ribo) <- colnames(mod_ribo)
    batch_cols <- grep("^batch", colnames(mod_ribo))
    group_cols <- grep("^group", colnames(mod_ribo))
    
    # Apply batch coefficients only to ribo data
    if (length(bselect) > 0) {
      for (i in seq_along(batch_cols)) {
        coeffs_ribo[bselect, batch_cols[i]] <- bcoeffs[bselect] * runif(length(bselect), 0.8, 1.2)
      }
    }
    coeffs_ribo[, group_cols] <- gcoeffs_ribo
    
    # RNA-seq data has no batch effect
    coeffs_rna <- matrix(0, nrow = nrow(counts_rna_filtered), ncol = ncol(mod_rna))
    colnames(coeffs_rna) <- colnames(mod_rna)
    group_cols <- grep("^group", colnames(mod_rna))
    coeffs_rna[, group_cols] <- gcoeffs_rna
    
  } else if (batch_scenario == "nested") {
    # Assign coefficients to the nested factor columns directly
    coeffs_ribo <- matrix(0, nrow = nrow(counts_ribo_filtered), ncol = ncol(mod))
    colnames(coeffs_ribo) <- colnames(mod)
    coeffs_rna <- matrix(0, nrow = nrow(counts_rna_filtered), ncol = ncol(mod))
    colnames(coeffs_rna) <- colnames(mod)
    
    for (col in colnames(mod)) {
      coeffs_ribo[, col] <- gcoeffs_ribo
      coeffs_rna[, col] <- gcoeffs_rna
    }
    # Note: The batch effect is integrated into the group effect in this scenario.
    # gcoeffs already contain the combined effect.
    
  } else {
    # This covers balanced, random, unbalanced, and confounded (which now only has `group`)
    coeffs_ribo <- matrix(0, nrow = nrow(counts_ribo_filtered), ncol = ncol(mod))
    colnames(coeffs_ribo) <- colnames(mod)
    
    coeffs_rna <- matrix(0, nrow = nrow(counts_rna_filtered), ncol = ncol(mod))
    colnames(coeffs_rna) <- colnames(mod)
    
    # Assign batch coefficients
    batch_cols <- grep("^batch", colnames(mod))
    if (length(bselect) > 0 && length(batch_cols) > 0) {
      # Assign coefficients per batch column
      for (i in seq_along(batch_cols)) {
        coeffs_ribo[bselect, batch_cols[i]] <- bcoeffs[bselect] * runif(length(bselect), 0.8, 1.2)
        coeffs_rna[bselect, batch_cols[i]] <- bcoeffs[bselect] * runif(length(bselect), 0.8, 1.2)
      }
    }
    
    # Assign group coefficients (DOT effect)
    group_cols <- grep("^group", colnames(mod))
    if (length(group_cols) > 0) {
      coeffs_ribo[, group_cols] <- gcoeffs_ribo
      coeffs_rna[, group_cols] <- gcoeffs_rna
    }
  }
  
  if (batch_scenario == "modality_specific") {
    sim_ribo_filtered <- polyester::create_read_numbers(params_ribo$mu, params_ribo$fit, params_ribo$p0, beta = coeffs_ribo, mod = mod_ribo, seed = seed)
    sim_rna_filtered <- polyester::create_read_numbers(params_rna$mu, params_rna$fit, params_rna$p0, beta = coeffs_rna, mod = mod_rna, seed = seed)
  } else {
    sim_ribo_filtered <- polyester::create_read_numbers(params_ribo$mu, params_ribo$fit, params_ribo$p0, beta = coeffs_ribo, mod = mod, seed = seed)
    sim_rna_filtered <- polyester::create_read_numbers(params_rna$mu, params_rna$fit, params_rna$p0, beta = coeffs_rna, mod = mod, seed = seed)
  }
  
  final_cols_ribo <- paste0("sample", 1:total_samples, ".condition", group, ".batch", batch, ".ribo")
  final_cols_rna <- paste0("sample", 1:total_samples, ".condition", group, ".batch", batch, ".rna")
  
  sim_ribo_full <- matrix(0, nrow = length(original_genes), ncol = total_samples, dimnames = list(original_genes, final_cols_ribo))
  sim_rna_full <- matrix(0, nrow = length(original_genes), ncol = total_samples, dimnames = list(original_genes, final_cols_rna))
  
  sim_ribo_full[common_filtered_genes, ] <- sim_ribo_filtered
  sim_rna_full[common_filtered_genes, ] <- sim_rna_filtered
  
  merged <- cbind(sim_ribo_full, sim_rna_full)
  replicate <- rep(rep(seq(1, num_samples), conditions), num_batches)
  
  # Ensure batch vector is correct for colData
  if (batch_scenario == "modality_specific") {
    batch_rna <- rep("none", total_samples)
    batch_ribo <- as.character(batch)
    batch_full <- c(batch_ribo, batch_rna)
  } else {
    batch_full <- as.character(rep(batch, times = 2))
  }
  
  coldata <- data.frame(
    run = colnames(merged),
    condition = factor(rep(group, 2)),
    replicate = factor(rep(replicate, 2)),
    strategy = factor(rep(c("ribo", "rna"), each = total_samples)),
    batch = factor(batch_full)
  )
  
  final_labels <- rep(0, length(original_genes))
  names(final_labels) <- original_genes
  if (!is.null(regulation_type)) {
    final_labels[common_filtered_genes] <- labels
  } else { # DTE
    final_labels[names(labels)] <- labels
  }
  
  final_change <- rep(0, length(original_genes))
  names(final_change) <- original_genes
  final_change[names(gcoeffs_ribo)] <- gcoeffs_ribo - gcoeffs_rna
  
  
  if (isTRUE(diagplot_ribo) | isTRUE(diagplot_rna)) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(merged),
      colData = coldata,
      design = ~ strategy
    )
    
    dds_ribo <- dds[, dds$strategy == "ribo"]
    dds_rna <- dds[, dds$strategy == "rna"]
    
    design(dds_ribo) <- ~ condition + batch
    design(dds_rna) <- ~ condition + batch
    
    vst_ribo <- DESeq2::vst(dds_ribo, blind = FALSE)
    pca_ribo <- prcomp(t(assay(vst_ribo)))
    
    vst_rna <- DESeq2::vst(dds_rna, blind = FALSE)
    pca_rna <- prcomp(t(assay(vst_rna)))
  }
  
  if (isTRUE(diagplot_ribo)) {
    tryCatch({
      # Ensure a graphics device is open
      if (dev.cur() == 1) dev.new()
      
      # Get current margin and device size
      current_mar <- par("mar")
      dev_dims <- dev.size("in")  # width, height in inches
      
      if (dev_dims[1] > 7) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 4))  # enough space for outside legend
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else if (dev_dims[1] > 5.5) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 2))  # medium space
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else if (dev_dims[1] > 4) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 0.5)) # small margins
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else {
        message("Legend may overlap with data. Try enlarging the plot area.")
        par(xpd = FALSE, mar = c(4, 4, 2, 2))  # tight margins, legend inside
        legend_inset <- 0
        legend_outside <- FALSE
      }
      
      # PCA data
      percentVar_ribo <- round(100 * pca_ribo$sdev^2 / sum(pca_ribo$sdev^2))
      colors_ribo <- as.numeric(as.factor(colData(vst_ribo)$condition))
      shapes_ribo <- as.numeric(as.factor(colData(vst_ribo)$batch))
      
      # PCA plot
      plot(
        x = pca_ribo$x[, 1],
        y = pca_ribo$x[, 2],
        col = colors_ribo,
        pch = shapes_ribo,
        xlab = paste0("PC1: ", percentVar_ribo[1], "% variance"),
        ylab = paste0("PC2: ", percentVar_ribo[2], "% variance"),
        main = paste0("PCA of Ribo-seq (gcoeff: ", gcoeff, ")\n", batch_scenario,
                      " (num_batches: ", num_batches, ", bcoeff: ", paste(bcoeff, collapse = ", "), ")")
      )
      
      # Condition legend
      legend(
        "topleft",
        inset = legend_inset,
        xpd = legend_outside,
        bty = "n",
        legend = unique(as.data.frame.array(colData(vst_ribo))$condition),
        col = unique(colors_ribo),
        pch = 16,
        title = "Condition"
      )
      
      # Batch legend
      legend(
        "bottomleft",
        inset = legend_inset,
        xpd = legend_outside,
        bty = "n",
        legend = unique(as.data.frame.array(colData(vst_ribo))$batch),
        pch = unique(shapes_ribo),
        col = "black",
        title = "Batch"
      )
      
      # Reset margins
      par(xpd = FALSE, mar = c(5, 4, 4, 2) + 0.1)
    }, error = function(e) {
      message("Skipping Ribo-seq PCA plot due to error: ", e$message)
    })
  }
  
  if (isTRUE(diagplot_rna)) {
    tryCatch({
      # Ensure a graphics device is open
      if (dev.cur() == 1) dev.new()
      
      # Get current margin and device size
      current_mar <- par("mar")
      dev_dims <- dev.size("in")  # width, height in inches
      
      # Adjust margins based on device width
      if (dev_dims[1] > 7) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 4))  # enough space for outside legend
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else if (dev_dims[1] > 5.5) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 2))  # medium space
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else if (dev_dims[1] > 4) {
        par(xpd = TRUE, mar = current_mar + c(0, 0, 0, 0.5)) # small margins
        legend_inset <- c(1.02, 0)
        legend_outside <- TRUE
      } else {
        message("Legend may overlap with data. Try enlarging the plot area.")
        par(xpd = FALSE, mar = c(4, 4, 2, 2))  # tight margins, legend inside
        legend_inset <- 0
        legend_outside <- FALSE
      }
      
      # PCA plot
      percentVar_rna <- round(100 * pca_rna$sdev^2 / sum(pca_rna$sdev^2))
      colors_rna <- as.numeric(as.factor(colData(vst_rna)$condition))
      shapes_rna <- as.numeric(as.factor(colData(vst_rna)$batch))
      
      plot(
        x = pca_rna$x[, 1],
        y = pca_rna$x[, 2],
        col = colors_rna,
        pch = shapes_rna,
        xlab = paste0("PC1: ", percentVar_rna[1], "% variance"),
        ylab = paste0("PC2: ", percentVar_rna[2], "% variance"),
        main = paste0("PCA of RNA-seq (gcoeff: ", gcoeff, ")\n", batch_scenario,
                      " (num_batches: ", num_batches, ", bcoeff: ", paste(bcoeff, collapse = ", "), ")")
      )
      
      # Condition legend
      legend(
        "topleft",
        inset = legend_inset,
        xpd = legend_outside,
        bty = "n",
        legend = unique(as.data.frame.array(colData(vst_rna))$condition),
        col = unique(colors_rna),
        pch = 16,
        title = "Condition"
      )
      
      # Batch legend
      legend(
        "bottomleft",
        inset = legend_inset,
        xpd = legend_outside,
        bty = "n",
        legend = unique(as.data.frame.array(colData(vst_rna))$batch),
        pch = unique(shapes_rna),
        col = "black",
        title = "Batch"
      )
      
      # Reset margins
      par(xpd = FALSE, mar = c(5, 4, 4, 2) + 0.1)
    }, error = function(e) {
      message("Skipping RNA-seq PCA plot due to error: ", e$message)
    })
  }
  
  
  return(list(simData = merged, colData = coldata, labels = final_labels, logFC = final_change))
}
