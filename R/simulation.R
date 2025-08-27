#' Simulate Differential Translation Efficiency
#'
#' This function simulates ribosome profiling and matched RNA-seq count matrices.
#'
#' @param ribo A matrix or data frame of ribosome profiling counts (genes x samples).
#' @param rna A matrix or data frame of RNA-seq counts (genes x samples).
#' @param te_genes Numeric (percentage). Proportion of genes to be assigned as differential translation efficiency (default: 10).
#' @param bgenes Numeric (percentage). Proportion of genes to carry a batch effect (default: 10).
#' @param num_samples Integer. Number of biological replicates per condition (default: 3).
#' @param conditions Integer. Number of experimental conditions (default: 2).
#' @param bcoeff Numeric. Magnitude of batch effect coefficient (default: 0.9).
#' @param num_batches Integer. Number of batches (default: 1).
#' @param batchScenario Character. One of "balanced", "confounded", "random", "unbalanced", "nested", "modalitySpecific".
#' @param seed An optional integer. Use to set the seed for reproducible results (default: NULL).
#'
#' @return A list with:
#' \item{simData}{Combined count matrix (genes x samples)}
#' \item{colData}{Sample-level metadata}
#' \item{labels}{Vector indicating true positives (1) and negatives (0)}
#' \item{log2FC}{Vector for log2 fold-change}
#'
simDTE <- function(
    ribo,
    rna,
    te_genes = 10,
    bgenes = 10,
    num_samples = 3,
    conditions = 2,
    bcoeff = 0.9,
    num_batches = 1,
    batchScenario = "balanced",
    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ribo <- as.matrix(ribo)
  rna <- as.matrix(rna)
  original_genes <- intersect(rownames(ribo), rownames(rna))
  ribo <- ribo[original_genes, , drop = FALSE]
  rna <- rna[original_genes, , drop = FALSE]
  
  ribo_filter <- apply(ribo, 1, function(x) length(x[x > 5]) >= 2)
  rna_filter <- apply(rna, 1, function(x) length(x[x > 5]) >= 2)
  common_filtered_genes <- intersect(names(ribo_filter[ribo_filter]), names(rna_filter[rna_filter]))
  
  counts_ribo_filtered <- ribo[common_filtered_genes, , drop = FALSE]
  counts_rna_filtered <- rna[common_filtered_genes, , drop = FALSE]
  
  total_samples <- num_samples * conditions * num_batches
  group <- rep(rep(seq(0, conditions - 1), each = num_samples), num_batches)
  
  # Batch assignment logic
  if (batchScenario == "balanced") {
    batch <- rep(seq(1, num_batches), each = num_samples * conditions)
    mod <- model.matrix(~ -1 + batch + group)
    
  } else if (batchScenario == "confounded") {
    if (num_batches>1) {
      batch <- rep(rep(seq(1, conditions), each = num_samples), num_batches)
      mod <- model.matrix(~ -1 + batch + group)
    } else {
      stop("Use num_batches >1 for confounded to generate simulated data")
    }
    
  } else if (batchScenario == "random") {
    if (num_batches>1) {
      batch <- sample(seq(1, num_batches), total_samples, replace = TRUE)
      mod <- model.matrix(~ -1 + batch + group)
    } else {
      stop("Use num_batches >1 for unbalanced to generate simulated data")
    }
    
  } else if (batchScenario == "unbalanced") {
    if (num_batches>1) {
      major <- rep(1, round(0.6 * total_samples))
      minor <- sample(seq(2, num_batches), total_samples - length(major), replace = TRUE)
      batch <- c(major, minor)
      mod <- model.matrix(~ -1 + batch + group)
    } else {
      stop("Use num_batches >1 for unbalanced to generate simulated data")
    }
    
  } else if (batchScenario == "nested") {
    batch_factor <- paste0("C", rep(seq(0, conditions - 1), each = num_samples * num_batches),
                           "_B", rep(rep(seq(1, num_batches), each = num_samples), conditions))
    
    # Create the model matrix using only the nested factor
    mod <- model.matrix(~ -1 + as.factor(batch_factor))
    batch <- batch_factor
    
  } else if (batchScenario == "modalitySpecific") {
    batch <- rep(seq(1, num_batches), each = num_samples * conditions)
    mod_rna <- model.matrix(~ -1 + factor(group))  # No batch effect for RNA
    
    if (num_batches > 1) {
      batch_factor <- factor(batch[1:total_samples])
      if (length(levels(batch_factor)) < 2) {
        stop("Batch factor must have at least two levels for modalitySpecific scenario.")
      }
      mod_ribo <- model.matrix(~ -1 + batch_factor + factor(group))
    } else {
      stop("num_batches must be >1 for modalitySpecific scenario.")
    }
  }
  
  params_ribo <- polyester::get_params(counts_ribo_filtered)
  params_rna <- polyester::get_params(counts_rna_filtered)
  
  dTE_genes <- round(te_genes * nrow(counts_ribo_filtered) / 100)
  bgenes <- round(bgenes * nrow(counts_ribo_filtered) / 100)
  
  gcoeffs <- rgamma(nrow(counts_ribo_filtered), shape = 0.6, scale = 0.5) * sample(c(-1, 1), nrow(counts_ribo_filtered), replace = TRUE)
  select <- sample(1:nrow(counts_ribo_filtered), dTE_genes)
  
  bcoeffs <- rep(0, nrow(counts_ribo_filtered))
  bselect <- if (bgenes <= length(select)) sample(select, bgenes) else c(select, sample(setdiff(1:nrow(counts_ribo_filtered), select), bgenes - length(select)))
  bcoeffs[bselect] <- bcoeff * sample(c(1, -1), bgenes, replace = TRUE)
  
  gcoeffs_ribo <- gcoeffs
  gcoeffs_ribo[select] <- gcoeffs[select] + sample(c(1.5, -1.5), length(select), replace = TRUE)
  
  if (batchScenario == "modalitySpecific") {
    # Match coeffs_ribo to mod_ribo
    coeffs_ribo <- matrix(0, nrow = nrow(counts_ribo_filtered), ncol = ncol(mod_ribo))
    colnames(coeffs_ribo) <- colnames(mod_ribo)
    batch_cols <- grep("^as.factor\\(batch\\)", colnames(mod_ribo))
    coeffs_ribo[, batch_cols] <- bcoeffs * runif(length(batch_cols), 0.8, 1.2)
    group_cols <- grep("^as.factor\\(group\\)", colnames(mod_ribo))
    coeffs_ribo[, group_cols] <- gcoeffs_ribo
    
    coeffs_rna <- matrix(0, nrow = nrow(counts_rna_filtered), ncol = ncol(mod_rna))
    colnames(coeffs_rna) <- colnames(mod_rna)
    group_cols <- grep("^as.factor\\(group\\)", colnames(mod_rna))
    coeffs_rna[, group_cols] <- gcoeffs
    
  } else if (batchScenario == "nested") {
    batch <- paste0("C", rep(seq(0, conditions - 1), each = num_samples * num_batches),
                    "_B", rep(rep(seq(1, num_batches), each = num_samples), conditions))
    mod <- model.matrix(~ -1 + as.factor(batch))
    
    # Dynamically build the coeffs matrices
    coeffs_ribo <- matrix(0, nrow = nrow(counts_ribo_filtered), ncol = ncol(mod))
    colnames(coeffs_ribo) <- colnames(mod)
    
    coeffs_rna <- matrix(0, nrow = nrow(counts_rna_filtered), ncol = ncol(mod))
    colnames(coeffs_rna) <- colnames(mod)
    
    # Assign baseline gcoeffs to all columns for RNA-seq
    for (col in colnames(mod)) {
      coeffs_rna[, col] <- gcoeffs
    }
    
    # Assign TE-modified gcoeffs_ribo to Ribo-seq for all columns
    for (col in colnames(mod)) {
      coeffs_ribo[, col] <- gcoeffs_ribo
    } 
    
  } else {
    coeffs_ribo <- cbind(bcoeffs, gcoeffs_ribo)
    coeffs_rna <- cbind(bcoeffs, gcoeffs)
  }
  
  if (batchScenario == "modalitySpecific") {
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
  
  # Construct colData
  # Batch assignment for colData
  if (batchScenario == "modalitySpecific") {
    # RNA samples: no batch effect (can be NA or a constant)
    batch_rna <- rep("none", total_samples)
    
    # Ribo samples: batch effect applied
    batch_ribo <- rep(seq_len(num_batches), each = num_samples * conditions)
    
    # Combine
    batch_full <- c(batch_ribo, batch_rna)
  } else {
    # For all other scenarios, batch should be defined earlier
    if (length(batch) == total_samples) {
      batch_full <- rep(batch, times = 2)
    } else if (length(batch) == 2 * total_samples) {
      batch_full <- batch
    } else {
      stop("Batch vector length mismatch")
    }
  }
  
  coldata <- data.frame(
    run = colnames(merged),
    condition = factor(rep(group, 2)),
    replicate = factor(rep(replicate, 2)),
    strategy = factor(rep(c("ribo", "rna"), each = total_samples)),
    batch = factor(batch_full)
  )
  
  labels <- rep(0, length(original_genes))
  names(labels) <- original_genes
  labels[common_filtered_genes[select]] <- 1
  
  change <- gcoeffs_ribo - gcoeffs
  names(change) <- common_filtered_genes
  
  return(list(simData = merged, colData = coldata, labels = labels, log2FC = change))
}

