#' Simulate Differentially Translated ORFs
#'
#' This function simulates matched ribosome profiling and RNA-seq count matrices.
#'
#' @param ribo A matrix or data frame of ribosome profiling counts (genes x samples).
#' @param rna A matrix or data frame of RNA-seq counts (genes x samples).
#' @param te_genes Numeric (percentage). Proportion of genes to be assigned as differential translation efficiency (default: 10).
#' @param bgenes Numeric (percentage). Proportion of genes to carry a batch effect (default: 10).
#' @param num_samples Integer. Number of biological replicates per condition (default: 4).
#' @param conditions Integer. Number of experimental conditions (default: 2).
#' @param bcoeff Numeric. Magnitude of batch effect coefficient (default: 0.9).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{simData} - A combined count matrix (genes x samples), 
#'     including both Ribo-seq and RNA-seq simulations.
#'   \item \code{colData} - A data frame of sample-level metadata, including 
#'     replicate, strategy ("ribo" or "rna"), and batch assignment.
#'  \item \code{labels} - A numeric vector indicating true positives and true negatives.
#' }
#'
#' @examples
#' # Simulate DOT with default settings
#' simRes <- simDOT(ribo_counts, rna_counts)
#' str(simRes$simData)
#' head(simRes$colData)
#'
#'
simDOT <- function(
    ribo,
    rna,
    te_genes = 10,
    bgenes = 10,
    num_samples = 4,
    conditions = 2,
    bcoeff = 0.9,
    num_batches = 1) {
  
  if (!(is.matrix(ribo) || is.data.frame(ribo)) || !(is.matrix(rna) || is.data.frame(rna))) {
    stop("Both 'ribo' and 'rna' must be matrices or data frames.")
  }
  
  # Coerce to matrix if needed
  ribo <- as.matrix(ribo)
  rna  <- as.matrix(rna)
  
  # Check for matching number of samples
  if (ncol(ribo) != ncol(rna)) {
    stop("Ribo and RNA matrices must have the same number of columns (samples).")
  }
  
  # Store original gene order and sample names for later use
  original_genes <- intersect(rownames(ribo), rownames(rna))
  ribo <- ribo[original_genes, , drop = FALSE]
  rna <- rna[original_genes, , drop = FALSE]
  
  # Filtering genes with sufficient counts
  ribo_filter <- apply(ribo, 1, function(x) length(x[x > 5]) >= 2)
  rna_filter <- apply(rna, 1, function(x) length(x[x > 5]) >= 2)
  common_filtered_genes <- intersect(names(ribo_filter[ribo_filter]), names(rna_filter[rna_filter]))
  
  counts_ribo_filtered <- as.matrix(ribo[common_filtered_genes, , drop = FALSE])
  counts_rna_filtered  <- as.matrix(rna[common_filtered_genes, , drop = FALSE])
  
  # Use the filtered data to get polyester's parameters.
  # Use a dummy matrix later on with the final number of samples as a template for simulation.
  
  # Create expanded design variables
  final_total_samples <- num_samples * conditions * num_batches
  group <- rep(rep(seq(0, conditions - 1), each = num_samples), num_batches)
  batch <- rep(seq(1, num_batches), each = num_samples * conditions)
  
  # Create the expanded model matrix
  mod <- model.matrix(~ -1 + batch + group)
  
  # Get polyester parameters from the original filtered data
  params_ribo <- get_params(counts_ribo_filtered)
  params_rna  <- get_params(counts_rna_filtered)
  
  # Define and select ORFs and batch-affected genes from the filtered set
  dTE_genes <- round(te_genes * nrow(counts_ribo_filtered) / 100, 0)
  bgenes <- round(bgenes * nrow(counts_ribo_filtered) / 100, 0)
  
  gcoeffs <- rgamma(nrow(counts_ribo_filtered), shape = 0.6, scale = 0.5) * sample(c(-1, 1), nrow(counts_ribo_filtered), replace = TRUE)
  select <- sample(1:nrow(counts_ribo_filtered), dTE_genes)
  
  bcoeffs <- rep(0, nrow(counts_ribo_filtered))
  if (bgenes <= length(select)) {
    bselect <- sample(select, bgenes)
  } else {
    bselect <- c(select, sample(setdiff(1:nrow(counts_ribo_filtered), select), bgenes - length(select)))
  }
  bcoeffs[bselect] <- bcoeff * sample(c(1, -1), bgenes, replace = TRUE)
  
  coeffs_ribo <- cbind(bcoeffs, gcoeffs)
  
  gcoeffs_rna <- gcoeffs
  gcoeffs_rna[select] <- gcoeffs[select] + sample(c(1.5, -1.5), length(select), replace = TRUE)
  coeffs_rna <- cbind(bcoeffs, gcoeffs_rna)
  
  # Run the simulation for the expanded number of samples
  sim_ribo_filtered <- create_read_numbers(params_ribo$mu, params_ribo$fit, params_ribo$p0,
                                           beta = coeffs_ribo, mod = mod, seed = 32332)
  sim_rna_filtered <- create_read_numbers(params_rna$mu, params_rna$fit, params_rna$p0,
                                          beta = coeffs_rna, mod = mod, seed = 32332)
  
  # Create final matrices with all original genes and new sample columns
  final_cols_ribo <- paste0("sample", 1:final_total_samples, ".condition", group, ".batch", batch, ".ribo")
  final_cols_rna <- paste0("sample", 1:final_total_samples, ".condition", group, ".batch", batch, ".rna")
  
  sim_ribo_full <- matrix(0, nrow = length(original_genes), ncol = final_total_samples,
                          dimnames = list(original_genes, final_cols_ribo))
  sim_rna_full  <- matrix(0, nrow = length(original_genes), ncol = final_total_samples,
                          dimnames = list(original_genes, final_cols_rna))
  
  # Overwrite simulated values for the filtered genes
  sim_ribo_full[common_filtered_genes, ] <- sim_ribo_filtered
  sim_rna_full[common_filtered_genes, ] <- sim_rna_filtered
  
  # Create final merged matrix and coldata with a replicate ID
  merged <- cbind(sim_ribo_full, sim_rna_full)
  
  # The replicate ID must be unique for each sample
  replicate <- rep(rep(seq(1, num_samples), conditions), num_batches)
  
  coldata <- data.frame(
    condition = rep(group, 2),
    replicate = rep(replicate, 2),
    strategy = rep(c("ribo", "rna"), each = final_total_samples),
    batch = rep(batch, 2),
    row.names = colnames(merged)
  )
  
  labels <- rep(0, length(original_genes))
  names(labels) <- original_genes
  labels[common_filtered_genes[select]] <- 1
  
  return(list(simData = merged, colData = coldata, labels = labels))
}

