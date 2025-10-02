#' Calculate Translational Efficiency and Shifts in ORF Usage 
#'
#' This function computes translational efficiency (TE) and shifts in ORF usage
#' for a set of normalized RNA-seq and Ribo-seq counts. TE is calculated
#' as the ratio of Ribo-seq counts to RNA-seq counts. Shifts in ORF usage
#' is calculated as the log2 ratio of Ribo-seq proportion to RNA-seq
#' proportion within each gene.
#'
#' @param norm_counts A numeric matrix or data frame of normalized counts,
#'   where rows correspond to ORFs and columns correspond to samples. 
#'   RNA-seq and Ribo-seq reads are distinguished by suffixes.
#' @param sample_delim _delimiter used in sample name (default: ".").
#' @param rna_suffix Character string indicating the suffix of RNA-seq
#'   samples in the column names (default is ".rna").
#' @param ribo_suffix Character string indicating the suffix of Ribo-seq
#'   samples in the column names (default is ".ribo").
#' @param pseudocount Numeric value added to counts to avoid division by zero
#'   (default is 1e-6).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{te}{A matrix of translational efficiency (Ribo / RNA).}
#'   \item{usage_shift}{A matrix of log2 ratios of Ribo-proportion / RNA-proportion within each gene.}
#' }
#'
#' @examples
#' result <- calculateTE(norm_counts)
#' head(result$te)
#' head(result$usage_shift)
#'
#' @export
calculateTE <- function(norm_counts, 
                        rna_suffix = ".rna", 
                        ribo_suffix = ".ribo", 
                        sample_delim = ".",
                        pseudocount = 1e-6) {
  
  # Identify RNA and Ribo columns (match anywhere in name)
  rna_cols <- grep(rna_suffix, colnames(norm_counts), value = TRUE)
  ribo_cols <- grep(ribo_suffix, colnames(norm_counts), value = TRUE)
  
  if (!is.null(sample_delim)) {
    # Escape special regex characters
    escaped_delim <- gsub("([\\^$.|?*+(){}])", "\\\\\\1", sample_delim)
    # Construct regex pattern using sample_delim
    pattern <- paste0("^[^", escaped_delim, "]+", escaped_delim)
    
    # Remove sample ID using the constructed pattern
    rna_prefixes <- sub(pattern, "", rna_cols)
    ribo_prefixes <- sub(pattern, "", ribo_cols)
  } else {
    # If no sample_delim, use full column name
    rna_prefixes <- rna_cols
    ribo_prefixes <- ribo_cols
  }
  
  # Extract sample _prefixes by removing everything from the suffix onward
  rna_prefixes <- sub(paste0(rna_suffix, ".*"), "", rna_prefixes)
  ribo_prefixes <- sub(paste0(ribo_suffix, ".*"), "", ribo_prefixes)
  
  # Find common sample _prefixes
  common_prefixes <- intersect(rna_prefixes, ribo_prefixes)

  if (length(common_prefixes) == 0) {
    stop("No matching RNA/Ribo sample _prefixes found.")
  }
  
  # Initialize TE matrix
  te <- matrix(NA, nrow = nrow(norm_counts), ncol = length(common_prefixes))
  rownames(te) <- rownames(norm_counts)
  colnames(te) <- paste0(common_prefixes, ".te")
  
  for (prefix in common_prefixes) {
    rnaCol <- grep(paste0(prefix, rna_suffix), colnames(norm_counts), value = TRUE)
    riboCol <- grep(paste0(prefix, ribo_suffix), colnames(norm_counts), value = TRUE)
    if (length(rnaCol) != 1 || length(riboCol) != 1) {
      warning("Ambiguous or missing match for prefix: ", prefix)
      next
    }
    rnaVals <- norm_counts[, rnaCol]
    riboVals <- norm_counts[, riboCol]
    te[, paste0(prefix, ".te")] <- (riboVals + pseudocount) / (rnaVals + pseudocount)
  }
  
  # Compute usage shift
  gene_ids <- sub(":O.*", "", rownames(norm_counts))
  ribo_mat <- norm_counts[, ribo_cols, drop = FALSE]
  rna_mat <- norm_counts[, rna_cols, drop = FALSE]
  ribo_sum <- rowsum(ribo_mat, group = gene_ids)
  rna_sum <- rowsum(rna_mat, group = gene_ids)
  ribo_prop <- ribo_mat / ribo_sum[gene_ids, ]
  rna_prop <- rna_mat / rna_sum[gene_ids, ]
  usage_shift <- log2((ribo_prop + pseudocount) / (rna_prop + pseudocount))
  
  return(list(
    te = te,
    usage_shift = usage_shift
  ))
}
