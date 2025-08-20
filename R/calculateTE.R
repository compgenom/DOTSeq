#' Calculate Translational Efficiency and Occupancy Shift
#'
#' This function computes translational efficiency (TE) and occupancy shift
#' for a set of normalized RNA-seq and Ribo-seq counts. TE is calculated
#' as the ratio of Ribo-seq counts to RNA-seq counts. Occupancy
#' shift is calculated as the log2 ratio of Ribo-seq proportion to RNA-seq
#' proportion within each gene.
#'
#' @param normCnts A numeric matrix or data frame of normalized counts,
#'   where rows correspond to ORFs and columns correspond to samples. 
#'   RNA-seq and Ribo-seq reads are distinguished by suffixes.
#' @param rnaSuffix Character string indicating the suffix of RNA-seq
#'   samples in the column names. Default is ".rna".
#' @param riboSuffix Character string indicating the suffix of Ribo-seq
#'   samples in the column names. Default is ".ribo".
#' @param pseudoCnt Numeric value added to counts to avoid division by zero.
#'   Default is 1e-6.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{absoluteTE}{A matrix of translational efficiency (Ribo / RNA).}
#'   \item{occupancyShift}{A matrix of log2 ratios of Ribo-proportion / RNA-proportion within each gene.}
#' }
#'
#' @examples
#' result <- calculateTE(normCounts)
#' head(result$absoluteTE)
#' head(result$occupancyShift)
#'
#' @export
calculateTE <- function(normCnts, rnaSuffix = ".rna", riboSuffix = ".ribo", pseudoCnt = 1e-6) {
  # Identify RNA and Ribo columns (match anywhere in name)
  rnaCols <- grep(rnaSuffix, colnames(normCnts), value = TRUE)
  riboCols <- grep(riboSuffix, colnames(normCnts), value = TRUE)
  
  # Extract sample prefixes by removing everything from the suffix onward
  rnaPrefixes <- sub(paste0(rnaSuffix, ".*"), "", rnaCols)
  riboPrefixes <- sub(paste0(riboSuffix, ".*"), "", riboCols)
  
  # Find common sample prefixes
  commonPrefixes <- intersect(rnaPrefixes, riboPrefixes)
  if (length(commonPrefixes) == 0) {
    stop("No matching RNA/Ribo sample prefixes found.")
  }
  
  # Initialize TE matrix
  te <- matrix(NA, nrow = nrow(normCnts), ncol = length(commonPrefixes))
  rownames(te) <- rownames(normCnts)
  colnames(te) <- paste0(commonPrefixes, ".te")
  
  for (prefix in commonPrefixes) {
    rnaCol <- grep(paste0("^", prefix, rnaSuffix), colnames(normCnts), value = TRUE)
    riboCol <- grep(paste0("^", prefix, riboSuffix), colnames(normCnts), value = TRUE)
    if (length(rnaCol) != 1 || length(riboCol) != 1) {
      warning("Ambiguous or missing match for prefix: ", prefix)
      next
    }
    rnaVals <- normCnts[, rnaCol]
    riboVals <- normCnts[, riboCol]
    te[, paste0(prefix, ".te")] <- (riboVals + pseudoCnt) / (rnaVals + pseudoCnt)
  }
  
  # Compute occupancy shift
  gene_ids <- sub(":O.*", "", rownames(normCnts))
  riboMat <- normCnts[, riboCols, drop = FALSE]
  rnaMat <- normCnts[, rnaCols, drop = FALSE]
  riboSum <- rowsum(riboMat, group = gene_ids)
  rnaSum <- rowsum(rnaMat, group = gene_ids)
  riboProportion <- riboMat / riboSum[gene_ids, ]
  rnaProportion <- rnaMat / rnaSum[gene_ids, ]
  occupancyShift <- log2((riboProportion + pseudoCnt) / (rnaProportion + pseudoCnt))
  
  return(list(
    absoluteTE = te,
    occupancyShift = occupancyShift
  ))
}
