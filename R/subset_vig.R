#' Create a Subsetted Count Table for DOTSeq
#'
#' @description
#' Generates a smaller, reproducible subset of ORFs from a large DOTSeq dataset.
#' This function selects a random subset of rows from a given results dataframe, extracts
#' the corresponding ORF IDs, and combines them with the original count table. 
#' It also randomly sets counts of a proportion of missing ORFs to zero to simulate sparsity.
#'
#' @param df A data frame from \code{fitDOT} containing ORF-level results.
#'   Must include a column named \code{Row.names} with ORF identifiers.
#' @param cnt A count table data frame containing all ORFs and their sample counts.
#'   Must include a column named \code{Geneid}.
#' @param seed Integer; random seed for reproducibility. Default is \code{123}.
#' @param subset_fraction Numeric; fraction of ORFs to sample from \code{df}. Default is \code{0.2}.
#' @param zero_fraction Numeric; fraction of missing ORFs to set counts to zero. Default is \code{0.8}.
#'
#' @return A data frame containing the combined counts:
#' \describe{
#'   \item{metadata columns}{Original metadata columns from \code{cnt}, e.g., Geneid, Chr, Start, End, Strand, Length.}
#'   \item{count columns}{Counts for each sample, with some missing ORFs randomly set to zero.}
#' }
#'
#' @details
#' This function is useful for vignette demonstrations or for testing DOTSeq on a smaller
#' subset of ORFs to reduce computation time. It ensures reproducibility by using a fixed random seed.
#'
#' @export
#' 
make_subset <- function(df, cnt, seed = 123) {
  # set reproducibility
  set.seed(seed)  
  
  # pick 20% of df rows
  # n_rows <- nrow(df)
  # sample_rows_20 <- sample(n_rows, size = floor(0.2 * n_rows))
  
  # subset and order
  # df_20 <- df[sample_rows_20, ]
  # df_20 <- df_20[order(df_20$Row.names), ]
  # rownames(df_20) <- NULL
  
  # extract gene IDs
  gene_ids <- unique(sub(":.*", "", df$Row.names))
  gene_ids_df <- data.frame(Geneid = gene_ids, stringsAsFactors = FALSE)
  
  # find missing ORFs
  missing_orfs <- cnt[!cnt$Geneid %in% gene_ids_df$Geneid, ]
  
  # metadata + counts
  metadata_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- grep("^SRR", colnames(missing_orfs), value = TRUE)
  cols_to_keep <- c(metadata_cols, count_cols)
  
  # subset missing ORFs
  missing_orfs_subset <- missing_orfs[, cols_to_keep]
  rows_to_zero <- sample(nrow(missing_orfs_subset), floor(0.8 * nrow(missing_orfs_subset)))
  missing_orfs_subset[rows_to_zero, count_cols] <- 0
  
  # merge kept ORFs with cnt
  pos <- merge(gene_ids_df, cnt, by = "Geneid", all.x = TRUE)
  
  # combine and order
  combined_df <- rbind(pos[, cols_to_keep], missing_orfs_subset)
  combined_df <- combined_df[order(combined_df$Geneid, combined_df$Chr, combined_df$Start, combined_df$End), ]
  rownames(combined_df) <- NULL
  
  return(combined_df)
}