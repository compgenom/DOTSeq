#' Calculate ORF usage across conditions from a SummarizedExperiment
#'
#' This function extracts post hoc results for a given gene ID and
#' computes condition-specific ORF usage using emmeans contrasts.
#'
#' @param sumExp A SummarizedExperiment object containing DOUResults
#' @param gene_id A character string specifying the gene ID of interest
#'
#' @return A data.frame with ORF IDs, conditions, and usage estimates
#'
#' @details
#' Usage is computed as the inverse logit of the difference between
#' ribo and RNA strategy emmeans per condition. Only valid emmGrid
#' objects are used. Strategy levels are assigned using a helper
#' function \code{\link{assign_strategy_levels}}.
#'
#' @importFrom emmeans emmeans
#' @importFrom boot inv.logit
#' 
#' @export
#' @examples
#' # Load test data
#' dir <- system.file("extdata", package = "DOTSeq")
#' cnt <- read.table(
#'     file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
#'     header = TRUE,
#'     comment.char = "#"
#' )
#' names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
#'
#' flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#'
#' m <- DOTSeqDataSet(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     bed = bed
#' )
#'
#' m$sumExp <- m$sumExp[
#'     SummarizedExperiment::rowData(m$sumExp)$is_kept == TRUE, ]
#' # Subset only one gene
#' m$sumExp <- m$sumExp[
#'     SummarizedExperiment::rowData(m$sumExp)$gene_id == "ENSG00000119402.18", ]
#'
#' m <- DOTSeq(dotseq_dataset = m, modules = "DOU")
#'
#' usage_df <- calculate_orf_usage(
#'     m$sumExp,
#'     gene_id = "ENSG00000119402.18"
#' )
#' print(usage_df)
#' 
calculate_orf_usage <- function(sumExp, gene_id) {
    # Filter by gene ID
    se_gene <- sumExp[rowData(sumExp)$gene_id == gene_id, ]
    orf_rows <- rowData(se_gene)
    
    # Extract DOUResults for those ORFs
    posthoc_results <- rowData(se_gene)[["DOUResults"]][rownames(orf_rows)]
    
    # Apply posthoc to each ORF
    emm_list <- lapply(posthoc_results, posthoc)
    
    # Filter out non-emmGrid objects
    emm_list <- Filter(function(x) inherits(x, "emmGrid"), emm_list)
    
    # Apply usage scale transformation
    usage_scale <- lapply(names(emm_list), function(orf_name) {
        emm <- emm_list[[orf_name]]
        emm_result <- summary(
            emmeans::emmeans(
                emm,
                specs = ~ condition * strategy,
                adjust = "none"
            )
        )
        
        strategy_levels <- assign_strategy_levels(
            emm_result, 
            strategy_col = "strategy"
        )
        
        # Create a named vector or list of usage values per condition
        usage_per_condition <- setNames(nm = levels(emm_result$condition))
        
        for (condition in levels(emm_result$condition)) {
            beta_ribo <- emm_result[emm_result$strategy == strategy_levels$ribo_level & emm_result$condition == condition, ]$emmean
            beta_rna  <- emm_result[emm_result$strategy == strategy_levels$rna_level  & emm_result$condition == condition, ]$emmean
            
            usage <- inv.logit(beta_ribo - beta_rna)
            usage_per_condition[[condition]] <- usage
        }
        
        return(usage_per_condition)
    })
    names(usage_scale) <- names(emm_list)
    
    # usage_scale is a named list of named vectors/lists
    # Each top-level name is an ORF, each sub-name is a condition
    usage_df <- do.call(rbind, lapply(names(usage_scale), function(orf) {
        data.frame(
            orf_id = orf,
            condition = names(usage_scale[[orf]]),
            usage = unlist(usage_scale[[orf]]),
            stringsAsFactors = FALSE
        )
    }))
    
    return(usage_df)
}



#' Calculate Translational Efficiency and Shifts in ORF Usage
#'
#' This function computes translational efficiency (TE) and shifts in ORF 
#' usage for a set of normalized RNA-seq and Ribo-seq counts. TE is 
#' calculated as the ratio of Ribo-seq counts to RNA-seq counts. Shifts 
#' in ORF usage are computed as the log2 ratio of Ribo-seq proportion to 
#' RNA-seq proportion within each gene.
#'
#' @param counts A numeric matrix or data frame of normalized counts, where 
#'     rows correspond to ORFs and columns correspond to samples. RNA-seq 
#'     and Ribo-seq reads are distinguished by suffixes.
#'
#' @param sample_delim Character. Delimiter used in sample names.
#'     Default is \code{"."}.
#'
#' @param rna_suffix Character string indicating the suffix of RNA-seq 
#'     samples in the column names. Default is \code{".rna"}.
#'
#' @param ribo_suffix Character string indicating the suffix of Ribo-seq 
#'     samples in the column names. Default is \code{".ribo"}.
#'
#' @param pseudocount Numeric value added to counts to avoid division by 
#'     zero. Default is \code{1e-6}.
#'
#' @return A list with two elements:
#'     \describe{
#'         \item{te}{
#'             A matrix of translational efficiency (Ribo / RNA).
#'         }
#'         \item{usage_shift}{
#'             A matrix of log2 ratios of Ribo-proportion / RNA-proportion
#'             within each gene.
#'         }
#'     }
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' result <- calculateTE(counts)
#' head(result$te)
#' }
calculateTE <- function(
        counts,
        rna_suffix = ".rna",
        ribo_suffix = ".ribo",
        sample_delim = ".",
        pseudocount = 1e-6
) {
    # Identify RNA and Ribo columns (match anywhere in name)
    rna_cols <- grep(rna_suffix, colnames(counts), value = TRUE)
    ribo_cols <- grep(ribo_suffix, colnames(counts), value = TRUE)
    
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
    te <- matrix(NA, nrow = nrow(counts), ncol = length(common_prefixes))
    rownames(te) <- rownames(counts)
    colnames(te) <- paste0(common_prefixes, ".te")
    
    for (prefix in common_prefixes) {
        rnaCol <- grep(
            paste0(prefix, rna_suffix), 
            colnames(counts), 
            value = TRUE
        )
        riboCol <- grep(
            paste0(prefix, ribo_suffix), 
            colnames(counts), 
            value = TRUE
        )
        if (length(rnaCol) != 1 || length(riboCol) != 1) {
            warning("Ambiguous or missing match for prefix: ", prefix)
            next
        }
        rnaVals <- counts[, rnaCol]
        riboVals <- counts[, riboCol]
        te[, paste0(prefix, ".te")] <- (riboVals + pseudocount) /
            (rnaVals + pseudocount)
    }
    
    # Compute usage shift
    gene_ids <- sub(":O.*", "", rownames(counts))
    ribo_mat <- counts[, ribo_cols, drop = FALSE]
    rna_mat <- counts[, rna_cols, drop = FALSE]
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


