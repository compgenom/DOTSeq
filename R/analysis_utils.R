#' Calculate Translational Efficiency and Shifts in ORF Usage 
#'
#' This function computes translational efficiency (TE) and shifts in ORF usage
#' for a set of normalized RNA-seq and Ribo-seq counts. TE is calculated
#' as the ratio of Ribo-seq counts to RNA-seq counts. Shifts in ORF usage
#' is calculated as the log2 ratio of Ribo-seq proportion to RNA-seq
#' proportion within each gene.
#'
#' @param counts A numeric matrix or data frame of normalized counts,
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
#'
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' result <- calculateTE(counts)
#' head(result$te)
#' }
calculateTE <- function(counts, 
                        rna_suffix = ".rna", 
                        ribo_suffix = ".ribo", 
                        sample_delim = ".",
                        pseudocount = 1e-6) {
  
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
    rnaCol <- grep(paste0(prefix, rna_suffix), colnames(counts), value = TRUE)
    riboCol <- grep(paste0(prefix, ribo_suffix), colnames(counts), value = TRUE)
    if (length(rnaCol) != 1 || length(riboCol) != 1) {
      warning("Ambiguous or missing match for prefix: ", prefix)
      next
    }
    rnaVals <- counts[, rnaCol]
    riboVals <- counts[, riboCol]
    te[, paste0(prefix, ".te")] <- (riboVals + pseudocount) / (rnaVals + pseudocount)
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


#' Plot Venn Diagram of DTE and DOU Significance Overlap
#'
#' @description
#' Generates a Venn diagram (Euler diagram) showing the overlap of significantly
#' differentially translated ORFs (DTE) and differentially used ORFs (DOU).
#' ORFs are classified as significant in DTE only, DOU only, or both.
#'
#' @param results A data frame containing DTE and DOU results. Must include row names
#'   corresponding to ORF identifiers, and columns for adjusted p-values from both tests.
#' @param dou_padj_col Character string specifying the column name for DOU significance values.
#'   Should correspond to local false sign rate (lfsr). Default is \code{"lfsr"}.
#' @param dte_padj_col Character string specifying the column name for DTE adjusted p-values.
#'   Default is \code{"padj"}.
#'
#' @return A Venn diagram (Euler plot) showing the number of ORFs significant in DTE only,
#'   DOU only, or both.
#'
#' @details
#' Significance is determined using a threshold of 0.05 on the adjusted p-values.
#' The plot uses a color-blind friendly palette and includes counts for each region.
#'
#' @importFrom eulerr euler
#' @importFrom graphics plot
#' 
#' @keywords internal
#' 
plot_venn <- function(results, dou_padj_col = "lfsr", dte_padj_col = "padj") {
  
  padj_sig <- !is.na(results[[dte_padj_col]]) & results[[dte_padj_col]] < 0.05
  fdr_sig <- !is.na(results[[dou_padj_col]]) & results[[dou_padj_col]] < 0.05
  
  padj_set <- rownames(results)[padj_sig]
  fdr_set <- rownames(results)[fdr_sig]
  
  # Create Euler fit
  fit <- euler(c(
    "DTE" = length(setdiff(padj_set, fdr_set)),
    "DOU" = length(setdiff(fdr_set, padj_set)),
    "DTE&DOU" = length(intersect(padj_set, fdr_set))
  ))
  
  # Define color-blind friendly palette
  venn_colors <- c("DTE" = "#0072B2", "DOU" = "#E69F00", "DTE&DOU" = "#CC79A7")
  
  # pdf("dte_dot_venn_cycling_arrest.pdf", 2.5, 2)
  plot(
    fit,
    fills = list(fill = venn_colors, alpha = 0.6),
    labels = list(font = 2),
    quantities = TRUE,
    main = "Differentially translated ORFs"
    )
  
}


#' Plot Composite Scatter and Marginal Histograms for DTE vs DOU
#'
#' @description
#' Visualizes the relationship between differential translation efficiency (DTE)
#' and differential ORF usage (DOU) using a scatter plot with marginal histograms.
#' ORFs are colored based on significance in DTE, DOU, or both tests. 
#' This plot helps assess the overlap and divergence between DTE and DOU signals.
#'
#' @param results A data frame containing results from DTE and DOU post hoc contrasts.
#'   Must include columns for DTE and DOU posterior estimates, and adjusted p-values
#'   (e.g., \code{padj} for DTE and \code{lfsr} for DOU).
#' @param dou_estimates_col Character string specifying the column name for DOU estimates.
#'   Default is \code{"PosteriorMean"}.
#' @param dou_padj_col Character string specifying the column name for DOU significance values.
#'   Should correspond to local false sign rate (lfsr). Default is \code{"lfsr"}.
#' @param dte_estimates_col Character string specifying the column name for DTE estimates.
#'   Default is \code{"log2FoldChange"}.
#' @param dte_padj_col Character string specifying the column name for DTE adjusted p-values.
#'   Default is \code{"padj"}.
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of the DOU estimates.
#'   Useful for aligning directionality with DTE. Default is \code{TRUE}.
#' @param lhist Integer; number of bins to use for the marginal histograms.
#'   Default is \code{20}.
#'
#' @return A composite plot with:
#' \describe{
#'   \item{Scatter plot}{Displays DTE (log2 fold-change) vs DOU (log-odds change) estimates,
#'     with points colored by significance category: DTE-only, DOU-only, both, or neither.}
#'   \item{Marginal histograms}{Show the distribution of DTE and DOU estimates along the top and right margins.}
#' }
#' Additionally, the function prints the Spearman correlation between DTE and DOU estimates to the console.
#' 
#' @details
#' The function uses base R graphics and a custom layout to combine scatter and histogram plots.
#' Significant ORFs are determined using a threshold of 0.05 on adjusted p-values.

#'
#' @importFrom graphics plot points legend par layout barplot plot.new hist text
#' @importFrom grDevices adjustcolor
#' @importFrom stats cor.test density
#' 
#' @keywords internal
#' 
plot_composite <- function(
    results,
    dou_estimates_col = "PosteriorMean",
    dou_padj_col = "lfsr",
    dte_estimates_col = "log2FoldChange",
    dte_padj_col = "padj",
    flip_sign = TRUE,
    lhist = 20
) {
  local({
    old_par <- par(no.readonly = TRUE)
    on.exit({
      par(old_par)
      layout(1)  # Reset layout to default
    }, add = TRUE)
    
    if (isTRUE(flip_sign)) {
      correlation_results <- cor.test(results[[dte_estimates_col]], results[[dou_estimates_col]]  * -1, method = "spearman")
      results[[dou_estimates_col]] <- results[[dou_estimates_col]] * -1
    } else {
      correlation_results <- cor.test(results[[dte_estimates_col]], results[[dou_estimates_col]], method = "spearman")
    }
    
    print(correlation_results)
    
    results <- na.omit(results)
    padj_sig <- !is.na(results[[dte_padj_col]]) & results[[dte_padj_col]] < 0.05
    fdr_sig <- !is.na(results[[dou_padj_col]]) & results[[dou_padj_col]] < 0.05
    both_sig <- padj_sig & fdr_sig
    padj_only <- padj_sig & !fdr_sig
    fdr_only <- fdr_sig & !padj_sig
    
    col_dte <- adjustcolor("#0072B2", alpha.f = 0.6)
    col_dou <- adjustcolor("#E69F00", alpha.f = 0.6)
    col_both <- adjustcolor("#CC79A7", alpha.f = 0.6)
    col_none <- adjustcolor("grey80", alpha.f = 0.6)
    
    # Layout setup
    layMat <- matrix(c(1, 4, 3, 2), ncol = 2)
    layout(layMat, widths = c(6/7, 1/7), heights = c(1/7, 6/7))
    ospc <- 0.5
    pext <- 4
    bspc <- 1
    par(mar = c(pext, pext, bspc, bspc), oma = rep(ospc, 4))
    
    xlim <- range(results[[dte_estimates_col]])
    ylim <- range(results[[dou_estimates_col]])
    
    # Top histogram
    xhist <- hist(results[[dte_estimates_col]], breaks = seq(xlim[1], xlim[2], length.out = lhist), plot = FALSE)
    par(mar = c(0, pext, 0, 0))
    barplot(xhist$density, axes = FALSE, ylim = c(0, max(xhist$density)), space = 0,
            col = "gray", border = "black")
    
    # Right histogram
    yhist <- hist(results[[dou_estimates_col]], breaks = seq(ylim[1], ylim[2], length.out = lhist), plot = FALSE)
    par(mar = c(pext, 0, 0, 0))
    barplot(yhist$density, axes = FALSE, xlim = c(0, max(yhist$density)), space = 0,
            col = "gray", border = "black", horiz = TRUE)
    
    # Placeholder
    par(mar = c(0, 0, 0, 0))
    plot.new()
    
    # Main scatter plot
    par(mar = c(pext, pext, 0, 0))
    plot(results[[dte_estimates_col]], results[[dou_estimates_col]],
         col = col_none, pch = 16,
         xlab = "log2 fold-change in ORF TE",
         ylab = "log-odds change in ORF usage")
    
    points(results[[dte_estimates_col]][padj_only], results[[dou_estimates_col]][padj_only], col = col_dte)
    points(results[[dte_estimates_col]][fdr_only], results[[dou_estimates_col]][fdr_only], col = col_dou)
    points(results[[dte_estimates_col]][both_sig], results[[dou_estimates_col]][both_sig], col = col_both)
    
    legend("bottomright", legend = c("DTE", "DOU", "Both"),
           col = c(col_dte, col_dou, col_both), pch = 1, bty = "n", inset = 0.05)
  })
  
}


#' Prepare Data and Metadata for DOU Heatmap
#'
#' @description
#' Prepares a matrix of differential ORF usage (DOU) estimates for heatmap visualization,
#' comparing mORFs with either uORFs or dORFs. Filters significant ORFs based on local false sign rate (lfsr),
#' clusters rows, and retrieves gene annotations from Ensembl for labeling.
#'
#' @param results A data frame from \code{testDOT} containing DOU estimates, gene labels, and group IDs.
#'   Must include columns for DOU effect sizes and lfsr values.
#' @param dou_estimates_col Character string specifying the column name for DOU effect size estimates.
#'   Default is \code{"PosteriorMean"}.
#' @param dou_padj_col Character string specifying the column name for DOU significance values.
#'   Should correspond to local false sign rate (lfsr). Default is \code{"lfsr"}.
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of DOU estimates to align directionality
#'   with translation efficiency. Default is \code{TRUE}.
#' @param orf_type Character string specifying the ORF type to compare with mORFs.
#'   Accepts \code{"uORF"} or \code{"dORF"}.
#' @param species_dataset Ensembl dataset name, e.g., \code{"hsapiens_gene_ensembl"}.
#' @param symbol_col Column name for gene symbols in Ensembl, e.g., \code{"hgnc_symbol"} or \code{"mgi_symbol"}.
#' @param threshold Numeric threshold for lfsr-based significance filtering. Default is \code{0.05}.
#'
#' @return A list containing:
#' \describe{
#'   \item{ordered_matrix}{Matrix of DOU estimates for significant ORFs.}
#'   \item{row_dend_clean}{Dendrogram object for row clustering.}
#'   \item{highlight_df}{Data frame indicating significant tiles to highlight.}
#'   \item{genes}{Full gene annotation data retrieved from Ensembl.}
#'   \item{gene_labels}{Vector of gene symbols for heatmap rows.}
#'   \item{color_palette}{Color palette for heatmap visualization.}
#'   \item{color_breaks}{Breaks used for color scaling.}
#'   \item{abs_max}{Maximum absolute value used for color scaling.}
#' }
#'
#' @details
#' This function is designed to support visualization of differential ribosome loading across ORFs
#' within genes. It uses empirical Bayes shrinkage via the \code{ashr} package to compute lfsr,
#' and filters ORFs with lfsr below the specified threshold. Gene annotations are retrieved from
#' Ensembl using the \code{biomaRt} package.
#'
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom stats setNames dist hclust as.dendrogram order.dendrogram complete.cases
#' @importFrom grDevices colorRampPalette
#' 
#' @keywords internal
#' 
cistronic_data <- function(results, dou_rowdata, dou_estimates_col = "PosteriorMean", dou_padj_col = "lfsr", flip_sign = TRUE, sorf_type = "uORF", 
                       species_dataset = "hsapiens_gene_ensembl", symbol_col = "hgnc_symbol", threshold = 0.05) {
  
  results$gene_id <- sub("\\..*", "", results$orf_id)
  
  sig_genes <- results[results[[dou_padj_col]]<threshold,]$gene_id
  sig_res <- results[results$gene_id %in% sig_genes, ]
  
  rowdata <- dou_rowdata
  rowdata$DOUResults <- NULL
  rowdata$orf_id <- rownames(rowdata)
  rowdata <- as.data.frame(rowdata)[, c("orf_id", "orf_type")]
  
  sig_res <- merge(sig_res, rowdata, by.x = "orf_id")[, c("gene_id", "orf_type", dou_estimates_col, dou_padj_col)]
  
  sorf_res <- sig_res[sig_res$orf_type == sorf_type, ][, c("gene_id", dou_estimates_col)]
  morf_res <- sig_res[sig_res$orf_type == "mORF", ][, c("gene_id", dou_estimates_col)]
  sig_mat <- merge(sorf_res, morf_res, by = "gene_id")
  names(sig_mat) <- c("gene_id", sorf_type, "mORF")
  sig_mat$delta <- abs(sig_mat$mORF - sig_mat[[sorf_type]])
  sig_mat <- sig_mat[order(sig_mat$gene_id, -sig_mat$delta), ]
  sig_mat <- sig_mat[!duplicated(sig_mat$gene_id), ]
  
  sig_mat <- sig_mat[!is.na(sig_mat$gene_id), ]
  rownames(sig_mat) <- sig_mat$gene_id
  
  sig_mat <- sig_mat[ , !(names(sig_mat) %in% c("gene_id", "delta"))]
  
  sig_mat[] <- lapply(sig_mat, as.numeric)
  sig_mat <- as.matrix(sig_mat)
  sig_mat_clean <- sig_mat[complete.cases(sig_mat), ]
  
  row_cluster_clean <- hclust(dist(sig_mat_clean))
  row_dend_clean <- as.dendrogram(row_cluster_clean)
  
  min_val <- min(sig_mat_clean, na.rm = TRUE)
  max_val <- max(sig_mat_clean, na.rm = TRUE)
  abs_max <- max(abs(min_val), abs(max_val))
  
  color_breaks <- seq(-abs_max, abs_max, length.out = 101)
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  row_order <- order.dendrogram(row_dend_clean)
  ordered_matrix <- sig_mat_clean[row_order, ]
  sig_res_filtered <- sig_res[
    sig_res$orf_type %in% colnames(ordered_matrix) &
      sig_res$gene_id %in% rownames(ordered_matrix),
  ]
  group_to_row <- setNames(seq_len(nrow(ordered_matrix)), rownames(ordered_matrix))
  label_to_col <- setNames(seq_len(ncol(ordered_matrix)), colnames(ordered_matrix))
  highlight_df <- data.frame(
    row = group_to_row[sig_res_filtered$gene_id],
    col = label_to_col[sig_res_filtered$orf_type],
    color = "black"
  )
  highlight_df <- highlight_df[complete.cases(highlight_df), ]
  
  ensembl_ids <- sub("\\..*", "", rownames(ordered_matrix))
  ensembl <- useEnsembl(biomart = "genes", dataset = species_dataset)
  genes <- getBM(
    attributes = c("ensembl_gene_id", symbol_col), #, "go_id", "name_1006", "namespace_1003"
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  genes_unique <- genes[!duplicated(genes[c("ensembl_gene_id", symbol_col)]), ]
  genes_unique_sorted <- genes_unique[match(ensembl_ids, genes_unique$ensembl_gene_id), ]
  
  if (isTRUE(flip_sign)) {
    ordered_matrix <- ordered_matrix * -1
  }
  
  return(list(
    ordered_matrix = ordered_matrix,
    row_dend_clean = row_dend_clean,
    highlight_df = highlight_df,
    genes = genes,
    gene_labels = genes_unique_sorted[[symbol_col]],
    color_palette = color_palette,
    color_breaks = color_breaks,
    abs_max = abs_max
  ))
}


#' Plot DOU Heatmap with Dendrogram and Highlights
#'
#' This internal function generates a heatmap of differential ORF usage (DOU) estimates,
#' with hierarchical clustering and optional tile highlighting for significant ORFs.
#' All required inputs are passed via a named list (`paired_data`).
#'
#' @param paired_data Output from cistronic_data. A named list containing the following elements:
#' \describe{
#'   \item{ordered_matrix}{Matrix of DOU estimates to plot.}
#'   \item{row_dend_clean}{Dendrogram object for row ordering.}
#'   \item{highlight_df}{Data frame specifying tiles to highlight, with columns \code{row}, \code{col}, and \code{color}.}
#'   \item{genes}{Full gene annotation data retrieved from Ensembl.}
#'   \item{gene_labels}{Character vector of gene symbols for row labels.}
#'   \item{color_palette}{Vector of colors for the heatmap.}
#'   \item{color_breaks}{Numeric vector of breaks for color scaling.}
#'   \item{abs_max}{Maximum absolute value for the color scale.}
#' }
#'
#' @importFrom graphics axis image mtext rect par layout plot
#'
#' @keywords internal
#'
#' @return No return value; used for plotting.
#' 
plot_heatmap <- function(paired_data) { # width = NULL, height = NULL, output_file = "dou.pdf"
  
  ordered_matrix <- paired_data$ordered_matrix
  row_dend_clean <- paired_data$row_dend_clean
  highlight_df <- paired_data$highlight_df
  gene_labels <- paired_data$gene_labels
  color_palette <- paired_data$color_palette
  color_breaks <- paired_data$color_breaks
  abs_max <- paired_data$abs_max
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  n_rows <- nrow(ordered_matrix)
  n_cols <- ncol(ordered_matrix)
  
  # Dynamic scaling
  cex_row <- ifelse(n_rows > 50, 0.5, 0.8)
  cex_col <- ifelse(n_cols > 50, 0.5, 0.8)
  right_margin <- max(5, min(5, n_rows / 5))
  layout_widths <- c(2, max(6, n_cols / 10))
  key_bottom <- ifelse(n_rows > 50, 0.13, 0.13)
  key_top <- key_bottom + 0.08
  
  # pdf(output_file, width = width, height = height)
  layout(matrix(c(1, 2), nrow = 1), widths = layout_widths)
  
  ## Panel 1: Dendrogram
  par(mar = c(5, 1.5, 4, 0))
  plot(row_dend_clean, horiz = TRUE, axes = FALSE, yaxs = "i")
  
  ## Panel 2: Heatmap
  par(mar = c(5, 0.5, 4, right_margin))
  image(x = 1:ncol(ordered_matrix),
        y = 1:nrow(ordered_matrix),
        z = t(ordered_matrix),
        col = color_palette,
        breaks = color_breaks,
        axes = FALSE,
        xlab = "")
  
  for (i in seq_len(nrow(highlight_df))) {
    rect(
      xleft = highlight_df$col[i] - 0.5,
      xright = highlight_df$col[i] + 0.5,
      ybottom = highlight_df$row[i] - 0.5,
      ytop = highlight_df$row[i] + 0.5,
      border = highlight_df$color[i],
      lwd = 1.5
    )
  }
  
  axis(1, at = 1:ncol(ordered_matrix), labels = colnames(ordered_matrix), las = 2, cex.axis = cex_col)
  axis(4, at = 1:nrow(ordered_matrix), labels = gene_labels, las = 2, cex.axis = cex_row)
  
  mtext("Differential ORF usage", side = 3, line = 1.5, cex = 1.2, adj = 0.5, font = 2)
  
  ## Color key
  par(fig = c(0.75, 0.9, key_bottom, key_top), new = TRUE, mar = c(0, 0, 2, 0))
  image(z = t(matrix(seq(-abs_max, abs_max, length.out = 100), nrow = 1)),
        col = color_palette, axes = FALSE)
  
  tick_vals <- pretty(c(-abs_max, abs_max), n = 5)
  tick_pos <- (tick_vals - (-abs_max)) / (2 * abs_max)
  axis(1, at = tick_pos, labels = tick_vals, las = 1, cex.axis = 0.7)
  mtext("log-odds", side = 1, line = 2, cex = 0.8)
  
  # dev.off()
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
}


#' Plot Volcano Plot for Differential ORF Usage (DOU)
#'
#' @description
#' Generates a volcano plot for visualizing differential ORF usage (DOU) results.
#' The plot displays log-odds changes (effect sizes) on the x-axis and the negative
#' log10-transformed local false sign rate (lfsr) on the y-axis. Points are colored
#' based on significance in DTE, DOU, or both.
#'
#' @param results A data frame containing DOU and DTE results. Must include columns for
#'   DOU estimates and lfsr, and optionally DTE significance indicators.
#' @param paired_data Output from cistronic_data. A named list containing the following elements:
#' \describe{
#'   \item{ordered_matrix}{Matrix of DOU estimates to plot.}
#'   \item{row_dend_clean}{Dendrogram object for row ordering.}
#'   \item{highlight_df}{Data frame specifying tiles to highlight, with columns \code{row}, \code{col}, and \code{color}.}
#'   \item{genes}{Full gene annotation data retrieved from Ensembl.}
#'   \item{gene_labels}{Character vector of gene symbols for row labels.}
#'   \item{color_palette}{Vector of colors for the heatmap.}
#'   \item{color_breaks}{Numeric vector of breaks for color scaling.}
#'   \item{abs_max}{Maximum absolute value for the color scale.}
#' }
#' @param symbol_col Column name for gene symbols in Ensembl, e.g., \code{"hgnc_symbol"} or \code{"mgi_symbol"}.
#' @param dou_estimates_col Character string specifying the column name for DOU estimates.
#'   Default is \code{"PosteriorMean"}.
#' @param dou_padj_col Character string specifying the column name for DOU significance values.
#'   Should correspond to local false sign rate (lfsr). Default is \code{"lfsr"}.
#' @param dte_estimates_col Character string specifying the column name for DTE estimates.
#'   Default is \code{"log2FoldChange"}.
#' @param dte_padj_col Character string specifying the column name for DTE adjusted p-values.
#'   Default is \code{"padj"}.
#' @param dou_estimates_threshold Numeric threshold for effect size significance. Default is \code{1}.
#' @param dou_padj_threshold Numeric threshold for lfsr significance. Default is \code{0.05}.
#' @param dou_estimates_threshold Numeric. Threshold for DOU effect size significance. Default is \code{NULL}.
#' @param extreme_threshold Numeric. Threshold for labeling extreme points (in -log10 LFSR). Default is \code{NULL}.
#' 
#' @return A volcano plot is displayed in a new graphics device.
#'
#' @importFrom graphics plot points abline legend
#'
#' @keywords internal
#'
plot_volcano <- function(
    results,
    paired_data,
    symbol_col = "hgnc_symbol",
    dou_estimates_col = "PosteriorMean",
    dou_padj_col = "lfsr",
    dte_estimates_col = "log2FoldChange",
    dte_padj_col = "padj",
    flip_sign = TRUE,
    dou_estimates_threshold = 1,
    dou_padj_threshold = 0.05,
    extreme_threshold = NULL
    ) {
  
  genes_unique <- paired_data$genes[!duplicated(paired_data$genes[c("ensembl_gene_id", symbol_col)]), ]
  results$ensembl_gene_id <- sub("\\..*", "", results$orf_id)
  results <- merge(results, genes_unique, by = "ensembl_gene_id", all.x = TRUE)
  
  # results <- na.omit(results)
  padj_sig <- !is.na(results[[dte_padj_col]]) & results[[dte_padj_col]] < 0.05
  fdr_sig <- !is.na(results[[dou_padj_col]]) & results[[dou_padj_col]] < 0.05
  both_sig <- padj_sig & fdr_sig
  padj_only <- padj_sig & !fdr_sig
  fdr_only <- fdr_sig & !padj_sig
  
  col_dte = "#0072B2"
  col_dou = "#E69F00"
  col_both = "#CC79A7"
  col_none = "grey80"
  
  if (isTRUE(flip_sign)) {
    estimates <- results[[dou_estimates_col]] * -1
  } else {
    estimates <- results[[dou_estimates_col]]
  }
  
  lfsr <- results[[dou_padj_col]]
  loglfsr <- -log10(lfsr)
  
  plot(estimates, loglfsr, pch = 20, col = col_none,
                 main = NULL, #"Volcano Plot for Differential ORF Usage",
                 xlab = "Log-odds change in ORF usage",
       ylab = expression(paste("-log"[10], "(LFSR)")))
  
  if (!is.null(extreme_threshold)) {
    extreme_idx <- which(loglfsr > extreme_threshold)
    text(
      x = estimates[extreme_idx],
      y = loglfsr[extreme_idx],
      labels = results$hgnc_symbol[extreme_idx],
      pos = 3, cex = 0.7
    )
  }
  
  points(estimates[padj_only], loglfsr[padj_only], col = col_dte)
  points(estimates[fdr_only], loglfsr[fdr_only], col = col_dou)
  points(estimates[both_sig], loglfsr[both_sig], col = col_both)
  
  abline(h = -log10(dou_padj_threshold), col = "black", lty = 2)
  abline(v = c(-dou_estimates_threshold, dou_estimates_threshold), col = "black", lty = 2)
  
  legend("topright", legend = c("DTE", "DOU", "Both"),
                   col = c(col_dte, col_dou, col_both),
                   pch = 1, bty = "n", inset = c(0.02, 0.05))
}


#' Generate Differential ORF Translation (DOT) Visualization Suite
#'
#' @description
#' A high-level wrapper that visualizes differential ORF usage (DOU) and translation efficiency (DTE)
#' relationships through Venn, volcano, composite, and heatmap plots. It integrates Ensembl gene
#' annotations and highlights significant ORFs based on empirical Bayes shrinkage (via `ashr`).
#'
#' @param results Data frame with DOU and DTE estimates and significance values.
#' @param dou_rowdata Data frame containing ORF-level metadata (e.g., gene IDs, ORF class).
#' @param dou_estimates_col, dou_padj_col Column names for DOU estimate and significance (default: `"PosteriorMean"`, `"lfsr"`).
#' @param dte_estimates_col, dte_padj_col Column names for DTE estimate and significance (default: `"log2FoldChange"`, `"padj"`).
#' @param dou_estimates_threshold Numeric threshold for effect size significance. Default is \code{1}.
#' @param dou_padj_threshold Numeric threshold for lfsr significance. Default is \code{0.05}.
#' @param extreme_threshold Numeric. Threshold for labeling extreme points (in -log10 LFSR). Default is \code{NULL}.
#' @param flip_sign Logical; if TRUE, flips the sign of DOU estimates for consistency with TE (default: TRUE).
#' @param sorf_type Character; `"uORF"` or `"dORF"`, specifying comparison type.
#' @param species_dataset,symbol_col Ensembl dataset and symbol column for annotation.
#'
#' @return Data frame of annotated gene–ORF pairs used in the heatmap visualization.
#' 
#' @export
#' 
plotDOT <- function(
    results, 
    dou_rowdata,
    dou_estimates_col = "PosteriorMean", 
    dou_padj_col = "lfsr", 
    dte_estimates_col = "log2FoldChange", 
    dte_padj_col = "padj", 
    dou_estimates_threshold = 1,
    dou_padj_threshold = 0.05,
    extreme_threshold = NULL,
    flip_sign = TRUE, 
    sorf_type = "uORF", 
    species_dataset = "hsapiens_gene_ensembl", 
    symbol_col = "hgnc_symbol"
    ) {
  
  old_par <- par(no.readonly = TRUE) 
  on.exit(par(old_par), add = TRUE)
  
  venn <- plot_venn(
    results = results, 
    dou_padj_col = dou_padj_col, 
    dte_padj_col = dte_padj_col
  )
  grid::grid.draw(venn)
  
  composite <- plot_composite(
    results = results, 
    dou_estimates_col = dou_estimates_col, 
    dou_padj_col = dou_padj_col, 
    dte_estimates_col = dte_estimates_col, 
    dte_padj_col = dte_padj_col, 
    flip_sign = flip_sign, 
    lhist=20
  )
  
  paired_data <- cistronic_data(
    results = results,
    dou_rowdata = dou_rowdata,
    dou_estimates_col = dou_estimates_col,
    dou_padj_col = dou_padj_col,
    flip_sign = flip_sign,
    sorf_type = sorf_type,
    species_dataset = species_dataset,
    symbol_col = symbol_col
    )

  volcano <- plot_volcano(
    results = results,
    paired_data = paired_data,
    dou_estimates_col = dou_estimates_col,
    dou_padj_col = dou_padj_col,
    dte_estimates_col = dte_estimates_col, 
    dte_padj_col = dte_padj_col, 
    dou_estimates_threshold = dou_estimates_threshold,
    dou_padj_threshold = dou_padj_threshold,
    extreme_threshold = extreme_threshold
  )
  
  tryCatch({
    heatmap <- plot_heatmap(paired_data)
  }, error = function(e) {
    message("Try enlarging the plot area: ", e$message)
  })
  
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, fig = c(0, 1, 0, 1))
  
  return(list(
    venn = venn,
    volcano = volcano,
    composite = composite,
    heatmap = heatmap,
    heatmap_data = paired_data$genes
  ))
}

