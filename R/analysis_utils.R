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


#' Plot DTE vs DOU Comparison with Marginal Histograms
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
#' @return Generates a composite plot with:
#' \describe{
#'   \item{Scatter plot}{Displays DTE (log2 fold-change) vs DOU (log-odds change) estimates,
#'     with points colored by significance category: DTE-only, DOU-only, both, or neither.}
#'   \item{Marginal histograms}{Show the distribution of DTE and DOU estimates along the top and right margins.}
#' }
#'
#' @details
#' The function uses base R graphics and a custom layout to combine scatter and histogram plots.
#' Significant ORFs are determined using a threshold of 0.05 on adjusted p-values.
#' The color scheme is designed to be accessible to color-blind viewers.
#'
#' @importFrom graphics plot points legend par layout barplot plot.new hist
#' @importFrom grDevices adjustcolor
#' @importFrom stats density
#' @importFrom eulerr euler
#' 
plotDOT <- function(results, dou_estimates_col = "PosteriorMean", dou_padj_col = "lfsr", dte_estimates_col = "log2FoldChange", dte_padj_col = "padj", flip_sign = TRUE, lhist=20) {
  if (isTRUE(flip_sign)) {
    res[[dou_estimates_col]] <- -res[[dou_estimates_col]]
  }
  
  padj_sig <- !is.na(res[[dte_padj_col]]) & res[[dte_padj_col]] < 0.05
  fdr_sig <- !is.na(res[[dou_padj_col]]) & res[[dou_padj_col]] < 0.05
  both_sig <- padj_sig & fdr_sig
  padj_only <- padj_sig & !fdr_sig
  fdr_only <- fdr_sig & !padj_sig
  
  col_dte <- adjustcolor("#0072B2", alpha.f = 0.6)
  col_dou <- adjustcolor("#E69F00", alpha.f = 0.6)
  col_both <- adjustcolor("#CC79A7", alpha.f = 0.6)
  col_none <- adjustcolor("grey80", alpha.f = 0.6)
  
  old_par <- par(no.readonly = TRUE) 
  
  withCallingHandlers({
    # https://stackoverflow.com/a/31750636
    ## set up layout and graphical parameters
    layMat <- matrix(c(1,4,3,2), ncol=2)
    layout(layMat, widths=c(6/7, 1/7), heights=c(1/7, 6/7))
    ospc <- 0.5                                                  # outer space
    pext <- 4                                                    # par extension down and to the left
    bspc <- 1                                                    # space between scatter plot and bar plots
    par. <- par(mar=c(pext, pext, bspc, bspc), oma=rep(ospc, 4)) # plot parameters
    
    xlim <- range(res[[dte_estimates_col]])
    ylim <- range(res[[dou_estimates_col]])
    
    # Top histogram
    xhist <- hist(res[[dte_estimates_col]], breaks=seq(xlim[1], xlim[2], length.out=lhist), plot=FALSE)
    par(mar=c(0, pext, 0, 0))
    barplot(xhist$density, axes = FALSE, ylim=c(0, max(xhist$density)), space=0,
            col = "gray", border = "black")
    
    # Right histogram (y-axis)
    yhist <- hist(res[[dou_estimates_col]], breaks=seq(ylim[1], ylim[2], length.out=lhist), plot=FALSE)
    par(mar=c(pext, 0, 0, 0))
    barplot(yhist$density, axes = FALSE, xlim=c(0, max(yhist$density)), space=0,
            col = "gray", border = "black", horiz = TRUE)
    
    # placeholder
    dx <- density(res[[dte_estimates_col]])
    dy <- density(res[[dou_estimates_col]])
    par(mar=c(0, 0, 0, 0))
    plot.new()
    
    # Main scatter plot
    par(mar=c(pext, pext, 0, 0))
    plot(res[[dte_estimates_col]], res[[dou_estimates_col]],
         col = col_none, pch = 16,
         xlab = "log2 fold-change in ORF TE",
         ylab = "log-odds change in ORF usage")
    
    points(res[[dte_estimates_col]][padj_only], res[[dou_estimates_col]][padj_only], col = col_dte)
    points(res[[dte_estimates_col]][fdr_only], res[[dou_estimates_col]][fdr_only], col = col_dou)
    points(res[[dte_estimates_col]][both_sig], res[[dou_estimates_col]][both_sig], col = col_both)
    
    legend("bottomright", legend = c("DTE", "DOU", "Both"), col = c(col_dte, col_dou, col_both), pch = 1, bty = "n", inset = 0.05)
    
  }, error = function(e) {
    par(old_par)
    stop(e)
  })
  
  par(old_par)
  
  padj_set <- rownames(res)[padj_sig]
  fdr_set <- rownames(res)[fdr_sig]
  
  # Create Euler fit
  fit <- eulerr::euler(c(
    "DTE" = length(setdiff(padj_set, fdr_set)),
    "DOU" = length(setdiff(fdr_set, padj_set)),
    "DTE&DOU" = length(intersect(padj_set, fdr_set))
  ))
  
  # Define color-blind friendly palette
  venn_colors <- c("DTE" = "#0072B2", "DOU" = "#E69F00", "DTE&DOU" = "#CC79A7")
  
  # pdf("dte_dot_venn_cycling_arrest.pdf", 2.5, 2)
  plot(fit,
       fills = list(fill = venn_colors, alpha = 0.6),
       labels = list(font = 2),
       quantities = TRUE,
       main = "Differentially translated ORFs")
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
#'   \item{orderedMatrix}{Matrix of DOU estimates for significant ORFs.}
#'   \item{rowDendClean}{Dendrogram object for row clustering.}
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
citronicData <- function(results, dou_estimates_col = "PosteriorMean", dou_padj_col = "lfsr", flip_sign = TRUE, orf_type = "uORF", 
                       species_dataset = "hsapiens_gene_ensembl", symbol_col = "hgnc_symbol", threshold = 0.05) {
  sigRes <- results[results[[dou_padj_col]]<threshold,]
  
  orfSig <- sigRes[sigRes$labels==orf_type, ][, c("groupID", dou_estimates_col)]
  mORFRes <- results[results$labels=="mORF", ][, c("groupID", dou_estimates_col)]
  sigMat <- merge(orfSig, mORFRes, by = "groupID")
  names(sigMat) <- c("groupID", orf_type, "mORF")
  sigMat$delta <- abs(sigMat$mORF - sigMat[[orf_type]])
  sigMat <- sigMat[order(sigMat$groupID, -sigMat$delta), ]
  sigMat <- sigMat[!duplicated(sigMat$groupID), ]
  
  sigMat <- sigMat[!is.na(sigMat$groupID), ]
  rownames(sigMat) <- sigMat$groupID
  
  sigMat <- sigMat[ , !(names(sigMat) %in% c("groupID", "delta"))]
  
  sigMat[] <- lapply(sigMat, as.numeric)
  sigMat <- as.matrix(sigMat)
  sigMatClean <- sigMat[complete.cases(sigMat), ]
  
  rowClusterClean <- hclust(dist(sigMatClean))
  rowDendClean <- as.dendrogram(rowClusterClean)
  
  minVal <- min(sigMatClean, na.rm = TRUE)
  maxVal <- max(sigMatClean, na.rm = TRUE)
  absMax <- max(abs(minVal), abs(maxVal))
  
  colorBreaks <- seq(-absMax, absMax, length.out = 101)
  colorPalette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  row_order <- order.dendrogram(rowDendClean)
  orderedMatrix <- sigMatClean[row_order, ]
  sigResFiltered <- sigRes[
    sigRes$labels %in% colnames(orderedMatrix) &
      sigRes$groupID %in% rownames(orderedMatrix),
  ]
  group_to_row <- setNames(seq_len(nrow(orderedMatrix)), rownames(orderedMatrix))
  label_to_col <- setNames(seq_len(ncol(orderedMatrix)), colnames(orderedMatrix))
  highlight_df <- data.frame(
    row = group_to_row[sigResFiltered$groupID],
    col = label_to_col[sigResFiltered$labels],
    color = "black"
  )
  highlight_df <- highlight_df[complete.cases(highlight_df), ]
  
  ensembIds <- sub("\\..*", "", rownames(orderedMatrix))
  ensembl <- useEnsembl(biomart = "genes", dataset = species_dataset)
  genes <- getBM(
    attributes = c("ensembl_gene_id", symbol_col, "go_id", "name_1006", "namespace_1003"),
    filters = "ensembl_gene_id",
    values = ensembIds,
    mart = ensembl
  )
  
  genesUnique <- genes[!duplicated(genes[c("ensembl_gene_id", symbol_col)]), ]
  genesUniqueSorted <- genesUnique[match(ensembIds, genesUnique$ensembl_gene_id), ]
  
  if (isTRUE(flip_sign)) {
    orderedMatrix <- orderedMatrix * -1
  }
  
  return(list(
    orderedMatrix = orderedMatrix,
    rowDendClean = rowDendClean,
    highlight_df = highlight_df,
    genes = genes,
    gene_labels = genesUniqueSorted[[symbol_col]],
    color_palette = colorPalette,
    color_breaks = colorBreaks,
    abs_max = absMax
  ))
}


#' An Internal Function for Plotting DOU Heatmap with Dendrogram and Highlights
#'
#' This function plots a heatmap of differential ORF usage (DOU) estimates,
#' with hierarchical clustering and optional tile highlighting for significant ORFs.
#'
#' @param orderedMatrix Matrix of DOU estimates to plot.
#' @param rowDendClean Dendrogram object for row ordering.
#' @param highlight_df Data frame specifying tiles to highlight (with \code{row}, \code{col}, and \code{color}).
#' @param gene_labels Vector of gene symbols for row labels.
#' @param color_palette Vector of colors for heatmap.
#' @param color_breaks Numeric vector of breaks for color scaling.
#' @param abs_max Maximum absolute value for color scale.
#' 
#' @importFrom graphics axis image mtext rect
#' 
.plotHeatmap_internal <- function(orderedMatrix,
                                  rowDendClean,
                                  highlight_df,
                                  gene_labels,
                                  color_palette,
                                  color_breaks,
                                  abs_max) { # width = NULL, height = NULL, output_file = "dou.pdf"
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  n_rows <- nrow(orderedMatrix)
  n_cols <- ncol(orderedMatrix)
  
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
  plot(rowDendClean, horiz = TRUE, axes = FALSE, yaxs = "i")
  
  ## Panel 2: Heatmap
  par(mar = c(5, 0.5, 4, right_margin))
  image(x = 1:ncol(orderedMatrix),
        y = 1:nrow(orderedMatrix),
        z = t(orderedMatrix),
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
  
  axis(1, at = 1:ncol(orderedMatrix), labels = colnames(orderedMatrix), las = 2, cex.axis = cex_col)
  axis(4, at = 1:nrow(orderedMatrix), labels = gene_labels, las = 2, cex.axis = cex_row)
  
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


#' Generate DOU Heatmap from Results
#'
#' @description
#' This wrapper function prepares data and plots a heatmap of differential ORF usage (DOU),
#' comparing mORFs with either uORFs or dORFs. It highlights significant ORFs based on
#' local false sign rate (lfsr) from empirical Bayes shrinkage (via the `ashr` package),
#' and retrieves gene annotations from Ensembl for labeling.
#'
#' @param results A data frame containing DOU estimates, ORF labels, and group IDs.
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
#'
#' @return A data frame containing full gene annotation data retrieved from Ensembl.
#'
#' @details
#' The function filters and reshapes DOU results for heatmap visualization, highlighting ORFs
#' with lfsr < 0.05 as significant. It uses Ensembl BioMart to annotate genes and supports
#' comparison of mORFs with upstream (uORFs) or downstream (dORFs) ORFs. The heatmap helps
#' visualize differential ribosome loading across ORFs within genes.
#'
#' @importFrom biomaRt useMart getBM
#' @importFrom stats setNames
#' @importFrom stats heatmap
#' 
#' @export
#' 
plotHeatmap <- function(results, dou_estimates_col = "PosteriorMean", dou_padj_col = "lfsr", flip_sign = TRUE, orf_type = "uORF", 
                        species_dataset = "hsapiens_gene_ensembl", symbol_col = "hgnc_symbol") {
  paired_data <- citronicData(
    results = results,
    dou_estimates_col = dou_estimates_col,
    dou_padj_col = dou_padj_col,
    flip_sign = flip_sign,
    orf_type = orf_type,
    species_dataset = species_dataset,
    symbol_col = symbol_col
  )
  
  tryCatch({
    .plotHeatmap_internal(
      orderedMatrix = paired_data$orderedMatrix,
      rowDendClean = paired_data$rowDendClean,
      highlight_df = paired_data$highlight_df,
      gene_labels = paired_data$gene_labels,
      color_palette = paired_data$color_palette,
      color_breaks = paired_data$color_breaks,
      abs_max = paired_data$abs_max
    )
  }, error = function(e) {
    message("Try enlarging the plot area: ", e$message)
  })
  
  return(paired_data$genes)
}

