#' Prepare Data and Metadata for DOU Heatmap
#'
#' This function prepares a matrix of DOU estimates for heatmap visualization,
#' comparing mORFs with either uORFs or dORFs. It filters significant ORFs,
#' clusters rows, and retrieves gene annotations from Ensembl.
#'
#' @param results A data frame from \code{testDOT} containing DOU estimates, gene labels, and group IDs.
#' @param padj_col Column name for adjusted p-values or local FDR (default: \code{"lfsr"}).
#' @param estimates_col Column name for DOU estimates (default: \code{"shrunkBeta"}).
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of DOU estimates (default: \code{TRUE}).
#' @param orf_type Character string specifying the ORF type to compare with mORFs. Accepts \code{"uORF"} or \code{"dORF"}.
#' @param species_dataset Ensembl dataset name, e.g., \code{"hsapiens_gene_ensembl"}.
#' @param symbol_col Column name for gene symbols in Ensembl, e.g., \code{"hgnc_symbol"} or \code{"mgi_symbol"}.
#' @param threshold Numeric threshold for significance filtering (default: \code{0.05}).
#'
#' @return A list containing:
#' \describe{
#'   \item{orderedMatrix}{Matrix of DOU estimates for significant ORFs.}
#'   \item{rowDendClean}{Dendrogram object for row clustering.}
#'   \item{highlight_df}{Data frame indicating significant tiles to highlight.}
#'   \item{genes}{Full gene annotation data from Ensembl.}
#'   \item{gene_labels}{Vector of gene symbols for heatmap rows.}
#'   \item{color_palette}{Color palette for heatmap.}
#'   \item{color_breaks}{Breaks used for color scaling.}
#'   \item{abs_max}{Maximum absolute value for color scaling.}
#' }
#'
#' @importFrom biomaRt useEnsembl getBM
pairedData <- function(results, padj_col = "lfsr", estimates_col = "shrunkBeta", flip_sign = TRUE, orf_type = "uORF", 
                       species_dataset = "hsapiens_gene_ensembl", symbol_col = "hgnc_symbol", threshold = 0.05) {
  sigRes <- results[results[[padj_col]]<threshold,]
  
  orfSig <- sigRes[sigRes$labels==orf_type, ][, c("groupID", estimates_col)]
  mORFRes <- results[results$labels=="mORF", ][, c("groupID", estimates_col)]
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
#' This wrapper function prepares data and plots a heatmap of differential ORF usage (DOU),
#' comparing mORFs with either uORFs or dORFs. It retrieves gene annotations from Ensembl
#' and highlights significant ORFs.
#'
#' @param results A data frame containing DOU estimates, labels, and group IDs.
#' @param padj_col Column name for adjusted p-values or local FDR (default: \code{"lfsr"}).
#' @param estimates_col Column name for DOU estimates (default: \code{"shrunkBeta"}).
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of DOU estimates (default: \code{TRUE}).
#' @param orf_type Character string specifying the ORF type to compare with mORFs. Accepts \code{"uORF"} or \code{"dORF"}.
#' @param species_dataset Ensembl dataset name, e.g., \code{"hsapiens_gene_ensembl"}.
#' @param symbol_col Column name for gene symbols in Ensembl, e.g., \code{"hgnc_symbol"} or \code{"mgi_symbol"}.
#' @param output_file Optional filename for PDF output (currently unused).
#'
#' @export
plotHeatmap <- function(results = df, padj_col = "lfsr", estimates_col = "shrunkBeta", flip_sign = TRUE, orf_type = "uORF", 
                       species_dataset = "hsapiens_gene_ensembl", symbol_col = "hgnc_symbol", output_file) {
  prep_out <- pairedData(
    results = results,
    estimates_col = estimates_col,
    flip_sign = flip_sign,
    padj_col = padj_col,
    orf_type = orf_type,
    species_dataset = species_dataset,
    symbol_col = symbol_col
  )
  
  tryCatch({
    .plotHeatmap_internal(
      orderedMatrix = prep_out$orderedMatrix,
      rowDendClean = prep_out$rowDendClean,
      highlight_df = prep_out$highlight_df,
      gene_labels = prep_out$gene_labels,
      color_palette = prep_out$color_palette,
      color_breaks = prep_out$color_breaks,
      abs_max = prep_out$abs_max
    )
  }, error = function(e) {
    message("Try enlarging the plot area: ", e$message)
  })
}
