#' Prepare Data and Plot DOT Heatmap
#'
#' This function generates a DOT (Differential ORF Translation) heatmap using 
#' estimates from a results data frame, comparing mORFs with either uORFs or dORFs.
#'
#' @param results An output data frame from testDOT containing estimates, gene labels, and groupID.
#' @param orf_type Character string specifying the ORF type to compare with mORF. 
#'        Accepts either "uORF" or "dORF".
#' @param species_dataset A string specifying the Ensembl dataset, e.g., "hsapiens_gene_ensembl".
#' @param symbol_col A string indicating the column name for gene symbols in Ensembl, 
#'        e.g., "hgnc_symbol" for human or "mgi_symbol" for mouse.
#' @param threshold Numeric value specifying the empirical FDR threshold (default = 0.05).
#'
#' @importFrom biomaRt useEnsembl getBM
#'
pairedData <- function(results, orf_type, species_dataset, symbol_col, threshold = 0.05) {
  sigRes <- results[results$empirical.fdr<threshold,]
  
  orfSig <- sigRes[sigRes$labels==orf_type, ][, c("groupID", "estimates")]
  mORFRes <- results[results$labels=="mORF", ][, c("groupID", "estimates")]
  sigMat <- merge(orfSig, mORFRes, by = "groupID")
  names(sigMat) <- c("groupID", orf_type, "mORF")
  sigMat$delta <- abs(sigMat$mORF - sigMat[[orf_type]])
  sigMat <- sigMat[ave(sigMat$delta, sigMat$groupID, FUN = function(x) x == max(x)) == 1, ]
  
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
#' Plot DOT heatmap
#'
#' @param orderedMatrix matrix of values to plot
#' @param rowDendClean dendrogram for row ordering
#' @param highlight_df dataframe for highlighting tiles
#' @param gene_labels vector of gene labels
#' @param color_palette vector of colors
#' @param color_breaks breaks for colors
#' @param abs_max maximum absolute value
#' @param output_file filename for PDF output
plotHeatmap <- function(orderedMatrix,
                        rowDendClean,
                        highlight_df,
                        gene_labels,
                        color_palette,
                        color_breaks,
                        abs_max,
                        width = NULL, height = NULL,
                        output_file = "dot.pdf") {

  n_rows <- nrow(orderedMatrix)
  n_cols <- ncol(orderedMatrix)
  
  # Dynamic scaling
  cex_row <- ifelse(n_rows > 50, 0.5, 0.8)
  cex_col <- ifelse(n_cols > 50, 0.5, 0.8)
  right_margin <- max(9.5, min(12, n_rows / 5))
  layout_widths <- c(2, max(6, n_cols / 10))
  key_bottom <- ifelse(n_rows > 50, 0.05, 0.13)
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
  
  mtext("Differential ORF translation", side = 3, line = 1.5, cex = 1.2, adj = 0.5, font = 2)
  
  ## Color key
  par(fig = c(0.55, 1, key_bottom, key_top), new = TRUE, mar = c(0, 0, 2, 4))
  image(z = t(matrix(seq(-abs_max, abs_max, length.out = 100), nrow = 1)),
        col = color_palette, axes = FALSE)
  
  tick_vals <- pretty(c(-abs_max, abs_max), n = 5)
  tick_pos <- (tick_vals - (-abs_max)) / (2 * abs_max)
  axis(1, at = tick_pos, labels = tick_vals, las = 1, cex.axis = 0.7)
  mtext("log-odds", side = 1, line = 2, cex = 0.8)
  
  # dev.off()
  par(mfrow = c(1, 1))
}

#' Run DOT heatmap wrapper
#'
#' @param results A data frame containing estimates, labels, and groupID
#' @param species_dataset A string specifying the Ensembl dataset, e.g., "hsapiens_gene_ensembl".
#' @param symbol_col A string indicating the column name for gene symbols in Ensembl, 
#'        e.g., "hgnc_symbol" for human or "mgi_symbol" for mouse.
#' @param output_file PDF file name
#'
#' @export
pairedDOTHeatmap <- function(results, orf_type, species_dataset, symbol_col, output_file, width = NULL, height = NULL) {
  prep_out <- pairedData(
    results = results,
    orf_type = orf_type,
    species_dataset = species_dataset,
    symbol_col = symbol_col
  )
  
  plotHeatmap(
    orderedMatrix = prep_out$orderedMatrix,
    rowDendClean = prep_out$rowDendClean,
    highlight_df = prep_out$highlight_df,
    gene_labels = prep_out$gene_labels,
    color_palette = prep_out$color_palette,
    color_breaks = prep_out$color_breaks,
    abs_max = prep_out$abs_max,
    output_file = output_file,
    width = width,
    height = height
  )
}
