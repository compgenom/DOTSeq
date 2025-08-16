#' Prepare data and plot DOT heatmap
#'
#' @param results A data frame containing estimates, labels, and groupID.
#' @param sigRes A data frame containing significant results.
#' @param output_file PDF file name for saving the heatmap.
#' @param species_dataset string, e.g. "hsapiens_gene_ensembl"
#' @param symbol_col Name of the gene symbol column in Ensembl for your species,
#'        e.g. "hgnc_symbol" for human, "mgi_symbol" for mouse
#'
#' @importFrom biomaRt useEnsembl getBM
#'
prepHeatMapData <- function(results, sigRes, species_dataset, symbol_col) {
  results$absEstimates <- abs(results$estimates)
  splitList <- split(results, results$groupID)

  maxList <- lapply(splitList, function(results) {
    sapply(unique(results$labels), function(label) {
      subdf <- results[results$labels == label, ]
      subdf$estimates[which.max(subdf$absEstimates)]
    })
  })

  allLabels <- unique(unlist(lapply(maxList, names)))
  maxDiff <- data.frame(groupID = names(maxList), stringsAsFactors = FALSE)
  for (label in allLabels) {
    maxDiff[[label]] <- sapply(maxList, function(x) if (label %in% names(x)) x[[label]] else NA)
  }
  
  allLabels <- unique(unlist(lapply(maxList, names)))
  maxDiff <- data.frame(groupID = names(maxList), stringsAsFactors = FALSE)
  for (label in allLabels) {
    maxDiff[[label]] <- sapply(maxList, function(x) if (label %in% names(x)) x[[label]] else NA)
  }

  uORFSig <- sigRes[sigRes$labels=="uORF", ]
  mORFSig <- sigRes[sigRes$labels=="mORF", ]
  dORFSig <- sigRes[sigRes$labels=="dORF", ]

  sigMaxDiff <- maxDiff[maxDiff$groupID %in% unique(uORFSig$groupID), ]
  sigMaxDiff$dORF <- NULL

  sigMaxDiffMat <- sigMaxDiff
  rownames(sigMaxDiffMat) <- sigMaxDiffMat$groupID
  sigMaxDiffMat$groupID <- NULL
  sigMaxDiffMat[] <- lapply(sigMaxDiffMat, as.numeric)
  sigMaxDiffMat <- as.matrix(sigMaxDiffMat)
  sigMaxDiffMat <- sigMaxDiffMat[, c("uORF", "mORF")]
  sigMaxDiffMatClean <- sigMaxDiffMat[complete.cases(sigMaxDiffMat), ]

  rowClusterClean <- hclust(dist(sigMaxDiffMatClean))
  rowDendClean <- as.dendrogram(rowClusterClean)

  minVal <- min(sigMaxDiffMatClean, na.rm = TRUE)
  maxVal <- max(sigMaxDiffMatClean, na.rm = TRUE)
  absMax <- max(abs(minVal), abs(maxVal))

  colorBreaks <- seq(-absMax, absMax, length.out = 101)
  colorPalette <- colorRampPalette(c("blue", "white", "red"))(100)

  ordered_matrix <- sigMaxDiffMatClean
  sigResFiltered <- sigRes[
    sigRes$labels %in% colnames(ordered_matrix) &
    sigRes$groupID %in% rownames(ordered_matrix),
    ]
  group_to_row <- setNames(seq_len(nrow(ordered_matrix)), rownames(ordered_matrix))
  label_to_col <- setNames(seq_len(ncol(ordered_matrix)), colnames(ordered_matrix))
  highlight_df <- data.frame(
    row = group_to_row[sigResFiltered$groupID],
    col = label_to_col[sigResFiltered$labels],
    color = "black"
  )
  highlight_df <- highlight_df[complete.cases(highlight_df), ]

  ensembIds <- sub("\\..*", "", rownames(ordered_matrix))
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
    sigMaxDiffMatClean = sigMaxDiffMatClean,
    rowDendClean = rowDendClean,
    highlight_df = highlight_df,
    gene_labels = genesUniqueSorted[[symbol_col]],
    color_palette = colorPalette,
    color_breaks = colorBreaks,
    abs_max = absMax
  ))
}
#' Plot DOT heatmap
#'
#' @param sigMaxDiffMatClean matrix of values to plot
#' @param rowDendClean dendrogram for row ordering
#' @param highlight_df dataframe for highlighting tiles
#' @param gene_labels vector of gene labels
#' @param color_palette vector of colors
#' @param color_breaks breaks for colors
#' @param abs_max maximum absolute value
#' @param output_file filename for PDF output
plotHeatmap <- function(sigMaxDiffMatClean,
                           rowDendClean,
                           highlight_df,
                           gene_labels,
                           color_palette,
                           color_breaks,
                           abs_max,
                           output_file = "dot.pdf") {

  n_rows <- nrow(sigMaxDiffMatClean)
  n_cols <- ncol(sigMaxDiffMatClean)
  
  # Dynamic scaling
  cex_row <- ifelse(n_rows > 50, 0.5, 0.8)
  cex_col <- ifelse(n_cols > 50, 0.5, 0.8)
  right_margin <- max(9.5, min(12, n_rows / 5))
  layout_widths <- c(2, max(6, n_cols / 10))
  key_bottom <- ifelse(n_rows > 50, 0.05, 0.13)
  key_top <- key_bottom + 0.12
  
  pdf(output_file, width = 4, height = 4.5)
  layout(matrix(c(1, 2), nrow = 1), widths = layout_widths)
  
  ## Panel 1: Dendrogram
  par(mar = c(5, 1.5, 4, 0))
  plot(rowDendClean, horiz = TRUE, axes = FALSE, yaxs = "i")
  
  ## Panel 2: Heatmap
  row_order <- order.dendrogram(rowDendClean)
  ordered_matrix <- sigMaxDiffMatClean[row_order, ]
  
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
  
  mtext("Differential TE", side = 3, line = 1.5, cex = 1.2, adj = 0.5)
  
  ## Color key
  par(fig = c(0.55, 1, key_bottom, key_top), new = TRUE, mar = c(0, 0, 2, 4))
  image(z = t(matrix(seq(-abs_max, abs_max, length.out = 100), nrow = 1)),
        col = color_palette, axes = FALSE)
  
  tick_vals <- pretty(c(-abs_max, abs_max), n = 5)
  tick_pos <- (tick_vals - (-abs_max)) / (2 * abs_max)
  axis(1, at = tick_pos, labels = tick_vals, las = 1, cex.axis = 0.7)
  mtext("Log-odds", side = 1, line = 2, cex = 0.8)
  
  dev.off()
  par(mfrow = c(1, 1))
}

#' Run DOT heatmap wrapper
#'
#' @param results A data frame containing estimates, labels, and groupID
#' @param sigRes A data frame containing significant results
#' @param species_dataset Ensembl dataset string, e.g. "hsapiens_gene_ensembl"
#' @param symbol_col Name of the gene symbol column in Ensembl for your species,
#'        e.g. "hgnc_symbol" for human, "mgi_symbol" for mouse
#' @param output_file PDF file name
#'
#' @export
runDotHeatmap <- function(results, sigRes, species_dataset, symbol_col, output_file) {
  prep_out <- prepHeatMapData(
    results = results,
    sigRes = sigRes,
    species_dataset = species_dataset,
    symbol_col = symbol_col
  )
  
  plotHeatmap(
    sigMaxDiffMatClean = prep_out$sigMaxDiffMatClean,
    rowDendClean = prep_out$rowDendClean,
    highlight_df = prep_out$highlight_df,
    gene_labels = prep_out$gene_labels,
    color_palette = prep_out$color_palette,
    color_breaks = prep_out$color_breaks,
    abs_max = prep_out$abs_max,
    output_file = output_file
  )
}
