#' Plot DTE vs DOU Comparison and Venn Diagram
#'
#' This function visualizes the relationship between differential translation efficiency (DTE)
#' and differential ORF usage (DOU) using a scatter plot and a Venn diagram.
#' It highlights ORFs that are significant in either or both tests and uses a color-blind friendly palette.
#'
#' @param dot A list containing two data frames: \code{dteRes} and \code{dotRes}, 
#'        each with rownames as ORF IDs and columns including \code{padj}, \code{log2FoldChange}, 
#'        \code{empirical.fdr}, and \code{estimates}.
#'
#' @return Generates two plots: a scatter plot of DTE vs DOU estimates, and a Venn diagram 
#'         showing the overlap of significant ORFs.
#'
#' @importFrom graphics plot points legend
#' @importFrom grDevices adjustcolor
#' @importFrom eulerr euler
#'
#' @export
plotDOT <- function(dot) {
  res <- merge(dot$dteRes, dot$douRes, by = "row.names", all = TRUE)
  
  padj_sig <- !is.na(res$padj) & res$padj < 0.05
  fdr_sig <- !is.na(res$empirical.fdr) & res$empirical.fdr < 0.05
  both_sig <- padj_sig & fdr_sig
  padj_only <- padj_sig & !fdr_sig
  fdr_only <- fdr_sig & !padj_sig
  
  col_dte <- adjustcolor("#0072B2", alpha.f = 0.6)
  col_dou <- adjustcolor("#E69F00", alpha.f = 0.6)
  col_both <- adjustcolor("#CC79A7", alpha.f = 0.6)
  
  plot(res$log2FoldChange, res$estimates,
       col = "grey", pch = ".",
       xlab = "log2 fold-change in ORF TE", ylab = "log-odds change in ORF usage",
       main = "Differential translation")
  
  points(res$log2FoldChange[padj_only], res$estimates[padj_only], col = col_dte, pch = ".")
  points(res$log2FoldChange[both_sig], res$estimates[both_sig], col = col_both, pch = ".")
  points(res$log2FoldChange[fdr_only], res$estimates[fdr_only], col = col_dou, pch = ".")
  
  legend("bottomright", inset = 0.06, legend = c("DTE", "DOU", "Both"),
         col = c("#0072B2", "#E69F00", "#CC79A7"), pch = 16, bty = "n")
  
  # Venn diagram
  padj_set <- rownames(res)[padj_sig]
  fdr_set <- rownames(res)[fdr_sig]
  
  fit <- eulerr::euler(c(
    "DTE" = length(setdiff(padj_set, fdr_set)),
    "DOU" = length(setdiff(fdr_set, padj_set)),
    "Both" = length(intersect(padj_set, fdr_set))
  ))
  
  venn_colors <- c("DTE" = "#0072B2", "DOU" = "#E69F00", "Both" = "#CC79A7")
  
  plot(fit,
       fills = list(fill = venn_colors, alpha = 0.6),
       labels = list(font = 2),
       quantities = TRUE,
       main = "Significant differentially translated ORFs")
}
