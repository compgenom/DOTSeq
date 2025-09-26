#' Plot DTE vs DOU Comparison and Venn Diagram
#'
#' This function visualizes the relationship between differential translation efficiency (DTE)
#' and differential ORF usage (DOU) using a scatter plot and a Venn diagram.
#' It highlights ORFs that are significant in either or both tests, using a color-blind friendly palette.
#'
#' @param res A data frame containing results from DTE and DOU analyses.
#'   Must include columns: \code{padj}, \code{lfsr}, and numeric columns for log2 fold-change and DOU estimates.
#' @param estimates_col Character string specifying the column name for DOU estimates (default: \code{"shrunkBeta"}).
#' @param lfc_col Character string specifying the column name for log2 fold-change values from DTE (default: \code{"log2FoldChange"}).
#' @param flip_sign Logical; if \code{TRUE}, flips the sign of the DOU estimates (default: \code{FALSE}).
#' @param lhist Integer; number of bins to use for the marginal histograms (default: \code{20}).
#'
#' @return Generates two plots:
#' \describe{
#'   \item{Scatter plot}{Shows DTE (log2 fold-change) vs DOU (log-odds change) estimates, with significance categories.}
#'   \item{Venn diagram}{(Optional, not currently implemented in this version) showing overlap of significant ORFs.}
#' }
#'
#' @importFrom graphics plot points legend par layout barplot plot.new hist
#' @importFrom grDevices adjustcolor
#' @importFrom stats density
#'
#' @export
plotDOT <- function(res, estimates_col = "shrunkBeta", lfc_col = "log2FoldChange", flip_sign = FALSE, lhist=20) {
  if (isTRUE(flip_sign)) {
    res[[estimates_col]] <- -res[[estimates_col]]
  }
  
  padj_sig <- !is.na(res$padj) & res$padj < 0.05
  fdr_sig <- !is.na(res$lfsr) & res$lfsr < 0.05
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
    
    xlim <- range(res[[lfc_col]])
    ylim <- range(res[[estimates_col]])
    
    # Top histogram
    xhist <- hist(res[[lfc_col]], breaks=seq(xlim[1], xlim[2], length.out=lhist), plot=FALSE)
    par(mar=c(0, pext, 0, 0))
    barplot(xhist$density, axes = FALSE, ylim=c(0, max(xhist$density)), space=0,
         col = "gray", border = "black")
    
    # Right histogram (y-axis)
    yhist <- hist(res[[estimates_col]], breaks=seq(ylim[1], ylim[2], length.out=lhist), plot=FALSE)
    par(mar=c(pext, 0, 0, 0))
    barplot(yhist$density, axes = FALSE, xlim=c(0, max(yhist$density)), space=0,
            col = "gray", border = "black", horiz = TRUE)
    
    # placeholder
    dx <- density(res[[lfc_col]])
    dy <- density(res[[estimates_col]])
    par(mar=c(0, 0, 0, 0))
    plot.new()
    
    # Main scatter plot
    par(mar=c(pext, pext, 0, 0))
    plot(res[[lfc_col]], res[[estimates_col]],
         col = col_none, pch = 16,
         xlab = "log2 fold-change in ORF TE",
         ylab = "log-odds change in ORF usage")
    
    points(res[[lfc_col]][padj_only], res[[estimates_col]][padj_only], col = col_dte)
    points(res[[lfc_col]][fdr_only], res[[estimates_col]][fdr_only], col = col_dou)
    points(res[[lfc_col]][both_sig], res[[estimates_col]][both_sig], col = col_both)
    
    legend("bottomright", legend = c("DTE", "DOU", "Both"), col = c(col_dte, col_dou, col_both), pch = 1, bty = "n", inset = 0.05)
    
  }, error = function(e) {
    par(old_par)
    stop(e)
  })
  
  par(old_par)
}
