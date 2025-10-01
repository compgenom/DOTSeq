#' Plot DTE vs DOU Comparison with Marginal Histograms
#'
#' @description
#' Visualizes the relationship between differential translation efficiency (DTE)
#' and differential ORF usage (DOU) using a scatter plot with marginal histograms.
#' ORFs are colored based on significance in DTE, DOU, or both tests. 
#' This plot helps assess the overlap and divergence between DTE and DOU signals.
#'
#' @param res A data frame containing results from DTE and DOU analyses.
#'   Must include columns for DTE log2 fold-change, DOU estimates, and adjusted p-values
#'   (e.g., \code{padj} for DTE and \code{tests.padj} for DOU).
#' @param dou_estimates_col Character string specifying the column name for DOU estimates.
#'   Default is \code{"PosteriorMean"}.
#' @param dou_padj_col Character string specifying the column name for DOU significance values.
#'   Should correspond to local false sign rate (lfsr). Default is \code{"lfsr"}.
#' @param dte_estimates_col Character string specifying the column name for DTE log2 fold-change.
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
plotDOT <- function(res, dou_estimates_col = "PosteriorMean", dou_padj_col = "lfsr", dte_estimates_col = "log2FoldChange", dte_padj_col = "padj", flip_sign = TRUE, lhist=20) {
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
