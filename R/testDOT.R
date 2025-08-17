#' Test Differential ORF Translation
#'
#' This function computes differential translation statistics for ORFs
#' using a fitted DOTSeq model (from \code{satuRn::fitDTU}) and contrasts.
#' It calculates log-odds estimates, standard errors, t-statistics, empirical p-values,
#' and adjusted FDR values. Optional diagnostic plots can visualize the z-score distribution
#' and empirical null fitting.
#'
#' @param sumExp A \code{SummarizedExperiment} object returned by \code{satuRn::fitDTU}.
#'   Must contain fitted models in \code{rowData(sumExp)[["fitDTUModels"]]}.
#' @param contrast Numeric vector specifying the contrast of coefficients to test.
#' @param df Numeric degrees of freedom for empirical null fitting (default: 7).
#' @param bre Numeric number of histogram breaks for empirical null estimation (default: 120).
#' @param diagplot1 Logical; if TRUE, produces the first diagnostic plot of z-scores (default: FALSE).
#' @param diagplot2 Logical; if TRUE, produces the second diagnostic plot of empirical z-scores (default: FALSE).
#' @param main Character string to set the title for diagnostic plots (default: NULL).
#' @param seed An optional integer. Use to set the seed for the random number
#'   generator to ensure reproducible results (default: NULL).
#'
#' @return A \code{data.frame} containing the following columns:
#' \describe{
#'   \item{estimates}{Log-odds estimates for the specified contrast.}
#'   \item{se}{Standard errors of the estimates.}
#'   \item{posterior}{Posterior degrees of freedom for each estimate.}
#'   \item{t.statistic}{t-statistic for each estimate.}
#'   \item{pval}{Two-sided p-value based on the t-statistic.}
#'   \item{empirical.pval}{P-values corrected using empirical null estimation.}
#'   \item{empirical.fdr}{FDR-adjusted p-values using Benjamini-Hochberg correction.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- testDOT(sumExp = mySummarizedExp, contrast = c(0,1,0), 
#'                    df = 7, bre = 120, diagplot1 = TRUE, diagplot2 = TRUE, main = "Sample ORFs", seed = 42)
#' head(results)
#' }
#'
#' @export
testDOT <- function(sumExp, contrast, df = 7, bre = 120, diagplot1 = FALSE, diagplot2 = FALSE, main = NULL, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- .GlobalEnv$.Random.seed
    on.exit({
      if (is.null(old_seed)) {
        rm(.Random.seed, envir = .GlobalEnv)
      } else {
        .GlobalEnv$.Random.seed <- old_seed
      }
    })
    
    # Set the seed
    set.seed(seed)
  }
  
  # Input validation checks
  if (!is(sumExp, "SummarizedExperiment")) {
    stop("Input `sumExp_filtered` must be a SummarizedExperiment object.")
  }
  if (is.null(rowData(sumExp)[["fitDTUModels"]])) {
    stop("The 'fitDTUModels' slot is empty. Did you run fitDTU first?")
  }
  if (!is.numeric(contrast) || !is.vector(contrast)) {
    stop("Input `contrast` must be a numeric vector.")
  }
  if (!is.numeric(df) || length(df) != 1 || df < 1) {
    stop("Input `df` must be a single positive numeric value.")
  }
  if (!is.numeric(bre) || length(bre) != 1 || bre < 1) {
    stop("Input `bre` must be a single positive numeric value.")
  }
  
  # Extract the necessary data from the sumExp object
  models <- rowData(sumExp)[["fitDTUModels"]]
  
  # A more robust check for contrast length
  first_model <- models[[which.max(vapply(models, function(x) length(satuRn:::getCoef(x)), numeric(1)))]]
  if(length(contrast) != length(satuRn:::getCoef(first_model))) {
    stop("The length of the contrast vector does not match the number of coefficients in the models.")
  }
  
  # Compute the usage estimate (log-odds)
  estimates <- vapply(models, satuRn:::getEstimates, contrast = contrast, FUN.VALUE = numeric(1))
  # Compute the standard error (SE) for each of the estimates
  se <- sqrt(vapply(models, satuRn:::varContrast, contrast = contrast, FUN.VALUE = numeric(1)))
  # Extract the posterior degrees of freedom for each model, which are necessary for the t-test
  posterior <- unlist(vapply(models, satuRn:::getDfPosterior, FUN.VALUE = numeric(1)))
  # Replace NULL or numeric(0) with NA
  posterior[!(vapply(posterior, length, numeric(1)))] <- NA
  
  # Calculate the t-statistic and convert to z-score
  # Compute the t-statistic for each transcript's contrast
  tstat <- estimates / se
  # Calculate the two-sided p-value for each t-statistic.
  pval <- pt(-abs(tstat), posterior) * 2
  # Convert the calculated p-values back into standard normal (z) scores.
  zvalues <- qnorm(pval / 2) * sign(tstat)
  
  # This code is adapted directly from satuRn's internal functions
  # Filter out extreme z-scores
  zvaluesMid <- zvalues[abs(zvalues) < 10]
  zvaluesMid <- zvaluesMid[!is.na(zvaluesMid)]
  
  par(mfrow = c(1, 2))
  if(diagplot1){
    plot <- locfdr(zz = zvaluesMid, df = df,
                   main = paste0("diagplot 1: ", main))
  }
  
  # Initial MLE estimation
  N <- length(zvaluesMid)
  b <- 4.3 * exp(-0.26 * log(N, 10))
  med <- median(zvaluesMid)
  sc <- diff(quantile(zvaluesMid)[c(2, 4)]) / (2 * qnorm(.75))
  mlests <- satuRn:::locfdr_locmle(zvaluesMid, xlim = c(med, b * sc))
  
  lo <- min(zvaluesMid)
  up <- max(zvaluesMid)
  
  breaks <- seq(lo, up, length = bre)
  zzz <- pmax(pmin(zvaluesMid, up), lo)
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  x <- (breaks[-1] + breaks[-length(breaks)]) / 2
  sw <- 0
  
  X <- cbind(1, poly(x, df = df))
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  y <- zh$counts
  glmfit <- glm(y ~ poly(x, df = df), poisson)
  f <- glmfit$fit
  # This is a proxy for the 'f(z) misfit' check
  # Issue a warning if the residual deviance is very high compared to the degrees of freedom
  if((diagplot1==FALSE) && (glmfit$deviance > glmfit$df.residual * 2)) {
    warning(paste0("GLM fit for empirical null is poor (residual deviance = ", 
                   round(glmfit$deviance, 2), 
                   "). Consider rerunning with a higher 'df' value."))
  }
  
  Cov.in <- list(x = x, X = X, f = f, sw = sw)
  ml.out <- satuRn:::locfdr_locmle(zvaluesMid,
                                   xlim = c(mlests[1], b * mlests[2]),
                                   d = mlests[1], s = mlests[2], Cov.in = Cov.in
  )
  mlests <- ml.out$mle
  
  # Use the derived mlests to correct the full z-score vector
  zval_empirical <- (zvalues - mlests[1]) / mlests[2]
  pval_empirical <- 2 * pnorm(-abs(zval_empirical), mean = 0, sd = 1)
  
  if (diagplot2) {
    zval_empirical_mid <- zval_empirical[abs(zval_empirical) < 10]
    zval_empirical_mid <- zval_empirical_mid[!is.na(zval_empirical_mid)]
    lo <- min(zval_empirical_mid)
    up <- max(zval_empirical_mid)
    
    lo <- min(lo, -1 * up) # to center the figure
    up <- max(up, -1 * lo)
    
    bre <- 120
    breaks <- seq(lo, up, length = bre)
    zzz <- pmax(pmin(zval_empirical_mid, up), lo)
    zh <- hist(zzz, breaks = breaks, plot = FALSE)
    yall <- zh$counts
    K <- length(yall)
    
    hist(zzz, breaks = breaks, xlab = "z-scores", 
         main = paste0("diagplot 2: ", main), freq = FALSE)
    xfit <- seq(min(zzz), max(zzz), length = 4000)
    yfit <- dnorm(xfit / mlests[3], mean = 0, sd = 1)
    lines(xfit, yfit, col = "darkgreen", lwd = 2)
  }
  par(mfrow = c(1, 1))
  
  FDR <- p.adjust(pval_empirical, method = "BH")
  
  # Define thresholds
  fdrThresh <- 0.05
  effectThresh <- 1
  
  # Calculate -log10(FDR)
  logFDR <- -log10(FDR)
  
  # Create logical vectors for significance and effect size
  sig <- FDR < fdrThresh
  posEffect <- estimates > effectThresh
  negEffect <- estimates < -effectThresh
  
  # Plot all points in gray
  plot(estimates, logFDR, pch = 20, col = "gray80",
       main = paste("Volcano Plot:", main),
       xlab = "log-odds change", ylab = "-log10(FDR)")
  points(estimates[sig & posEffect], logFDR[sig & posEffect], pch = 20, col = "red")
  points(estimates[sig & negEffect], logFDR[sig & negEffect], pch = 20, col = "blue")
  # Add threshold lines
  abline(h = -log10(fdrThresh), col = "black", lty = 2)
  abline(v = c(-effectThresh, effectThresh), col = "black", lty = 2)
  
  
  results <- data.frame(estimates = estimates, se = se, posterior = posterior,
                        t.statistic = tstat, pval = pval, 
                        empirical.pval = pval_empirical, empirical.fdr = FDR)
  results <- na.omit(results)
  
  return(results)
}
