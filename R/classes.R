#' @title The PostHoc class for DOTSeq
#'
#' @description
#' The \code{PostHoc} class represents post hoc summaries derived from a beta-binomial GLM/GLMM 
#' fitted to ORF-level data. It is also used to store diagnostics and metadata for each ORF.
#'
#' Objects of this class are typically created by the user-level function \code{fitDOU()},
#' or manually using the \code{PostHoc()} constructor. In the DOTSeq pipeline,
#' each ORF is assigned a \code{PostHoc} object, which is stored in a \code{DataFrame}
#' and embedded in the \code{rowData} slot of a \code{SummarizedExperiment}.
#'
#' @slot type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#'   If the model is successfully fitted, the type is typically \code{"glmmTMB"}.
#' @slot results A \code{list} containing model results, parameters, and test statistics.
#' @slot posthoc An object of class \code{ANY} storing post hoc summary objects (e.g., from \code{emmeans}).
#'
#' @export
#'
#' @examples
#' ## Create a dummy PostHoc object
#' PostHocRes <- PostHoc(
#'   type = "glmmTMB",
#'   results = list(
#'     model_fit = list(aic = 96.03),
#'     tests = list(pvalue_best = 0.355)
#'   )
#' )
#' PostHocRes
.PostHoc <- setClass("PostHoc", 
                     slots = c(
                       type = "character",
                       results = "list",
                       posthoc = "ANY"
                       )
                     )


#' @title Construct a PostHoc Object
#'
#' @description
#' This function constructs a new \code{PostHoc} object, which stores the results
#' of a statistical model fitted to ORF-level data. It is typically used internally
#' by the DOTSeq pipeline or manually for testing and diagnostics.
#'
#' @param type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#' @param results A \code{list} containing model results, parameters, and test statistics.
#' @param posthoc An optional object storing post hoc summary objects (e.g., from \code{emmeans}). Default is \code{NA_real_}.
#'
#' @return A \code{PostHoc} S4 object.
#'
#' @examples
#' ## Create a dummy PostHoc object
#' PostHocRes <- PostHoc(
#'   type = "glmmTMB",
#'   results = list(x = 3, y = 7, b = 4)
#' )
#' PostHocRes
#'
#' @importFrom methods new
#' @export
PostHoc <- function(type = "fitError",
                      results = list(),
                      posthoc = NA_real_) {
  out <- new("PostHoc")
  out@type <- type
  out@results <- results
  out@posthoc <- posthoc
  return(out)
}