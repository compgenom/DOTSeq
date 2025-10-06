#' @title The StatModel class for DOTSeq
#'
#' @description
#' The \code{StatModel} class represents a statistical model fitted to ORF-level data
#' using a beta-binomial GLMM via \code{glmmTMB}. It is used to store model results,
#' diagnostics, and metadata for each ORF.
#'
#' Objects of this class are typically created by the user-level function \code{fitDOU()},
#' or manually using the \code{StatModel()} constructor. In the DOTSeq pipeline,
#' each ORF is assigned a \code{StatModel} object, which is stored in a \code{DataFrame}
#' and embedded in the \code{rowData} slot of a \code{SummarizedExperiment}.
#'
#' @slot type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#'   If the model is successfully fitted, the type is typically \code{"glmmTMB"}.
#' @slot results A \code{list} containing model results, parameters, and test statistics.
#' @slot model An object of class \code{ANY} storing the fitted model object (e.g., from \code{glmmTMB}).
#'
#' @export
#'
#' @examples
#' ## Create a dummy StatModel object
#' myModel <- StatModel(
#'   type = "glmmTMB",
#'   results = list(
#'     model_fit = list(aic = 96.03),
#'     tests = list(pvalue_best = 0.355)
#'   )
#' )
#' myModel
.StatModel <- setClass("StatModel",
                       slots = c(
                         type = "character",
                         results = "list",
                         model = "ANY"
                       )
)


#' @title Construct a StatModel Object
#'
#' @description
#' This function constructs a new \code{StatModel} object, which stores the results
#' of a statistical model fitted to ORF-level data. It is typically used internally
#' by the DOTSeq pipeline or manually for testing and diagnostics.
#'
#' @param type A \code{character(1)} string indicating the model type. Default is \code{"fitError"}.
#' @param results A \code{list} containing model results, parameters, and test statistics.
#' @param model An optional fitted model object (e.g., from \code{glmmTMB}). Default is \code{NA_real_}.
#'
#' @return A \code{StatModel} S4 object.
#'
#' @examples
#' ## Create a dummy StatModel object
#' myModel <- StatModel(
#'   type = "glmmTMB",
#'   results = list(x = 3, y = 7, b = 4)
#' )
#' myModel
#'
#' @importFrom methods new
#' @export
StatModel <- function(type = "fitError",
                      results = list(),
                      model = NA_real_) {
  out <- new("StatModel")
  out@type <- type
  out@results <- results
  out@model <- model
  return(out)
}