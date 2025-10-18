#' @title The PostHoc class for DOTSeq
#'
#' @description
#' The \code{PostHoc} class represents post hoc summaries derived from a
#' beta-binomial GLM/GLMM fitted to ORF-level data. It is also used to store
#' diagnostics and metadata for each ORF.
#'
#' Objects of this class are typically created by the user-level function
#' \code{fitDOU()}, or manually using the \code{PostHoc()} constructor. In the
#' DOTSeq pipeline, each ORF is assigned a \code{PostHoc} object, which is
#' stored in a \code{DataFrame} and embedded in the \code{rowData} slot of a
#' \code{SummarizedExperiment}.
#'
#' @slot type A \code{character(1)} string indicating the model type.
#'     Default is \code{"fitError"}. If the model is successfully fitted,
#'     the type is typically \code{"glmmTMB"}.
#'
#' @slot results A \code{list} containing model results, parameters, and
#'     test statistics.
#'
#' @slot posthoc An object of class \code{ANY} storing post hoc summary
#'     objects (e.g., from \code{emmeans}).
#'
#' @export
#'
#' @examples
#' ## Create a dummy PostHoc object
#' PostHocRes <- PostHoc(
#'     type = "glmmTMB",
#'     results = list(
#'         model_fit = list(aic = 96.03),
#'         tests = list(pvalue_best = 0.355)
#'     )
#' )
#' PostHocRes
#'
.PostHoc <- setClass(
    "PostHoc",
    slots = c(
        type = "character",
        results = "list",
        posthoc = "ANY"
    )
)


#' @title Construct a PostHoc Object
#'
#' @description
#' This function constructs a new \code{PostHoc} object, which stores the
#' results of a statistical model fitted to ORF-level data. It is typically
#' used internally by the DOTSeq pipeline or manually for testing and
#' diagnostics.
#'
#' @param type A \code{character(1)} string indicating the model type.
#'     Default is \code{"fitError"}.
#'
#' @param results A \code{list} containing model results, parameters, and
#'     test statistics.
#'
#' @param posthoc An optional object storing post hoc summary objects
#'     (e.g., from \code{emmeans}). Default is \code{NA_real_}.
#'
#' @return A \code{PostHoc} S4 object.
#'
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#' ## Create a dummy PostHoc object
#' PostHocRes <- PostHoc(
#'     type = "glmmTMB",
#'     results = list(x = 3, y = 7, b = 4)
#' )
#' PostHocRes
#'
PostHoc <- function(
    type = "fitError",
    results = list(),
    posthoc = NA_real_
) {
    
    out <- new("PostHoc")
    out@type <- type
    out@results <- results
    out@posthoc <- posthoc
    return(out)
}


#' @title Access the model type from a PostHoc object
#' @description Retrieves the model type string from a \code{PostHoc} object.
#' @param object A \code{PostHoc} object.
#' @return A \code{character(1)} string indicating the model type.
#' @rdname PostHoc-accessors
#' @export
#' @examples
#' ph <- PostHoc(type = "glmmTMB")
#' type(ph)
setGeneric("type", function(object) standardGeneric("type"))

#' @describeIn PostHoc-accessors Access the model type from a PostHoc object.
#' @export
setMethod("type", "PostHoc", function(object) object@type)


#' @title Access the results list from a PostHoc object
#' @description Retrieves model results, parameters, and diagnostics.
#' @param object A \code{PostHoc} object.
#' @return A \code{list} containing model results and diagnostics.
#' @rdname PostHoc-accessors
#' @export
#' @examples 
#' ph <- PostHoc(results = list(aic = 100))
#' results(ph)
setGeneric("results", function(object) standardGeneric("results"))

#' @describeIn PostHoc-accessors Access the results list.
#' @export
setMethod("results", "PostHoc", function(object) object@results)


#' @title Access the post hoc summary from a PostHoc object
#' @description Retrieves the post hoc summary object (e.g. from \code{emmeans}).
#' @param object A \code{PostHoc} object.
#' @return A post hoc summary object.
#' @rdname PostHoc-accessors
#' @export
#' @examples
#' ph <- PostHoc(posthoc = "dummy_emmeans")
#' posthoc(ph)
setGeneric("posthoc", function(object) standardGeneric("posthoc"))

#' @describeIn PostHoc-accessors Access the post hoc summary.
#' @export
setMethod("posthoc", "PostHoc", function(object) object@posthoc)
