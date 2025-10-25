#' DOTSeqDataSet-class
#'
#' A \code{RangedSummarizedExperiment} object extended to store modeling 
#' components for Differential ORF Usage (DOU) analysis in the DOTSeq 
#' framework.
#'
#' This class contains raw counts, sample metadata, and additional slots 
#' for the model formula, \code{\link[emmeans]{emmeans}} specifications, 
#' and contrast results. It supports flexible modeling of 
#' translation-specific effects using GLM / GLMM and post hoc contrasts.
#'
#' @slot formula A \code{formula} object specifying the model design 
#' (e.g., ~ condition * strategy).
#' @slot specs A \code{formula} object used to generate 
#' \code{\link[emmeans]{emmeans}} specifications for post hoc contrasts.
#' @slot interactionResults A \code{DFrame} or \code{data.frame} containing 
#' results interaction-specific contrasts
#' @slot strategyResults A \code{DFrame} or \code{data.frame} containing 
#' results strategy-specific contrasts
#'
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}}, 
#' \code{\link[glmmTMB]{glmmTMB}}, \code{\link[emmeans]{emmeans}}
#' 
#' @name DOTSeqDataSet-class
#' @rdname DOTSeqDataSet
#' @export
#' 
setClass(
    "DOTSeqDataSet",
    contains = "RangedSummarizedExperiment",
    slots = list(
        formula = "formula",
        specs = "formula",
        interactionResults = "DFrame",
        strategyResults = "DFrame"
    )
)

#' Validity method for DOTSeqDataSet
#'
#' Checks that the slots in a DOTSeqDataSet object are correctly populated.
#'
#' @param object A \code{DOTSeqDataSet} object.
#' @return TRUE if valid, or a character string describing the error.
#' @name DOTSeqDataSet-validity
#' @rdname DOTSeqDataSet
#' @importFrom SummarizedExperiment assayNames assay
setValidity("DOTSeqDataSet", function(object) {
    
    # Check counts
    if (!("counts" %in% assayNames(object))) {
        return("The assays slot must contain a matrix named 'counts'.")
    }
    
    cnt <- assay(object)
    if (!is.numeric(cnt)) return("The count data must be numeric.")
    if (any(is.na(cnt))) return("NA values are not allowed in the count matrix.")
    if (any(cnt < 0)) return("The count data contains negative values.")
    
    fmla <- conditionalFormula(object)
    specs <- emmSpecs(object)
    interaction_df <- interactionResults(object)
    strategy_df <- strategyResults(object)
    
    if (!is(fmla, "formula")) return("The 'formula' slot must be a formula object.")
    if (!is(specs, "formula")) return("The 'specs' slot must be a formula object.")
    if (!(is(interaction_df, "DFrame")) || !(is(strategy_df, "DFrame"))) {
        return("The 'interactionResults' and 'strategyResults' slots must be a DFrame object.")
    }
    
    # Validate formula variables
    fmla_vars <- all.vars(fmla)
    if (!all(fmla_vars %in% names(colData(object)))) {
        return("All variables in the formula must be columns in colData.")
    }
    
    var_classes <- vapply(
        fmla_vars, 
        function(v) class(colData(object)[[v]])[1],
        character(1)
    )
    if (any(var_classes == "character")) {
        return(
            "Variables in the formula are character vectors. ", 
            "Convert them to factors before use."
        )
    }
    
    factor_vars <- fmla_vars[var_classes == "factor"]
    has_duplicate_levels <- vapply(
        factor_vars, 
        function(v) {
            lvls <- levels(colData(object)[[v]])
            any(duplicated(make.names(lvls)))
        }, 
        logical(1)
    )
    if (any(has_duplicate_levels)) {
        return(
            "Factor levels in the formula have non-unique names ", 
            "after applying make.names()."
        )
    }
    
    has_unsafe_levels <- vapply(
        factor_vars, 
        function(v) {
            lvls <- levels(colData(object)[[v]])
            any(!grepl("^[A-Za-z0-9_.]+$", lvls))
        }, 
        logical(1)
    )
    if (any(has_unsafe_levels)) {
        message(
            "Factor levels contain characters other than letters, numbers, ", 
            "'_' or '.'. These may cause issues in downstream modeling."
        )
    }
    
    TRUE
})


#' Access the model formula used in DOTSeq analysis
#'
#' Retrieves the conditional formula stored in a \code{DOTSeqDataSet} 
#' or \code{DESeqDataSet} object. This formula defines the fixed effects 
#' structure used in GLM / GLMM modeling.
#'
#' @param object A \code{DOTSeqDataSet} or \code{DESeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname conditionalFormula
#' @export
setGeneric("conditionalFormula", function(object) standardGeneric("conditionalFormula"))

#' @rdname conditionalFormula
#' @export
setMethod("conditionalFormula", "DOTSeqDataSet", function(object) object@formula)

#' @rdname conditionalFormula
#' @export
setMethod("conditionalFormula", "DESeqDataSet", function(object) {
    if (!"formula" %in% names(metadata(object))) {
        stop("No 'formula' found in metadata.")
    }
    metadata(object)$formula
})


#' Access the emmeans specification formula
#'
#' Retrieves the \code{\link[emmeans]{emmeans}} specification formula 
#' used for post hoc contrast generation.
#'
#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname emmSpecs
#' @export
setGeneric("emmSpecs", function(object) standardGeneric("emmSpecs"))

#' @rdname emmSpecs
#' @export
setMethod("emmSpecs", "DOTSeqDataSet", function(object) object@specs)

#' @rdname emmSpecs
#' @export
setMethod("emmSpecs", "DESeqDataSet", function(object) {
    if (!"specs" %in% names(metadata(object))) {
        stop("No 'specs' found in metadata.")
    }
    metadata(object)$specs
})


#' Access `interactionResults``` from post hoc analysis
#'
#' Retrieves the contrast results stored in a \code{DOTSeqDataSet} object.
#' These results may include shrinkage-adjusted estimates and contrast 
#' statistics
#'
#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname interactionResults
#' @export
setGeneric("interactionResults", function(object) standardGeneric("interactionResults"))

#' @rdname interactionResults
#' @export
setMethod("interactionResults", "DOTSeqDataSet", function(object) object@interactionResults)

#' @rdname interactionResults
#' @export
setMethod("interactionResults", "DESeqDataSet", function(object) {
    if (!"interaction_results" %in% names(metadata(object))) {
        stop("No 'interaction_results' found in metadata.")
    }
    metadata(object)$interaction_results
})

#' @rdname interactionResults
#' @export
setGeneric("interactionResults<-", function(object, value) standardGeneric("interactionResults<-"))

#' @importFrom methods validObject
#' @rdname interactionResults
#' @export
setReplaceMethod("interactionResults", "DOTSeqDataSet", function(object, value) {
    object@interactionResults <- value
    validObject(object)
    return(object)
})

#' @rdname interactionResults
#' @export
setReplaceMethod("interactionResults", "DESeqDataSet", function(object, value) {
    metadata(object)$interaction_results <- value
    return(object)
})


#' Access `strategyResults``` from post hoc analysis
#'
#' Retrieves the contrast results stored in a \code{DOTSeqDataSet} object.
#' These results may include shrinkage-adjusted estimates and contrast 
#' statistics
#'
#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname strategyResults
#' @export
setGeneric("strategyResults", function(object) standardGeneric("strategyResults"))

#' @rdname strategyResults
#' @export
setMethod("strategyResults", "DOTSeqDataSet", function(object) object@strategyResults)

#' @rdname strategyResults
#' @export
setMethod("strategyResults", "DESeqDataSet", function(object) {
    if (!"strategy_results" %in% names(metadata(object))) {
        stop("No 'strategy_results' found in metadata.")
    }
    metadata(object)$strategy_results
})

#' @rdname strategyResults
#' @export
setGeneric("strategyResults<-", function(object, value) standardGeneric("strategyResults<-"))

#' @importFrom methods validObject
#' @rdname strategyResults
#' @export
setReplaceMethod("strategyResults", "DOTSeqDataSet", function(object, value) {
    object@strategyResults <- value
    validObject(object)
    return(object)
})

#' @rdname strategyResults
#' @export
setReplaceMethod("strategyResults", "DESeqDataSet", function(object, value) {
    metadata(object)$strategy_results <- value
    return(object)
})


#' DOTSeqResults-class
#'
#' A wrapper class to store both DOU and DTE results from DOTSeq analysis.
#'
#' @slot DOU A \code{DOTSeqDataSet} object.
#' @slot DTE A \code{DESeqDataSet} object.
#' @rdname DOTSeqResults
#' @export
setClass(
    "DOTSeqResults",
    slots = list(
        DOU = "DOTSeqDataSet",
        DTE = "DESeqDataSet"
    )
)

#' @rdname DOTSeqResults
#' @export
setGeneric("getDOU", function(object) standardGeneric("getDOU"))

#' @rdname DOTSeqResults
#' @export
setMethod("getDOU", "DOTSeqResults", function(object) object@DOU)

#' @rdname DOTSeqResults
#' @export
setGeneric("getDOU<-", function(object, value) standardGeneric("getDOU<-"))

#' @importFrom methods validObject
#' @rdname DOTSeqResults
#' @export
setReplaceMethod("getDOU", "DOTSeqResults", function(object, value) {
    object@DOU <- value
    validObject(object)
    return(object)
})

#' @rdname DOTSeqResults
#' @export
setGeneric("getDTE", function(object) standardGeneric("getDTE"))

#' @rdname DOTSeqResults
#' @export
setMethod("getDTE", "DOTSeqResults", function(object) object@DTE)

#' @rdname DOTSeqResults
#' @export
setGeneric("getDTE<-", function(object, value) standardGeneric("getDTE<-"))

#' @importFrom methods validObject
#' @rdname DOTSeqResults
#' @export
setReplaceMethod("getDTE", "DOTSeqResults", function(object, value) {
    object@DTE <- value
    validObject(object)
    return(object)
})

#' @rdname DOTSeqResults
#' @export
setMethod("show", "DOTSeqResults", function(object) {
    message("DOTSeqResults")
    message("  DOU: ", class(object@DOU))
    message("  DTE: ", class(object@DTE))
    message("Use getDOU(), getDTE(), or interactionResults() to access contents.")
})

#' @rdname DOTSeqResults
#' @export
setMethod("interactionResults", "DOTSeqResults", function(object) {
    list(
        DOU = interactionResults(object@DOU),
        DTE = interactionResults(object@DTE)
    )
})


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
#'     objects (e.g., from \code{\link[emmeans]{emmeans}}).
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
#'     (e.g., from \code{\link[emmeans]{emmeans}}). Default is 
#'     \code{NA_real_}.
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
#' @description 
#' Retrieves the model type string from a \code{PostHoc} object.
#' @param object A \code{PostHoc} object.
#' @return A \code{character(1)} string indicating the model type.
#' @rdname PostHoc-accessors
#' @export
#' @examples
#' ph <- PostHoc(type = "glmmTMB")
#' modelType(ph)
setGeneric("modelType", function(object) standardGeneric("modelType"))

#' @describeIn PostHoc-accessors Access the model type from a PostHoc object.
#' @export
setMethod("modelType", "PostHoc", function(object) object@type)


#' @title Access the results list from a PostHoc object
#' @description 
#' Retrieves model results, parameters, and diagnostics.
#' @param object A \code{PostHoc} object.
#' @return A \code{list} containing model results and diagnostics.
#' @rdname PostHoc-accessors
#' @export
#' @examples 
#' ph <- PostHoc(results = list(aic = 100))
#' fitResults(ph)
setGeneric("fitResults", function(object) standardGeneric("fitResults"))

#' @describeIn PostHoc-accessors Access the results list.
#' @export
setMethod("fitResults", "PostHoc", function(object) object@results)


#' @title Access the post hoc summary from a PostHoc object
#' @description 
#' Retrieves the post hoc summary object (e.g. from \code{\link[emmeans]{emmeans}}).
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
