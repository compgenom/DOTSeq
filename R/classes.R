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
#' @param object A \code{DOTSeqDataSet} object.
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
#' @examples
#' # Read in count matrix, condition table, and annotation files
#' dir <- system.file("extdata", package = "DOTSeq")
#'
#' cnt <- read.table(
#'     file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
#'     header = TRUE,
#'     comment.char = "#"
#' )
#' names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
#'
#' flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#'
#' # Create a DOTSeqObjects object
#' dot <- DOTSeqDataSet(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     bed = bed
#' )
#' 
#' getDOU(dot)
#' 
#' getDTE(dot)
#' 
#' fmla(getDOU(dot))
#' 
#' specs(getDOU(dot))
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
    
    fmla <- fmla(object)
    specs <- specs(object)
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
#' @rdname DOTSeqDataSet
#' @export
setGeneric("fmla", function(object) standardGeneric("fmla"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("fmla", "DOTSeqDataSet", function(object) object@formula)

#' @param object A \code{DESeqDataSet} object.
#' @return A \code{formula} object.
#' @importFrom DESeq2 design
#' @rdname DOTSeqDataSet
#' @export
setMethod("fmla", "DESeqDataSet", function(object) design(object))


#' Access the emmeans specification formula
#'
#' Retrieves the \code{\link[emmeans]{emmeans}} specification formula 
#' used for post hoc contrast generation.
#'
#' @param object A \code{DOTSeqDataSet} or \code{DOTSeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname DOTSeqDataSet
#' @export
setGeneric("specs", function(object) standardGeneric("specs"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("specs", "DOTSeqDataSet", function(object) object@specs)

#' @param object A \code{DESeqDataSet} object.
#' @return A \code{formula} object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("specs", "DESeqDataSet", function(object) {
    if (!"specs" %in% names(metadata(object))) {
        stop("No 'specs' found in metadata.")
    }
    return(metadata(object)$specs)
})


#' Access `interactionResults` from post hoc analysis
#'
#' Retrieves the contrast results stored in a \code{DOTSeqDataSet} object.
#' These results may include shrinkage-adjusted estimates and contrast 
#' statistics for interaction-specific contrasts.
#'
#' @param object A \code{DOTSeqDataSet} or \code{DESeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setGeneric("interactionResults", function(object) standardGeneric("interactionResults"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("interactionResults", "DOTSeqDataSet", function(object) object@interactionResults)

#' @param object A \code{DESeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("interactionResults", "DESeqDataSet", function(object) {
    if (!"interaction_results" %in% names(metadata(object))) {
        stop("No 'interaction_results' found in metadata.")
    }
    metadata(object)$interaction_results
})

#' @param object A \code{DOTSeqDataSet} or \code{DESeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @param value A replacement object (e.g., a \code{DOTSeqDataSet}).
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setGeneric("interactionResults<-", function(object, value) standardGeneric("interactionResults<-"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @param value A replacement object (e.g., a \code{DOTSeqDataSet}).
#' @importFrom methods validObject
#' @rdname DOTSeqDataSet
#' @export
setReplaceMethod("interactionResults", "DOTSeqDataSet", function(object, value) {
    object@interactionResults <- value
    validObject(object)
    return(object)
})

#' @param object A \code{DESeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @param value A replacement object (e.g., a \code{DESeqDataSet}).
#' @rdname DOTSeqDataSet
#' @export
setReplaceMethod("interactionResults", "DESeqDataSet", function(object, value) {
    metadata(object)$interaction_results <- value
    return(object)
})


#' Access `strategyResults``` from post hoc analysis
#'
#' Retrieves the contrast results stored in a \code{DOTSeqDataSet} object.
#' These results may include shrinkage-adjusted estimates and contrast 
#' statistics for strategy-specific contrasts.
#'
#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setGeneric("strategyResults", function(object) standardGeneric("strategyResults"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setMethod("strategyResults", "DOTSeqDataSet", function(object) object@strategyResults)

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @rdname DOTSeqDataSet
#' @export
setGeneric("strategyResults<-", function(object, value) standardGeneric("strategyResults<-"))

#' @param object A \code{DOTSeqDataSet} object.
#' @return A \code{DFrame} or \code{data.frame} containing shrinkage-
#' adjusted estimates and contrast statistics. It is not row-aligned and 
#' will not be subset when subsetting rows of the object.
#' @importFrom methods validObject
#' @rdname DOTSeqDataSet
#' @export
setReplaceMethod("strategyResults", "DOTSeqDataSet", function(object, value) {
    object@strategyResults <- value
    validObject(object)
    return(object)
})


#' DOTSeqObjects-class
#'
#' A wrapper class to store both DOU and DTE results from DOTSeq analysis.
#'
#' @slot DOU A \code{DOTSeqDataSet} object.
#' @slot DTE A \code{DESeqDataSet} object.
#' @rdname DOTSeqObjects
#' @export
#' @examples
#' # Create dummy objects
#' dou <- new("DOTSeqDataSet")
#' dte <- new("DESeqDataSet")
#' res <- new("DOTSeqObjects", DOU = dou, DTE = dte)
#' getDOU(res)
#' getDTE(res)
setClass(
    "DOTSeqObjects",
    slots = list(
        DOU = "DOTSeqDataSet",
        DTE = "DESeqDataSet"
    )
)

#' @param object A \code{DOTSeqObjects} object.
#' @rdname DOTSeqObjects
#' @export
setGeneric("getDOU", function(object) standardGeneric("getDOU"))

#' @param object A \code{DOTSeqObjects} object.
#' @rdname DOTSeqObjects
#' @return A \code{DOTSeqDataSet} object containing DOU results.
#' @export
setMethod("getDOU", "DOTSeqObjects", function(object) object@DOU)

#' @param object A \code{DOTSeqObjects} object.
#' @param value A replacement object (e.g., a \code{DOTSeqDataSet}).
#' @rdname DOTSeqObjects
#' @export
setGeneric("getDOU<-", function(object, value) standardGeneric("getDOU<-"))

#' @param object A \code{DOTSeqObjects} object.
#' @param value A replacement object (e.g., a \code{DOTSeqDataSet}).
#' @importFrom methods validObject
#' @rdname DOTSeqObjects
#' @export
setReplaceMethod("getDOU", "DOTSeqObjects", function(object, value) {
    object@DOU <- value
    validObject(object)
    return(object)
})

#' @param object A \code{DOTSeqObjects} object.
#' @rdname DOTSeqObjects
#' @export
setGeneric("getDTE", function(object) standardGeneric("getDTE"))

#' @param object A \code{DOTSeqObjects} object.
#' @param value A replacement object (e.g., a \code{DESeqDataSet}).
#' @rdname DOTSeqObjects
#' @export
setMethod("getDTE", "DOTSeqObjects", function(object) object@DTE)

#' @param object A \code{DOTSeqObjects} object.
#' @param value A replacement object (e.g., a \code{DESeqDataSet}).
#' @rdname DOTSeqObjects
#' @export
setGeneric("getDTE<-", function(object, value) standardGeneric("getDTE<-"))

#' @param object A \code{DOTSeqObjects} object.
#' @param value A replacement object (e.g., a \code{DOTSeqDataSet}).
#' @importFrom methods validObject
#' @rdname DOTSeqObjects
#' @export
setReplaceMethod("getDTE", "DOTSeqObjects", function(object, value) {
    object@DTE <- value
    validObject(object)
    return(object)
})

#' @param object A \code{DOTSeqObjects} object.
#' @rdname DOTSeqObjects
#' @export
setMethod("show", "DOTSeqObjects", function(object) {
    message("DOTSeqObjects")
    message("  DOU: ", class(object@DOU))
    message("  DTE: ", class(object@DTE))
    message("Use getDOU(), getDTE(), or interactionResults() to access contents.")
})

#' @param object A \code{DOTSeqObjects} object.
#' @rdname DOTSeqObjects
#' @export
setMethod("interactionResults", "DOTSeqObjects", function(object) {
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
