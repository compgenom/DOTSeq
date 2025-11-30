#' DOUData-class
#'
#' A \code{RangedSummarizedExperiment} object extended to store modeling 
#' components for Differential ORF Usage (DOU) analysis in the 
#' \pkg{DOTSeq} framework.
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
#' @slot interaction A \code{DFrame} containing results 
#' from interaction-specific contrasts.
#' @slot strategy A \code{DFrame} containing results 
#' from strategy-specific contrasts.
#'
#' @seealso \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}, 
#' \code{\link[glmmTMB]{glmmTMB}}, \code{\link[emmeans]{emmeans}}
#'
#' @name DOUData-class
#' @rdname DOUData-class
#' @aliases DOUData-class
#' @exportClass DOUData
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
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setClass(
    "DOUData",
    contains = "RangedSummarizedExperiment",
    slots = list(
        formula = "formula",
        specs = "formula",
        interaction = "DFrame",
        strategy = "DFrame"
    )
)


#' Validity method for DOUData
#'
#' Checks that the slots in a DOUData object are correctly populated.
#'
#' @param object A \code{\link{DOUData-class}} object.
#' @return TRUE if valid, or a character string describing the error.
#' @rdname DOUData-validity
#' @name DOUData-validity
#' @import SummarizedExperiment
setValidity("DOUData", function(object) {
    
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
    interaction_df <- getContrasts(object, type = "interaction")
    strategy_df <- getContrasts(object, type = "strategy")
    
    if (!is(fmla, "formula")) return("The 'formula' slot must be a formula object.")
    if (!is(specs, "formula")) return("The 'specs' slot must be a formula object.")
    if (!(is(interaction_df, "DFrame")) || !(is(strategy_df, "DFrame"))) {
        return("The 'interaction' and 'strategy' slots must be a DFrame object.")
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


#' DTEData-class
#'
#' A \code{DESeqDataSet} object extended to store modeling 
#' components for Differential Translation Efficiency (DTE) analysis in 
#' the \pkg{DOTSeq} framework.
#' 
#' This class contains raw counts, sample metadata, and additional slots 
#' for the \code{\link[emmeans]{emmeans}} specifications, and contrast 
#' results. 
#' 
#' Inherits all slots from \code{DESeqDataSet}, and adds:
#' - specs: A \code{formula} object used to generate 
#' \code{\link[emmeans]{emmeans}} specifications for post hoc contrasts.
#' - interaction: A \code{DFrame} or \code{data.frame} containing 
#' results interaction-specific contrasts
#' - strategy: A \code{DFrame} or \code{data.frame} containing 
#' results strategy-specific contrasts
#' @seealso \code{\link[DESeq2]{DESeqDataSet-class}}
#' @name DTEData-class
#' @rdname DTEData-class
#' @aliases DTEData-class
#' @exportClass DTEData
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
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setClass(
    "DTEData",
    contains = "DESeqDataSet",
    slots = list(
        formula = "formula",
        specs = "formula",
        interaction = "DFrame",
        strategy = "DFrame"
    )
)


#' @importFrom DESeq2 DESeqDataSet
setClassUnion("DESeqOrDTE", c("DESeqDataSet", "DTEData"))

#' DOTSeqDataSets-class
#'
#' A wrapper class to store both DOU and DTE results from \pkg{DOTSeq} 
#' analysis.
#'
#' @slot DOU A \code{\link{DOUData-class}} object.
#' @slot DTE A \code{\link{DTEData-class}} object.
#' @return A \code{\link{DOTSeqDataSets-class}} S4 object containing DOU and DTE results.
#' @name DOTSeqDataSets-class
#' @rdname DOTSeqDataSets-class
#' @exportClass DOTSeqDataSets
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
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setClass(
    "DOTSeqDataSets",
    slots = list(
        DOU = "DOUData",
        DTE = "DESeqOrDTE"
    )
)

#' Show method for DOTSeqDataSets objects
#'
#' Displays a summary of the DOTSeqDataSets object.
#'
#' @param object A \code{\link{DOTSeqDataSets-class}} object.
#' @name DOTSeqDataSets-class
#' @rdname DOTSeqDataSets-class
#' @aliases show,DOTSeqDataSets-method
#' @export
setMethod("show", "DOTSeqDataSets", function(object) {
    message("DOTSeqDataSets")
    message("  DOU: ", class(object@DOU))
    message("  DTE: ", class(object@DTE))
    message("Use getDOU() or getDTE() to access contents.")
    message("Use getContrasts() if post hoc inference has been performed.")
})


#' Accessor and replacement methods for DOU slot
#'
#' These methods allow access to and replacement of the \code{\link{DOUData-class}} 
#' object stored within a \code{\link{DOTSeqDataSets-class}} container.
#'
#' @param object A \code{\link{DOTSeqDataSets-class}} object.
#' @param value A replacement object (e.g., a \code{\link{DOUData-class}}).
#' @return For the accessor, a \code{\link{DOUData-class}} object. For the replacement, 
#' an updated \code{\link{DOTSeqDataSets-class}} object.
#' @rdname getDOU-method
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
#' gtf <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
#' bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
#'
#' meta <- read.table(file.path(dir, "metadata.txt.gz"))
#' names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
#' cond <- meta[meta$treatment == "chx", ]
#' cond$treatment <- NULL
#'
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = gtf,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setGeneric("getDOU", function(object) standardGeneric("getDOU"))

#' @rdname getDOU-method
#' @export
setMethod("getDOU", "DOTSeqDataSets", function(object) object@DOU)

#' @rdname getDOU-method
#' @export
setGeneric("getDOU<-", function(object, value) standardGeneric("getDOU<-"))

#' @rdname getDOU-method
#' @importFrom methods validObject
#' @export
setReplaceMethod("getDOU", "DOTSeqDataSets", function(object, value) {
    object@DOU <- value
    validObject(object)
    return(object)
})



#' Accessor and replacement methods for DTE slot
#'
#' These methods allow access to and replacement of the \code{\link{DTEData-class}} 
#' object stored within a \code{\link{DOTSeqDataSets-class}} container.
#'
#' @param object A \code{\link{DOTSeqDataSets-class}} object.
#' @param value A replacement object (e.g., a \code{\link{DTEData-class}}).
#' @return For the accessor, a \code{\link{DTEData-class}} object. For the replacement, 
#' an updated \code{\link{DOTSeqDataSets-class}} object.
#' @rdname getDTE-method
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
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setGeneric("getDTE", function(object) standardGeneric("getDTE"))

#' @rdname getDTE-method
#' @export
setMethod("getDTE", "DOTSeqDataSets", function(object) object@DTE)

#' @rdname getDTE-method
#' @export
setGeneric("getDTE<-", function(object, value) standardGeneric("getDTE<-"))

#' @rdname getDTE-method
#' @importFrom methods validObject
#' @export
setReplaceMethod("getDTE", "DOTSeqDataSets", function(object, value) {
    object@DTE <- value
    validObject(object)
    return(object)
})


# setValidity("DTEData", function(object) {
#   if (!inherits(design(object), "formula")) {
#     return("The design must be a formula object.")
#   }
#   TRUE
# })


#' Access the model formula used in \pkg{DOTSeq} analysis
#'
#' Retrieves the conditional formula stored in a 
#' \code{\link{DOUData-class}} or \code{\link{DTEData-class}} object. 
#' This formula defines the fixed effects structure used in GLM / 
#' GLMM modeling.
#'
#' @param object A \code{\link{DOUData-class}} or 
#' \code{\link{DTEData-class}} object.
#' @return A \code{formula} object.
#' @keywords internal
#' @rdname fmla-method
setGeneric("fmla", function(object) standardGeneric("fmla"))

#' @rdname fmla-method
setMethod("fmla", "DOUData", function(object) object@formula)

#' @rdname fmla-method
#' @importFrom DESeq2 design
setMethod("fmla", "DTEData", function(object) design(object))


#' Access the emmeans specification formula
#'
#' Retrieves the \code{\link[emmeans]{emmeans}} specification formula 
#' used for post hoc contrast generation.
#'
#' @param object A \code{\link{DOUData-class}} or \code{\link{DTEData-class}} object.
#' @return A \code{formula} object.
#' @keywords internal
#' @rdname specs-method
setGeneric("specs", function(object) standardGeneric("specs"))

#' @rdname specs-method
setMethod("specs", "DOUData", function(object) object@specs)

#' @rdname specs-method
setMethod("specs", "DTEData", function(object) object@specs)


#' Access and replace contrast results from post hoc analysis
#'
#' These methods retrieve or replace either interaction-specific or 
#' strategy-specific contrast results from a \code{\link{DOUData-class}}, 
#' \code{\link{DTEData-class}}, or \code{\link{DOTSeqDataSets-class}} object.
#'
#' @param object A \code{\link{DOUData-class}} or \code{\link{DTEData-class}} object.
#' @param type A character string, either \code{"interaction"} or 
#' \code{"strategy"}.
#' @param value A \code{DFrame} or \code{data.frame} to replace the 
#' contrast results.
#' @return For the accessor, a \code{DFrame} or \code{data.frame} 
#' containing contrast results. For the replacement, an updated 
#' \code{\link{DOUData-class}}, \code{\link{DTEData-class}} or \code{\link{DOTSeqDataSets-class}} object.
#' @rdname getContrasts-method
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
#' # Create a DOTSeqDataSets object
#' d <- DOTSeqDataSetsFromFeatureCounts(
#'     count_table = cnt,
#'     condition_table = cond,
#'     flattened_gtf = flat,
#'     flattened_bed = bed
#' )
#' 
#' getDOU(d)
#' 
#' getDTE(d)
#' 
setGeneric("getContrasts", function(object, type = c("interaction", "strategy")) standardGeneric("getContrasts"))

#' @rdname getContrasts-method
#' @export
setMethod("getContrasts", "DOUData", function(object, type = c("interaction", "strategy")) {
    type <- match.arg(type)
    if (type == "interaction") object@interaction else object@strategy
})

#' @rdname getContrasts-method
#' @export
setMethod("getContrasts", "DTEData", function(object, type = c("interaction", "strategy")) {
    type <- match.arg(type)
    if (type == "interaction") object@interaction else object@strategy
})

#' @rdname getContrasts-method
#' @export
setMethod("getContrasts", "DOTSeqDataSets", function(object, type = c("interaction", "strategy")) {
    type <- match.arg(type)
    
    dou_res <- getContrasts(object@DOU, type = type)
    dte_res <- getContrasts(object@DTE, type = type)
    
    if (is.null(dou_res)) warning("No DOU contrast results found.")
    if (is.null(dte_res)) warning("No DTE contrast results found.")
    
    return(list(DOU = dou_res, DTE = dte_res))
})

#' @rdname getContrasts-method
#' @export
setGeneric("getContrasts<-", function(object, type = c("interaction", "strategy"), value) standardGeneric("getContrasts<-"))

#' @rdname getContrasts-method
#' @importFrom methods validObject
#' @export
setReplaceMethod("getContrasts", "DOUData", function(object, type = c("interaction", "strategy"), value) {
    type <- match.arg(type)
    if (!is(value, "DFrame")) stop("Replacement value must be a DFrame.")
    if (type == "interaction") object@interaction <- value else object@strategy <- value
    validObject(object)
    return(object)
})

#' @rdname getContrasts-method
#' @importFrom methods validObject
#' @export
setReplaceMethod("getContrasts", "DTEData", function(object, type = c("interaction", "strategy"), value) {
    type <- match.arg(type)
    if (!is(value, "DFrame")) stop("Replacement value must be a DFrame.")
    if (type == "interaction") object@interaction <- value else object@strategy <- value
    validObject(object)
    return(object)
})



#' @title The PostHoc class for DOTSeq
#'
#' @description
#' The \code{PostHoc} class represents post hoc summaries derived from a
#' beta-binomial GLM/GLMM fitted to ORF-level data. It is also used to store
#' diagnostics and metadata for each ORF.
#'
#' Objects of this class are typically created by the user-level function
#' \code{fitDOU()}, or manually using the \code{PostHoc()} constructor. 
#' In the \pkg{DOTSeq} workflow, each ORF is assigned a \code{PostHoc} 
#' object, which is stored in a \code{DataFrame} and embedded in the 
#' \code{rowData} slot of a \code{SummarizedExperiment}.
#'
#' @slot type A \code{character(1)} string indicating the model type.
#' Default is \code{"fitError"}. If the model is successfully fitted,
#' the type is typically \code{"glmmTMB"}.
#'
#' @slot results A \code{list} containing model results, parameters, and
#' test statistics.
#'
#' @slot posthoc An object of class \code{ANY} storing post hoc summary
#' objects (e.g., from \code{\link[emmeans]{emmeans}}).
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
#' used internally by the \pkg{DOTSeq} workflow or manually for testing and
#' diagnostics.
#'
#' @param type A \code{character(1)} string indicating the model type.
#' Default is \code{"fitError"}.
#'
#' @param results A \code{list} containing model results, parameters, and
#' test statistics.
#'
#' @param posthoc An optional object storing post hoc summary objects
#' (e.g., from \code{\link[emmeans]{emmeans}}). Default is 
#' \code{NA_real_}.
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
