#' Format elapsed time between two timestamps
#'
#' @description
#' Computes the elapsed time between a start and end timestamp, returning a list
#' with minutes and seconds. If the duration is less than one minute, only seconds
#' are returned. Useful for logging runtimes in a human-readable format.
#'
#' @param end_time POSIXct. The end time of the process.
#' @param start_time POSIXct. The start time of the process.
#' @param units Character. Units for `difftime()` calculation. Default is `"secs"`.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{mins} (optional) — integer number of minutes
#'   \item \code{secs} — integer number of seconds
#' }
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' start <- Sys.time()
#' Sys.sleep(3.5)
#' end <- Sys.time()
#' runtime(end, start)
#' }
runtime <- function(end_time, start_time, units = "secs") {
  total_seconds <- as.numeric(difftime(end_time, start_time, units = units))
  
  mins <- floor(total_seconds / 60)
  secs <- total_seconds - (mins * 60)
  
  if (mins > 0) {
    return(list(mins = mins, secs = secs))
  } else {
    return(list(secs = secs))
  }
}


#' Extract Model Results into a Unified Data Frame
#'
#' @description
#' Extracts scalar results from a named list of model objects (typically from `glmmTMB` fits)
#' into a tidy data frame. Handles nested lists, missing or `NULL` models, and ensures
#' consistent columns across all entries.
#'
#' @param sumExp A SummarizedExperiment object containing fitted model objects 
#' stored in `rowData(sumExp)[['DOUResults']]`.
#' @param verbose Logical. If \code{TRUE}, messages will be printed for skipped or failed models.
#'
#' @return A data frame where each row corresponds to an ORF and its associated model results.
#'   Columns include scalar parameters extracted from the model, \code{orf_id}, and \code{model_type}.
#'
#' @details
#' The function uses a recursive helper (`flatten_scalars`) to extract scalar values
#' from nested lists. It supports models of type \code{"glmmTMB"} and \code{"glmmTMB_joint"},
#' and gracefully handles \code{NULL} or unsupported model types by returning minimal rows.
#'
#' Missing columns across models are filled with \code{NA} to ensure a consistent structure.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom stats setNames
#' 
#' @keywords internal
#' 
extract_results <- function(sumExp, verbose = TRUE) {
  
  models_list <- rowData(sumExp)[["DOUResults"]]
  # Helper to recursively flatten and extract scalar values
  flatten_scalars <- function(x, prefix = NULL) {
    if (is.atomic(x) && length(x) == 1) {
      name <- if (is.null(prefix)) "value" else prefix
      return(setNames(list(x), name))
    } else if (is.list(x)) {
      result <- list()
      # Use `names(x)` and `seq_along` to handle unnamed elements gracefully
      list_names <- names(x)
      if (is.null(list_names)) {
        list_names <- paste0("elem", seq_along(x))
      }
      for (i in seq_along(x)) {
        n <- list_names[i]
        sub_prefix <- if (is.null(prefix)) n else paste0(prefix, ".", n)
        result <- c(result, flatten_scalars(x[[i]], sub_prefix))
      }
      return(result)
    } else {
      return(list()) # Return an empty list for non-atomic, non-list objects
    }
  }
  
  # Process a single model
  process_model <- function(model_obj, name) {
    # Check for NULL objects
    if (is.null(model_obj)) {
      if (verbose) message("Skipping orf_id ", name, " due to NULL model object.")
      return(data.frame(
        orf_id = name,
        model_type = "NULL", # Or "Failed", etc.
        stringsAsFactors = FALSE
      ))
    }
    
    model_type <- model_obj@type
    
    if ((model_type == "glmmTMB_joint") | (model_type == "glmmTMB")) {
      params <- model_obj@results
      param_values <- flatten_scalars(params)
      
      param_values$orf_id <- name
      param_values$model_type <- model_type
      
      df <- as.data.frame(param_values, stringsAsFactors = FALSE)
      
    } else {
      # Return minimal row with NA for failed or single-ORF models
      df <- data.frame(
        orf_id = name,
        model_type = model_type,
        stringsAsFactors = FALSE
      )
    }
    
    return(df)
  }
  
  # Apply to all models
  results_df_list <- lapply(names(models_list), function(name) {
    process_model(models_list[[name]], name)
  })
  
  # Get all unique column names
  all_cols <- unique(unlist(lapply(results_df_list, names)))
  
  # Fill missing columns with NA
  results_df_list_filled <- lapply(results_df_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]
  })
  
  # Combine all rows
  results_df <- do.call(rbind, results_df_list_filled)
  rownames(results_df) <- NULL
  
  return(results_df)
}


