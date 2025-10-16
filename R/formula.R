#' @title Remove random effects from a formula
#' @description This function takes an R formula object that may contain random effect terms and returns a new formula with all such terms removed.
#' @param formula An R formula object. Random effect terms are identified by the pattern `(1 | group)`.
#' @return A new R formula object containing only the fixed effect terms.
#'
#' @importFrom stats as.formula
#'
remove_random_effects <- function(formula) {
  # Convert the formula to a character string and remove the '~'
  formula_str <- as.character(formula)[2]
  # Split the string by the '+' sign to get individual terms
  terms <- unlist(strsplit(formula_str, "\\s*\\+\\s*"))
  # Identify and remove any terms that contain a random effect pattern
  clean_terms <- terms[!grepl("\\(1\\s*\\|.*?\\)", terms)]
  # Collapse the remaining terms back into a single formula string
  clean_formula_str <- paste(clean_terms, collapse = " + ")
  # Convert the string back to a formula object
  as.formula(paste("~", clean_formula_str))
}


#' Reduce a formula by removing terms with insufficient variation in the data
#'
#' This function takes a formula and a data frame, and returns a simplified formula
#' that includes only terms corresponding to variables with at least two levels.
#' Interaction terms using `*` or `:` are retained only if all involved variables are valid.
#' If the formula is reduced, a message is printed to inform the user.
#'
#' @param formula_input A formula object (e.g., `~ batch + condition * strategy`) or a character string representing a formula.
#' @param data A data frame containing the variables referenced in the formula.
#'
#' @return A reduced formula object that excludes terms with only one level in the data.
#' If no valid terms remain, the function will throw an error.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   condition = factor(c(0, 0, 1, 1)),
#'   strategy = factor(c("ribo", "ribo", "rna", "rna")),
#'   batch = factor(rep(1, 4)) # single-level factor
#' )
#' formula_input <- ~ batch + condition * strategy
#' reduce_formula(formula_input, df)
#' # Returns: ~ condition * strategy, with a message about reduction
#' }
reduce_formula <- function(formula_input, data) {
  # Convert to character if it's a formula
  formula_str <- if (inherits(formula_input, "formula")) {
    deparse(formula_input)
  } else if (is.character(formula_input)) {
    formula_input
  } else {
    stop("Input must be a formula or a character string.")
  }
  
  # Extract RHS of formula
  rhs <- strsplit(formula_str, "~")[[1]][2]
  terms <- strsplit(rhs, "\\+")[[1]]
  terms <- trimws(terms)
  
  matched_terms <- terms[grepl("\\*|:", terms)]
  
  if (length(matched_terms) > 1) {
    stop("Expected one interaction term (e.g., condition * strategy), but got multiple: ", 
         paste(matched_terms, collapse = ", "))
  } else if (length(matched_terms) == 0) {
    stop("Please specify an interaction term using '*' or ':' (e.g., condition * strategy).")
  }
  
  # Identify valid terms
  valid_terms <- c()

  for (term in terms) {
    if (grepl("\\*|:", term)) {
      emm_specs <- as.formula(paste0("~",term))
      components <- trimws(unlist(strsplit(term, "\\*|:")))
      missing_vars <- setdiff(components, colnames(data))
      if (length(missing_vars) > 0) {
        stop(paste0(
          "Invalid formula: interaction term '", term, "' includes missing variable(s): ",
          paste(missing_vars, collapse = ", "), "."
        ))
      }
      valid_terms <- c(valid_terms, term)
    } else if (term %in% colnames(data)) {
      valid_terms <- c(valid_terms, term)
    }
  }
  
  if (length(valid_terms) == 0) {
    stop("No valid terms left in formula after reduction.")
  }
  
  new_formula_str <- paste("~", paste(valid_terms, collapse = " + "))
  reduced_formula <- as.formula(new_formula_str)
  
  if (!isTRUE(all.equal(reduced_formula, formula_input, ignore.environment = TRUE))) {
    message(
      "formula has been reduced due to missing variables or terms with only one level\n",
      "   original formula: ", deparse(formula_input), "\n",
      "   reduced formula:  ", deparse(reduced_formula)
    )
  }
  
  return(list(reduced_formula = reduced_formula, 
              emm_specs = emm_specs))
}
