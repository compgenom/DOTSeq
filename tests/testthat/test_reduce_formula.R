test_that("reduce_formula removes single-level terms and retains valid ones", {
    df <- data.frame(
        condition = factor(c("A", "A", "B", "B")),
        strategy = factor(c("ribo", "ribo", "rna", "rna")),
        batch = factor(rep(1, 4)) # single-level factor
    )
    
    formula_input <- ~ batch + condition * strategy
    
    result <- reduce_formula(formula_input, df)
    
    expect_s3_class(result$reduced_formula, "formula")
    expect_true(
        all.equal(
            result$reduced_formula,
            ~ condition * strategy,
            ignore.environment = TRUE
        )
    )
    
    expect_true("emm_specs" %in% names(result))
    expect_true(
        all.equal(
            result$emm_specs,
            ~ condition * strategy,
            ignore.environment = TRUE
        )
    )
})

test_that("reduce_formula throws error when no valid terms remain", {
    df <- data.frame(
        batch = factor(rep(1, 4)) # only one variable, single-level
    )
    
    formula_input <- ~ batch
    
    expect_error(
        reduce_formula(formula_input, df),
        "Please specify an interaction term using '\\*' or ':'"
    )
    
    formula_input <- ~ condition * strategy
    
    expect_error(
        reduce_formula(formula_input, df),
        "Invalid formula: interaction term"
    )
})

test_that("reduce_formula throws error for multiple interaction terms", {
    df <- data.frame(
        condition = factor(c("A", "A", "B", "B")),
        strategy = factor(c("ribo", "ribo", "rna", "rna")),
        group = factor(c("X", "X", "Y", "Y"))
    )
    
    formula_input <- ~ condition * strategy + condition:group
    
    expect_error(
        reduce_formula(formula_input, df),
        "Expected one interaction term"
    )
})