test_that("plotDOT generates all plot types without error", {
    set.seed(42)
    
    # Helper to generate random strings
    random_string <- function(n, length = 12) {
        replicate(n, paste0(sample(c(LETTERS, 0:9), length, replace = TRUE), collapse = ""))
    }
    
    # Generate 50 random ORF IDs
    gene_ids <- random_string(50, 12)
    orf_suffixes <- random_string(50, 4)
    orf_ids <- paste0(gene_ids, ":", orf_suffixes)
    
    # Create results_df
    results_df <- data.frame(
        orf_id = orf_ids,
        lfsr = runif(50, 0, 0.1),
        padj = runif(50, 0, 0.1),
        posterior = rnorm(50),
        log2FoldChange = rnorm(50)
    )
    
    # Plot Venn diagram and composite plots
    expect_message(
        plotDOT(
            results = results_df,
            plot_type = "composite",
            verbose = FALSE,
            force_new_device = FALSE
        ),
        regexp = "Spearman"
    )
    
    testthat::skip_if_not_installed("eulerr", minimum_version = NULL)
    expect_message(
        plotDOT(
            results = results_df,
            plot_type = "venn",
            verbose = FALSE,
            force_new_device = TRUE
        ),
        regexp = "resetting graphics device"
    )
})