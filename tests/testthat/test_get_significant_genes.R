test_that("get_significant_genes returns correct gene IDs below threshold", {
    results_df <- data.frame(
        orf_id = c("ENSG00000139618.19:O001", "ENSG00000139618.19:O002", "ENSG00000157764.15:O003"),
        lfsr = c(0.01, 0.2, 0.03)
    )

    sig_genes <- get_significant_genes(results_df, padj_col = "lfsr", padj_threshold = 0.05)

    expect_type(sig_genes, "character")
    expect_true(all(sig_genes %in% c("ENSG00000139618", "ENSG00000157764")))
    expect_equal(length(sig_genes), 2)
})
