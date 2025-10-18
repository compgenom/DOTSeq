test_that("plotDOT generates a Venn diagram without error", {
    # Example ORF-level results
    results_df <- data.frame(
        orf_id = c("ENSG00000139618.19:O001", "ENSG00000139618.19:O002", "ENSG00000157764.15:O003"),
        lfsr = c(0.01, 0.2, 0.03),
        padj = c(0.02, 0.01, 0.1)
    )

    # Example rowData
    rowdata_df <- data.frame(
        orf_id = c("ENSG00000139618.19:O001", "ENSG00000139618.19:O002", "ENSG00000157764.15:O003"),
        gene_id = c("ENSG00000139618.19", "ENSG00000139618.19", "ENSG00000157764.15"),
        orf_type = c("mORF", "dORF", "mORF")
    )

    # Suppress plotting to avoid opening a device during tests
    expect_silent(
        plotDOT(
            results = results_df,
            rowdata = rowdata_df,
            plot_types = "venn",
            verbose = FALSE,
            force_new_device = FALSE
        )
    )
})
