test_that("mapIDs returns gene symbols", {
    fake_df <- data.frame(
        ensembl_gene_id = "ENSG00000139618",
        external_gene_name = "BRCA2"
    )
    
    with_mocked_bindings(
        bm_use_ensembl = function(biomart, dataset, host = NULL) "fake_mart",
        bm_get = function(attributes, filters, values, mart) fake_df,
        {
            res <- mapIDs(
                "ENSG00000139618",
                dataset = "hsapiens_gene_ensembl",
                mart_source = "ensembl"
            )
            
            expect_equal(res$external_gene_name, "BRCA2")
        }
    )
})
