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


test_that("mapIDs falls back to description when symbol_col is NA", {
    fake_df_na <- data.frame(
        ensembl_gene_id = "ENSG00000139618",
        external_gene_name = NA,
        description = "BRCA2 DNA repair associated"
    )
    
    
    with_mocked_bindings(
        bm_use_ensembl = function(biomart, dataset, host = NULL) "fake_mart",
        bm_get = function(attributes, filters, values, mart) {
            if ("description" %in% attributes) {
                data.frame(
                    ensembl_gene_id = "ENSG00000139618",
                    description = "BRCA2 DNA repair associated"
                )
            } else {
                data.frame(
                    ensembl_gene_id = "ENSG00000139618",
                    external_gene_name = NA
                )
            }
        },
        {
            res <- mapIDs("ENSG00000139618", dataset = "hsapiens_gene_ensembl")
            expect_equal(res$external_gene_name, "BRCA2 DNA repair associated")
        }
    )
})


test_that("mapIDs validates inputs", {
    expect_error(mapIDs(NULL, dataset = "hsapiens_gene_ensembl"), "must be a non-empty character vector")
    expect_error(mapIDs("ENSG00000139618", dataset = NULL), "must be a string")
    expect_error(mapIDs("ENSG00000139618", dataset = "hsapiens_gene_ensembl", include_go = "yes"), "must be TRUE or FALSE")
    expect_error(mapIDs("ENSG00000139618", dataset = "hsapiens_gene_ensembl", mart_source = "invalid"), "should be one of")
})


test_that("mapIDs includes GO terms when include_go = TRUE", {
    fake_df_go <- data.frame(
        ensembl_gene_id = "ENSG00000139618",
        external_gene_name = "BRCA2",
        go_id = "GO:0006281",
        name_1006 = "DNA repair",
        namespace_1003 = "biological_process"
    )
    
    with_mocked_bindings(
        bm_use_ensembl = function(biomart, dataset, host = NULL) "fake_mart",
        bm_get = function(attributes, filters, values, mart) fake_df_go,
        {
            res <- mapIDs(
                "ENSG00000139618",
                dataset = "hsapiens_gene_ensembl",
                include_go = TRUE
            )
            
            expect_true(all(c("go_id", "name_1006", "namespace_1003") %in% colnames(res)))
        }
    )
})


test_that("mapIDs passes custom host to bm_use_ensembl", {
    host_used <- NULL
    with_mocked_bindings(
        bm_use_ensembl = function(biomart, dataset, host = NULL) {
            host_used <<- host
            "fake_mart"
        },
        bm_get = function(attributes, filters, values, mart) data.frame(ensembl_gene_id = "ENSG00000139618", external_gene_name = "BRCA2"),
        {
            mapIDs("ENSG00000139618", dataset = "hsapiens_gene_ensembl", host = "https://archive.ensembl.org")
            expect_equal(host_used, "https://archive.ensembl.org")
        }
    )
})
