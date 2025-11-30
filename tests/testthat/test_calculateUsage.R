test_that("calculateUsage validates inputs", {
    expect_error(calculateUsage(NULL, "ENSG00000119402"), "Please provide gene_id")
    expect_error(calculateUsage("not_a_DOUData", "ENSG00000119402"), "'data' must be a DOTSeqDataSets object")
    expect_error(calculateUsage(new("DOUData"), NULL), "Please provide gene_id")
    expect_error(calculateUsage(new("DOUData"), 123), "'gene_id' must be a character vector")
})
