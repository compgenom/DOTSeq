test_that("calculateTE computes TE correctly for simple input", {
    counts <- matrix(
        c(10, 20, 30, 40), 
        nrow = 2, ncol = 2,
        dimnames = list(c("ORF1", "ORF2"), c("sample1.rna", "sample1.ribo"))
    )
    
    te <- calculateTE(counts, sample_delim = NULL)  # Disable prefix stripping
    
    expect_equal(dim(te), c(2, 1))
    expect_equal(colnames(te), "sample1.te")
    
    expected_te <- (counts[, "sample1.ribo"] + 1e-6) / (counts[, "sample1.rna"] + 1e-6)
    expect_equal(te[, "sample1.te"], expected_te)
})


test_that("calculateTE handles multiple samples", {
    counts <- matrix(
        c(10, 20, 30, 40, 50, 60), 
        nrow = 3, ncol = 2,
        dimnames = list(c("ORF1", "ORF2", "ORF3"), c("s1.rna", "s1.ribo"))
    )
    
    te <- calculateTE(counts, sample_delim = NULL)
    
    expect_equal(dim(te), c(3, 1))
    expect_true(all(grepl("\\.te$", colnames(te))))
})


test_that("calculateTE stops when no matching prefixes", {
    counts <- matrix(
        c(10, 20), 
        nrow = 1, ncol = 2,
        dimnames = list("ORF1", c("sample1.rna", "sample2.ribo"))
    )
    
    expect_error(calculateTE(counts, sample_delim = NULL), "No matching RNA/Ribo sample _prefixes found")
})


test_that("calculateTE warns on ambiguous matches", {
    counts <- matrix(
        c(10, 20, 30, 40, 50, 60, 70, 80, 90),
        nrow = 3, ncol = 3,
        dimnames = list(c("ORF1", "ORF2", "ORF3"), 
                        c("sample1.rna", "sample1.rna.extra", "sample1.ribo"))
    )
    
    expect_warning(calculateTE(counts, sample_delim = NULL), "Ambiguous or missing match")
})


test_that("calculateTE respects pseudocount", {
    counts <- matrix(
        c(0, 0), 
        nrow = 1, ncol = 2,
        dimnames = list("ORF1", c("sample1.rna", "sample1.ribo"))
    )
    
    te <- calculateTE(counts, pseudocount = 1, sample_delim = NULL)
    expect_equal(te[1, 1], 1) # (0+1)/(0+1)
})
