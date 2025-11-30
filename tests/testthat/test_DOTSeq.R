test_that("DOTSeq wrapper returns expected structure", {
    # Load test data from extdata
    dir <- system.file("extdata", package = "DOTSeq")

    cnt <- read.table(
        file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
        header = TRUE,
        comment.char = "#"
    )
    names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))

    gtf <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
    bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")

    meta <- read.table(file.path(dir, "metadata.txt.gz"))
    names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
    cond <- meta[meta$treatment == "chx", ]
    cond$treatment <- NULL

    # Run DOTSeq wrapper
    datasets <- DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = gtf,
        flattened_bed = bed
    )
    d <- DOTSeq(datasets = datasets, modules = "DTE")
    expect_message(show(d), regexp = "DOTSeqDataSets")
    
    # Check structure
    expect_type(d, "S4")
    expect_s4_class(d, "DOTSeqDataSets")
    expect_s4_class(getDTE(d), "DTEData")
    expect_s4_class(getDOU(d), "DOUData")
    expect_type(specs(getDTE(d)), "language")
    expect_type(fmla(getDTE(d)), "language")
    
    res <- getContrasts(d, type = "interaction")
    expect_type(res, "list")
    
    # Do replacement
    dou <- getDOU(d)
    getDOU(d) <- dou
    
    dte <- getDTE(d)
    getDTE(d) <- dte
    
    # Check structure again
    expect_type(d, "S4")
    expect_s4_class(d, "DOTSeqDataSets")
    expect_s4_class(getDTE(d), "DTEData")
    expect_s4_class(getDOU(d), "DOUData")
    
    te <- getContrasts(dte, type = "interaction")
    expect_type(te, "S4")
    
    getContrasts(dte, type = "interaction") <- te
    expect_s4_class(dte, "DTEData")
    
    # Negative tests
    expect_error(DOTSeq(datasets = NULL), regexp = "must be either a DOTSeqDataSets object")
    
    expect_error(DOTSeq(datasets = NULL, modules = NULL), regexp = "Please specify at least one module")
})


# Minimal unit test for assign_strategy_levels
test_that("assign_strategy_levels basic functionality", {
    df <- data.frame(strategy = c("RNA", "Ribo"))
    levels <- assign_strategy_levels(df)
    expect_equal(levels$rna_level, "RNA")
    expect_equal(levels$ribo_level, "Ribo")

    df_multi_rna <- data.frame(strategy = c("RNA", "rna", "Ribo"))
    expect_warning(assign_strategy_levels(df_multi_rna), "Multiple RNA-seq levels found")

    df_no_rna <- data.frame(strategy = c("Ribo", "Ribo2"))
    expect_error(assign_strategy_levels(df_no_rna), "No RNA-seq level found")
})


# Integration test for DOTSeqDataSetsFromFeatureCounts
test_that("DOTSeqDataSetsFromFeatureCounts handles invalid condition_table", {
    dir <- system.file("extdata", package = "DOTSeq")
    cnt <- read.table(file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
                      header = TRUE, comment.char = "#")
    names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))

    flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
    bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")

    meta <- read.table(file.path(dir, "metadata.txt.gz"))
    names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
    cond <- meta[meta$treatment == "chx", ]
    cond$treatment <- NULL

    # Missing columns
    bad_cond1 <- cond
    bad_cond1$strategy <- NULL
    expect_error(DOTSeqDataSetsFromFeatureCounts(cnt, bad_cond1, flat, bed), "Invalid formula")

    bad_cond2 <- cond
    bad_cond2$condition <- NULL
    expect_error(DOTSeqDataSetsFromFeatureCounts(cnt, bad_cond2, flat, bed), "Invalid formula")

    # Empty table
    bad_cond3 <- cond[0, ]
    expect_error(DOTSeqDataSetsFromFeatureCounts(cnt, bad_cond3, flat, bed), "Only 0 runs are matched")

    # Invalid strategy values
    bad_cond4 <- cond
    bad_cond4$strategy <- rep("invalid", nrow(bad_cond4))
    expect_error(DOTSeqDataSetsFromFeatureCounts(cnt, bad_cond4, flat, bed), "No RNA-seq level found")

    # Mixed numeric representation (should succeed)
    bad_cond5 <- cond
    bad_cond5$strategy <- c("0", rep("Ribo", nrow(bad_cond5) - 1))
    expect_error(
        DOTSeqDataSetsFromFeatureCounts(cnt, bad_cond5, flat, bed),
        "model matrix is not full rank"
    )
    
})
