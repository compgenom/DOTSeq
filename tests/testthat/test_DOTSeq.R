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
    dot <- DOTSeq(datasets = datasets, modules = "DTE")
    
    # Check structure
    expect_type(dot, "S4")
    expect_s4_class(dot, "DOTSeqDataSets")
    expect_s4_class(getDTE(dot), "DTEData")
    expect_s4_class(getDOU(dot), "DOUData")
})
