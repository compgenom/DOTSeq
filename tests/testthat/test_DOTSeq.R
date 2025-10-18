test_that("DOTSeq wrapper returns expected structure", {
    # Load test data from extdata
    dir <- system.file("extdata", package = "DOTSeq")

    cnt <- read.table(
        file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
        header = TRUE,
        comment.char = "#"
    )
    names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))

    flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
    bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")

    meta <- read.table(file.path(dir, "metadata.txt.gz"))
    names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
    cond <- meta[meta$treatment == "chx", ]
    cond$treatment <- NULL

    # Run DOTSeq wrapper
    m <- DOTSeq(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        bed = bed,
        modules = "DTE"
    )

    # Check structure
    expect_type(m, "list")
    expect_true("dds" %in% names(m))
    expect_true("sumExp" %in% names(m))
    expect_s4_class(m$dds, "DESeqDataSet")
    expect_s4_class(m$sumExp, "SummarizedExperiment")
})
