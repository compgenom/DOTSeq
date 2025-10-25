test_that("DOTSeqDataSet returns valid SummarizedExperiment object", {
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

    dot <- DOTSeqDataSet(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        bed = bed
    )

    expect_type(dot, "S4")
    expect_s4_class(dot, "DOTSeqObjects")
    
    dou <- getDOU(dot)

    # Check metadata
    expect_true(nrow(dou) > 0)
    expect_true(ncol(dou) > 0)
    expect_true("strategy" %in% colnames(colData(dou)))
    expect_true("gene_id" %in% colnames(rowData(dou)))
})
