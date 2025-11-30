test_that("simDOT simulates DOTSeqDataSets correctly and handles invalid inputs", {
    dir <- system.file("extdata", package = "DOTSeq")
    
    # Load example data
    cnt <- read.table(file.path(dir, "featureCounts.cell_cycle_subset.txt.gz"),
                      header = TRUE, comment.char = "#")
    names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
    
    flat <- file.path(dir, "gencode.v47.orf_flattened_subset.gtf.gz")
    bed <- file.path(dir, "gencode.v47.orf_flattened_subset.bed.gz")
    
    meta <- read.table(file.path(dir, "metadata.txt.gz"))
    names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
    cond <- meta[meta$treatment == "chx", ]
    cond$treatment <- NULL
    
    # Create DOTSeqDataSets object
    d <- DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = bed
    )
    
    raw_counts <- assay(getDOU(d))
    raw_counts <- raw_counts[, grep("Cycling|Interphase", colnames(raw_counts))]
    ribo <- raw_counts[, grep("ribo", colnames(raw_counts))]
    rna <- raw_counts[, grep("rna", colnames(raw_counts))]
    rowranges <- rowRanges(getDOU(d))
    
    # simulate DOT
    sim <- simDOT(
        ribo = ribo,
        rna = rna,
        annotation = rowranges,
        regulation_type = "uORF_up_mORF_down",
        gcoeff = 1.5,
        num_samples = 1,
        num_batches = 2
    )
    
    expect_s4_class(sim, "DOTSeqDataSets")
    expect_true(nrow(getDOU(sim)) > 0)
    expect_true(ncol(getDOU(sim)) > 0)
    expect_true("strategy" %in% colnames(colData(getDOU(sim))))
    expect_true("gene_id" %in% colnames(rowData(getDOU(sim))))
    
    # Negative tests for invalid inputs
    # Missing ribo argument
    expect_error(simDOT(rna = rna, annotation = rowranges),
                 regexp = "argument \"ribo\" is missing")
    
    # Missing rna argument
    expect_error(simDOT(ribo = ribo, annotation = rowranges),
                 regexp = "argument \"rna\" is missing")
    
    # Invalid annotation type
    expect_error(simDOT(ribo = ribo, rna = rna, annotation = "not_GRanges"),
                 regexp = "object 'orfs' not found")
    
    # Invalid regulation_type
    expect_error(simDOT(ribo = ribo, rna = rna, annotation = rowranges,
                        regulation_type = "invalid_type"),
                 regexp = "Invalid scenario specified")
    
    # Invalid num_samples (negative)
    expect_error(simDOT(ribo = ribo, rna = rna, annotation = rowranges,
                        num_samples = -1),
                 regexp = "invalid 'each' argument")
    
    # Invalid conditions (zero)
    expect_error(simDOT(ribo = ribo, rna = rna, annotation = rowranges,
                        conditions = 0),
                 regexp = "variable lengths differ")
    
    # Invalid batch_scenario
    expect_error(simDOT(ribo = ribo, rna = rna, annotation = rowranges,
                        batch_scenario = "wrong"),
                 regexp = "object 'mod' not found")
    
    # Test with diagnostic plots enabled
    sim <- suppressMessages(simDOT(
        ribo = ribo,
        rna = rna,
        annotation = rowranges,
        regulation_type = "uORF_up_mORF_down",
        gcoeff = 1.5,
        num_samples = 1,
        num_batches = 2,
        diagplot_ribo = TRUE,
        diagplot_rna = TRUE
    ))
    expect_s4_class(sim, "DOTSeqDataSets")
    expect_true(nrow(getDOU(sim)) > 0)
    expect_true(ncol(getDOU(sim)) > 0)
})
