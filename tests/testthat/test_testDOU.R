test_that("testDOU adds post hoc results to fitted ORFs", {
    # Load test data
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

    m <- DOTSeqDataSet(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        bed = bed
    )

    m$sumExp <- m$sumExp[rowRanges(m$sumExp)$is_kept == TRUE, ]
    set.seed(42)
    m$sumExp <- m$sumExp[sample(seq_len(nrow(m$sumExp)), size = 25), ]

    # Fit model
    suppressMessages(
        suppressWarnings({
            rowData(m$sumExp)[["DOUResults"]] <- fitDOU(
                count_table = assay(m$sumExp),
                rowdata = rowData(m$sumExp),
                coldata = colData(m$sumExp),
                formula = ~ condition * strategy,
                emm_specs = ~ condition * strategy,
                dispersion_modeling = "auto",
                lrt = FALSE,
                optimizers = FALSE,
                diagnostic = FALSE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = FALSE
            )
        })
    )
    # Run post hoc tests
    m$sumExp <- testDOU(m$sumExp, verbose = FALSE)

    # Check metadata for contrast results
    meta_data <- metadata(m$sumExp)

    expect_true("interaction_results" %in% names(meta_data))
    expect_true("strategy_results" %in% names(meta_data))

    interaction_df <- meta_data$interaction_results
    strategy_df <- meta_data$strategy_results

    expect_s4_class(interaction_df, "DataFrame")
    expect_s4_class(strategy_df, "DataFrame")

    # Check expected columns in interaction_results
    expect_true(all(c("betahat", "sebetahat", "PosteriorMean", "qvalue", "lfsr", "lfdr", "contrast") %in% colnames(interaction_df)))

    # Check expected columns in strategy_results
    expect_true(all(c("betahat", "sebetahat", "PosteriorMean", "qvalue", "lfsr", "contrast", "strategy") %in% colnames(strategy_df)))

    # Check that contrasts are Rle-compressed
    expect_s4_class(interaction_df$contrast, "Rle")
    expect_s4_class(strategy_df$contrast, "Rle")
    expect_s4_class(strategy_df$strategy, "Rle")
})
