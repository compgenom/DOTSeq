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

    dot <- DOTSeqDataSet(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        bed = bed
    )
    
    dou <- getDOU(dot)

    dou <- dou[rowRanges(dou)$is_kept == TRUE, ]
    set.seed(42)
    dou <- dou[sample(seq_len(nrow(dou)), size = 25), ]

    # Fit model
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                dou = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
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
    dou <- testDOU(dou, verbose = FALSE)

    # Check post hoc results
    interaction_df <- interactionResults(dou)
    strategy_df <- strategyResults(dou)

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
