test_that("testDOU adds post hoc results to fitted ORFs", {
    testthat::skip_if_not_installed("DHARMa", minimum_version = NULL)
    
    # Load test data
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
    
    dot <- DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = gtf,
        flattened_bed = bed
    )
    
    dou <- getDOU(dot)

    dou <- dou[rowRanges(dou)$is_kept == TRUE, ]
    set.seed(42)
    dou <- dou[sample(seq_len(nrow(dou)), size = 100), ]
    
    # Store DOUData subset
    getDOU(dot) <- dou

    # Test different parameters and error handing
    # lrt optimizers diagnostic verbose TRUE
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "auto",
                lrt = TRUE,
                optimizers = TRUE,
                diagnostic = TRUE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = TRUE
            )
        })
    )
    expect_s4_class(dou, "DOUData")
    
    # lrt optimizers diagnostic verbose TRUE, shared dispersion_modeling
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "shared",
                lrt = TRUE,
                optimizers = TRUE,
                diagnostic = TRUE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = TRUE
            )
        })
    )
    expect_s4_class(dou, "DOUData")
    
    # lrt optimizers diagnostic verbose TRUE, shared dispersion_modeling
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "custom",
                dispformula = ~1,
                lrt = TRUE,
                optimizers = TRUE,
                diagnostic = TRUE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = TRUE
            )
        })
    )
    expect_s4_class(dou, "DOUData")
    
    # shared dispersion_modeling
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "shared",
                lrt = FALSE,
                optimizers = FALSE,
                diagnostic = FALSE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = FALSE
            )
        })
    )
    expect_s4_class(dou, "DOUData")
    
    # custom dispersion_modeling
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "custom",
                dispformula = ~1,
                lrt = FALSE,
                optimizers = FALSE,
                diagnostic = FALSE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = FALSE
            )
        })
    )
    expect_s4_class(dou, "DOUData")
    
    # auto dispersion_modeling
    suppressMessages(
        suppressWarnings({
            rowData(dou)[["DOUResults"]] <- fitDOU(
                data = dou,
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
    expect_s4_class(dou, "DOUData")
    
    # Run post hoc tests
    dou <- testDOU(dou, verbose = FALSE)

    # Check post hoc results
    interaction_df <- getContrasts(dou, type = "interaction")
    strategy_df <- getContrasts(dou, type = "strategy")

    expect_s4_class(interaction_df, "DataFrame")
    expect_s4_class(strategy_df, "DataFrame")

    # Check expected columns in interaction_results
    required_cols <- c("betahat", "sebetahat", "posterior", "qvalue", "lfsr", "lfdr", "contrast")
    expect_true(all(vapply(required_cols, function(x) any(grepl(x, colnames(interaction_df))), logical(1))))

    # Check expected columns in strategy_results
    expect_true(all(vapply(required_cols, function(x) any(grepl(x, colnames(strategy_df))), logical(1))))

    # Check that contrasts are Rle-compressed
    expect_s4_class(interaction_df$contrast, "Rle")
    expect_s4_class(strategy_df$contrast, "Rle")
    expect_s4_class(strategy_df$strategy, "Rle")
    
    res <- extract_results(sumExp = dou)
    
    expect_s3_class(res, "data.frame")
    
    # Check that essential columns exist
    expect_true(all(c("orf_id", "model_type") %in% colnames(res)))
    
    # Check number of rows matches number of ORFs processed
    expect_equal(nrow(res), length(rowRanges(dou)))
    
    # Ensure model_type includes "glmmTMB" for fitted models
    expect_true(any(grepl("glmm", res$model_type)))
    
    # Check that scalar columns (e.g., logLik or AIC) exist if models were fitted
    expect_true(any(grepl("aic", colnames(res))))
    
    # Ensure no row is completely NA except orf_id/model_type
    expect_false(any(rowSums(is.na(res)) == ncol(res)))
    
    # Optional: Check that orf_id values match rownames of dou
    expect_true(all(res$orf_id %in% rownames(dou)))
    
    getContrasts(dou, type = "interaction") <- interaction_df
    getContrasts(dou, type = "strategy") <- strategy_df
    expect_s4_class(dou, "DOUData")
    
    # Test calculateUsage
    gene_id <- strsplit(interaction_df$orf_id[1], ":")[[1]][1]
    gene_id <- strsplit(gene_id, "\\.")[[1]][1]
    usage <- calculateUsage(dou, gene_id) 
    expect_s3_class(usage, "data.frame")
    expect_true(nrow(usage) > 1)
})
