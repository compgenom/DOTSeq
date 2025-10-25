test_that("fitDOU returns expected structure for successfully fitted ORFs", {
    # Skip test if DHARMa not installed
    testthat::skip_if_not_installed("DHARMa")
    
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
    
    # Subset for speed
    dou <- getDOU(dot)
    
    dou <- dou[rowRanges(dou)$is_kept == TRUE, ]
    set.seed(42)
    dou <- dou[sample(seq_len(nrow(dou)), size = 25), ]
    
    # Fit models (with diagnostics enabled because DHARMa is available)
    suppressMessages(
        suppressWarnings({
            results <- fitDOU(
                dou = dou,
                formula = ~ condition * strategy,
                specs = ~ condition * strategy,
                dispersion_modeling = "auto",
                lrt = FALSE,
                optimizers = FALSE,
                diagnostic = TRUE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = FALSE
            )
        })
    )
    
    # Filter successfully fitted models
    fitted_orfs <- Filter(function(x) modelType(x) == "glmmTMB", results)
    expect_gt(length(fitted_orfs), 0)
    
    for (res in fitted_orfs) {
        keys <- names(fitResults(res))
        
        expect_true("model_fit" %in% keys)
        expect_true("estimates" %in% keys)
        expect_true("dispersion" %in% keys)
        
        if ("tests" %in% keys)
            expect_type(fitResults(res)$tests, "list")
        
        expect_true("diagnostics" %in% keys)
        expect_type(fitResults(res)$diagnostics, "list")
    }
})
