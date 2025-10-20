test_that("fitDOU returns expected structure for successfully fitted ORFs", {
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
    
    m <- DOTSeqDataSet(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        bed = bed
    )
    
    m$sumExp <- m$sumExp[rowRanges(m$sumExp)$is_kept == TRUE, ]
    set.seed(42)
    m$sumExp <- m$sumExp[sample(seq_len(nrow(m$sumExp)), size = 50), ]
    
    suppressMessages(
        suppressWarnings({
            results <- fitDOU(
                count_table = assay(m$sumExp),
                rowdata = rowData(m$sumExp),
                anno = colData(m$sumExp),
                formula = ~ condition * strategy,
                emm_specs = ~ condition * strategy,
                dispersion_modeling = "auto",
                lrt = FALSE,
                optimizers = FALSE,
                diagnostic = TRUE,
                parallel = list(n = 2L, autopar = TRUE),
                verbose = FALSE
            )
        })
    )
    
    fitted_orfs <- Filter(function(x) model_type(x) == "glmmTMB", results)
    expect_gt(length(fitted_orfs), 0)
    
    for (res in fitted_orfs) {
        keys <- names(fit_results(res))
        
        expect_true("model_fit" %in% keys)
        expect_true("estimates" %in% keys)
        expect_true("dispersion" %in% keys)
        
        if ("tests" %in% keys) {
            expect_type(fit_results(res)$tests, "list")
        }
        if (has_dharma) {
            expect_true("diagnostics" %in% keys)
            expect_type(fit_results(res)$diagnostics, "list")
        } else {
            expect_false("diagnostics" %in% keys)
        }
    }
})
