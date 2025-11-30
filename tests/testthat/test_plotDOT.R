test_that("plotDOT error handling", {
    testthat::skip_if_not_installed("eulerr", minimum_version = NULL)
    testthat::skip_if_not_installed("ggplot2", minimum_version = NULL)
    testthat::skip_if_not_installed("ggsignif", minimum_version = NULL)
    
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
    # set.seed(42)
    # dou <- dou[sample(seq_len(nrow(dou)), size = 100), ]
    
    # Get the first 100 rows
    dou <- dou[seq_len(100), ]
    
    # Store DOUData subset
    getDOU(dot) <- dou
    
    # Test plots
    suppressMessages(
        suppressWarnings({
            d <- DOTSeq(dot)
        })
    )
    
    expect_s4_class(d, "DOTSeqDataSets")
    
    results <- getContrasts(d, type = "interaction")
    ou <- results$DOU[results$DOU$contrast == "Mitotic_Cycling - Interphase", ]
    te <- results$DTE[results$DTE$contrast == "Mitotic_Cycling - Interphase", ]
    results <- merge(ou, te, by = c("orf_id", "contrast"), all = TRUE)
    results$lfsr <- results$lfsr * 0.0001
    
    # Test plots using ID mapping
    # gene_id <- strsplit(results$orf_id[1], ":")[[1]][1]
    # gene_id <- strsplit(gene_id, "\\.")[[1]][1]
    gene_ids <- vapply(strsplit(results$orf_id, ":"), function(x) x[1], character(1))
    mapping <- data.frame(gene_ids, results$orf_id, gene_ids, gene_ids)
    names(mapping) <- c("gene_id", "orf_id", "ensembl_gene_id", "external_gene_name")
    
    # Successful tests
    expect_message(
        plotDOT(
            plot_type = "composite", 
            results = results, 
            plot_params = list(color_by = "significance", legend_position = "bottomright"),
            force_new_device = FALSE
        ), regexp = "correlation"
    )
    
    expect_message(
        plotDOT(
            plot_type = "composite", 
            results = results, 
            data = getDOU(d),
            plot_params = list(color_by = "orf_type", legend_position = "bottomright"),
            force_new_device = FALSE
        ), regexp = "correlation",
        fixed = FALSE
    )
    
    expect_message(
        plotDOT(
            plot_type = "composite",
            results = results,
            plot_params = list(color_by = "significance", legend_position = "bottomright"),
            force_new_device = FALSE
        ), regexp = "correlation"
    )
    
    expect_message(
        plotDOT(
            plot_type = "heatmap", 
            results = results, 
            data = getDOU(d), 
            id_mapping = mapping, 
            plot_params = list(rank_by = "score", top_hits = 20),
            force_new_device = TRUE
        ), regexp = "heatmap",
        fixed = FALSE
    )
    
    # Plotting failed
    orderby <- c("Mitotic_Cycling", "Mitotic_Arrest", "Interphase")
    expect_message(
        suppressWarnings({
            plotDOT(
                plot_type = "usage",
                data = getDOU(d), 
                gene_id = "gene_id", 
                id_mapping = mapping, 
                plot_params = list(order_by = orderby),
                force_new_device = TRUE
            )
        }), regexp = "failed",
        fixed = FALSE
    )
    
    # Test plots without ID mapping
    # Error handling
    set.seed(42)
    res <- results[sample(seq_len(nrow(results)), size = 5), ]
    expect_warning(
        plotDOT(
            plot_type = "composite", 
            results = res, 
            plot_params = list(color_by = "significance", legend_position = "bottomright"),
            force_new_device = FALSE
        ), regexp = "no non-missing arguments"
    )
    
    expect_warning(
        plotDOT(
            plot_type = "composite", 
            results = res, 
            data = getDOU(d),
            plot_params = list(color_by = "orf_type", legend_position = "bottomright"),
            force_new_device = FALSE
        ), regexp = "no non-missing arguments",
        fixed = FALSE
    )
    
    expect_message(
        plotDOT(
            plot_type = "volcano", 
            results = results,
            data = getDOU(d),
            id_mapping = FALSE,
            plot_params = list(color_by = "significance", legend_position = "topright"),
            force_new_device = FALSE
        ), regexp = "volcano",
        fixed = FALSE
    )
    
    expect_message(
        plotDOT(
            plot_type = "volcano", 
            results = results,
            data = getDOU(d),
            id_mapping = FALSE,
            plot_params = list(color_by = "orf_type", legend_position = "topright"),
            force_new_device = FALSE
        ), regexp = "volcano",
        fixed = FALSE
    )
    
    expect_message(
        plotDOT(
            plot_type = "heatmap", 
            results = results, 
            data = getDOU(d), 
            id_mapping = FALSE, 
            plot_params = list(rank_by = "score", top_hits = 20),
            force_new_device = FALSE
        ), regexp = "heatmap",
        fixed = FALSE
    )
    
    id <- "ENSG00000016864.19" # "ENSG00000112245.14"
    if (id %in% gene_ids) {
        expect_message(
            plotDOT(
                plot_type = "usage",
                data = getDOU(d),
                gene_id = id,
                id_mapping = FALSE,
                plot_params = list(order_by = orderby),
                force_new_device = FALSE
            ), regexp = "usage",
            fixed = FALSE
        )
    }
})


test_that("plotDOT generates all plot types without error", {
    set.seed(42)
    
    # Helper to generate random strings
    random_string <- function(n, length = 12) {
        replicate(n, paste0(sample(c(LETTERS, 0:9), length, replace = TRUE), collapse = ""))
    }
    
    # Generate 50 random ORF IDs
    gene_ids <- random_string(50, 12)
    orf_suffixes <- random_string(50, 4)
    orf_ids <- paste0(gene_ids, ":", orf_suffixes)
    
    # Create results_df
    results_df <- data.frame(
        orf_id = orf_ids,
        lfsr = runif(50, 0, 0.1),
        padj = runif(50, 0, 0.1),
        posterior = rnorm(50),
        log2FoldChange = rnorm(50)
    )
    
    # Plot Venn diagram and composite plots
    expect_message(
        plotDOT(
            results = results_df,
            plot_type = "composite",
            verbose = FALSE,
            force_new_device = FALSE
        ),
        regexp = "Spearman"
    )
    
    testthat::skip_if_not_installed("eulerr", minimum_version = NULL)
    expect_message(
        plotDOT(
            results = results_df,
            plot_type = "venn",
            verbose = FALSE,
            force_new_device = TRUE
        ),
        regexp = "resetting graphics device"
    )
})
