test_that("DOTSeqDataSetsFromFeatureCounts returns valid DOTSeqDataSets object and handles invalid inputs", {
    dir <- system.file("extdata", package = "DOTSeq")
    
    # Load valid input files
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
    
    dot <- DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = bed
    )
    
    expect_s4_class(dot, "DOTSeqDataSets")
    dou <- getDOU(dot)
    expect_true(nrow(dou) > 0)
    expect_true(ncol(dou) > 0)
    expect_true("strategy" %in% colnames(colData(dou)))
    expect_true("gene_id" %in% colnames(rowData(dou)))
    
    # Negative tests for invalid inputs
    # Missing count_table
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = NULL,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = bed
    ), regexp = "Missing required arguments: count_table")
    
    # Missing condition_table
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = NULL,
        flattened_gtf = flat,
        flattened_bed = bed
    ), regexp = "Missing required arguments: condition_table")
    
    # Missing annotation files
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = NULL,
        flattened_bed = bed
    ), regexp = "Missing required arguments: flattened_gtf")
    
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = NULL
    ), regexp = "Missing required arguments: flattened_bed")
    
    # Condition table missing columns → formula error
    bad_cond <- data.frame(run = cond$run) # missing replicate, condition, strategy
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = bad_cond,
        flattened_gtf = flat,
        flattened_bed = bed
    ), regexp = "Invalid formula.*missing variable")
    
    # Non-integer counts → rowname mismatch occurs first
    bad_cnt <- cnt
    bad_cnt[] <- runif(length(bad_cnt)) # convert to numeric
    expect_error(DOTSeqDataSetsFromFeatureCounts(
        count_table = bad_cnt,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = bed
    ), regexp = "identical\\(rownames\\(dcounts\\), names\\(gff_granges\\)\\) is not TRUE")
    
    # Invalid file paths → rtracklayer/read.table error
    expect_error(suppressWarnings(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = "nonexistent.gtf",
        flattened_bed = bed
    )), regexp = "cannot open the connection")
    
    expect_error(suppressWarnings(DOTSeqDataSetsFromFeatureCounts(
        count_table = cnt,
        condition_table = cond,
        flattened_gtf = flat,
        flattened_bed = "nonexistent.bed"
    )), regexp = "cannot open the connection")
})


test_that("DOTSeqDataSetsFromSummarizeOverlaps returns valid DOTSeqDataSets object", {
    testthat::skip_if_not_installed("withr")
    testthat::skip_if_not_installed("pasillaBamSubset")
    testthat::skip_if_not_installed("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
    
    library(withr)
    library(pasillaBamSubset)
    library(GenomeInfoDb)
    library(AnnotationDbi)
    library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    
    # Save a subset of TxDb as an sqlite file
    txdb_chr4 <- keepSeqlevels(
        TxDb.Dmelanogaster.UCSC.dm3.ensGene, 
        "chr4", 
        pruning.mode = "coarse"
    )
    
    temp_dir <- tempdir()
    
    txdb_path <- file.path(temp_dir, "dm3_chr4.sqlite")
    saveDb(txdb_chr4, file = txdb_path)
    withr::defer(unlink(txdb_path))  # cleanup
    
    # Create a GRanges object with a link to the TxDb sqlite file, 
    # which is required for getExonicReads()
    gr <- genes(txdb_chr4)
    metadata(gr)$txdb <- txdb_path
    
    # Prepare GRanges for testing
    bam_list <- c(untreated1_chr4(), untreated3_chr4())
    
    # Set the chromosome names according to the BAM files
    seqlevelsStyle(gr) = "UCSC"
    
    # Keep only reads mapped to the exons of coding genes
    getExonicReads(gr = gr, bam_files = bam_list, bam_output_dir = temp_dir, coding_genes_only = TRUE)
    
    bam_list <- list.files(temp_dir, pattern = "exonic", full.names = TRUE)
    
    # Run countReads
    cnt <- countReads(gr = gr, bam_files = bam_list, verbose = FALSE)
    withr::defer(unlink(bam_list))

    # Expectations
    expect_s3_class(cnt, "data.frame")  # Output should be a data frame
    expect_true(nrow(cnt) == length(gr))  # Rows should match number of features
    expect_true(ncol(cnt) == length(bam_list))  # Columns should match number of BAM files
    expect_true(all(colnames(cnt) %in% basename(bam_list)))  # Column names should correspond to BAM files
    expect_true(all(cnt >= 0))  # Counts should be non-negative integers
    
    set.seed(42)
    
    # Create count_table
    # Create two replicates for each condition with random scaling
    rna_treated_reps <- do.call(cbind, replicate(2, cnt[["untreated3_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), min = 0.5, max = 2), simplify = FALSE))
    rna_control_reps <- do.call(cbind, replicate(2, cnt[["untreated3_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), min = 0.1, max = 0.5), simplify = FALSE))
    ribo_treated_reps <- do.call(cbind, replicate(2, cnt[["untreated1_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), min = 0.5, max = 2), simplify = FALSE))
    ribo_control_reps <- do.call(cbind, replicate(2, cnt[["untreated1_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), min = 0.1, max = 0.5), simplify = FALSE))
    
    # Combine and name columns
    colnames(rna_treated_reps) <- paste0("rna_treated", 1:2)
    colnames(rna_control_reps) <- paste0("rna_control", 1:2)
    colnames(ribo_treated_reps) <- paste0("ribo_treated", 1:2)
    colnames(ribo_control_reps) <- paste0("ribo_control", 1:2)
    
    cnt_expanded <- cbind(
        rna_treated_reps, 
        rna_control_reps, 
        ribo_treated_reps, 
        ribo_control_reps
    )
    
    # Convert numbers to integer
    cnt_expanded <- round(cnt_expanded)
    storage.mode(cnt_expanded) <- "integer"
    
    rownames(cnt_expanded) <- rownames(cnt)
    cnt_expanded <- as.data.frame(cnt_expanded)
    
    
    # Create condition_table
    # Sample names from cnt_expanded
    sample_names <- colnames(cnt_expanded)
    
    # Define condition and strategy for each sample
    condition <- c(
        rep("treated", 2), 
        rep("control", 2), 
        rep("treated", 2), 
        rep("control", 2)
    )
    strategy <- c(rep("RNA", 4), rep("Ribo", 4))
    
    cond <- data.frame(
        run = sample_names,
        replicate = c(1,2),
        condition = factor(condition, levels = c("control", "treated")),
        strategy = factor(strategy, levels = c("RNA", "Ribo"))
    )
    
    # Create a DOTSeqDataSets object
    d <- DOTSeqDataSetsFromSummarizeOverlaps(
        count_table = cnt_expanded, 
        condition_table = cond, 
        annotation = gr
    )
    
    expect_type(d, "S4")
    expect_s4_class(d, "DOTSeqDataSets")
    
    dou <- getDOU(d)
    
    # Check metadata
    expect_true(nrow(dou) > 0)
    expect_true(ncol(dou) > 0)
    expect_true("strategy" %in% colnames(colData(dou)))
    expect_true("gene_id" %in% colnames(rowData(dou)))
})


test_that("DOTSeqDataSetsFromSummarizeOverlaps returns valid DOTSeqDataSets object and handles invalid inputs", {
    testthat::skip_if_not_installed("withr")
    testthat::skip_if_not_installed("pasillaBamSubset")
    testthat::skip_if_not_installed("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
    
    library(pasillaBamSubset)
    library(GenomeInfoDb)
    library(AnnotationDbi)
    library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    
    set.seed(42)
    
    # Save a subset of TxDb as an sqlite file
    txdb_chr4 <- keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "chr4", pruning.mode = "coarse")
    txdb_path <- tempfile(fileext = ".sqlite")
    saveDb(txdb_chr4, file = txdb_path)
    withr::defer(unlink(txdb_path))  # cleanup
    
    # Create GRanges linked to TxDb
    gr <- genes(txdb_chr4)
    gr <- gr[sample(seq_along(gr), size = 50)]
    metadata(gr)$txdb <- txdb_path
    seqlevelsStyle(gr) <- "UCSC"
    
    # Prepare BAM files
    bam_list <- c(untreated1_chr4(), untreated3_chr4())
    
    # Filter exonic reads
    temp_dir <- tempdir()
    getExonicReads(gr = gr, bam_files = bam_list, bam_output_dir = temp_dir, coding_genes_only = TRUE)
    
    # Get filtered BAM files
    bam_list <- list.files(path = temp_dir, pattern = "*exonic.*", full.names = TRUE)
    withr::defer(unlink(bam_list))  # cleanup
    
    # Count reads
    cnt <- countReads(gr = gr, bam_files = bam_list, verbose = FALSE)
    
    # Expectations for happy path
    expect_s3_class(cnt, "data.frame")
    expect_equal(nrow(cnt), length(gr))
    expect_equal(ncol(cnt), length(bam_list))
    expect_true(all(colnames(cnt) %in% basename(bam_list)))
    expect_true(all(cnt >= 0))
    
    # Expand counts for multiple replicates
    rna_treated_reps <- do.call(cbind, replicate(2, cnt[["untreated3_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), 0.5, 2), simplify = FALSE))
    rna_control_reps <- do.call(cbind, replicate(2, cnt[["untreated3_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), 0.1, 0.5), simplify = FALSE))
    ribo_treated_reps <- do.call(cbind, replicate(2, cnt[["untreated1_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), 0.5, 2), simplify = FALSE))
    ribo_control_reps <- do.call(cbind, replicate(2, cnt[["untreated1_chr4.exonic.sorted.bam"]] * runif(nrow(cnt), 0.1, 0.5), simplify = FALSE))
    
    colnames(rna_treated_reps) <- paste0("rna_treated", 1:2)
    colnames(rna_control_reps) <- paste0("rna_control", 1:2)
    colnames(ribo_treated_reps) <- paste0("ribo_treated", 1:2)
    colnames(ribo_control_reps) <- paste0("ribo_control", 1:2)
    
    cnt_expanded <- cbind(rna_treated_reps, rna_control_reps, ribo_treated_reps, ribo_control_reps)
    cnt_expanded <- round(cnt_expanded)
    storage.mode(cnt_expanded) <- "integer"
    rownames(cnt_expanded) <- rownames(cnt)
    cnt_expanded <- as.data.frame(cnt_expanded)
    
    # Condition table
    sample_names <- colnames(cnt_expanded)
    condition <- c(rep("treated", 2), rep("control", 2), rep("treated", 2), rep("control", 2))
    strategy <- c(rep("RNA", 4), rep("Ribo", 4))
    
    cond <- data.frame(
        run = sample_names,
        replicate = rep(1:2, 4),
        condition = factor(condition, levels = c("control", "treated")),
        strategy = factor(strategy, levels = c("RNA", "Ribo"))
    )
    
    # Create DOTSeqDataSets object
    d <- DOTSeqDataSetsFromSummarizeOverlaps(cnt_expanded, cond, gr)
    expect_s4_class(d, "DOTSeqDataSets")
    
    dou <- getDOU(d)
    expect_true(nrow(dou) > 0)
    expect_true(ncol(dou) > 0)
    expect_true("strategy" %in% colnames(colData(dou)))
    expect_true("gene_id" %in% colnames(rowData(dou)))

    # Negative tests for invalid inputs
    minimal_cnt <- data.frame(sample1 = 1:10, sample2 = 11:20)
    valid_cond <- data.frame(
        run = c("sample1", "sample2"),
        replicate = c(1, 2),
        condition = factor(c("control", "treated")),
        strategy = factor(c("RNA", "Ribo"))
    )
    valid_gr <- GRanges(seqnames = "chr4", ranges = IRanges(1, 100))
    
    # Missing arguments
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(NULL, valid_cond, valid_gr),
                 "Missing required arguments: count_table")
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(minimal_cnt, NULL, valid_gr),
                 "Missing required arguments: condition_table")
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(minimal_cnt, valid_cond, NULL),
                 "Missing required arguments: annotation")
    
    # Condition table missing columns triggers formula error
    bad_cond <- data.frame(run = c("sample1", "sample2"))
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(minimal_cnt, bad_cond, valid_gr),
                 "Invalid formula: interaction term 'condition \\* strategy'")
    
    # Non-integer counts triggers run ID mismatch first
    bad_cnt <- data.frame(sample1 = runif(10), sample2 = runif(10))
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(bad_cnt, valid_cond, valid_gr),
                 "Run ID mismatch")
    
    # Empty GRanges triggers run ID mismatch first
    expect_error(DOTSeqDataSetsFromSummarizeOverlaps(minimal_cnt, valid_cond, GRanges()),
                 "Run ID mismatch")
})
