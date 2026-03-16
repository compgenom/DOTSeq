test_that("countReadsSingleCell: tiny BAM with CB/UB tags yields expected sparse matrix", {
    skip_on_cran()
    skip_if_not_installed("Rsamtools")
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("BiocParallel")
    
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(Rsamtools)
        library(Matrix)
        library(BiocParallel)
    })
    
    # Use a temporary directory for all files and clean up automatically
    td <- withr::local_tempdir()
    
    # Two ORF-like features on chr1
    gr <- GRanges("chr1", IRanges(c(90, 140), width = 40))
    names(gr) <- c("orf1", "orf2")
    
    # Tiny SAM with CB/UB tags: two cells (cellA, cellB), 4 reads total
    # r1 and r2: same (cellA, umi1), overlapping orf1 -> dedup to 1 when dedup=TRUE
    # r3: cellA / umi2 overlaps orf2
    # r4: cellB / umi1 overlaps orf2
    sam <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        "r1\t0\tchr1\t100\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tFFFFFFFFFFFFFFFFFFFF\tCB:Z:cellA\tUB:Z:umi1",
        "r2\t0\tchr1\t110\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tFFFFFFFFFFFFFFFFFFFF\tCB:Z:cellA\tUB:Z:umi1",
        "r3\t0\tchr1\t150\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tFFFFFFFFFFFFFFFFFFFF\tCB:Z:cellA\tUB:Z:umi2",
        "r4\t0\tchr1\t150\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tFFFFFFFFFFFFFFFFFFFF\tCB:Z:cellB\tUB:Z:umi1"
    )
    sam_path <- file.path(td, "toy.sam")
    writeLines(sam, sam_path)
    
    dest_base <- file.path(td, "toy")
    bam_path <- Rsamtools::asBam(sam_path, destination = dest_base, overwrite = TRUE)
    bam_path <- paste0(dest_base, ".bam")
    Rsamtools::indexBam(bam_path)
    
    # Count with UMI de-duplication enabled
    mat_dedup <- countReadsSingleCell(
        gr = gr,
        bam = bam_path,
        seqlevels_style = "UCSC",
        tags = list(cell = "CB", umi = "UB"),
        mapq = 10,
        dedup = TRUE,
        BPPARAM = BiocParallel::SerialParam() # portable across OS
    )
    
    expect_s4_class(mat_dedup, "dgCMatrix")
    expect_identical(rownames(mat_dedup), c("orf1", "orf2"))
    expect_true(all(c("cellA", "cellB") %in% colnames(mat_dedup)))
    expect_identical(dim(mat_dedup), c(2L, length(colnames(mat_dedup))))
    
    m1 <- as.matrix(mat_dedup)
    # With dedup=TRUE:
    # - orf1 in cellA is 1 (r1+r2 share CB/UB, collapsed)
    # - orf2 in cellA is 1 (r3)
    # - orf2 in cellB is 1 (r4)
    # - orf1 in cellB is 0 (no overlap)
    expect_equal(m1["orf1", "cellA"], 1)
    expect_equal(m1["orf2", "cellA"], 1)
    expect_equal(m1["orf2", "cellB"], 1)
    expect_equal(m1["orf1", "cellB"], 0)
    
    # Count again with UMI de-duplication disabled (orf1/cellA should be 2)
    mat_nodedup <- countReadsSingleCell(
        gr = gr,
        bam = bam_path,
        seqlevels_style = "UCSC",
        tags = list(cell = "CB", umi = "UB"),
        mapq = 10,
        dedup = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
    m2 <- as.matrix(mat_nodedup)
    expect_equal(m2["orf1", "cellA"], 2)  # r1 and r2 now both counted
    expect_equal(m2["orf2", "cellA"], 1)
    expect_equal(m2["orf2", "cellB"], 1)
})