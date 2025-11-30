library(rtracklayer)
library(GenomicRanges)
library(S4Vectors)
library(withr)

# Helper to create a valid GTF file
create_valid_gtf <- function(path) {
    gtf_lines <- c(
        "chr1\tENSEMBL\ttranscript\t100\t150\t.\t+\t.\thgnc_id \"HGNC:1\"; protein_id \"P1\"; ccdsid \"CCDS1\";",
        "chr1\tENSEMBL\ttranscript\t200\t250\t.\t-\t.\thgnc_id \"HGNC:2\"; protein_id \"P2\"; ccdsid \"CCDS2\";",
        "chr2\tGENCODE\ttranscript\t300\t350\t.\t+\t.\thgnc_id \"HGNC:3\"; protein_id \"P3\"; ccdsid \"CCDS3\";"
    )
    writeLines(gtf_lines, path)
}


test_that("filter_gtf filters by required IDs", {
    gtf_path <- tempfile(fileext = ".gtf")
    create_valid_gtf(gtf_path)
    withr::defer(unlink(gtf_path))  # cleanup
    
    gr <- filter_gtf(gtf_path, verbose = FALSE)
    expect_s4_class(gr, "GRanges")
    expect_true(all(!is.na(mcols(gr)$hgnc_id)))
    expect_true(all(!is.na(mcols(gr)$protein_id)))
    expect_true(all(!is.na(mcols(gr)$ccdsid)))
})


test_that("filter_gtf applies source filter", {
    gtf_path <- tempfile(fileext = ".gtf")
    create_valid_gtf(gtf_path)
    withr::defer(unlink(gtf_path))
    
    gr <- filter_gtf(gtf_path, source_filter = "ENSEMBL", verbose = FALSE)
    
    expect_true(all(mcols(gr)$source == "ENSEMBL"))
})


test_that("filter_gtf returns empty GRanges if no match", {
    gtf_path <- tempfile(fileext = ".gtf")
    create_valid_gtf(gtf_path)
    withr::defer(unlink(gtf_path))
    
    gr <- filter_gtf(gtf_path, source_filter = "NON_EXISTENT", verbose = FALSE)
    expect_equal(length(gr), 0)
})


test_that("filter_gtf warns if required column missing", {
    gtf_path <- tempfile(fileext = ".gtf")
    create_valid_gtf(gtf_path)
    withr::defer(unlink(gtf_path))
    
    gr <- import(gtf_path, feature.type = "transcript")
    mcols(gr)$hgnc_id <- NULL
    export(gr, gtf_path)
    
    expect_warning(filter_gtf(gtf_path, verbose = FALSE), "Column hgnc_id not found")
})
